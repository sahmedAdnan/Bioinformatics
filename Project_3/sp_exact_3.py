#!/usr/bin/env python3
import os
import sys
import numpy as np
import numba
from rich import print
from rich.console import Console
from rich.prompt import Prompt, Confirm
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich import box

# Import the test function from msa_sp_score_3k.py
from msa_sp_score_3k import compute_sp_score

console = Console()

# ------------------------------------------------------------------
# Sanitize sequence (replace all non-ACGT characters with 'A')
# ------------------------------------------------------------------
def sanitize_sequence(seq):
    seq = seq.upper()
    sanitized = "".join(ch if ch in "ACGT" else "A" for ch in seq)
    if sanitized != seq:
        console.print("[yellow]Warning: Ambiguous nucleotide(s) found; substituting with 'A'.[/yellow]")
    return sanitized

# ---------------------------
# Exceptions and Safe Input
# ---------------------------
class BackException(Exception):
    """Raised when the user wants to go back (by typing 'esc')."""
    pass

def safe_prompt(prompt_text):
    response = Prompt.ask(prompt_text)
    if response.strip().lower() == "esc":
        raise BackException()
    if response.strip().lower() in ("terminate", "quit"):
        console.print("[bold red]Terminating the program...[/bold red]")
        sys.exit(0)
    return response

# ---------------------------
# FASTA File Handling Functions
# ---------------------------
def read_fasta_file(filepath):
    if not os.path.isfile(filepath):
        console.print(f"[red]Error: File '{filepath}' does not exist.[/red]")
        return None

    sequences = []
    header = None
    seq_lines = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if header and seq_lines:
                        raw_seq = "".join(seq_lines)
                        sanitized = sanitize_sequence(raw_seq)
                        sequences.append((header, sanitized))
                        seq_lines = []
                    header = line[1:].strip()
                else:
                    seq_lines.append(line)
            if header and seq_lines:
                raw_seq = "".join(seq_lines)
                sanitized = sanitize_sequence(raw_seq)
                sequences.append((header, sanitized))
    except Exception as e:
        console.print(f"[red]Error reading file: {e}[/red]")
        return None

    if not sequences:
        console.print("[red]No sequences were found in the file.[/red]")
        return None

    return sequences

def display_fasta_options(sequences):
    table = Table(title="Available Sequences", header_style="bold blue")
    table.add_column("Index", justify="center")
    table.add_column("Header", justify="left")
    table.add_column("Sequence Preview", justify="left")
    for idx, (header, seq) in enumerate(sequences, start=1):
        preview = seq[:30] + ("..." if len(seq) > 30 else "")
        table.add_row(str(idx), header, preview)
    console.print(table)

def choose_n_sequences(sequences, n):
    while True:
        try:
            display_fasta_options(sequences)
            choice = safe_prompt(f"[green]Enter {n} comma-separated indices (or type 'esc' to go back):[/green]")
            indices = [int(x.strip()) for x in choice.split(",")]
            if len(indices) != n:
                console.print(f"[red]Please enter exactly {n} indices.[/red]")
                continue
            if any(idx < 1 or idx > len(sequences) for idx in indices):
                console.print("[red]One or more indices are out of range. Please try again.[/red]")
                continue
            return [sequences[idx-1][1] for idx in indices]
        except ValueError:
            console.print("[red]Invalid input. Please enter integers separated by commas.[/red]")
        except BackException:
            raise

def get_manual_sequence(prompt_text):
    while True:
        try:
            seq = safe_prompt(prompt_text).strip()
            if not seq:
                console.print("[red]Sequence cannot be empty. Please try again.[/red]")
                continue
            seq = sanitize_sequence(seq)  # Convert invalid characters to 'A'
            return seq
        except BackException:
            raise

def get_missing_sequences(missing_count):
    additional = []
    while missing_count > 0:
        try:
            use_file = Confirm.ask(f"[bold green]Would you like to load {missing_count} missing sequence(s) from a FASTA file?[/bold green]", default=True)
            if use_file:
                fasta_file = safe_prompt("[green]Enter the path to the FASTA file:[/green]").strip()
                new_seqs = read_fasta_file(fasta_file)
                if new_seqs is None or len(new_seqs) == 0:
                    console.print("[red]No sequences found. Please try again.[/red]")
                    continue
                if len(new_seqs) < missing_count:
                    console.print(f"[yellow]Only {len(new_seqs)} sequence(s) found. They will be used, and you'll be asked to enter the remaining {missing_count - len(new_seqs)} sequence(s) manually.[/yellow]")
                    additional.extend([t[1] for t in new_seqs])
                    missing_count -= len(new_seqs)
                elif len(new_seqs) == missing_count:
                    additional.extend([t[1] for t in new_seqs])
                    missing_count = 0
                else:
                    chosen = choose_n_sequences(new_seqs, missing_count)
                    additional.extend(chosen)
                    missing_count = 0
            else:
                for i in range(missing_count):
                    seq = get_manual_sequence(f"Enter missing sequence {i+1}:")
                    additional.append(seq)
                missing_count = 0
        except BackException:
            console.print("[yellow]Returning to previous step for missing sequences.[/yellow]")
            raise
    return additional

# ---------------------------
# Substitution Matrix Functions
# ---------------------------
default_subst = np.array([
    [0, 5, 2, 5],
    [5, 0, 5, 2],
    [2, 5, 0, 5],
    [5, 2, 5, 0]
], dtype=np.int32)

def display_substitution_matrix(matrix):
    table = Table(title="Substitution Matrix", header_style="bold blue", box=box.DOUBLE_EDGE)
    table.add_column("", justify="center")
    for col in ["A", "C", "G", "T"]:
        table.add_column(col, justify="center")
    rows = ["A", "C", "G", "T"]
    for i, row_label in enumerate(rows):
        row_values = [str(matrix[i, j]) for j in range(4)]
        table.add_row(row_label, *row_values)
    console.print(table)

def show_format_example():
    example = "0 5 2 5\n5 0 5 2\n2 5 0 5\n5 2 5 0"
    panel = Panel(example, title="[bold yellow]Accepted Format (PHYLIP-like)[/bold yellow]", border_style="yellow", box=box.ROUNDED)
    console.print(panel)

def save_default_matrix_file():
    filename = "default_substitution_matrix.txt"
    try:
        with open(filename, "w") as f:
            f.write("0 5 2 5\n5 0 5 2\n2 5 0 5\n5 2 5 0\n")
        console.print(f"[cyan]Default substitution matrix saved to [bold]{filename}[/bold].[/cyan]")
    except Exception as e:
        console.print(f"[red]Failed to save file: {e}[/red]")

def load_substitution_matrix(filename):
    try:
        matrix = np.loadtxt(filename)
        if matrix.shape != (4, 4):
            console.print("[red]Error: The substitution matrix file must be 4x4 (for A, C, G, T).[/red]")
            show_format_example()
            if Confirm.ask("[bold yellow]Would you like to save the default substitution matrix to see the accepted format?[/bold yellow]", default=True):
                save_default_matrix_file()
            return None
        return matrix.astype(np.int32)
    except Exception as e:
        console.print(f"[red]Error loading substitution matrix file: {e}[/red]")
        show_format_example()
        if Confirm.ask("[bold yellow]Would you like to save the default substitution matrix to see the accepted format?[/bold yellow]", default=True):
            save_default_matrix_file()
        return None

def get_substitution_matrix():
    while True:
        try:
            console.print("\n[bold yellow]Substitution Matrix Input Options:[/bold yellow]")
            console.print("1. Use default substitution matrix")
            console.print("2. Load substitution matrix from a file")
            console.print("3. Input substitution matrix manually")
            option = safe_prompt("[green]Choose option (1, 2, or 3):[/green]").strip()
            if option == "1":
                console.print("[cyan]Using default substitution matrix.[/cyan]")
                display_substitution_matrix(default_subst)
                return default_subst
            elif option == "2":
                while True:
                    try:
                        matrix_file = safe_prompt("[green]Enter the path to the 4x4 substitution matrix file:[/green]").strip()
                        if not matrix_file:
                            console.print("[red]You must enter a file path.[/red]")
                            continue
                        matrix = load_substitution_matrix(matrix_file)
                        if matrix is not None:
                            console.print(f"[cyan]Substitution matrix loaded from [bold]{matrix_file}[/bold].[/cyan]")
                            display_substitution_matrix(matrix)
                            return matrix
                        else:
                            if not Confirm.ask("[bold yellow]Do you want to try loading another file?[/bold yellow]", default=True):
                                raise BackException()
                    except BackException:
                        raise
            elif option == "3":
                console.print("[yellow]Please enter the substitution matrix manually.")
                console.print("Enter 4 rows, each with 4 integers separated by spaces.")
                console.print("Example row: [bold]0 5 2 5[/bold]")
                rows = []
                for i in range(4):
                    while True:
                        try:
                            row_str = safe_prompt(f"[green]Enter row {i+1}:[/green]").strip()
                            row = [int(x) for x in row_str.split()]
                            if len(row) != 4:
                                console.print("[red]Each row must have exactly 4 integers.[/red]")
                                continue
                            rows.append(row)
                            break
                        except ValueError:
                            console.print("[red]Invalid input. Please enter integers separated by spaces.[/red]")
                        except BackException:
                            raise
                matrix = np.array(rows, dtype=np.int32)
                console.print("[cyan]Substitution matrix successfully inputted.[/cyan]")
                display_substitution_matrix(matrix)
                return matrix
            else:
                console.print("[red]Invalid option. Please choose 1, 2, or 3.[/red]")
        except BackException:
            console.print("[yellow]Returning to previous menu in substitution matrix selection.[/yellow]")
            continue

def build_full_cost_matrix(subst_matrix, gap_penalty):
    full_matrix = np.empty((5, 5), dtype=np.int32)
    full_matrix[:4, :4] = subst_matrix
    full_matrix[4, :4] = gap_penalty
    full_matrix[:4, 4] = gap_penalty
    full_matrix[4, 4] = 0
    return full_matrix

# ---------------------------
# Alignment Functions
# ---------------------------
def letter_to_index(letter):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    try:
        return mapping[letter]
    except KeyError:
        # If we somehow get an unexpected character, treat it as 'A' by default:
        return 0

def seq_to_array(seq):
    arr = np.empty(len(seq), dtype=np.int32)
    for i, letter in enumerate(seq):
        arr[i] = letter_to_index(letter)
    return arr

# The seven possible moves for the 3D dynamic programming.
moves = np.array([
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 1, 0),
    (1, 0, 1),
    (0, 1, 1),
    (1, 1, 1)
], dtype=np.int32)

INF = 10**9

@numba.njit
def sp_exact_3_dp(seq1, seq2, seq3, cost_matrix):
    l1 = len(seq1)
    l2 = len(seq2)
    l3 = len(seq3)
    dp = np.full((l1+1, l2+1, l3+1), INF, dtype=np.int32)
    back = np.full((l1+1, l2+1, l3+1), -1, dtype=np.int32)
    dp[0, 0, 0] = 0
    for i in range(l1+1):
        for j in range(l2+1):
            for k in range(l3+1):
                if dp[i, j, k] == INF:
                    continue
                for m in range(7):
                    di = moves[m, 0]
                    dj = moves[m, 1]
                    dk = moves[m, 2]
                    ni = i + di
                    nj = j + dj
                    nk = k + dk
                    if ni > l1 or nj > l2 or nk > l3:
                        continue
                    a = seq1[i] if di == 1 else 4
                    b = seq2[j] if dj == 1 else 4
                    c = seq3[k] if dk == 1 else 4
                    col_cost = cost_matrix[a, b] + cost_matrix[a, c] + cost_matrix[b, c]
                    new_score = dp[i, j, k] + col_cost
                    if new_score < dp[ni, nj, nk]:
                        dp[ni, nj, nk] = new_score
                        back[ni, nj, nk] = m
    return dp, back

def sp_exact_3_align(seq1_str, seq2_str, seq3_str, cost_matrix):
    seq1 = seq_to_array(seq1_str)
    seq2 = seq_to_array(seq2_str)
    seq3 = seq_to_array(seq3_str)
    dp, back = sp_exact_3_dp(seq1, seq2, seq3, cost_matrix)
    l1, l2, l3 = len(seq1), len(seq2), len(seq3)
    score = dp[l1, l2, l3]
    aligned1, aligned2, aligned3 = [], [], []
    i, j, k = l1, l2, l3
    while i > 0 or j > 0 or k > 0:
        m = back[i, j, k]
        di = moves[m, 0]
        dj = moves[m, 1]
        dk = moves[m, 2]
        aligned1.append(seq1_str[i-1] if di == 1 else '-')
        aligned2.append(seq2_str[j-1] if dj == 1 else '-')
        aligned3.append(seq3_str[k-1] if dk == 1 else '-')
        i -= di
        j -= dj
        k -= dk
    return ''.join(aligned1[::-1]), ''.join(aligned2[::-1]), ''.join(aligned3[::-1]), score

def save_alignment_fasta(aln1, aln2, aln3, filename):
    if not filename.lower().endswith(".fasta"):
        filename += ".fasta"
    try:
        with open(filename, "w") as f:
            f.write(">seq1\n")
            f.write(aln1 + "\n")
            f.write(">seq2\n")
            f.write(aln2 + "\n")
            f.write(">seq3\n")
            f.write(aln3 + "\n")
        console.print(f"[cyan]Alignment successfully saved to [bold]{filename}[/bold].[/cyan]")
        return filename
    except Exception as e:
        console.print(f"[red]Error saving file: {e}[/red]")
        return None

# ---------------------------
# Main Interactive Function
# ---------------------------
def main():
    console.clear()
    # Create a designer-level header panel
    header = Text("SP_EXACT_3\n", style="bold white on dark_blue")
    header.append("Exact Multiple Sequence Alignment for 3 Sequences", style="bold bright_yellow")
    header_panel = Panel(header, box=box.DOUBLE_EDGE, border_style="bright_blue", padding=(1, 4), expand=True)
    console.print(header_panel, justify="center")
    console.print("[italic]At any prompt, type 'esc' to go back or 'terminate' (or 'quit') to exit.[/italic]\n", style="dim", justify="center")
    
    # Substitution Matrix Input
    while True:
        try:
            subst_matrix = get_substitution_matrix()
            break
        except BackException:
            console.print("[yellow]Returning to previous menu in substitution matrix selection.[/yellow]")
            continue

    # Gap Penalty Input
    while True:
        try:
            gap_input = safe_prompt("[green]Enter the gap penalty (an integer):[/green]").strip()
            gap_penalty = int(gap_input)
            break
        except ValueError:
            console.print("[red]Invalid input. Please enter a valid integer.[/red]")
        except BackException:
            console.print("[yellow]Returning to substitution matrix menu.[/yellow]")
            subst_matrix = get_substitution_matrix()
    cost_matrix = build_full_cost_matrix(subst_matrix, gap_penalty)
    
    # Sequence Input Mode Selection
    console.print("\n[bold yellow]Sequence Input Options:[/bold yellow]")
    console.print("1. Manual Entry")
    console.print("2. Load from FASTA file")
    while True:
        try:
            mode = safe_prompt("[green]Choose input mode (1 or 2):[/green]").strip()
            if mode not in ("1", "2"):
                console.print("[red]Invalid choice. Please enter 1 or 2.[/red]")
                continue
            break
        except BackException:
            console.print("[yellow]Returning to previous menu.[/yellow]")
    
    sequences = []
    if mode == "1":
        console.print("\n[bold cyan]Please enter three sequences manually (any invalid chars will be replaced with 'A'):[/bold cyan]")
        try:
            seq1 = get_manual_sequence("Sequence 1:")
            seq2 = get_manual_sequence("Sequence 2:")
            seq3 = get_manual_sequence("Sequence 3:")
            sequences = [seq1, seq2, seq3]
        except BackException:
            console.print("[yellow]Returning to sequence input mode selection.[/yellow]")
            return main()  # Restart main if back is chosen.
    else:
        while True:
            try:
                fasta_file = safe_prompt("[green]Enter the path to the FASTA file:[/green]").strip()
                if not fasta_file:
                    console.print("[red]You must enter a file path.[/red]")
                    continue
                fasta_sequences = read_fasta_file(fasta_file)
                if fasta_sequences is None or len(fasta_sequences) == 0:
                    console.print("[red]No sequences found. Please try another file.[/red]")
                    continue
                break
            except BackException:
                console.print("[yellow]Returning to sequence input mode selection.[/yellow]")
                return main()
        if len(fasta_sequences) == 3:
            console.print("[cyan]Exactly 3 sequences found. They will be used for alignment.[/cyan]")
            sequences = [t[1] for t in fasta_sequences]
        elif len(fasta_sequences) > 3:
            console.print(f"\n[cyan]{len(fasta_sequences)} sequences found in the FASTA file.[/cyan]")
            try:
                sequences = choose_n_sequences(fasta_sequences, 3)
            except BackException:
                console.print("[yellow]Returning to FASTA file selection menu.[/yellow]")
                return main()
        else:
            console.print(f"[yellow]Only {len(fasta_sequences)} sequence(s) found in the FASTA file.[/yellow]")
            sequences = [t[1] for t in fasta_sequences]
            missing_count = 3 - len(sequences)
            console.print(f"[yellow]You need {missing_count} more sequence(s).[/yellow]")
            try:
                additional = get_missing_sequences(missing_count)
            except BackException:
                console.print("[yellow]Returning to FASTA file selection menu.[/yellow]")
                return main()
            sequences.extend(additional)
    
    if len(sequences) != 3:
        console.print("[red]Error: Exactly three sequences are required for the alignment.[/red]")
        return

    console.print("\n[bold yellow]Computing the optimal alignment using exact dynamic programming...[/bold yellow]\n")
    aln1, aln2, aln3, score = sp_exact_3_align(sequences[0], sequences[1], sequences[2], cost_matrix)
    
    # Alignment Display in a monospaced panel with labels
    aligned_text = Text()
    aligned_text.append("Alignment Score: ", style="bold bright_magenta")
    aligned_text.append(f"{score}\n\n", style="bold bright_magenta")
    aligned_text.append("Sequence 1:\n", style="bold green")
    aligned_text.append(aln1 + "\n", style="monospace")
    aligned_text.append("Sequence 2:\n", style="bold green")
    aligned_text.append(aln2 + "\n", style="monospace")
    aligned_text.append("Sequence 3:\n", style="bold green")
    aligned_text.append(aln3, style="monospace")
    
    result_panel = Panel(
        aligned_text,
        title="[bold green]Optimal Alignment Result[/bold green]",
        border_style="bright_blue",
        box=box.HEAVY,
        padding=(1, 2)
    )
    console.print(result_panel)
    console.rule("[bold bright_magenta]End of Alignment[/bold bright_magenta]")
    
    # Ask the user if they want to save the alignment to a FASTA file.
    try:
        save_choice = safe_prompt("[green]Would you like to save the alignment to a FASTA file? (yes/no):[/green]").strip().lower()
        if save_choice in ("yes", "y"):
            while True:
                try:
                    filename = safe_prompt("[green]Enter the file name to save the alignment:[/green]").strip()
                    if not filename:
                        console.print("[red]File name cannot be empty. Please try again.[/red]")
                        continue
                    saved_filename = save_alignment_fasta(aln1, aln2, aln3, filename)
                    break
                except BackException:
                    console.print("[yellow]Returning to file name prompt.[/yellow]")
    except BackException:
        console.print("[yellow]Skipping save alignment step.[/yellow]")
        saved_filename = None
    
    # Ask the user if they want to verify the alignment using the test function.
    try:
        verify_choice = safe_prompt("[green]Would you like to verify the alignment using the test program? (yes/no):[/green]").strip().lower()
        if verify_choice in ("yes", "y"):
            if not saved_filename:
                console.print("[red]No alignment file was saved. Unable to verify alignment.[/red]")
            else:
                test_score = compute_sp_score(saved_filename)
                if test_score == score:
                    console.print("[bold green]Alignment is verified and passed![/bold green]")
                else:
                    console.print(f"[bold red]Alignment verification failed! Test score: {test_score}, Expected: {score}[/bold red]")
    except BackException:
        console.print("[yellow]Skipping alignment verification step.[/yellow]")

if __name__ == '__main__':
    main()
