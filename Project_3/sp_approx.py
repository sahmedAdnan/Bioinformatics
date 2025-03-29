#!/usr/bin/env python3
import os
import sys
import numpy as np
from rich import print
from rich.console import Console
from rich.prompt import Prompt, Confirm
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich import box
from concurrent.futures import ProcessPoolExecutor, as_completed

# Import the test function for verifying the alignment
from msa_sp_score_3k import compute_sp_score

console = Console()

# Global gap cost and score matrix definitions.
GAP_COST = 5
score_matrix = {
    'A': {'A':0, 'C':5, 'G':2, 'T':5, '-':5},
    'C': {'A':5, 'C':0, 'G':5, 'T':2, '-':5},
    'G': {'A':2, 'C':5, 'G':0, 'T':5, '-':5},
    'T': {'A':5, 'C':2, 'G':5, 'T':0, '-':5},
    '-': {'A':5, 'C':5, 'G':5, 'T':5, '-':0}
}

# Utility: Sanitize sequence (replace ambiguous nucleotides with 'A')
def sanitize_sequence(seq):
    sanitized = "".join(c if c in "ACGT" else "A" for c in seq.upper())
    if sanitized != seq.upper():
        console.print("[yellow]Warning: Ambiguous nucleotide(s) found; substituting with 'A'.[/yellow]")
    return sanitized

# --------------------------------------------------
# Exceptions and Safe Input
# --------------------------------------------------
class BackException(Exception):
    """Raised when the user types 'esc' to go back."""
    pass

def safe_prompt(prompt_text):
    response = Prompt.ask(prompt_text)
    lower = response.strip().lower()
    if lower == "esc":
        raise BackException()
    if lower in ("terminate", "quit"):
        console.print("[bold red]Terminating the program...[/bold red]")
        sys.exit(0)
    return response

# --------------------------------------------------
# FASTA File Handling Functions
# --------------------------------------------------
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
                        sequences.append((header, sanitize_sequence("".join(seq_lines))))
                        seq_lines = []
                    header = line[1:].strip()
                else:
                    seq_lines.append(line)
            if header and seq_lines:
                sequences.append((header, sanitize_sequence("".join(seq_lines))))
    except Exception as e:
        console.print(f"[red]Error reading file: {e}[/red]")
        return None
    if not sequences:
        console.print("[red]No sequences were found in the file.[/red]")
        return None
    return sequences

def display_fasta_overview(sequences):
    table = Table(title="Available Sequences", header_style="bold blue", box=box.DOUBLE_EDGE)
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
            #display_fasta_overview(sequences)
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

def get_sequences_manually():
    console.print("\n[bold yellow]Enter sequences manually (one per line). Type 'done' when finished:[/bold yellow]")
    seq_list = []
    while True:
        try:
            inp = safe_prompt("[green]Enter sequence (or 'done' to finish):[/green]").strip()
            if inp.lower() == "done":
                break
            if not inp:
                console.print("[red]Empty sequence. Please enter a valid sequence.[/red]")
                continue
            seq_list.append(sanitize_sequence(inp))
        except BackException:
            console.print("[yellow]Returning to previous menu.[/yellow]")
            raise
    if len(seq_list) < 2:
        console.print("[red]At least two sequences are required for MSA.[/red]")
        return None
    headers = [f"Sequence {i+1}" for i in range(len(seq_list))]
    return list(zip(headers, seq_list))

# --------------------------------------------------
# Pairwise Alignment Functions (Needleman–Wunsch style)
# --------------------------------------------------
def pairwise_score(s1, s2):
    m, n = len(s1), len(s2)
    dp = [[0]*(n+1) for _ in range(m+1)]
    for i in range(m+1):
        dp[i][0] = i * GAP_COST
    for j in range(n+1):
        dp[0][j] = j * GAP_COST
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = dp[i-1][j-1] + score_matrix[s1[i-1]][s2[j-1]]
            delete = dp[i-1][j] + GAP_COST
            insert = dp[i][j-1] + GAP_COST
            dp[i][j] = min(match, delete, insert)
    return dp[m][n]

def pairwise_align(ref, seq):
    m, n = len(ref), len(seq)
    dp = [[0]*(n+1) for _ in range(m+1)]
    for i in range(m+1):
        dp[i][0] = i * GAP_COST
    for j in range(n+1):
        dp[0][j] = j * GAP_COST
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = dp[i-1][j-1] + score_matrix[ref[i-1]][seq[j-1]]
            delete = dp[i-1][j] + GAP_COST
            insert = dp[i][j-1] + GAP_COST
            dp[i][j] = min(match, delete, insert)
    i, j = m, n
    ref_aligned, seq_aligned = [], []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + score_matrix[ref[i-1]][seq[j-1]]:
            ref_aligned.append(ref[i-1])
            seq_aligned.append(seq[j-1])
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j] + GAP_COST:
            ref_aligned.append(ref[i-1])
            seq_aligned.append('-')
            i -= 1
        else:
            ref_aligned.append('-')
            seq_aligned.append(seq[j-1])
            j -= 1
    return (''.join(reversed(ref_aligned)), ''.join(reversed(seq_aligned)))

def merge_alignments(a1, a2):
    merged = []
    i = j = 0
    while i < len(a1) or j < len(a2):
        if i < len(a1) and j < len(a2) and a1[i] == a2[j]:
            merged.append(a1[i])
            i += 1
            j += 1
        else:
            merged.append('-')
            if i < len(a1) and a1[i] == '-':
                i += 1
            elif j < len(a2) and a2[j] == '-':
                j += 1
            else:
                if i < len(a1):
                    merged[-1] = a1[i]
                    i += 1
                elif j < len(a2):
                    merged[-1] = a2[j]
                    j += 1
    return ''.join(merged)

# --------------------------------------------------
# Parallelization of Pairwise Score Computation
# --------------------------------------------------
from concurrent.futures import ProcessPoolExecutor, as_completed

def compute_pairwise_score_for_pair(i, j, sequences):
    score = pairwise_score(sequences[i], sequences[j])
    return (i, j, score)

def parallel_pairwise_scores(sequences):
    n = len(sequences)
    scores = np.zeros((n, n), dtype=np.int32)
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(compute_pairwise_score_for_pair, i, j, sequences): (i, j) for i, j in pairs}
        for future in as_completed(futures):
            i, j, sc = future.result()
            scores[i, j] = sc
            scores[j, i] = sc
    return scores

# --------------------------------------------------
# sp_approx Implementation (Center Method) with Parallelization
# --------------------------------------------------
def sp_approx_main(headers, sequences):
    if len(sequences) < 2:
        console.print("[red]At least two sequences are required for MSA.[/red]")
        sys.exit(1)
    n = len(sequences)
    # Compute pairwise scores in parallel.
    scores_matrix = parallel_pairwise_scores(sequences)
    total_scores = [scores_matrix[i].sum() for i in range(n)]
    center_idx = int(np.argmin(total_scores))
    center_header = headers[center_idx] if headers else f"Sequence {center_idx+1}"
    center_seq = sequences[center_idx]
    console.print(f"[bold green]Center selected:[/bold green] {center_header}")

    # Align every sequence to the center.
    aligned_pairs = []
    for i in range(n):
        if i == center_idx:
            aligned_pairs.append((center_seq, center_seq))
        else:
            try:
                aligned_ref, aligned_seq = pairwise_align(center_seq, sequences[i])
                aligned_pairs.append((aligned_ref, aligned_seq))
            except Exception as e:
                console.print(f"[red]Error aligning sequence {i}: {e}[/red]")
                sys.exit(1)
    
    # Build consensus by merging all aligned center strings.
    consensus = aligned_pairs[0][0]
    for i in range(1, len(aligned_pairs)):
        consensus = merge_alignments(consensus, aligned_pairs[i][0])
    
    # Build final multiple alignment by adjusting each alignment to match the consensus.
    final_alignment = []
    for (aligned_center, aligned_seq) in aligned_pairs:
        adjusted = ""
        ptr = 0
        for c in consensus:
            if ptr < len(aligned_center) and aligned_center[ptr] == c:
                adjusted += aligned_seq[ptr]
                ptr += 1
            else:
                adjusted += '-'
        final_alignment.append(adjusted)
    
    # Compute SP score.
    sp_score = 0
    for col in zip(*final_alignment):
        for i in range(len(col)):
            for j in range(i+1, len(col)):
                sp_score += score_matrix[col[i]][col[j]]
    return sp_score, headers, final_alignment, center_header, center_seq

# --------------------------------------------------
# Saving Alignment as FASTA (Multiple Sequences)
# --------------------------------------------------
def save_alignment_fasta_multiple(alignments, headers, filename):
    if not filename.lower().endswith(".fasta"):
        filename += ".fasta"
    try:
        with open(filename, "w") as f:
            for i, aln in enumerate(alignments):
                head = headers[i] if headers else f"seq{i+1}"
                f.write(f">{head}\n")
                f.write(f"{aln}\n")
        console.print(f"[cyan]Alignment successfully saved to [bold]{filename}[/bold].[/cyan]")
        return filename
    except Exception as e:
        console.print(f"[red]Error saving file: {e}[/red]")
        return None

# --------------------------------------------------
# Main Interactive Function for sp_approx
# --------------------------------------------------
def main():
    try:
        console.clear()
    except Exception as e:
        print(f"Error clearing console: {e}")
    
    header_text = Text("SP_APPROX\n", style="bold white on dark_green")
    header_text.append("2-Approximation MSA using Center-Star Method", style="bold bright_yellow")
    header_panel = Panel(header_text, box=box.DOUBLE_EDGE, border_style="bright_green", padding=(1,4), expand=True)
    console.print(header_panel, justify="center")
    console.print("[italic]At any prompt, type 'esc' to go back or 'terminate' to exit.[/italic]\n", style="dim", justify="center")
    
    console.print("[bold yellow]Sequence Input Options:[/bold yellow]")
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
            sys.exit(1)
    
    headers = None
    sequences = []
    if mode == "1":
        manual = get_sequences_manually()
        if manual is None:
            sys.exit(1)
        headers, seq_list = zip(*manual)
        headers = list(headers)
        sequences = list(seq_list)
    else:
        while True:
            try:
                fasta_path = safe_prompt("[green]Enter the path to the FASTA file:[/green]").strip()
                if not fasta_path:
                    console.print("[red]File path cannot be empty.[/red]")
                    continue
                fasta_seqs = read_fasta_file(fasta_path)
                if not fasta_seqs:
                    console.print("[red]No sequences found. Please try another file.[/red]")
                    continue
                # Display overview once.
                display_fasta_overview(fasta_seqs)
                break
            except BackException:
                console.print("[yellow]Returning to input mode selection.[/yellow]")
                sys.exit(1)
        if len(fasta_seqs) < 2:
            console.print("[red]At least two sequences are required for MSA.[/red]")
            sys.exit(1)
        elif len(fasta_seqs) == 2:
            headers, seqs = zip(*fasta_seqs)
            headers = list(headers)
            sequences = list(seqs)
        else:
            console.print(f"\n[cyan]{len(fasta_seqs)} sequences found in the FASTA file.[/cyan]")
            valid = False
            while not valid:
                try:
                    count_choice = safe_prompt("[green]Enter the number of sequences to align (min 2):[/green]").strip()
                    count = int(count_choice)
                    if count < 2 or count > len(fasta_seqs):
                        console.print("[red]Invalid number of sequences. Please try again.[/red]")
                        continue
                    valid = True
                except ValueError:
                    console.print("[red]Invalid input. Please enter a numeric value.[/red]")
            try:
                chosen = choose_n_sequences(fasta_seqs, count)
                selected = [(h, s) for h, s in fasta_seqs if s in set(chosen)]
                if len(selected) != count:
                    console.print("[red]Error in selection. Exiting.[/red]")
                    sys.exit(1)
                headers, seqs = zip(*selected)
                headers = list(headers)
                sequences = list(seqs)
            except (ValueError, BackException):
                console.print("[yellow]Returning to FASTA file selection menu.[/yellow]")
                sys.exit(1)
    
    sp_score, headers, final_alignment, center_header, center_seq = sp_approx_main(headers, sequences)
    
    result_text = Text()
    result_text.append("SP Score: ", style="bold bright_magenta")
    result_text.append(f"{sp_score}\n\n", style="bold bright_magenta")
    for i, aln in enumerate(final_alignment):
        label = headers[i] if headers else f"Sequence {i+1}"
        if headers and headers[i] == center_header:
            label += " (Center)"
        else:
            label += " (Aligned)"
        result_text.append(label + ":\n", style="bold green")
        result_text.append(aln + "\n", style="monospace")
    result_panel = Panel(result_text, title="[bold green]Approximate MSA Result[/bold green]",
                           border_style="bright_green", box=box.HEAVY, padding=(1,2))
    console.print(result_panel)
    console.rule("[bold bright_magenta]End of Alignment[/bold bright_magenta]")
    
    saved_filename = None
    try:
        save_choice = safe_prompt("[green]Would you like to save the alignment to a FASTA file? (yes/no):[/green]").strip().lower()
        if save_choice in ("yes", "y"):
            while True:
                try:
                    fname = safe_prompt("[green]Enter the file name to save the alignment:[/green]").strip()
                    if not fname:
                        console.print("[red]File name cannot be empty. Please try again.[/red]")
                        continue
                    saved_filename = save_alignment_fasta_multiple(final_alignment, headers, fname)
                    break
                except BackException:
                    console.print("[yellow]Returning to file name prompt.[/yellow]")
    except BackException:
        console.print("[yellow]Skipping save alignment step.[/yellow]")
    
    try:
        verify_choice = safe_prompt("[green]Would you like to verify the alignment using the test program? (yes/no):[/green]").strip().lower()
        if verify_choice in ("yes", "y"):
            if saved_filename is None:
                console.print("[red]No alignment file was saved. Unable to verify alignment.[/red]")
            else:
                test_score = compute_sp_score(saved_filename)
                if test_score == sp_score:
                    console.print("[bold green]Alignment is verified and passed![/bold green]")
                else:
                    console.print(f"[bold red]Alignment verification failed! Test score: {test_score}, Expected: {sp_score}[/bold red]")
    except BackException:
        console.print("[yellow]Skipping alignment verification step.[/yellow]")

if __name__ == "__main__":
    main()
