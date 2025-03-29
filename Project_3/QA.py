#!/usr/bin/env python3
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import time
import tracemalloc
from rich import print
from rich.console import Console
from rich.prompt import Prompt
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from rich.rule import Rule
from rich import box

# Import our alignment modules (assumed implemented)
from sp_exact_3 import sp_exact_3_align, build_full_cost_matrix, default_subst  # Returns (aln1, aln2, aln3, score)
from sp_approx_gpu import sp_approx_main     # Returns (score, headers, final_alignment, center_header, center_seq)

console = Console()
# Build cost matrix using sp_exact_3's method.
cost_matrix = build_full_cost_matrix(default_subst, gap_penalty=5)

# ------------------------------------------------------------
# Utility: Sanitize sequence (replace ambiguous bases with 'A')
# ------------------------------------------------------------
def sanitize_sequence(seq):
    # Replace any character not in A, C, G, T with 'A'
    return "".join(c if c in "ACGT" else "A" for c in seq.upper())

# ------------------------------------------------------------
# Function to read FASTA file
# Returns a list of (header, sanitized sequence) tuples.
# ------------------------------------------------------------
def read_fasta(filepath):
    if not os.path.exists(filepath):
        console.print(f"[red]Error: File '{filepath}' does not exist.[/red]")
        sys.exit(1)
    sequences = []
    header = None
    seq_lines = []
    try:
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line: 
                    continue
                if line.startswith(">"):
                    if header and seq_lines:
                        sequences.append((header, sanitize_sequence("".join(seq_lines))))
                        seq_lines = []
                    header = line.strip()
                else:
                    seq_lines.append(line)
            if header and seq_lines:
                sequences.append((header, sanitize_sequence("".join(seq_lines))))
    except Exception as e:
        console.print(f"[red]Error reading file: {e}[/red]")
        sys.exit(1)
    return sequences

# ------------------------------------------------------------
# Function to run sp_exact_3 on first 3 sequences of brca1-testseqs.fasta
# ------------------------------------------------------------
def run_exact_on_brca1():
    fasta_file = "brca1-testseqs.fasta"
    sequences = read_fasta(fasta_file)
    if len(sequences) < 3:
        console.print("[red]Error: Less than 3 sequences found in brca1-testseqs.fasta[/red]")
        sys.exit(1)
    selected = sequences[:3]
    headers = [h for h, s in selected]
    seqs = [s for h, s in selected]
    
    tracemalloc.start()
    start_time = time.time()
    aln1, aln2, aln3, score = sp_exact_3_align(seqs[0], seqs[1], seqs[2], cost_matrix)
    end_time = time.time()
    current_mem, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    exec_time = end_time - start_time
    mem_usage_mb = peak_mem / (1024 * 1024)
    
    result_panel = Panel.fit(
        f"[bold green]sp_exact_3 took {exec_time:.3f} seconds[/bold green]\n"
        f"[bold green]Peak Memory: {mem_usage_mb:.3f} MB[/bold green]\n\n"
        f"[bold blue]Exact Alignment Score:[/bold blue] {score}\n\n"
        f"[bold blue]{headers[0]}:[/bold blue]\n{aln1}\n\n"
        f"[bold blue]{headers[1]}:[/bold blue]\n{aln2}\n\n"
        f"[bold blue]{headers[2]}:[/bold blue]\n{aln3}",
        title="[bold underline]Exact Alignment (sp_exact_3)[/bold underline]",
        border_style="bright_blue", box=box.HEAVY)
    console.print(result_panel)
    return score, (aln1, aln2, aln3), exec_time, mem_usage_mb

# ------------------------------------------------------------
# Function to run sp_approx on first 5 sequences of brca1-testseqs.fasta
# ------------------------------------------------------------
def run_approx_on_brca1():
    fasta_file = "brca1-testseqs.fasta"
    sequences = read_fasta(fasta_file)
    if len(sequences) < 5:
        console.print("[red]Error: Less than 5 sequences found in brca1-testseqs.fasta[/red]")
        sys.exit(1)
    selected = sequences[:5]
    headers = [h for h, s in selected]
    seqs = [s for h, s in selected]
    
    tracemalloc.start()
    start_time = time.time()
    score, approx_headers, final_alignment, center_header, center_seq = sp_approx_main(headers, seqs)
    end_time = time.time()
    current_mem, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    exec_time = end_time - start_time
    mem_usage_mb = peak_mem / (1024 * 1024)
    
    panel_content = (f"[bold green]sp_approx took {exec_time:.3f} seconds[/bold green]\n"
                     f"[bold green]Peak Memory: {mem_usage_mb:.3f} MB[/bold green]\n\n"
                     f"[bold blue]Approximate Alignment Score:[/bold blue] {score}\n"
                     f"[bold blue]Center Sequence Chosen:[/bold blue] {center_header}\n\n")
    for h, aln in zip(approx_headers, final_alignment):
        panel_content += f"[bold blue]{h}:[/bold blue]\n{aln}\n\n"
    
    result_panel = Panel.fit(panel_content, title="[bold underline]Approximate Alignment (sp_approx)[/bold underline]",
                               border_style="bright_magenta", box=box.HEAVY)
    console.print(result_panel)
    return score, center_header, exec_time, mem_usage_mb

# ------------------------------------------------------------
# Experiment: Compare sp_exact_3 and sp_approx for 3 sequences (from testseqs files)
# ------------------------------------------------------------
def experiment():
    results = []  # (Seq Length, Exact Score, Approx Score, Ratio, Time_exact, Mem_exact, Time_approx, Mem_approx)
    
    for L in range(10, 210, 10):
        filename = f"testseqs/testseqs_{L}_3.fasta"
        if not os.path.exists(filename):
            console.print(f"[red]File {filename} not found.[/red]")
            continue
        sequences = read_fasta(filename)
        if len(sequences) < 3:
            continue
        seqs = [s for h, s in sequences[:3]]
        # Exact alignment measurement
        try:
            tracemalloc.start()
            start_time = time.time()
            _, _, _, exact_score = sp_exact_3_align(seqs[0], seqs[1], seqs[2], cost_matrix)
            end_time = time.time()
            current_mem, peak_mem = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            time_exact = end_time - start_time
            mem_exact = peak_mem / (1024 * 1024)
        except Exception as e:
            console.print(f"[red]Error computing exact alignment for length {L}: {e}[/red]")
            continue
        
        # Approximate alignment measurement
        try:
            tracemalloc.start()
            start_time = time.time()
            approx_result = sp_approx_main([">s1", ">s2", ">s3"], seqs)
            end_time = time.time()
            current_mem, peak_mem = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            time_approx = end_time - start_time
            mem_approx = peak_mem / (1024 * 1024)
            approx_score = approx_result[0]
        except Exception as e:
            console.print(f"[red]Error computing approximate alignment for length {L}: {e}[/red]")
            continue
        
        ratio = approx_score / exact_score if exact_score != 0 else 0
        results.append((L, exact_score, approx_score, ratio, time_exact, mem_exact, time_approx, mem_approx))
        console.print(f"[cyan]Length: {L}, Exact: {exact_score}, Approx: {approx_score}, Ratio: {ratio:.3f}[/cyan]")
    
    # Build and display a comparison table using Rich.
    table = Table(title="Comparison: Exact vs. Approximate Alignment", header_style="bold magenta", box=box.DOUBLE_EDGE)
    table.add_column("Seq Length", justify="center")
    table.add_column("Exact Score", justify="center")
    table.add_column("Approx Score", justify="center")
    table.add_column("Ratio", justify="center")
    table.add_column("Time (Exact)", justify="center")
    table.add_column("Mem (Exact) MB", justify="center")
    table.add_column("Time (Approx)", justify="center")
    table.add_column("Mem (Approx) MB", justify="center")
    
    for res in results:
        L, ex, ap, r, t_ex, m_ex, t_ap, m_ap = res
        table.add_row(str(L), str(ex), str(ap), f"{r:.3f}",
                      f"{t_ex:.3f}", f"{m_ex:.3f}", f"{t_ap:.3f}", f"{m_ap:.3f}")
    
    console.print(table)
    
    # Plot the approximation ratio vs sequence length.
    lengths = [r[0] for r in results]
    ratios = [r[3] for r in results]
    plt.figure(figsize=(8, 5))
    plt.plot(lengths, ratios, marker="o", linestyle="-", color="darkorange")
    plt.xlabel("Sequence Length", fontsize=12)
    plt.ylabel("Approx Score / Exact Score", fontsize=12)
    plt.title("Approximation Ratio vs. Sequence Length (3 sequences)", fontsize=14)
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.show()

# ------------------------------------------------------------
# Main Program: Answer Questions and Run Experiment
# ------------------------------------------------------------
def main():
    splash = Panel.fit(
        "[bold bright_yellow]Welcome to the MSA QA Program![/bold bright_yellow]\n"
        "[bold cyan]Aligning sequences using sp_exact_3 and sp_approx methods[/bold cyan]",
        title="[bold underline]MSA Project QA[/bold underline]", 
        border_style="bright_yellow", box=box.ROUNDED)
    console.print(splash, justify="center")
    console.print(Rule(style="bright_yellow"))
    
    console.print("[bold underline]Question 1: Optimal Alignment (sp_exact_3)[/bold underline]\n")
    run_exact_on_brca1()
    
    console.print(Rule(style="bright_cyan"))
    
    console.print("[bold underline]Question 2: Approximate Alignment (sp_approx)[/bold underline]\n")
    run_approx_on_brca1()
    
    console.print(Rule(style="bright_magenta"))
    
    console.print("[bold underline]Question 3: Experiment to Validate Approximation Ratio[/bold underline]\n")
    experiment()

if __name__ == "__main__":
    main()
