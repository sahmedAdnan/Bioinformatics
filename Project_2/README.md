# DNA Sequence Alignment Tool

This project provides a suite of tools for performing global DNA sequence alignments using both linear and affine gap cost models. It is designed for both interactive use and batch evaluation, and it includes performance analysis capabilities.

## Features

- **Global Linear Alignment:** Implements the Needleman-Wunsch algorithm using a linear gap penalty model (located in the `linear` folder).
- **Global Affine Alignment:** Implements an affine gap penalty model with separate gap opening and extension penalties (located in the `affine` folder).
- **Interactive Command-Line Interface:** Users can input sequences manually or load sequences from FASTA files.
- **Evaluation Module:** Computes alignment scores between multiple sequences and presents the results as a score matrix.
- **Performance Analysis:** Generates scatter plots to analyze and compare the running times of the alignment algorithms.

**Installation** 
Clone the Repository:

```bash
git clone https://github.com/yourusername/DNA-Sequence-Alignment-Tool.git
```

## Usage
**Interactive Alignment (main.py)**
**Run the command-line tool for interactive DNA sequence alignment:**

```bash
python main.py
```
- **Sequence Input:** Choose to input sequences manually or via FASTA files.
- **Substitution Matrix:** The tool automatically loads the substitution matrix from experiments/eval_score_matrix.txt.
- **Gap Cost Methods:** Select between global linear or affine alignment.
- **Display and Save Alignments:** View possible alignments and optionally save one in FASTA format.

## Batch Evaluation (evaluation.py)
**The evaluation script computes alignment scores for multiple sequences:**

```bash
python evaluation.py
```
- Reads sequences from FASTA files (e.g., `experiments/eval_seq1.fasta` to `eval_seq5.fasta`).
- Calculates scores using both linear and affine alignment methods.
- Outputs score matrices for further analysis.
- Performance Analysis (illustrations.py)

## To analyze the running time performance of the alignment algorithms, run:

```bash
python illustrations.py <max_sequence_length> <number_of_replications>
```
**For example:**
```bash
python illustrations.py 50 10
```

- **Sequence Lengths:** The script generates sequences with lengths ranging from a minimum (e.g., 10) to the specified maximum.
- **Replications:** The experiments are repeated a specified number of times.
- **Plots:** Two scatter plots are generated to compare the running times for linear and affine alignments.
- **Sequence Lengths:** The script generates sequences with lengths ranging from a minimum (e.g., 10) to the specified maximum.
Replications: Experiments are repeated for a specified number of replications.
Plots: Two scatter plots are displayed comparing the running times for linear and affine alignments.
Contributing
Contributions to improve or extend the functionality of this project are welcome! Please submit issues or pull requests with detailed descriptions of your changes.
