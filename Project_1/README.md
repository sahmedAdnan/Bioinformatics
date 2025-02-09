# Project 1: Pairwise Sequence Alignment

## Description
This project implements **global pairwise sequence alignment** using **dynamic programming** and **Biopython**. It computes the **optimal alignment score** and **all possible optimal alignments** for two given sequences based on a predefined substitution matrix and gap cost.

## Features
- **Dynamic Programming Table** for alignment computation.
- **Backtracking** to find **all possible optimal alignments**.
- **Biopython Implementation** for sequence alignment scoring.
- **Visualization** of the DP table and alignment paths.
- **FASTA File Support** to process external sequence data.

## Requirements
Ensure you have the following dependencies installed:

```bash
pip install biopython numpy matplotlib nbconvert[webpdf] playwright
```

Additionally, install Playwright for `nbconvert` PDF conversion:

```bash
playwright install
```

## Usage
Run the Jupyter Notebook to:
1. Compute the **dynamic programming table**.
2. Find **all optimal alignments**.
3. Generate **visual representations** of the alignment process.

To convert the notebook to PDF:
```bash
jupyter nbconvert --to webpdf Project_1.ipynb
```

## Project Structure
```
Project_1/
│── data/
│   ├── seq1.fasta
│   ├── seq2.fasta
│── Project_1.ipynb
│── README.md
│── results/
│   ├── dp_table.png
│   ├── alignments.txt
```

## Access All Projects
You can find all project folders here: [Project Repository](your_github_repo_link)
