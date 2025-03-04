from common.io import read_fasta, load_score_matrix
from linear.alignment import global_linear_alignment
from affine.alignment import global_affine_alignment
import pandas as pd

# read sequences
all_sequences = {f"seq{i}": read_fasta(f"experiments/eval_seq{i}.fasta") for i in range(1, 6)}

# read score matrix
score_matrix = load_score_matrix("experiments\eval_score_matrix.txt")

### ANSWER THE QUESTIONS ###
# question 1
q1_score, q1_alignment = global_linear_alignment(all_sequences['seq1'], all_sequences['seq2'], score_matrix, gap_penalty=5)
print(f"\nAnswer to question 1: {q1_score}")

# question 2
q2_score, q2_alignment = global_affine_alignment(all_sequences['seq1'], all_sequences['seq2'], score_matrix, gap_extend=5, gap_open=5)
print(f"\nAnswer to question 2: {q2_score}")

# question 3
score_table_q3 = {seq:{} for seq in all_sequences}

for seqi in all_sequences:
    for seqj in all_sequences:
        score, alignment = global_linear_alignment(all_sequences[seqi], all_sequences[seqj], score_matrix, gap_penalty=5)
        score_table_q3[seqi][seqj] = score

# for printing prettier
st3 = pd.DataFrame.from_dict(score_table_q3)
print("\nAnswer to question 3:")
print(st3)
    

# question 4
score_table_q4 = {seq:{} for seq in all_sequences}

for seqi in all_sequences:
    for seqj in all_sequences:
        score, alignment = global_affine_alignment(all_sequences[seqi], all_sequences[seqj], score_matrix, gap_extend=5, gap_open=5)
        score_table_q4[seqi][seqj] = score

# for printing prettier
st4 = pd.DataFrame.from_dict(score_table_q4)
print("\nAnswer to question 4:")
print(st4)