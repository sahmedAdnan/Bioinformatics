# run_experiments.py

import time
import random
import matplotlib.pyplot as plt
from linear.alignment import global_linear_alignment
from affine.alignment import global_affine_alignment
from common.io import load_score_matrix

def generate_random_sequence(length):
    return ''.join(random.choice("ACGT") for _ in range(length))

def main():
    score_matrix_file = "eval_score_matrix.txt"
    score_matrix = load_score_matrix(score_matrix_file)

    gap_linear = 5   # for linear
    gap_open = 5        # gap opening (for affine)
    gap_extend = 5         # gap extension (for affine)

    lengths = [100, 200, 400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800]

    linear_times = []
    affine_times = []

    for n in lengths:
        seqA = generate_random_sequence(n)
        seqB = generate_random_sequence(n)

        # --- Measure linear alignment
        start = time.time()
        score_lin = global_linear_alignment(seqA, seqB, score_matrix, gap_linear, backtracking=False)
        end = time.time()
        linear_times.append(end - start)

        # --- Measure affine alignment
        start = time.time()
        score_aff = global_affine_alignment(seqA, seqB, score_matrix, gap_open, gap_extend, backtracking=False)
        end = time.time()
        affine_times.append(end - start)

        print(f"n={n}, Linear Time={linear_times[-1]:.4f}s, Affine Time={affine_times[-1]:.4f}s")

    # Plot
    plt.plot(lengths, linear_times, label='Linear Gap')
    plt.plot(lengths, affine_times, label='Affine Gap')
    plt.xlabel("Sequence length (n)")
    plt.ylabel("Running time (seconds)")
    plt.legend()
    plt.title("Comparison of Linear vs. Affine Gap Alignment")
    plt.show()

if __name__ == "__main__":
    main()
