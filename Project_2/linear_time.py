import random
import string
import numpy as np
import time
from common.io import load_score_matrix
import matplotlib.pyplot as plt
from linear.alignment import global_linear_alignment

def generate_input(n):
    nucleotides = ['A', 'C', 'T', 'G']
    sequence1 = ''.join(random.choices(nucleotides, k=n))
    sequence2 = ''.join(random.choices(nucleotides, k=n))
    substitution_matrix = load_score_matrix("eval_score_matrix.txt")
    gap_cost = 5 
    return sequence1, sequence2, substitution_matrix, gap_cost

np.random.seed(42)

input_sizes = [100, 200, 400, 500, 600, 800, 1000, 1200, 1400, 1600, 1800]

measured_times = []

for n in input_sizes:
    sequence1, sequence2, substitution_matrix, gap_cost = generate_input(n)
    start_time = time.time()
    global_linear_alignment(sequence1, sequence2, substitution_matrix, gap_cost, backtracking=False)
    end_time = time.time()
    running_time = (end_time - start_time)
    measured_times.append(running_time)

# Calculate the theoretical running times
theoretical_times = [n**2 for n in input_sizes]

fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot measured times on the left y-axis.
ax1.set_xlabel('Input Size', fontsize=12)
ax1.set_ylabel('Measured Time (s)', color='tab:red', fontsize=12)
line1, = ax1.plot(input_sizes, measured_times, color='tab:red', marker='o', 
                  label='Measured Time')
ax1.tick_params(axis='y', labelcolor='tab:red')
ax1.grid(True, linestyle='--', alpha=0.5)


ax2 = ax1.twinx()  
ax2.set_ylabel('Theoretical Time (s)', color='tab:blue', fontsize=12)
line2, = ax2.plot(input_sizes, theoretical_times, color='tab:blue', linestyle='--', 
                  marker='s', label='Theoretical Time')
ax2.tick_params(axis='y', labelcolor='tab:blue')

lines = [line1, line2]
labels = [line.get_label() for line in lines]
ax1.legend(lines, labels, loc='upper left', fontsize=12)

plt.title('Measured (linear) vs. Theoretical Runtime', fontsize=14)
fig.tight_layout()  
plt.show()