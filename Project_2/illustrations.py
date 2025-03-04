# illustrations.py max_sequence_length nr_replicates

### Analyse running time
import sys
from time import time
import matplotlib.pyplot as plt
from experiments.testing_data import testing_data
from common.io import load_score_matrix
from linear.alignment import global_linear_alignment
from affine.alignment import global_affine_alignment

# score matrix and gap scores
score_matrix = load_score_matrix("experiments\eval_score_matrix.txt")
gap, ge, go = 5, 2, 5
start, end, number = 10, int(sys.argv[1]), 6
nrep = int(sys.argv[2])

replication = []
# affine and linear alignments
running_time_linear = []
running_time_affine = []

for rep in range(nrep):
    # generate testing data
    lengths = [int(start + i * (end - start) / (number - 1)) for i in range(number)]
    sequence_pairs = {p:testing_data(p) for p in lengths}

    for pair in sequence_pairs:
        replication.append(rep)

        seq1 = sequence_pairs[pair][0]
        seq2 = sequence_pairs[pair][1]

        # for linear gap cost
        start_time = time()
        score, alignments = global_linear_alignment(seq1, seq2, score_matrix, gap)
        end_time = time()

        running_time_linear.append(end_time - start_time)
        
        # for affine gap cost
        start_time = time()
        score, alignments = global_affine_alignment(seq1, seq2, score_matrix, ge, go)
        end_time = time()

        running_time_affine.append(end_time - start_time)

# illustrations

# plotting the points 
plt.subplot(1, 2, 1)
plt.scatter(lengths*nrep, running_time_linear, marker='o', label='Linear', alpha=0.3)
plt.xlabel('Sequence length')
plt.ylabel('Running time (sec)')
plt.legend()  

plt.subplot(1, 2, 2)
plt.scatter(lengths*nrep, running_time_affine, marker='s', label='Affine', alpha=0.3)
plt.xlabel('Sequence length')
plt.ylabel('Running time (sec)')
plt.legend()  

plt.suptitle(f"Number of replications: {nrep}; max length of sequence: {end}")
plt.tight_layout()
plt.show()
