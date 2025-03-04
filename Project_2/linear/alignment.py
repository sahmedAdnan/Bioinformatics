# linear/alignment.py
import numpy as np

def global_linear_alignment(seq1, seq2, score_matrix, gap_penalty):
    m, n = len(seq1), len(seq2)
    dp = np.zeros((m+1, n+1), dtype=int)
    traceback = [[[] for j in range(n+1)] for i in range(m+1)]
    
    # Initialize first row and column.
    for i in range(1, m+1):
        dp[i][0] = i * gap_penalty
        traceback[i][0].append('up')
    for j in range(1, n+1):
        dp[0][j] = j * gap_penalty
        traceback[0][j].append('left')
    
    # Fill DP table and record optimal moves.
    for i in range(1, m+1):
        for j in range(1, n+1):
            diag = dp[i-1][j-1] + score_matrix[seq1[i-1]][seq2[j-1]]
            up = dp[i-1][j] + gap_penalty
            left = dp[i][j-1] + gap_penalty
            min_cost = min(diag, up, left)
            dp[i][j] = min_cost
            if diag == min_cost:
                traceback[i][j].append('diag')
            if up == min_cost:
                traceback[i][j].append('up')
            if left == min_cost:
                traceback[i][j].append('left')
    
    # Backtracking to retrieve all alignments.
    alignments = []
    def backtrack(i, j, aligned1, aligned2):
        if i == 0 and j == 0:
            alignments.append((aligned1[::-1], aligned2[::-1]))
            return
        for direction in traceback[i][j]:
            if direction == 'diag':
                backtrack(i-1, j-1, aligned1 + seq1[i-1], aligned2 + seq2[j-1])
            elif direction == 'up':
                backtrack(i-1, j, aligned1 + seq1[i-1], aligned2 + '-')
            elif direction == 'left':
                backtrack(i, j-1, aligned1 + '-', aligned2 + seq2[j-1])
    
    backtrack(m, n, "", "")
    return dp[m][n], alignments
