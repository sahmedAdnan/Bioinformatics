# affine/alignment.py
import numpy as np

def global_affine_alignment(seq1, seq2, score_matrix, gap_open, gap_extend, show_all=False):
    m, n = len(seq1), len(seq2)
    
    # Initialize matrices with infinity.
    M = np.full((m+1, n+1), np.inf)
    X = np.full((m+1, n+1), np.inf)
    Y = np.full((m+1, n+1), np.inf)
    
    # Initialize traceback structure.
    trace = np.empty((m+1, n+1), dtype=object)
    for i in range(m+1):
        for j in range(n+1):
            trace[i][j] = {'M': [], 'X': [], 'Y': []}
    
    # Base case.
    M[0][0] = 0
    
    # Initialize first column.
    for i in range(1, m+1):
        X[i][0] = gap_open + (i-1) * gap_extend
        trace[i][0]['X'].append(('X', i-1, 0))
    
    # Initialize first row.
    for j in range(1, n+1):
        Y[0][j] = gap_open + (j-1) * gap_extend
        trace[0][j]['Y'].append(('Y', 0, j-1))
    
    # Fill matrices.
    for i in range(1, m+1):
        for j in range(1, n+1):
            # Update M matrix.
            match_cost = score_matrix[seq1[i-1]][seq2[j-1]]
            candidates = [M[i-1][j-1], X[i-1][j-1], Y[i-1][j-1]]
            M[i][j] = match_cost + min(candidates)
            for idx, val in enumerate(candidates):
                if M[i][j] == match_cost + val:
                    trace[i][j]['M'].append((['M', 'X', 'Y'][idx], i-1, j-1))
            
            # Update X matrix.
            x_options = [
                M[i-1][j] + gap_open + gap_extend,
                X[i-1][j] + gap_extend
            ]
            X[i][j] = min(x_options)
            if X[i][j] == x_options[0]:
                trace[i][j]['X'].append(('M', i-1, j))
            if X[i][j] == x_options[1]:
                trace[i][j]['X'].append(('X', i-1, j))
            
            # Update Y matrix.
            y_options = [
                M[i][j-1] + gap_open + gap_extend,
                Y[i][j-1] + gap_extend
            ]
            Y[i][j] = min(y_options)
            if Y[i][j] == y_options[0]:
                trace[i][j]['Y'].append(('M', i, j-1))
            if Y[i][j] == y_options[1]:
                trace[i][j]['Y'].append(('Y', i, j-1))
    
    # Determine optimal score.
    final_score = min(M[m][n], X[m][n], Y[m][n])
    
    # Backtracking to recover alignments.
    all_alignments = []
    
    def backtrack(i, j, a1, a2, state):
        if i == 0 and j == 0:
            all_alignments.append((a1[::-1], a2[::-1]))
            return
        pointers = trace[i][j][state]
        for ptr in pointers:
            new_state, pi, pj = ptr
            if state == 'M':
                new_a1 = a1 + seq1[i-1]
                new_a2 = a2 + seq2[j-1]
                backtrack(pi, pj, new_a1, new_a2, new_state)
            elif state == 'X':
                new_a1 = a1 + seq1[i-1]
                new_a2 = a2 + '-'
                backtrack(pi, pj, new_a1, new_a2, new_state)
            elif state == 'Y':
                new_a1 = a1 + '-'
                new_a2 = a2 + seq2[j-1]
                backtrack(pi, pj, new_a1, new_a2, new_state)
    
    if M[m][n] == final_score:
        backtrack(m, n, '', '', 'M')
    if X[m][n] == final_score:
        backtrack(m, n, '', '', 'X')
    if Y[m][n] == final_score:
        backtrack(m, n, '', '', 'Y')
    
    return final_score, all_alignments if show_all else all_alignments[:1]
