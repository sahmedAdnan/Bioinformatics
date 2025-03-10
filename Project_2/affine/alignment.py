import numpy as np

def global_affine_alignment(seq1, seq2, score_matrix, gap_open, gap_extend, backtracking=True, show_all=True):
    m, n = len(seq1), len(seq2)

    # Initialize DP tables
    M = np.full((m+1, n+1), np.inf)
    X = np.full((m+1, n+1), np.inf)
    Y = np.full((m+1, n+1), np.inf)

    if backtracking:
        trace = np.empty((m+1, n+1), dtype=object)
        for i in range(m+1):
            for j in range(n+1):
                trace[i][j] = {'M': [], 'X': [], 'Y': []}
    else:
        trace = None  

    M[0][0] = 0

    for i in range(1, m+1):
        X[i][0] = gap_open + i * gap_extend
        if backtracking:
            trace[i][0]['X'].append(('X', i-1, 0))

    for j in range(1, n+1):
        Y[0][j] = gap_open + j * gap_extend
        if backtracking:
            trace[0][j]['Y'].append(('Y', 0, j-1))

    for i in range(1, m+1):
        for j in range(1, n+1):
            match_cost = score_matrix[seq1[i-1]][seq2[j-1]]
            candidates = [
                M[i-1][j-1],  
                X[i-1][j-1],  
                Y[i-1][j-1]   
            ]
            best_prev = min(candidates)
            M[i][j] = match_cost + best_prev

            if backtracking:
                for idx, val in enumerate(candidates):
                    if best_prev == val:
                        prev_state = ['M', 'X', 'Y'][idx]
                        trace[i][j]['M'].append((prev_state, i-1, j-1))

            x_options = [
                M[i-1][j] + gap_open + gap_extend,  
                X[i-1][j] + gap_extend            
            ]
            best_x = min(x_options)
            X[i][j] = best_x

            if backtracking:
                if x_options[0] == best_x:
                    trace[i][j]['X'].append(('M', i-1, j))
                if x_options[1] == best_x:
                    trace[i][j]['X'].append(('X', i-1, j))

            y_options = [
                M[i][j-1] + gap_open + gap_extend, 
                Y[i][j-1] + gap_extend            
            ]
            best_y = min(y_options)
            Y[i][j] = best_y

            if backtracking:
                if y_options[0] == best_y:
                    trace[i][j]['Y'].append(('M', i, j-1))
                if y_options[1] == best_y:
                    trace[i][j]['Y'].append(('Y', i, j-1))

    final_score = min(M[m][n], X[m][n], Y[m][n])

    if not backtracking:
        return final_score
    
    all_alignments = []

    def backtrack(i, j, alignedA, alignedB, state):
        if i == 0 and j == 0:
            all_alignments.append((alignedA[::-1], alignedB[::-1]))
            return
        pointers = trace[i][j][state]
        for (prev_state, pi, pj) in pointers:
            if state == 'M':
                # matched or mismatched characters
                newA = alignedA + seq1[i-1]
                newB = alignedB + seq2[j-1]
                backtrack(pi, pj, newA, newB, prev_state)
            elif state == 'X':
                # gap in seq2
                newA = alignedA + seq1[i-1]
                newB = alignedB + '-'
                backtrack(pi, pj, newA, newB, prev_state)
            elif state == 'Y':
                # gap in seq1
                newA = alignedA + '-'
                newB = alignedB + seq2[j-1]
                backtrack(pi, pj, newA, newB, prev_state)

    if M[m][n] == final_score:
        backtrack(m, n, '', '', 'M')
    if X[m][n] == final_score:
        backtrack(m, n, '', '', 'X')
    if Y[m][n] == final_score:
        backtrack(m, n, '', '', 'Y')

    # If show_all=False, return only the first alignment found.
    if show_all:
        return final_score, all_alignments
    else:
        return final_score, all_alignments[:1]
