from random import choices

def testing_data(n):

    bases = ["A", "G", "C", "T"] 

    seq1 = "".join(choices(bases, k=n))
    seq2 = "".join(choices(bases, k=n))

    return seq1, seq2