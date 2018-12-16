import numpy as np
from Bio.SubsMat import MatrixInfo
def needleman_wunsch(seq1, seq2, gap):
    n = len(seq1) + 1  # Dimension of matrix column
    m = len(seq2) + 1  # Dimension of matrix row
    align_matrix = np.zeros((m, n), dtype=int)  # Create a matrix of dimensions (n,m)
    matrix = MatrixInfo.blosum62
    # Fill the first row elements of the matrix with multiple of gap penalty.
    for i in range(m):
        align_matrix[i][0] = gap * i

    # Fill the first column elements of the matrix with multiple of gap penalty.
    for j in range(n):
        align_matrix[0][j] = gap * j

    # Fill the rest of the matrix.
    for i in range(1, m):
        for j in range(1, n):
            try:
                m = matrix[seq1[j-1], seq2[i-1]]
            except KeyError:
                m = matrix[seq2[i-1], seq1[j-1]]
            di = align_matrix[i-1][j-1] + m
            ho = align_matrix[i][j-1] + gap
            ve = align_matrix[i-1][j] + gap
            align_matrix[i][j] = max(di, ho, ve)

    i = len(seq2)
    j = len(seq1)
    seq2_t = '>'
    seq1_s = '>'

    while i > 0 and j > 0:
        try:
            m = matrix[seq1[j-1], seq2[i-1]]
        except KeyError:
            m = matrix[seq2[i-1], seq1[j-1]]

        if align_matrix[i, j] - m == align_matrix[i-1, j-1]:
            seq2_t = seq2_t + seq2[i-1]
            seq1_s = seq1_s + seq1[j-1]
            i = i-1
            j = j-1

        elif align_matrix[i, j] - gap == align_matrix[i, j-1]:
            seq1_s = seq1_s + seq1[j-1]
            seq2_t = seq2_t + '_'
            j = j-1

        elif align_matrix[i, j] - gap == align_matrix[i-1, j]:
            seq1_s = seq1_s + '_'
            seq2_t = seq2_t + seq2[i-1]
            i = i-1

    if i > 0:
        while i > 0:
            seq1_s = seq1_s + '_'
            seq2_t = seq2_t + seq2[i-1]
            i = i -1

    elif j > 0:
        while j > 0:
            seq1_s = seq1_s + seq1[j-1]
            seq2_t = seq2_t + '_'
            j = j -1

    print(seq2_t[::-1])
    print(seq1_s[::-1])
needleman_wunsch(seq1="MGGETFA", seq2 ="GGVTTF", gap=-4)
