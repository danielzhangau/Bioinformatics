def baseCounts(sequence, ):
    print('Base counts:')
    print('Adenine: ', sequence.count('A'))
    print('Cytosine: ', sequence.count('C'))
    print('Guanine: ', sequence.count('G'))
    print('Thymine: ', sequence.count('T'))

    print('DNA: ', sequence.count('ACGU'))

exercise1_seq = "AAAACCTCTCTGTTCAGCACTTCCTCTCTCTTGGTCTGGTCTCAACGGTCACCATGGCGAGACCCTTGGAGGAGGCCCTGGATGTAATAGTGTCCACCTTCCACAAATACTCAGGCAACGAGGGTGACAAGTTCAAGCTGAACAAGACAGAGCTCAAGGAGCTACTGACCAGGGAGCTGCCTAGCTTCCTGGGGAGAAGGACAGACGAAGCTGCATTCCA"
baseCounts(exercise1_seq)
