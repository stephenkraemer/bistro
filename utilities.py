def reverse_complement_seq(seq):
    return seq.translate(str.maketrans('CGHD', 'GCDH')).reverse()


def complement_seq(seq):
    return seq.translate(str.maketrans('CGHD', 'GCDH'))
