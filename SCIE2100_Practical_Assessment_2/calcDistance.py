def calcDistances(self, measure, a=1.0):
    """ Calculate the evolutionary distance between all pairs of sequences
    in this alignment, using the given measure. Measure can be one of
    'fractional', 'poisson' or 'gamma'.
    Definitions of each distance metric are found in Zvelebil and Baum p268-276.
    These are mostly intended for DNA, but adapted for protein (as below).
    Note however that there are alternative distance matrices for proteins (p276).
    """
    measure = measure.lower()
    if not measure in ['fractional', 'poisson', 'gamma']:
        raise RuntimeError('Unsupported evolutionary distance measure: %s' % measure)
    if len(self.alphabet) not in [4, 20]:
        raise RuntimeError('Invalid sequence alphabet: %s' % str(self.alphabet))
    distmat = numpy.zeros((len(self.seqs), len(self.seqs)))
    # Loop through each pair of sequences
    for i in range(len(self.seqs)):
        for j in range(i + 1, len(self.seqs)):
            seqA = self.seqs[i]
            seqB = self.seqs[j]
            # Calculate the fractional distance (p) first
            # The two sequences of interest are in seqA and seqB
            L = 0
            D = 0
            for k in range(self.alignlen):
                # For every non-gapped column, put to L
                # For every non-gapped column where the sequences are
                # different, put to D
                if seqA[k] != '-' and seqB[k] != '-':
                    L += 1
                    if seqA[k] != seqB[k]:
                        D += 1
            p = float(D) / L
            # Now calculate the specified measure based on p
            if measure == 'fractional':
                dist = p
            elif measure == 'poisson':
                dist = -math.log(1 - p)
            elif measure == 'gamma':
                dist = a * ((1 - p) ** (-1 / a) - 1)
            distmat[i, j] = distmat[j, i] = dist
    return distmat