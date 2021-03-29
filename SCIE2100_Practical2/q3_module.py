def calcSubstMatrix(self, background=None):
    """ Return a substitutionMatrix whose foreground (fg) is based on this un-gapped
    multiple sequence alignment. Scores are given in half-bits. """
    # Get a list of the amino acids
    aminoAcids = self.alphabet.symbols
    columns = self.alignlen  # Length of sequences in alignment
    numSeqs = len(self.seqs)  # Number of sequences in alignment

    # seqPairs is the total number of possible combinations of two amino acid alignments
    # in all given sequences it can be written as
    #    numSeqs!       1
    # -------------- * ---
    # (numSeqs - 2)!    2
    #
    # aaPairs is the number of possible pair combinations in a column multiplied
    # by the number of columns
    seqPairs = math.factorial(numSeqs) / (2 * math.factorial(numSeqs - 2))
    # Number of pairs of sequences in ungapped alignment
    # i.e. how many ways can you pair off 20 sequences?
    aaPairs = columns * seqPairs  # Number of pairs of amino acids in ungapped alignment
    # i.e. how many pairs of aa's would you get across all columns?

    # For each pair of amino acids in this alignment,
    # calculate the proportion of all aligned amino acids made up of that pair
    # (i.e., q[ab] = fab / aaPairs, where fab is the number of times
    #  a and b are aligned in this alignment)
    # See page 122 in Understanding Bioinformatics.
    q = {}
    for i in range(len(aminoAcids)):
        a = aminoAcids[i]
        for j in range(i, len(aminoAcids)):
            b = aminoAcids[j]
            # Count the number of times each pair of amino acids is aligned
            fab = 0
            for column in range(columns):
                # Count number of each amino acid in each column
                col = [seq[column] for seq in self.seqs]
                # Hint: if you don't understand what a and b represent at this
                # point in the code, go back and check. Remember: we want to explore
                # every amino acid pair
                if a == b:
                    # Number of ways of pairing up n occurrences of amino
                    # acid a is n*(n-1)/2
                    cnt = col.count(a)
                    fab += cnt * (cnt - 1) / 2
                else:
                    # Number of ways of pairing up n & m occurrences of
                    # amino acids a & b is n*m
                    fab += col.count(a) * col.count(b)
            # Calculate proportion of all aligned pairs of amino acids
            q[a + b] = q[b + a] = float(fab) / aaPairs
            if q[a + b] == 0:  # This is so we don't end up doing log(0)
                q[a + b] = q[b + a] = 0.001
    # Background frequency calculation if required
    p = background or self.calcBackground()
    # Calculate log-odds ratio for each pair of amino acids
    s = SubstMatrix(self.alphabet)
    for a in aminoAcids:
        for b in aminoAcids:
            # Calculate random chance probability (eab)
            if a == b:
                eab = p[a] * p[a]
                # eaa = pa ^ 2
            else:
                eab = 2 * p[a] * p[b]
                # eab = 2*pa*pb
            if eab == 0:
                eab = 0.001
            # Calculate final score to be set in the substitution matrix
            odds = q[a + b] / eab
            sab = math.log(odds, 2)  # log_2 transform
            sab = sab * 2  # units in half bits
            s.set(a, b, int(round(sab)))
    return s
