"""
Module *** sequence ***

This module depends on the following modules

sym -- defines an alphabet
prob -- defines structures to hold probabilities (prob also depends on sym) 

sym and prob should be in the same directory as this module so it can function

This module incorporates classes for

Sequence -- names and defines a sequence of symbols; computes various transformations and pairwise alignments
Alignment -- defines a multiple sequence alignment; computes stats for use in substitution matrices
SubstMatrix -- substitution matrix class to support alignment methods

Incorporates methods for loading and saving files relevant to the above (e.g. FASTA, ALN, substitution matrices)
and methods for retrieving relevant data from web services 
"""

import numpy

from prob import *
from webservice import *


# Sequence ------------------

class Sequence(object):
    """ A biological sequence. Stores the sequence itself (as a compact array), 
    the alphabet (i.e., type of sequence it is), and optionally a name and further 
    information. """

    sequence = None  # The array of symbols that make up the sequence
    alphabet = None  # The alphabet from which symbols come
    name = None  # The name (identifier) of a sequence
    info = None  # Other information (free text; e.g. annotations)
    length = None  # The number of symbols that the sequence is composed of
    gappy = None  # True if the sequence has "gaps", i.e. positions that represent deletions relative another sequence

    def __init__(self, sequence, alphabet=None, name='', info='', gappy=False):
        """ Create a sequence with the sequence data. Specifying the alphabet,
        name and other information about the sequence are all optional.
        The sequence data is immutable (stored as a string).
        Example:
        >>> myseq = Sequence('MVSAKKVPAIAMSFGVSF')
        will create a sequence with no name, and assign one of the predefined
        alphabets on the basis of what symbols were used.
        >>> myseq.alphabet.symbols
        will output the standard protein alphabet:
        ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
        'R', 'S', 'T', 'V', 'W', 'Y'] """

        self.sequence = sequence

        # Assign an alphabet
        # If no alphabet is provided, attempts to identify the alphabet from sequence
        self.alphabet = None
        if not alphabet is None:
            for sym in self.sequence:
                if not sym in alphabet and (sym != '-' or not gappy):  # error check: bail out
                    raise RuntimeError('Invalid symbol: %c in sequence %s' % (sym, name))
            self.alphabet = alphabet
        else:
            for alphaName in preferredOrder:
                alpha = predefAlphabets[alphaName]
                valid = True
                for sym in self.sequence:
                    if not sym in alpha and (sym != '-' or not gappy):
                        valid = False
                        break
                if valid:
                    self.alphabet = alpha
                    break
            if self.alphabet is None:
                raise RuntimeError('Could not identify alphabet from sequence: %s' % name)

        # Store other information
        self.name = name
        self.info = info
        self.length = len(self.sequence)
        self.gappy = gappy

    def __len__(self):
        """ Defines what the "len" operator returns for an instance of Sequence, e.g.
        >>> seq = Sequence('ACGGTAGGA', DNA_Alphabet)
        >>> print (len(seq))
        9
        """
        return len(self.sequence)

    def __str__(self):
        """ Defines what should be printed when the print statement is used on a Sequence instance """
        str = self.name + ': '
        for sym in self:
            str += sym
        return str

    def __iter__(self):
        """ Defines how a Sequence should be "iterated", i.e. what its elements are, e.g.
        >>> seq = Sequence('AGGAT', DNA_Alphabet)
        >>> for sym in seq:
                print (sym)
        will print A, G, G, A, T (each on a separate row)
        """
        tsyms = tuple(self.sequence)
        return tsyms.__iter__()

    def __contains__(self, item):
        """ Defines what is returned when the "in" operator is used on a Sequence, e.g.
        >>> seq = Sequence('ACGGTAGGA', DNA_Alphabet)
        >>> print ('T' in seq)
        True
            which is equivalent to 
        >>> print (seq.__contains__('T'))
        True
        >>> print ('X' in seq)
        False
        """
        for sym in self.sequence:
            if sym == item:
                return True
        return False

    def __getitem__(self, ndx):
        """ Retrieve a specified index (or a "slice" of indices) of the sequence data.
            Calling self.__getitem__(3) is equivalent to self[3] 
        """
        if type(ndx) is slice:
            return ''.join(self.sequence[ndx])
        else:
            return self.sequence[ndx]

    def writeFasta(self):
        """ Write one sequence in FASTA format to a string and return it. """
        fasta = '>' + self.name + ' ' + self.info + '\n'
        data = ''.join(self.sequence)
        nlines = int(math.ceil((len(self.sequence) - 1) / 60 + 1))
        for i in range(nlines):
            lineofseq = ''.join(data[i * 60: (i + 1) * 60]) + '\n'
            fasta += lineofseq
        return fasta

    def count(self, findme=None):
        """ Get the number of occurrences of specified symbol findme OR
            if findme = None, return a dictionary of counts of all symbols in alphabet """
        if findme != None:
            cnt = 0
            for sym in self.sequence:
                if findme == sym:
                    cnt = cnt + 1
            return cnt
        else:
            symbolCounts = {}
            for symbol in self.alphabet:
                symbolCounts[symbol] = self.count(symbol)
            return symbolCounts

    def find(self, findme):
        """ Find the position of the specified symbol or sub-sequence """
        return ''.join(self.sequence).find(findme)


"""
Below are some useful methods for loading data from strings and files.
Recognize the FASTA format (nothing fancy). 
"""


def readFasta(string, alphabet=None):
    """ Read the given string as FASTA formatted data and return the list of
    sequences contained within it. """
    seqlist = []  # list of sequences contained in the string
    seqname = None  # name of *current* sequence 
    seqinfo = None
    seqdata = []  # sequence data for *current* sequence
    for line in string.splitlines():  # read every line
        if len(line) == 0:  # ignore empty lines
            continue
        if line[0] == '>':  # start of new sequence            
            if seqname:  # check if we've got one current
                try:
                    current = Sequence(seqdata, alphabet, seqname, seqinfo)
                    seqlist.append(current)
                except RuntimeError as errmsg:
                    print("Error for sequence %s: %s" % (seqname, errmsg))
            # now collect data about the new sequence
            seqinfo = line[1:].split()  # skip first char
            if len(seqinfo) > 0:
                parsed = parseDefline(seqinfo[0])
                seqname = parsed[0]
                seqinfo = line[1:]
            else:
                seqname = ''
                seqinfo = ''
            seqdata = []
        else:  # we assume this is (more) data for current
            cleanline = line.split()
            for thisline in cleanline:
                seqdata.extend(tuple(thisline.strip('*')))
    # we're done reading the file, but the last sequence remains
    if seqname:
        try:
            lastseq = Sequence(seqdata, alphabet, seqname, seqinfo)
            seqlist.append(lastseq)
        except RuntimeError as errmsg:
            print("Error for sequence %s: %s" % (seqname, errmsg))
    return seqlist


def parseDefline(string):
    """ Parse the FASTA defline (see http://en.wikipedia.org/wiki/FASTA_format)
        GenBank, EMBL, etc                gi|gi-number|gb|accession|locus
        SWISS-PROT, TrEMBL                sp|accession|name
        ...
        Return a tuple with 
        [0] primary search key, e.g. UniProt accession, Genbank GI
        [1] secondary search key, e.g. UniProt name, Genbank accession 
        [2] source, e.g. 'sp' (SwissProt/UniProt), 'tr' (TrEMBL), 'gb' (Genbank)
    """
    if len(string) == 0: return ('', '', '', '')
    s = string.split()[0]
    if re.match("sp\|[A-Z][A-Z0-9]{5}\|\S+", s):
        arg = s.split('|');
        return (arg[1], arg[2], arg[0], '')
    elif re.match("tr\|[A-Z][A-Z0-9]{5}\|\S+", s):
        arg = s.split('|');
        return (arg[1], arg[2], arg[0], '')
    elif re.match("gi\|[0-9]*\|gb|emb|dbj\|\S+\|\S+", s):
        arg = s.split('|');
        return (arg[1], arg[3], arg[2], arg[4])
    elif re.match("refseq\|\S+\|\S+", s):
        arg = s.split('|');
        return (arg[1], arg[2], arg[0], '')
    else:
        return (s, '', '', '')


def readFastaFile(filename, alphabet=None):
    """ Read the given FASTA formatted file and return the list of sequences 
    contained within it. Note that if no alphabet is specified, it will take a 
    separate guess for each sequence. """
    fh = open(filename)
    seqlist = []
    batch = ''  # a batch of rows including one or more complete FASTA entries
    rowcnt = 0
    for row in fh:
        row = row.strip()
        if len(row) > 0:
            if row.startswith('>') and rowcnt > 0:
                more = readFasta(batch, alphabet)
                if len(more) > 0:
                    seqlist.extend(more)
                batch = ''
                rowcnt = 0
            batch += row + '\n'
            rowcnt += 1
    if len(batch) > 0:
        more = readFasta(batch, alphabet)
        if len(more) > 0:
            seqlist.extend(more)
    fh.close()
    return seqlist


def writeFastaFile(filename, seqs):
    """ Write the specified sequences to a FASTA file. """
    fh = open(filename, 'w')
    for seq in seqs:
        fh.write(seq.writeFasta())
    fh.close()


def getMarkov(seqs, order=0):
    """ Retrieve the Markov stats for a set of sequences. """
    myseqs = seqs
    if seqs is Sequence:
        myseqs = list([seqs])
    myalpha = None
    for seq in myseqs:
        if myalpha == None:
            myalpha = seq.alphabet
        else:
            if seq.alphabet != myalpha:
                raise RuntimeError('Sequence ' + seq.name + ' uses an invalid alphabet ')
    jp = Joint([myalpha for _ in range(order + 1)])
    for seq in myseqs:
        for i in range(len(seq) - order):
            sub = seq[i:i + order + 1]
            jp.observe(sub)
    return jp


def getCount(seqs, findme=None):
    if findme != None:
        cnt = 0
        for seq in seqs:
            cnt += seq.count(findme)
        return cnt
    else:
        if len(seqs) > 0:
            alpha = seqs[0].alphabet
            patcnt = {}
            for a in alpha:
                patcnt[a] = getCount(seqs, a)
        return patcnt


# Alignment ------------------

class Alignment:
    """ A sequence alignment class. Stores two or more sequences of equal length where
    one symbol is gap '-' 
    Example usage:
    >>> seqs = [Sequence('THIS-LI-NE-', Protein_Alphabet, gappy = True), Sequence('--ISALIGNED', Protein_Alphabet, gappy = True)]
    >>> print (Alignment(seqs))
     THIS-LI-NE-
     --ISALIGNED """

    alignlen = None
    seqs = None
    alphabet = None

    def __init__(self, seqs):
        self.alignlen = -1
        self.seqs = seqs
        self.alphabet = None
        for s in seqs:
            if self.alignlen == -1:
                self.alignlen = len(s)
            elif self.alignlen != len(s):
                raise RuntimeError("Alignment invalid: different lengths")
            if self.alphabet != None and self.alphabet != s.alphabet:
                raise RuntimeError("Alignment invalid: different alphabets")
            self.alphabet = s.alphabet

    def getnamelen(self):
        namelen = 0
        for seq in self.seqs:
            namelen = max(len(seq.name), namelen)
        return namelen

    def __str__(self):
        string = ''
        namelen = self.getnamelen()
        for seq in self.seqs:
            string += seq.name.ljust(namelen + 1)
            for sym in seq:
                string += sym
            string += '\n'
        return string

    def writeClustal(self, filename=None):
        """ Write the alignment to a string or file using the Clustal file format. """
        symbolsPerLine = 60
        maxNameLength = self.getnamelen() + 1
        string = ''
        wholeRows = self.alignlen / symbolsPerLine
        for i in range(wholeRows):
            for j in range(len(self.seqs)):
                string += self.seqs[j].name.ljust(maxNameLength) + ' '
                string += self.seqs[j][i * symbolsPerLine:(i + 1) * symbolsPerLine] + '\n'
            string += '\n'
        # Possible last row
        lastRowLength = self.alignlen - wholeRows * symbolsPerLine
        if lastRowLength > 0:
            for j in range(len(self.seqs)):
                if maxNameLength > 0:
                    string += self.seqs[j].name.ljust(maxNameLength) + ' '
                string += self.seqs[j][-lastRowLength:] + '\n'
        if filename != None:
            fh = open(filename, 'w')
            fh.write('CLUSTAL W (1.83) multiple sequence alignment\n\n\n')  # fake header so that clustal believes it
            fh.write(string)
            fh.close()
            return
        return string

    def getProfile(self, pseudo=0.0):
        """ Determine the probability matrix from the alignment, assuming
        that each position is independent of all others. """
        p = IndepJoint([self.alphabet for _ in range(self.alignlen)], pseudo)
        for seq in self.seqs:
            p.observe(seq)
        return p

    def calcBackground(self):
        """ Count the proportion of each amino acid's occurrence in the
            alignment, and return as a probability distribution. """
        p = Distrib(self.alphabet)
        for seq in self.seqs:
            for sym in seq:
                if sym in self.alphabet:  # ignore "gaps"
                    p.observe(sym)
        return p

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
        seqPairs = math.factorial(numSeqs) / (2 * math.factorial(numSeqs - 2))
        # Number of pairs of sequences in ungapped alignment
        # i.e. how many ways can you pair off 20 sequences?
        # aaPairs is the number of possible pair combinations in a column multiplied
        # by the number of columns
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

    def writeHTML(self, filename=None):
        """ Generate HTML that displays the alignment in color.
            Requires that the alphabet is annotated with the label 'html-color' (see Sequence.annotateSym)
            and that each symbol maps to a text string naming the color, e.g. 'blue'
        """
        html = '''<html><head><meta content="text/html; charset=ISO-8859-1" http-equiv="Content-Type">\n<title>Sequence Alignment</title>\n</head><body><pre>\n'''
        maxNameLength = self.getnamelen()
        html += ''.ljust(maxNameLength) + ' '
        for i in range(self.alignlen - 1):
            if (i + 1) % 10 == 0:
                html += str(i / 10 + 1)[0]
            else:
                html += ' '
        html += '%s\n' % (self.alignlen)

        if self.alignlen > 10:
            html += ''.ljust(maxNameLength) + ' '
            for i in range(self.alignlen - 1):
                if (i + 1) % 10 == 0:
                    index = len(str(i / 10 + 1).split('.')[0])
                    html += str(i / 10 + 1).split('.')[0][(index * -1) + 1] if (
                            len(str(i / 10 + 1).split('.')[0]) > 1) else '0'
                else:
                    html += ' '
            html += '\n'

        if self.alignlen > 100:
            html += ''.ljust(maxNameLength) + ' '
            for i in range(self.alignlen - 1):
                if (i + 1) % 10 == 0 and i >= 99:
                    index = len(str(i / 10 + 1).split('.')[0])
                    html += str(i / 10 + 1).split('.')[0][-1] if (len(str(i / 10 + 1).split('.')[0]) > 2) else '0'

                else:
                    html += ' '
            html += '\n'

        if self.alignlen > 1000:
            html += ''.ljust(maxNameLength) + ' '
            for i in range(self.alignlen - 1):
                if (i + 1) % 10 == 0:
                    html += '0' if (len(str(i / 10 + 1).split('.')[0]) > 2) else ' '

                else:
                    html += ' '
            html += '\n'
        for seq in self.seqs:
            html += seq.name.ljust(maxNameLength) + ' '
            for sym in seq:
                color = self.alphabet.getAnnotation('html-color', sym)
                if not color:
                    color = 'white'
                html += '<font style="BACKGROUND-COLOR: %s">%s</font>' % (color, sym)
            html += '\n'
        html += '</pre></body></html>'
        if filename:
            fh = open(filename, 'w')
            fh.write(html)
            fh.close()
        return html


def alignGlobal(seqA, seqB, substMatrix, gap=-1):
    """ Align seqA with seqB using the Needleman-Wunsch
    (global) algorithm. subsMatrix is the substitution matrix to use and
    gap is the linear gap penalty to use. """
    lenA, lenB = len(seqA), len(seqB)
    # Create the scoring matrix (S)
    S = numpy.zeros((lenA + 1, lenB + 1))
    # Fill the first row and column of S with multiples of the gap penalty
    for i in range(lenA + 1):
        S[i, 0] = i * gap
    for j in range(lenB + 1):
        S[0, j] = j * gap
    # Calculate the optimum score at each location in the matrix S
    # (where the score represents the best possible score for an alignment
    #  that ends at sequence indices i and j, for A and B, resp.)
    for i in range(1, lenA + 1):
        for j in range(1, lenB + 1):
            match = S[i - 1, j - 1] + substMatrix.get(seqA[i - 1], seqB[j - 1])
            fromTop = S[i - 1, j] + gap
            fromLeft = S[i, j - 1] + gap
            S[i, j] = max([match, fromTop, fromLeft])
    #######################################################################
    # Exercise 3: Complete the following block of code
    # Traceback the optimal alignment
    alignA = ''  # a string for sequence A when aligned (e.g. 'THIS-LI-NE-', initially empty).
    alignB = ''  # a string for sequence B when aligned (e.g. '--ISALIGNED', initially empty).
    # Start at the end (bottom-right corner of S)
    i = lenA
    j = lenB
    # Stop when we hit the beginning of at least one sequence
    while i > 0 and j > 0:
        if S[i, j] == S[i - 1, j] + gap:
            # got here by a gap in sequence B (go up)
            alignA = seqA[i - 1] + alignA
            alignB = '-' + alignB  # fill with gap
            i -= 1
        elif S[i, j] == S[i, j - 1] + gap:
            # got here by a gap in sequence A (go left)
            alignA = '-' + alignA
            alignB = seqB[j - 1] + alignB
            j -= 1
        else:
            # got here by aligning the bases (go diagonally)
            alignA = seqA[i - 1] + alignA
            alignB = seqB[j - 1] + alignB
            i -= 1
            j -= 1
    # Fill in the rest of the alignment if it begins with gaps
    # (i.e., traceback all the way to S[0, 0])
    while i > 0:
        # Go up
        alignA = seqA[i - 1] + alignA
        alignB = '-' + alignB
        i -= 1
    while j > 0:
        # Go left
        alignA = '-' + alignA
        alignB = seqB[j - 1] + alignB
        j -= 1
    return Alignment([Sequence(alignA, seqA.alphabet, seqA.name, gappy=True),
                      Sequence(alignB, seqB.alphabet, seqB.name, gappy=True)])


def alignLocal(seqA, seqB, substMatrix, gap=-1):
    """ Align seqA with seqB using the Smith-Waterman
    (local) algorithm. subsMatrix is the substitution matrix to use and
    gap is the linear gap penalty to use. """
    lenA, lenB = len(seqA), len(seqB)
    # Create the scoring matrix (S)
    S = numpy.zeros((lenA + 1, lenB + 1))
    # Fill the first row and column of S with multiples of the gap penalty
    for i in range(lenA + 1):
        S[i, 0] = 0  # Local: init 0
    for j in range(lenB + 1):
        S[0, j] = 0  # Local: init 0
    # Calculate the optimum score at each location in the matrix S
    # (where the score represents the best possible score for an alignment
    #  that ends at sequence indices i and j, for A and B, resp.)
    for i in range(1, lenA + 1):
        for j in range(1, lenB + 1):
            match = S[i - 1, j - 1] + substMatrix.get(seqA[i - 1], seqB[j - 1])
            fromTop = S[i - 1, j] + gap
            fromLeft = S[i, j - 1] + gap
            S[i, j] = max([match, fromTop, fromLeft, 0])  # Local: add option that we re-start alignment from "0"
    #######################################################################
    # Traceback the optimal alignment
    alignA = ''
    alignB = ''
    # Local: start at the cell which has the highest score; find it
    i = 0
    j = 0
    for ii in range(1, lenA + 1):
        for jj in range(1, lenB + 1):
            if S[ii, jj] > S[i, j]:
                i = ii
                j = jj
    # Stop when we hit the beginning of at least one sequence
    # Local: also stop when we hit a score 0
    while i > 0 and j > 0 and S[i, j] > 0:
        if S[i, j] == S[i - 1, j] + gap:
            # Got here by a gap in sequence B (go up)
            alignA = seqA[i - 1] + alignA
            alignB = '-' + alignB
            i -= 1
        elif S[i, j] == S[i, j - 1] + gap:
            # Got here by a gap in sequence A (go left)
            alignA = '-' + alignA
            alignB = seqB[j - 1] + alignB
            j -= 1
        else:
            # Got here by aligning the bases (go diagonally)
            alignA = seqA[i - 1] + alignA
            alignB = seqB[j - 1] + alignB
            i -= 1
            j -= 1

    return Alignment([Sequence(alignA, seqA.alphabet, seqA.name, gappy=True),
                      Sequence(alignB, seqB.alphabet, seqB.name, gappy=True)])


def readClustal(string, alphabet):
    """ Read a ClustalW2 alignment in the given string and return as an
    Alignment object. """
    seqs = {}  # sequence data
    for line in string.splitlines():
        if line.startswith('CLUSTAL') or line.startswith('STOCKHOLM') \
                or line.startswith('#'):
            continue
        if len(line.strip()) == 0:
            continue
        if line[0] == ' ' or '*' in line or ':' in line:
            continue
        sections = line.split()
        name, seqstr = sections[0:2]
        if name in seqs:
            seqs[name] += seqstr
        else:
            seqs[name] = seqstr
    sequences = []
    for name, seqstr in list(seqs.items()):
        sequences.append(Sequence(seqstr, alphabet, name, gappy=True))
    return Alignment(sequences)


def readClustalFile(filename, alphabet):
    """ Read a ClustalW2 alignment file and return an Alignment object
    containing the alignment. """
    fh = open(filename)
    data = fh.read()
    fh.close()
    aln = readClustal(data, alphabet)
    return aln


# Substitution Matrix ------------------

class SubstMatrix:

    def __init__(self, alphabet):
        self.scoremat = {}  # set as an empty dictionary to be filled by set function
        self.alphabet = alphabet

    def _getkey(self, sym1, sym2):
        """ Construct canonical (unordered) key for two symbols """
        if sym1 <= sym2:
            return tuple([sym1, sym2])
        else:
            return tuple([sym2, sym1])

    def set(self, sym1, sym2, score):
        """ Add a score to the substitution matrix """
        self.scoremat[self._getkey(sym1, sym2)] = score

    def get(self, sym1, sym2):
        return self.scoremat[self._getkey(sym1, sym2)]

    def __str__(self):
        symbols = self.alphabet.symbols  # what symbols are in the alphabet
        i = len(symbols)
        string = ''
        for a in symbols:
            string += a + ' '
            for b in symbols[:len(symbols) - i + 1]:
                score = self.scoremat[self._getkey(a, b)]
                if score != None:
                    string += str(score).rjust(3) + ' '
                else:
                    string += "?".rjust(3) + ' '
            string += '\n'
            i -= 1
        string += '    ' + '   '.join(self.alphabet.symbols)
        return string

    def writeFile(self, filename):
        """ Write this substitution matrix to the given file. """
        fh = open(filename, 'w')
        file = ''
        for key in self.scoremat:
            file += ''.join(key) + ': ' + str(self.scoremat[key]) + '\n'
        fh.write(file)
        fh.close()


def readSubstMatrix(filename, alphabet):
    """ Read in the substitution matrix stored in the given file. """
    mat = SubstMatrix(alphabet)
    fh = open(filename, 'r')
    data = fh.read()
    fh.close()
    lines = data.splitlines()
    for line in lines:
        if len(line.strip()) == 0:
            continue
        symbols, score = line.split(':')
        score = int(score)
        mat.set(symbols[0], symbols[1], score)
    return mat


# import os
# os.chdir('/Users/mikael/workspace/binf/data/')  # set to the directory where you keep your files
# BLOSUM62 = readSubstMatrix('blosum62.matrix', Protein_Alphabet)

# Web Service Functions -------------------

def getSequence(id, database='uniprotkb', start=None, end=None):
    """ Get the sequence identified by the given ID from the given database
    (e.g. 'uniprotkb', 'refseqn' or 'refseqp'), and return it as a Sequence
    object. An error is caused if the sequence ID is not found. If start and
    end are given, then only that section of the sequence is returned. 
    Note: more flexible search options are supported by using webservice.fetch
    directly."""

    MAX_TRY = 5

    for i in range(MAX_TRY):
        try:
            fastaData = fetch(id, database)
            seq = readFasta(fastaData)[0]
            break
        except:
            from time import sleep
            print('Failed on {i}th try for id {id}'.format(i=i, id=id))
            sleep(0.1)
    try:
        return Sequence(seq[start:end], seq.alphabet, seq.name, seq.info)
    except:
        raise RuntimeError(
            'An error occurred while retrieving the specified sequence: %s (maybe the ID doesn\'t exist)' % id)


def searchSequences(query, database='uniprot'):
    """ Search for sequences matching the given query in the given database
    (must be 'uniprot'), and return a list of sequence IDs. """
    ids = search(query, limit=None)
    return ids


def runClustal(sequences, method='slow'):
    """ Run a ClustalOmega alignment of the given list of Sequence objects.
    Return an Alignment object. Method should be one of 'fast' or 'slow'. """
    alpha = None
    for seq in sequences:
        if alpha == None:
            alpha = seq.alphabet
        elif alpha != seq.alphabet:
            raise RuntimeError("Invalid alphabet: " + str(seq.alphabet) + ". Not compatible with " + str(alpha))
    serviceName = 'clustalo'
    resultType = 'aln-clustal'
    fastaSeqs = ''.join([seq.writeFasta() for seq in sequences])
    params = {'alignment': method.lower(), 'sequence': fastaSeqs}
    service = EBI(serviceName)
    result = service.submit(params, resultType)
    alignment = readClustal(result, alpha)
    return alignment


def createTree(alignment, type):
    """ Run a ClustalW 2 phylogeny tree creation of either a 'Neighbour-joining'
    or 'UPGMA' type tree from the given multiple sequence Alignment object. """
    if not type in ['Neighbour-joining', 'UPGMA']:
        raise RuntimeError('type must be either \'Neighbour-joining\' or \'UPGMA\'.')
    serviceName = 'clustalw2_phylogeny'
    resultType = 'tree'
    output = 'dist'
    clustalAln = alignment.writeClustal()
    params = {'tree': output, 'sequence': clustalAln, 'clustering': type, 'tossgaps': 'true'}
    service = EBI(serviceName)
    tree = service.submit(params, resultType)
    return tree


def runBLAST(sequence, program='blastp', database='uniprotkb', exp='1e-1'):
    """ Run a BLAST search of nucleotide mouse databases using the given
    sequence as a query. Return a list of matched sequence IDs, in descending
    order of similarity to query sequence. 
    program: either blastn (nucleotide) or blastp (protein)
    database: many available, e.g. uniprotkb, pdb (protein); em_rel, nrnl1 (EMBL nucleotide, non-redundant resp)
        (for protein see http://www.ebi.ac.uk/Tools/sss/ncbiblast/help/index-protein.html#database)
        (for nucleotide see http://www.ebi.ac.uk/Tools/sss/ncbiblast/help/index-nucleotide.html#database)
    exp: E-value threshold (select only hits that have a better E-value than this)
    """
    if sequence.alphabet == predefAlphabets['DNA']:
        stype = 'dna'
    elif sequence.alphabet == predefAlphabets['RNA']:
        stype = 'rna'
    else:
        stype = 'protein'
    serviceName = 'ncbiblast'
    resultTypes = ['ids', 'out']  # request
    fastaSeq = sequence.writeFasta()
    databases = [database]
    params = {'program': program, 'database': databases, 'sequence': fastaSeq,
              'stype': stype, 'exp': exp}
    service = EBI(serviceName)
    idsData, output = service.submit(params, resultTypes)
    ids = []
    for id in idsData.splitlines():
        if len(id) > 0:
            ids.append(id.split(':')[1])
    return ids
