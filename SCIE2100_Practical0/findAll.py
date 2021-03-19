def findAll(self, findme):
    """ Find all occurrences of the specified symbol or sub-sequence.
        Uses a sliding window approach - very common when searching for motifs
        or contamination.
        Then, calculate the percentage of sequence made up of 'findme'.
        If the percentage is above 20%, replace all occurrences of 'findme' with Ns. 
    """
    word = len(findme)
    store = []
    count = 0
    for s in range(0, len(sequence)):
        cur = self.sequence[s:s+word]
        if cur == findme:
            store.append(s)
            count += 1
            print(cur, "Match at position ", s)
        else:
            print(cur)
    print("---------------")
    percentage = count*word / len(self.sequence)
    print("Percentage of findme in sequence: ", percentage)
    if percentage > 20:
        print("Removing contamination")
        newSeq = self.sequence.replace(findme, "N")
    return store, percentage, newSeq
