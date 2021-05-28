def getScoresWithForLoop(seq, index=0):
    """ Same functionality as getScores, but using a for loop instead of a list comprehension.
    """
    scoreList = []
    for s in seq:
        scoreList.append(cf_dict[s.upper()][index])

    return scoreList