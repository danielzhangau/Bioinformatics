def getConsensus(self):
    seq = ""
    for index in range(self.getSize()):
        seq += self.getConsensusForColumn(index)
    return seq