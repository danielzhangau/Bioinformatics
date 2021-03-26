from sequence import *

sigseqs = searchSequences('"signal+peptide"+AND+organism:3702+AND+length:[700+TO+*]')
lipseqs = searchSequences('"lipid+metabolism"+AND+organism:3702+AND+length:[700+TO+*]')

ids1 = set(sigseqs)
ids2 = set(lipseqs)

common_ids = list(ids1.intersection(ids2))  # the clever method
print(common_ids)
print(len(common_ids))
