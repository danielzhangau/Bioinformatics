from sequence import *
from phylo import *
import csv

yeasts = dict()
yeasts_list = []
with open('sugars.csv', 'rt') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if row[0] != "Yeast": #Ignore the header
            yeasts_list.append(row[0])
            #Update dictionary where key is species (row[0])
            #Value is list of booleans
            yeasts[row[0]] = [y == 'True' for y in row[1:]]
print(yeasts)
print("\nyeasts_list: ",yeasts_list)

# read fasta file
seqs = readFastaFile('MalS.fa')

# remove unmatched seq
for item in seqs:
    if item.name not in yeasts_list:
        seqs.remove(item)
print(len(seqs))

for seq in seqs:
    # print(seq.name, len(seq), len(seq.alphabet))
    seq.name = seq.info.replace(" ", "_")
    print(seq.name)

# save filtered seqs to a new fasta file
writeFastaFile('select.fa', seqs)