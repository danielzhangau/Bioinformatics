from sequence import *
from sym import *
from sstruct import *

# read both protein and structure sequences into lists
prot2 = readFastaFile("prot2.fa", Protein_Alphabet)
sstr3 = readFastaFile("sstr3.fa", DSSP3_Alphabet)

# store protein and structure sequences in dictionary for easy retrieval
protein_map = {}
structure_map = {}

for protein in prot2:
    protein_map.update({protein.name: protein})

for structure in sstr3:
    structure_map.update({structure.name: structure})

protein = protein_map.get("1EVH")
structure = structure_map.get("1EVH")
print("protein: ", protein)
print("structure: ", structure)

print("Steps 1 to 4")
# Step 2: Alpha-helix
alpha = getScores(protein, 0) # values from column 0
beta = getScores(protein, 1) # values from column 1
calls_a1 = markCountAbove(alpha, width=6, call_cnt=4)
print(makesstr(calls_a1, 'H'))

calls_a2 = extendDownstream(alpha, calls_a1, width=4)
calls_a3 = extendUpstream(alpha, calls_a2, width=4)

calls_b1 = markCountAbove(beta, width=5, call_cnt=3)
print(makesstr(calls_b1, 'E'))
calls_b2 = extendDownstream(beta, calls_b1, width=4)
calls_b3 = extendUpstream(beta, calls_b2, width=4)

avg_a = calcRegionAverage(alpha, calls_a3)
avg_b = calcRegionAverage(beta, calls_b3)

diff_a = [avg_a[i] - avg_b[i] for i in range(len(avg_a))]
diff_b = [avg_b[i] - avg_a[i] for i in range(len(avg_a))]

print("Step 5")
calls_a4 = checkSupport(calls_a3, diff_a)
calls_b4 = checkSupport(calls_b3, diff_b)
print(makesstr(calls_a4, 'H'))
print(makesstr(calls_b4, 'E'))

print("Final prediction")
# Combine calls and replace remaining '-' symbols to C (coil)
prediction = predicted_sequence = ''.join(['H' if alpha else ('E' if beta else 'C') for (alpha, beta) in zip(calls_a4, calls_b4)])
print(prediction)

position = 0
tp = 0  # number of true positives (correctly identified calls)
tn = 0  # number of true negatives (correctly missed no-calls)
fp = 0  # number of false positives (incorrectly identified no-calls)
fn = 0  # number of false negatives (incorrectly missed calls)
for element in structure.sequence:
    if element == "H" and calls_a4[position]:
        tp += 1 # matched alpha helix
    if element != "H" and not(calls_a4[position]):
        tn += 1 # matched beta sheet
    if element != "H" and calls_a4[position]:
        fp += 1 # matched beta sheet
    if element == "H" and not(calls_a4[position]):
        fn += 1 # matched beta sheet
    position += 1

print("alpha-helix Accuracy %.2f%%" % ((float((tp+tn) /  (tp+tn+fp+fn)) * 100)))

position = 0
tp = 0  # number of true positives (correctly identified calls)
tn = 0  # number of true negatives (correctly missed no-calls)
fp = 0  # number of false positives (incorrectly identified no-calls)
fn = 0  # number of false negatives (incorrectly missed calls)
for element in structure.sequence:
    if element == "E" and calls_b4[position]:
        tp += 1 # matched alpha helix
    if element != "E" and not(calls_b4[position]):
        tn += 1 # matched beta sheet
    if element != "E" and calls_b4[position]:
        fp += 1 # matched beta sheet
    if element == "E" and not(calls_b4[position]):
        fn += 1 # matched beta sheet
    position += 1

print("beta-strand Accuracy %.2f%%" % ((float((tp+tn) /  (tp+tn+fp+fn)) * 100)))