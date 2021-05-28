from sequence import *
from sym import *
from sstruct import *

atp = 0  # number of true positives (correctly identified calls)
atn = 0  # number of true negatives (correctly missed no-calls)
afp = 0  # number of false positives (incorrectly identified no-calls)
afn = 0  # number of false negatives (incorrectly missed calls)

btp = 0
btn = 0
bfp = 0
bfn = 0

ctp = 0
ctn = 0
cfp = 0
cfn = 0
tp = 0
tn = 0
fp = 0
fn = 0

# read both protein and structure sequences into lists
prot2 = readFastaFile("prot2.fa", Protein_Alphabet)
sstr3 = readFastaFile("sstr3.fa", DSSP3_Alphabet)

for index in range(5):
    protein = prot2[index]
    structure = sstr3[index]
    #     print("protein: ", protein)
    #     print("structure: ", structure)

    #     print("Steps 1 to 4")
    # Step 2: Alpha-helix
    alpha = getScores(protein, 0)  # values from column 0
    beta = getScores(protein, 1)  # values from column 1
    calls_a1 = markCountAbove(alpha, width=6, call_cnt=4)
    #     print(makesstr(calls_a1, 'H'))

    calls_a2 = extendDownstream(alpha, calls_a1, width=4)
    calls_a3 = extendUpstream(alpha, calls_a2, width=4)

    calls_b1 = markCountAbove(beta, width=5, call_cnt=3)
    #     print(makesstr(calls_b1, 'E'))
    calls_b2 = extendDownstream(beta, calls_b1, width=4)
    calls_b3 = extendUpstream(beta, calls_b2, width=4)

    avg_a = calcRegionAverage(alpha, calls_a3)
    avg_b = calcRegionAverage(beta, calls_b3)

    diff_a = [avg_a[i] - avg_b[i] for i in range(len(avg_a))]
    diff_b = [avg_b[i] - avg_a[i] for i in range(len(avg_a))]

    #     print("Step 5")
    calls_a4 = checkSupport(calls_a3, diff_a)
    calls_b4 = checkSupport(calls_b3, diff_b)
    #     print(makesstr(calls_a4, 'H'))
    #     print(makesstr(calls_b4, 'E'))

    #     print("Final prediction")
    # Combine calls and replace remaining '-' symbols to C (coil)
    prediction = predicted_sequence = ''.join(
        ['H' if alpha else ('E' if beta else 'C') for (alpha, beta) in zip(calls_a4, calls_b4)])
    #     print(prediction)

    position = 0
    for element in structure.sequence:
        if element == "H" and calls_a4[position]:
            atp += 1  # matched alpha helix
        if element != "H" and not (calls_a4[position]):
            atn += 1  # matched beta sheet
        if element != "H" and calls_a4[position]:
            afp += 1  # matched beta sheet
        if element == "H" and not (calls_a4[position]):
            afn += 1  # matched beta sheet
        position += 1

    position = 0
    for element in structure.sequence:
        if element == "E" and calls_b4[position]:
            btp += 1  # matched alpha helix
        if element != "E" and not (calls_b4[position]):
            btn += 1  # matched beta sheet
        if element != "E" and calls_b4[position]:
            bfp += 1  # matched beta sheet
        if element == "E" and not (calls_b4[position]):
            bfn += 1  # matched beta sheet
        position += 1

    position = 0
    for element in structure.sequence:
        if element == "C" and not (calls_a4[position]) and not (calls_b4[position]):
            ctp += 1  # matched alpha helix
        if element != "C" and ((calls_a4[position]) or (calls_b4[position])):
            ctn += 1  # matched beta sheet
        if element != "C" and not (calls_a4[position]) and not (calls_b4[position]):
            cfp += 1  # matched beta sheet
        if element == "C" and ((calls_a4[position]) or (calls_b4[position])):
            cfn += 1  # matched beta sheet
        position += 1

    tp = atp + btp + ctp
    tn = atn + btn + ctn
    fp = afp + bfp + cfp
    fn = afn + bfn + cfn

####### Accuracy calculations #######
print("alpha-helices True  Positive = %d" % atp)
print("alpha-helices True  Negative = %d" % atn)
print("alpha-helices False Positive = %d" % afp)
print("alpha-helices False Negative = %d" % afn)
print("alpha-helix Accuracy %.2f%%" % ((float((atp + atn) / (atp + atn + afp + afn)) * 100)))

print("beta-strands True  Positive = %d" % btp)
print("beta-strands True  Negative = %d" % btn)
print("beta-strands False Positive = %d" % bfp)
print("beta-strands False Negative = %d" % bfn)
print("beta-strands Accuracy %.2f%%" % ((float((btp + btn) / (btp + btn + bfp + bfn)) * 100)))

print("coil True  Positive = %d" % ctp)
print("coil True  Negative = %d" % ctn)
print("coil False Positive = %d" % cfp)
print("coil False Negative = %d" % cfn)
print("coil Accuracy %.2f%%" % ((float((ctp + ctn) / (ctp + ctn + cfp + cfn)) * 100)))

print("True  Positive = %d" % tp)
print("True  Negative = %d" % tn)
print("False Positive = %d" % fp)
print("False Negative = %d" % fn)
print("Accuracy %.2f%%" % ((float((tp + tn) / (tp + tn + fp + fn)) * 100)))