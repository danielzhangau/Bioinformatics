from sequence import *
from prob import *

cf_dict = {
# Chou-Fasman table
#     P(H), P(E), P(T),    f(i), f(i+1), f(i+2), f(i+3)
'A': ( 142,   83,   66,   0.060,  0.076,  0.035,  0.058 ),    # Alanine
'R': (  98,   93,   95,   0.070,  0.106,  0.099,  0.085 ),    # Arginine
'N': ( 101,   54,  146,   0.147,  0.110,  0.179,  0.081 ),    # Aspartic Acid
'D': (  67,   89,  156,   0.161,  0.083,  0.191,  0.091 ),    # Asparagine
'C': (  70,  119,  119,   0.149,  0.050,  0.117,  0.128 ),    # Cysteine
'E': ( 151,   37,   74,   0.056,  0.060,  0.077,  0.064 ),    # Glutamic Acid
'Q': ( 111,  110,   98,   0.074,  0.098,  0.037,  0.098 ),    # Glutamine
'G': (  57,   75,  156,   0.102,  0.085,  0.190,  0.152 ),    # Glycine
'H': ( 100,   87,   95,   0.140,  0.047,  0.093,  0.054 ),    # Histidine
'I': ( 108,  160,   47,   0.043,  0.034,  0.013,  0.056 ),    # Isoleucine
'L': ( 121,  130,   59,   0.061,  0.025,  0.036,  0.070 ),    # Leucine
'K': ( 114,   74,  101,   0.055,  0.115,  0.072,  0.095 ),    # Lysine
'M': ( 145,  105,   60,   0.068,  0.082,  0.014,  0.055 ),    # Methionine
'F': ( 113,  138,   60,   0.059,  0.041,  0.065,  0.065 ),    # Phenylalanine
'P': (  57,   55,  152,   0.102,  0.301,  0.034,  0.068 ),    # Proline
'S': (  77,   75,  143,   0.120,  0.139,  0.125,  0.106 ),    # Serine
'T': (  83,  119,   96,   0.086,  0.108,  0.065,  0.079 ),    # Threonine
'W': ( 108,  137,   96,   0.077,  0.013,  0.064,  0.167 ),    # Tryptophan
'Y': (  69,  147,  114,   0.082,  0.065,  0.114,  0.125 ),    # Tyrosine
'V': ( 106,  170,   50,   0.062,  0.048,  0.028,  0.053 )}    # Valine

seq = 'PPDFVYYFKGMCYFTNGTERVRLVTRYIYNREEYARFDSDVGVYRAVTPLGPPAAEYWNSQKEVLERTRAELDTVCRHNYQLELRT'

# predicting the secondary structure
predicted_seq = ''
myprot = Sequence(seq, Protein_Alphabet, 'seq')
for char in myprot:
    if cf_dict[char][0] >= 100 or cf_dict[char][1] >= 100:
        if cf_dict[char][0] > cf_dict[char][1]:
            predicted_seq += 'H'
        else:
            predicted_seq += 'E'
    else:
        predicted_seq += 'C'
print("predicted_seq: ", predicted_seq)

# This is the actual (3-class) secondary structure for the protein above:
ss3 = 'CCCCEEEEEEEEEEECCCCEEEEEEEEEECCEEEEEEECCCCCEEECCCCCHHHHHHHHCCHHHHHHHHHHHHHCHHHHHHHHHHC'

# class E
tp = 0
tn = 0
fp = 0
fn = 0
position = 0
for element in ss3:
    if element == "E" and predicted_seq[position] == "E":
        tp += 1 # matched alpha helix
    if element != "E" and predicted_seq[position] != "E":
        tn += 1 # matched beta sheet
    if element != "E" and predicted_seq[position] == "E":
        fp += 1 # matched beta sheet
    if element == "E" and predicted_seq[position] != "E":
        fn += 1 # matched beta sheet
    position += 1
Se = tp/(tp + fn)
print("sensitivity: ", Se)
Sp = tn/(tn + fp)
print("specificity: ", Sp)