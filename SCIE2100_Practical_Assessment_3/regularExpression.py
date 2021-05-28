from sequence import *
from prob import *

# Convert the JASPAR motif into a PWM
distribs = readMultiCount('FHL1.jaspar')
FHL1_pwm = PWM(distribs)
FHL1_pwm.display()

# read fasta file
yeast_promoters = readFastaFile('yeast_promoters.fa')
len(yeast_promoters)

# Create a new empty list yeast_subseqs
yeast_subseqs = []
# identify which yeast promoter sequences have a motif match with a score greater than 2.8
score_dict = {}
for s in yeast_promoters: # yeast_prom is an array of sequences
    if FHL1_pwm.maxscore(s)[0] > 2.8:
        score_dict[s.name] = FHL1_pwm.maxscore(s)[1] # save index only
        # Given each of these sequences, extract the subsequence (with length equal to the PWM)
        # located at the position reported by maxscore and store it in yeast_subseqs
        yeast_subseqs.append((''.join(s.sequence))[FHL1_pwm.maxscore(s)[1]:FHL1_pwm.maxscore(s)[1]+8])
print(score_dict)
print(yeast_subseqs)

# construct a regular expression to represent the FHL1 transcription factor binding site.
regular_expression = ""
char_list = []
for index in range(0,8):
    for item in yeast_subseqs:
        if item[index] not in char_list:
            char_list.append(item[index])
    if len(char_list) > 1:
        regular_expression = regular_expression + '['+ ''.join(char_list) + ']'
    else:
        regular_expression = regular_expression + ''.join(char_list)
    char_list = []
print("regular_expression: ", regular_expression)

tf = re.compile('[CG]ACGC[TA][GACT][GACT]')
for seq in yeast_subseqs:
    m = tf.match(seq)
    if m:
        print(seq, m)