from sequence import *
from phylo import *

my_tree = readNewick('cyp1a1.tree') # Read the newick file
aln = readClustalFile('cyp1a1.aln', Protein_Alphabet) # Load the alignment

my_tree.putAlignment(aln)   # Associate the tree with the alignment
# Extract the nodes for each clade
left_root_children = my_tree.root.left.getDescendants(transitive = True)
right_root_children = my_tree.root.right.getDescendants(transitive = True)

# Get the corresponding sequences
print("Left child")
left_list = []
right_list = []
for ele in left_root_children:
    if ele.getSequence() is not None:
        left_list.append(ele.sequence)
        print(ele.getSequence())
print("\nRight child")
for ele in right_root_children:
    if ele.getSequence() is not None:
        right_list.append(ele.sequence)
        print(ele.getSequence())

# Build two new alignments for each clade
left_aln = Alignment(left_list)
# Calculate consensus sequences
print("Mammals consensus:", left_aln.getConsensus())
right_aln = Alignment(right_list)
print("Fish consensus:", right_aln.getConsensus())