from sequence import *
from phylo import *

my_tree = readNewick('tyrosine_protein_kinase_fer.newick') # Read the newick file
mynode1 = my_tree.findLabel("A0A1A8AU08")
print(mynode1) # Look at selected node after loading the tree
mynode2 = my_tree.findLabel("G3PS36")
print(mynode2)

# print all the ancestor of the node A0A1A8AU08 with distance
mynode1 = my_tree.findLabel("A0A1A8AU08") # Get the node with the same name as the sequence
a_node1 = my_tree.getAncestorsOf(mynode1, transitive = True) # Get the direct ancestor for the node
for i in a_node1:
    print('has-distance', i.dist, 'from ancestor', i.label) # Print the evolutionary distance

# print all the ancestor of the node G3PS36 with distance
mynode2 = my_tree.findLabel("G3PS36") # Get the node with the same name as the sequence
a_node2 = my_tree.getAncestorsOf(mynode2, transitive = True) # Get the direct ancestor for the node
for i in a_node2:
    print('has-distance', i.dist, 'from ancestor', i.label) # Print the evolutionary distance

# find all common ancestors
common_ancestors = []
for x in a_node1:
    for y in a_node2:
        if x.label == y.label:
            common_ancestors.append(x)
print(common_ancestors)
# find the nearest one
nearest_dis = 100
nearest = ''
for i in common_ancestors:
    if i.dist < nearest_dis:
        nearest = i
print("nearest: ", nearest.label)