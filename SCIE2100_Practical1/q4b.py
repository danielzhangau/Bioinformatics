from sequence import *
## do again on Lipid metabolism
spat = searchSequences('"lipid+metabolism"+AND+organism:3702+AND+length:[700+TO+*]')
print(len(spat))