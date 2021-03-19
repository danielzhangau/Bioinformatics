from translate import* 
from guide import Alphabet  

DNA_Alphabet = Alphabet('ACGT')
RNA_Alphabet = Alphabet('ACGU')
Protein_Alphabet = Alphabet('ACDEFGHIKLMNPQRSTVWY')
Protein_wX = Alphabet('ACDEFGHIKLMNPQRSTVWYX')
Hydrophobic_Alphabet=Alphabet('VILFWCAYHSPGT')

def baseCounts(sequence,Alphabet):
    seq=sequence.upper()
    if Alphabet==DNA_Alphabet:
        DNA={}
        for nuc in DNA_Alphabet:
            DNA[nuc]=seq.count(nuc)
        return DNA
    elif Alphabet==RNA_Alphabet:
        RNA={}
        for nuc in RNA_Alphabet:
            RNA[nuc]=seq.count(nuc)
        return RNA
    elif Alphabet==Protein_Alphabet:
        amino={}
        for prot in Protein_Alphabet:
            amino[prot]=seq.count(prot)
        return amino 
        
exercise1_seq = "AAAACCTCTCTGTTCAGCACTTCCTCTCTCTTGGTCTGGTCTCAACGGTCACCATGGCGAGACCCTTGGAGGAGGCCCTGGATGTAATAGTGTCCACCTTCCACAAATACTCAGGCAACGAGGGTGACAAGTTCAAGCTGAACAAGACAGAGCTCAAGGAGCTACTGACCAGGGAGCTGCCTAGCTTCCTGGGGAGAAGGACAGACGAAGCTGCATTCCA"
result = baseCounts(exercise1_seq, DNA_Alphabet)
print(result)

class MySequence():
    def __init__(self,sequence,alphabet):
        """Myseq=MySequence('ACGTACGT'DNA_Alphabet) there is also Protein_Alphabet,RNA_Alphabet"""
        self.alphabet=alphabet
        self.sequence=sequence

    def getCounts(self):
        """Gets the count of AA and nucleotides"""
        seq=self.sequence 
        if self.alphabet==DNA_Alphabet:
            DNA={}
            for nuc in DNA_Alphabet:
                DNA[nuc]=seq.count(nuc) 
            return DNA
        elif self.alphabet==RNA_Alphabet:
            RNA={}
            for nuc in RNA_Alphabet:
                RNA[nuc]=round(float(seq.count(nuc))/float(len(seq)) ,2)
            return RNA
        elif self.alphabet==Protein_Alphabet:
            amino={}
            for prot in Protein_Alphabet:
                amino[prot]=round(float(seq.count(prot))/float(len(seq)),2)
            return amino
        elif self.alphabet==Hydrophobic_Alphabet:
            hydro_amino={}
            counts=[] 
            for prot in Hydrophobic_Alphabet:
                hydro_amino[prot]=seq.count(prot)
            for prot in Hydrophobic_Alphabet:
                counts.append(hydro_amino[prot])
            return sum(counts)  #returns sum of counts list   
 
        
    def reverseComplement(self):
        """ Return a new sequence: the reverse complement of this sequence.Causes an error if this isn't DNA or RNA. """
        self.sequence.upper()
        newseq=''
        symbols={'A':'T','C':'G','T':'A','G':'C'} # reverse complement dictionary
        if self.alphabet==DNA_Alphabet: # Checks if the sequence is DNA
            for symbol in self.sequence[::-1]:
                newsymbol=symbols[symbol] # uses the reverse complement symbols in dictionary
                newseq+=newsymbol
            return newseq  # returns RC sequences
        else:
            print("ERROR! this is not DNA") # error if not DNA

    def translateDNA(self, readingFrame=0, forwards=True):
        """ Return a new sequence: the amino acid sequence translated from this DNA sequence. Causes an error if this isn't DNA. Optionally the reading frame and direction of translation may be specified. """
        newseq=''
        codons=[]
        translated=[]
        if self.alphabet==DNA_Alphabet: # Checks if the sequence is DNA
            if readingFrame==0:
                i=0
                while i < len(self.sequence):
                    if len(self.sequence[i:i+3])==3:
                        codons.append(self.sequence[i:i+3])
                    i+=3
                for codon in codons:
                    translated.append(standardTranslation[codon])
                for amino in translated:
                    newseq+=amino
                return newseq
            elif readingFrame==1:
                i=1
                while i < len(self.sequence):
                    if len(self.sequence[i:i+3])==3:
                        codons.append(self.sequence[i:i+3])
                    i+=3
                for codon in codons:
                    translated.append(standardTranslation[codon])
                for amino in translated:
                    newseq+=amino
                return newseq
            else:
                i=2
                while i < len(self.sequence):
                    if len(self.sequence[i:i+3])==3:
                        codons.append(self.sequence[i:i+3])
                    i+=3
                for codon in codons:
                    translated.append(standardTranslation[codon])
                for amino in translated:
                    newseq+=amino
                return newseq               
        else:
            print("ERROR! this is not DNA") # error if not DNA 

        
        
    
        


            
            
        
