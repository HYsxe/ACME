import random
import numpy as np

def protein_scanning(proteome, global_args):
    '''
    Randomly sample peptides from the proteome
    Args:
        1. proteome: A dictionary of the human proteome.
        Output of the function read_proteome
    Return values:
        1. peptides: Sampled peptides.
    '''
    [blosum_matrix, aa, main_dir, output_path] = global_args
    #randomly generate peptides from the proteome
    #Encoded peptides
    peptides = []
    #peptide seuqences
    pep_seq = []
    #Sample 10000 peptides
    while len(peptides) < 10000:
        #Randomly choose a protein
        protein = random.choice(proteome.values())
        if len(protein) < 9:
            continue
        pep_start = random.randint(0, len(protein) - 9)
        pep = protein[pep_start:pep_start + 9]
        if pep in pep_seq:
            continue
        pep_blosum = []
        try:
            pep_blosum = []#Encoded peptide seuqence
            for residue_index in range(12):
                #Encode the peptide sequence in the 1-12 columns, with the N-terminal aligned to the left end
                #If the peptide is shorter than 12 residues, the remaining positions on
                #the rightare filled will zero-padding
                if residue_index < len(pep):
                    pep_blosum.append(blosum_matrix[aa[pep[residue_index]]])
                else:
                    pep_blosum.append(np.zeros(20))
            for residue_index in range(12):
                #Encode the peptide sequence in the 13-24 columns, with the C-terminal aligned to the right end
                #If the peptide is shorter than 12 residues, the remaining positions on
                #the left are filled will zero-padding
                if 12 - residue_index > len(pep):
                    pep_blosum.append(np.zeros(20)) 
                else:
                    pep_blosum.append(blosum_matrix[aa[pep[len(pep) - 12 + residue_index]]])
        except KeyError:
            continue
        pep_seq.append(pep)
        peptides.append(pep_blosum)
    return peptides