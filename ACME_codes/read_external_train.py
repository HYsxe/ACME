import re
import numpy as np
from math import log

def read_external_train(path_external,pseq_dict, global_args):
    '''
    read binding data for peptide-MHC pairs, which are used to train the network
    Data obtained from "MHC class I associated peptides derive from selective 
    regions of the human genome"
    Args:
        1. path_external: path to the data file external_training_set.txt
        2. pseq_dict: the output of pseudo_seq()
    Return values:
        1. external_dict, format: {MHC_allele_name:
            [encoded pep sequence, encoded MHC pseudo sequence, len of pep, affinity]}
    '''
    [blosum_matrix, aa, main_dir, output_path] = global_args
    external_dict = {}
    f = open(path_external)
    invalid = 0
    valid = 0
    for line in f:
        info = re.split("\t", line[:-1])
        pep = info[0]#Sequence of the peptide in the form of a string, like "AAVFPPLEP"
        affinity = 1-log(float(info[2]))/log(50000)
        allele = info[3][4:8]+":"+info[3][-2:]
        if allele in pseq_dict.keys():
            valid += 1
            pep_blosum = []#Encoded peptide seuqence
            #Encode the peptide sequence in the 1-12 columns, with the N-terminal aligned to the left end
            #If the peptide is shorter than 12 residues, the remaining positions on
            #the rightare filled will zero-padding
            for residue_index in range(12):
                if residue_index < len(pep):
                    pep_blosum.append(blosum_matrix[aa[pep[residue_index]]])
                else:
                    pep_blosum.append(np.zeros(20))
            #Encode the peptide sequence in the 13-24 columns, with the C-terminal aligned to the right end
            #If the peptide is shorter than 12 residues, the remaining positions on
            #the left are filled will zero-padding
            for residue_index in range(12):
                if 12 - residue_index > len(pep):
                    pep_blosum.append(np.zeros(20)) 
                else:
                    pep_blosum.append(blosum_matrix[aa[pep[len(pep) - 12 + residue_index]]])
            new_data = [pep_blosum, pseq_dict[allele], affinity, len(pep), pep]
            #new_data = [encoded pep sequence, encoded MHC pseudo sequence, affinity, len of pep, peptide_sequence]
            if allele not in external_dict.keys():
                external_dict[allele] = [new_data]
            else:
                external_dict[allele].append(new_data)
        else:
            invalid += 1
    print "valid", valid, "invalid", invalid
    return external_dict