def pseudo_seq(seq_dict,global_args):
    '''
    Generate the pseudo sequence for each allele
    Args:
        1. seq_dict: the output of allele_seq()
    Return values:
        1. A dictionary whose keys are the name of MHC alleles and the corresponding dict
        values are the pseudo-sequences of those alleles.    
    '''
    [blosum_matrix, aa, main_dir, output_path] = global_args

    pseq_dict = {}#pseudo sequence dictionary

    #First exclude the alleles that cannot align to the majority of the alleles
    #Most of the alleles can be directly aligned to each other, with a few exceptions having particularly short
    #Sequences and low homogeneity with other alleles
            
    #Remove the sequence of the signal peptide
    for allele in seq_dict.keys():
        seq_dict[allele] = seq_dict[allele][24:]
    
    #Actual indices of selected residues
    residue_indices = [7,9,24,45,59,62,63,66,67,69,70,73,74,76,77,80,81,84,95,97,99,114,116,118,143,147,150,
                       152,156,158,159,163,167,171]
    
    #Indices of these residues in a python list
    residue_indices = [i - 1 for i in residue_indices]
    
    #Now encode the MHC sequences into pseudo-sequences.
    for allele in seq_dict.keys():
        new_pseq = []
        pseq =""
        for index in residue_indices:
            pseq+= seq_dict[allele][index]
            new_pseq.append(blosum_matrix[aa[seq_dict[allele][index]]]) 
        pseq_dict[allele] = new_pseq
    
    return pseq_dict
