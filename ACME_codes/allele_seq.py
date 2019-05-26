import numpy as np
import re

def allele_seq(path):
    '''
    Read the sequece of each allele. If there are multiple subtypes, choose the first one.
    Args:
        1. path: path to the data file binding_data_train.txt
    Return values:
        1. A dictionary whose keys are the name of MHC alleles and the corresponding dict
        values are amino acid sequences of those alleles.
    '''
    seq_dict = {}
    f = open(path,"r")
    allele = None
    for line in f:
        if line[0] == ">":#A new allele
            match = re.search("(\w\*\d+:\d+)",line)#The standard allele id are like "A*01:01:..."
            #For alleles with the same two digits, like A*01:01:01 and A*01:01:02, we take the first one as the representative
            allele = None#While reading the sequence of the same alleles from different
            #lines of the file, allele is not None, so that the sequences of each line
            #will be added to the end of the correspondong sequence
            #Some times we run into alleles with incorrect names, so that allele is set to None
            #and match == None so allele will not be reset, then the following lines of sequences 
            #will not be recorded
            if match != None:#If the current allele has a name with the correct format
                if match.groups()[0] not in seq_dict.keys():
                    allele = match.groups()[0]#A new allele
                    seq_dict[allele] = ""#And its sequence
        elif allele != None:
            seq_dict[allele] = seq_dict[allele] + line[:-1]
            #Each line contains only 60 redidues, so add the sequence of the current line 
            #to the end of the corresponding sequence
            
    for allele in list(seq_dict.keys()):
        if len(seq_dict[allele]) < len(seq_dict['B*07:02']):
            #Some sequences lack certain parts like the leader peptide, and cannot
            #be aligned to other sequences well. The ones longer than B*07:02 can 
            #be aligned well with the majority of HLA A and B alleles, (see this in uniprot)
            seq_dict.pop(allele)      
            
    return seq_dict
