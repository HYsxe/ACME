import re
import numpy as np
from math import log

def read_external_test(path_external,pseq_dict,global_args):
    '''
    read binding data for peptide-MHC pairs, which are used as external
    test sets to evaluate the generalizability of our model.
    The datasets are downloaded from the IEDB weekly benchmark website
    http://tools.iedb.org/auto_bench/mhci/weekly/
    Args:
        1. path_external: path to the test data file
        2. pseq_dict: the output of pseudo_seq()
    Return values:
        1. external_dict: A dictionary containing the binding data o different
                        test datasets
    '''
    [blosum_matrix, aa, main_dir, output_path] = global_args
    external_dict = {}
    f = open(path_external)
    for line in f:
        info = re.split("\t", line[:-1])
        if info[0] == "" or info[0] == "Date":
            continue
        pep = info[5]
        allele = info[2][4:]
        dataset = info[1]
        pep_len = len(pep)
        #Datasets are categorized accoding their id, allele and peptide lengths
        if dataset not in external_dict:
            external_dict[dataset] = {}
        if allele not in external_dict[dataset].keys():
            external_dict[dataset][allele] = {}
        if pep_len not in external_dict[dataset][allele].keys():
            #[[sequence], [encoded sequence], [allele], [affinity]]
            external_dict[dataset][allele][pep_len] = [[],[],[],[]]
        if pep not in external_dict[dataset][allele][pep_len][0] and allele in pseq_dict.keys():# :
            external_dict[dataset][allele][pep_len][0].append(pep)
            pep_blosum = []
            for residue_index in range(12):
                if residue_index < len(pep):
                    pep_blosum.append(blosum_matrix[aa[pep[residue_index]]])
                else:
                    pep_blosum.append(np.zeros(20))
            for residue_index in range(12):
                if 12 - residue_index > len(pep):
                    pep_blosum.append(np.zeros(20)) 
                else:
                    pep_blosum.append(blosum_matrix[aa[pep[len(pep) - 12 + residue_index]]])
            external_dict[dataset][allele][pep_len][1].append(pep_blosum)
            external_dict[dataset][allele][pep_len][2].append(pseq_dict[allele])
            #For different kinds of affinity measurements, we transform them into values between 0 and 1
            if info[4] == "ic50":
                external_dict[dataset][allele][pep_len][3].append(1-log(float(info[6]))/log(50000))
            elif info[4] == "binary":
                external_dict[dataset][allele][pep_len][3].append(float(info[6]))
    return external_dict