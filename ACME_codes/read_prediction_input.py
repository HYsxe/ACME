import re
import numpy as np
def read_prediction_input(path):
    '''
    Read the peptide sequences and MHC alleles from a file.
    These are used as the model input to make binding predictions.
    Args:
        1. path: Path to the file containing input data
        input format:
            peptide_seq_1 MHC_allele_1 (separated using tab)
            peptide_seq_2 MHC_allele_2
            ......
    Return values:
        1. input_data:[peptides, MHC_alleles]
    '''
    input_data = [[],[]]
    f = open(path, "r")
    for line in f:
        if len(line) > 10:
            info = re.split("\t", line[:-1])
            input_data[0].append(info[0])
            input_data[1].append(info[1])
    return input_data