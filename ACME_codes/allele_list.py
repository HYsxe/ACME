import re

def allele_list(path):
    '''
    Read the names of all the alleles whose binding data are available
    Args:
        1. path: path to the data file binding_data_train.txt
    Return values:
        1. allele_list: A list of alleles names.
    '''
    allele_list = []
    f = open(path,"r")
    for line in f:
        #See if the current allele name has the right format. example: HLA-A*01:01
        match = re.search("HLA-(\w.+?:.+?)\t",line)
        if match != None and match.groups()[0] not in allele_list:
            allele_list.append(match.groups()[0])
    return allele_list