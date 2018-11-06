import re

def read_blosum(path):
    '''
    Read the blosum matrix from the file blosum50.txt
    Args:
        1. path: path to the file blosum50.txt
    Return values:
        1. The blosum50 matrix
    '''
    f = open(path,"r")
    blosum = []
    for line in f:
        blosum.append([(float(i))/10 for i in re.split("\t",line)])
        #The values are rescaled by a factor of 1/10 to facilitate training
    f.close()
    return blosum