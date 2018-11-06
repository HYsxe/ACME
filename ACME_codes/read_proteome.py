import re

def read_proteome(path):
    '''
    Read the sequences of the proteins in the human proteome. 
    Sequence data are stored in a fasta file     
    Args:
        1. path: The input file containing sequence data of the proteome
            downloaded from the ensemble biomart FTP:
            ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/.
    Return values:
        1. proteome: A dictionary whose keys are protein ensembl IDs
                    and values are protein sequences
    '''   
    f = open( path, "r" )
    proteome = {}
    seq = ""
    for tmp in f:
        if tmp[0]==">":
            #Start reading the sequence for a new protein
            id = re.search("transcript:(ENST\d+)",tmp)
            seq = ""
        else:
            #Update the sequence. Use len(tmp)-1 because the last char in each line is \n
            seq += tmp[0:len(tmp)-1]
            #update the sequence in the dictionary
            proteome[id.groups()[0]] = seq
    f.close()
    return proteome 
    


    