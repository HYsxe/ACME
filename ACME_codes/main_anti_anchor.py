import numpy as np
import random
import copy
from allele_seq import *
from pseudo_seq import *
from scoring import *
from foutput import *
from keras.models import model_from_json
import scipy.stats as ss
from read_proteome import *
from protein_scanning import *

def main_anti_anchor(global_args):
    '''
    Test whether or not a certain residue would be detrimental to binding.
    A set of random peptides are modifies in two different ways:
        (1) Replacing their original residue at the position in question with the anti-anchor residue.
        (2) Replacing their original residue at the same position with a random residue.
    See if the first type of modification leads to lower predicted scores
    Args:
        1.global arguments
    Return values:
        1. Write the results to a file: 
            (1)average/median score after two different kinds of modifications.
            (2)p-value of Kruskal test.
    '''
    [blosum_matrix, aa, main_dir, output_path] = global_args
    allele = "B*35:01"
    anti_anchor_aa = "K"#The amino acid type to be tested
    anti_anchor_pos = 7#The position of the anti-anchor
    foutput("anti-anchor: "+ allele + " " + anti_anchor_aa + " " + str(anti_anchor_pos), output_path)
    path_train = main_dir+ "binding_data/binding_data.txt"
    path_seq = main_dir+ "HLA_A_B.txt"
    seq_dict = allele_seq(path_seq)
    pseq_dict = pseudo_seq(seq_dict, global_args)   
    
    #Load the models trained previously
    models = []
    for i in range(25):
        json_f = open(main_dir + "models/model_"+str(i)+".json", 'r')
        loaded_model_json = json_f.read()
        json_f.close()
        loaded_model = model_from_json(loaded_model_json)
        loaded_model.load_weights((main_dir + "models/model_"+str(i)+".h5"))
        models.append(loaded_model)      

    #Randomly sample peptides from the proteome
    proteome_path = main_dir + "Homo_sapiens.GRCh38.pep.all.fa"
    proteome = read_proteome(proteome_path)    
    peptides = protein_scanning(proteome, global_args)
    
    #Modify the peptides in two ways
    peptides_1 = peptides
    peptides_2 = copy.deepcopy(peptides)
    #Replacing the original residue with the anti-anchor residue
    for pep in peptides_1:
        pep[anti_anchor_pos - 1] = blosum_matrix[aa[anti_anchor_aa]]
    #Replacing the original residue with a random residue
    for pep in peptides_2:
        pep[anti_anchor_pos - 1] = blosum_matrix[random.randint(0,19)]
        
    scores_1 = scoring(models, [np.array(peptides_1), \
                                np.array([pseq_dict[allele] for i in range(len(peptides_1))])])
    scores_2 = scoring(models, [np.array(peptides_2), \
                                np.array([pseq_dict[allele] for i in range(len(peptides_2))])])
    
    foutput("anti-anchor: " + str(np.median(scores_1)) + " random replacement: " + str(np.median(scores_2)) +"\n"\
            + str(ss.kruskal(scores_1, scores_2)), output_path)
    return