from allele_seq import *
from pseudo_seq import *
from scoring import *
from foutput import *
from read_prediction_input import *
from keras.models import model_from_json
from math import log
from sklearn.metrics import roc_curve, auc
import scipy.stats as ss
import numpy as np

def main_binding_prediction(global_args):
    [blosum_matrix, aa, main_dir, output_path] = global_args
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
 
    #Read the MHC alleles and peptide sequences
    [peptides, alleles] = read_prediction_input(main_dir + "binding_prediction/prediction_input.txt")

    #Encode the peptides
    input_pep = []
    for pep in peptides:
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
        input_pep.append(pep_blosum)

    #Encode the MHC alleles
    input_mhc = [pseq_dict[allele] for allele in alleles]
    
    #Making predictions
    scores = scoring(models, [np.array(input_pep),np.array(input_mhc)])   
    
    #Output
    for i in range(len(scores)):
        foutput(peptides[i] + "\t" + alleles[i] + "\t" + str(scores[i]), output_path)
    