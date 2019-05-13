from allele_seq import *
from pseudo_seq import *
from read_binding_data import *
from foutput import *
from read_external_train import *
from allele_list import *
import numpy as np
import distance
  
def main_MHC_clustering(global_args):
    #Read MHC sequence data and MHC-pep binding data
    [blosum_matrix, aa, main_dir, output_path] = global_args
    path_train = main_dir+ "binding_data/binding_data.txt"
    path_seq = main_dir+ "HLA_A_B.txt"
    seq_dict = allele_seq(path_seq)
    pseq_dict = pseudo_seq(seq_dict, global_args)
    data_dict = read_binding_data(path_train,pseq_dict,global_args)
    path_external = main_dir+ "binding_data/external_training_set.txt"#This is the one used in NetMHCpan4
    #Published in "MHC class I associated peptides derive from selective regions of the human genome"
    external_dict = read_external_train(path_external, pseq_dict,global_args)

    #Remove the redundant data between the training sets
    for allele in sorted(data_dict.keys()):
        if allele in external_dict.keys():
            print allele
            data_dict[allele] = data_dict[allele] + external_dict[allele]
            unique_seq = []
            unique_data = []
            for dt in data_dict[allele]:
                if dt[4] not in unique_seq:
                    unique_data.append(dt)
                    unique_seq.append(dt[4])
            print "unique", len(unique_data)
            data_dict[allele] = unique_data

    #Count the number of samples for each allele
    dataset_size = [len(data_dict[allele]) for allele in sorted(data_dict.keys())]
    alleles = sorted(data_dict.keys())
    n_alleles = len(alleles)
    similarity_matrix = np.zeros((n_alleles, n_alleles))
    for i in range(n_alleles):
        for j in range(n_alleles):
            similarity_matrix[i][j] = distance.levenshtein(pseq_dict[alleles[i]],\
                             pseq_dict[alleles[j]])
    
    #Measuring how much information one allele can obtain from the data of other alleles
    #We assume the information obtained is proportional to the size of the dataset and
    #inversely prportional to the editing distance between the MHC seuqences of the two alleles
    #To avoid division by zero, we add a smoothing term
    reference_info = []
    for i in range(n_alleles):
        sum_info = 0
        for j in range(n_alleles):
            if i != j:#Exclude the allele in question itself
                sum_info += dataset_size[j]*pow((similarity_matrix[i][j])+1, -2)
        reference_info.append(sum_info)
    for i in range(n_alleles):
        foutput(alleles[i] + "\t" + str(reference_info[i]), output_path)