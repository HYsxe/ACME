from allele_seq import *
from pseudo_seq import *
from read_binding_data import *
from read_external_train import *
from read_validation_data import *
from preparing_data import *
from model_training import *
from foutput import *
from simplified_panANN_training import *
from simplified_panCNN_training import *
import copy

def main_leave_one_out(global_args):
    [blosum_matrix, aa, main_dir, output_path_1, output_path_2] = global_args
    global_args = [blosum_matrix, aa, main_dir, output_path_2]
    path_train = main_dir+ "binding_data/binding_data.txt"
    path_seq = main_dir+ "HLA_A_B.txt"
    seq_dict = allele_seq(path_seq)
    pseq_dict = pseudo_seq(seq_dict, global_args)
    data_dict = read_binding_data(path_train,pseq_dict,global_args)
    n_splits = 5
    path_external = main_dir+ "binding_data/external_training_set.txt"#This is the one used in NetMHCpan4
    #Published in "MHC class I associated peptides derive from selective regions of the human genome"
    external_dict = read_external_train(path_external, pseq_dict,global_args)
    path_val = main_dir+ "binding_data/binding_data_val.txt"
    validation_data, validation_target = read_validation_data(path_val,pseq_dict,global_args)   

    alleles = ["A*02:01"]
    
    #Remove the redundant data between the training sets
    for allele in alleles:
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

    for allele in sorted(data_dict.keys()):
        if len(data_dict[allele]) > 100:
            print allele,len(data_dict[allele])

    #Leave-one-out training
    foutput("allele\tPCC\tAUROC\tACC", output_path_1)
    for left_out_allele in  ["B*48:01"]:
        foutput(left_out_allele, output_path_1)
        
        #Training data and test data
        temp_data_dict = copy.deepcopy(data_dict)
        [test_pep, test_mhc, test_target] = [[i[j] for i in data_dict[left_out_allele]] for j in range(3)] 
        print(np.shape(test_pep), np.shape(test_mhc), np.shape(test_target))
        #Remove the data of the left out allele
        temp_data_dict.pop(left_out_allele, None)
        
        #Model training
        models = []
        for cross_validation in range(1):            
            training_data, test_dicts = preparing_data(temp_data_dict, n_splits)    
            for partition in range(5):
                training_data_partition = training_data[partition]
                models.extend(model_training(np.array(training_data_partition), \
                        np.array(validation_data),np.array(validation_target),global_args, n_estimators = 5))
        #Test the performance of the models
        pcc, roc_auc, max_acc = model_eval(models,[np.array(test_pep),np.array(test_mhc)], np.array(test_target))
        foutput(left_out_allele + "\t" + str(pcc)+"\t"+str(roc_auc)+"\t"+str(max_acc), output_path_1)
    

