from allele_seq import *
from pseudo_seq import *
from read_binding_data import *
from read_validation_data import *
from preparing_data import *
from model_training import *
from foutput import *
from model_performance import *
from cross_validation_training_attention_fc import *
from redundancy_removal import *
      
def main_cross_validation_attention_fc(global_args):
    #Reading sequence data and peptide-MHC binding data
    [blosum_matrix, aa, main_dir, output_path] = global_args
    path_seq = main_dir+ "HLA_A_B.txt"
    seq_dict = allele_seq(path_seq)
    pseq_dict = pseudo_seq(seq_dict, global_args)
    path_train = main_dir+ "binding_data/binding_data_train.txt"
    data_dict = read_binding_data(path_train,pseq_dict,global_args)    
    data_dict = redundancy_removal(data_dict)
    path_val = main_dir+ "binding_data/binding_data_val.txt"
    validation_data, validation_target = read_validation_data(path_val,pseq_dict,global_args)
    
    #Data partition for cross-validation
    n_splits = 5
    training_data, test_dicts = preparing_data(data_dict, n_splits, test_len = 11)  
    print "Finished data loading"
    print "shape of training data", np.shape(training_data)       

    #Cross-validation
    performance_dicts = []
    for split in range(n_splits):
        performance_dict = cross_validation_training_attention_fc(np.array(training_data[split]), test_dicts[split], 
                                          validation_data, validation_target, global_args)
        performance_dicts.append(performance_dict)
        
    #Output the results
    foutput("allele\tPCC\tAUROC\tACC", output_path)
    for allele in sorted(performance_dicts[1].keys()):
        try:
            performances = [perf_dict[allele] for perf_dict in performance_dicts]
            mean_performance = [np.mean(metric) for metric in zip(*performances)]
            foutput(allele+"\t"+str(mean_performance[0])+"\t"+str(mean_performance[1])+"\t"+str(mean_performance[2]), output_path)
        except KeyError:
            pass
