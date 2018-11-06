from allele_seq import *
from pseudo_seq import *
from read_binding_data import *
from read_external_train import *
from read_validation_data import *
from preparing_data import *
from model_training import *
from foutput import *

def main_model_training(global_args):
    [blosum_matrix, aa, main_dir, output_path] = global_args
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
    
    #Train the models
    models = []
    for cross_validation in range(5):            
        training_data, test_dicts = preparing_data(data_dict, n_splits)   
        for partition in range(5):
            training_data_partition = training_data[partition]
            models.extend(model_training(np.array(training_data_partition), np.array(validation_data),\
                                         np.array(validation_target),global_args, n_estimators = 1))

    #Save model and weights to file
    for i in range(len(models)):
        model = models[i]
        model_json = model.to_json()
        with open(main_dir+ "models/model_"+str(i)+".json", "w") as json_file:
            json_file.write(model_json)
        model.save_weights(main_dir+ "models/model_"+str(i)+".h5")
