from allele_seq import *
from pseudo_seq import *
from scoring import *
from read_external_test import *
from foutput import *
from keras.models import model_from_json
from math import log
import scipy.stats as ss
from sklearn.metrics import roc_curve, auc

def main_external_testing(global_args):
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
 
    dataset_dates = ["20180511","20170901","20170323","20161209","20160503","20160219","20150807","20150731","20150717","20150619", "20150515"]

    foutput("dataste_date \t dataset \t allele \t len \t SRCC \t AUC\t", output_path)    
    for dataset_date in dataset_dates:
        #Read the data of external dataset
        path_external = main_dir + "IEDB_benchmarking_datasets/"+dataset_date+".txt"
        external_dict = read_external_test(path_external,pseq_dict,global_args)
        #Test the model on these datasets
        for dataset in external_dict.keys():
            for allele in external_dict[dataset].keys():
                #Peptides of different lengths are sorted into different datasets
                for len_pep in external_dict[dataset][allele].keys():
                    external_dataset = external_dict[dataset][allele][len_pep]
                    #Exclude the datasets with less than five samples
                    if len(external_dataset[0]) < 5:
                        continue
                    #Use the model to make prediction scores for the peptides
                    val_scores = scoring(models, [np.array(external_dataset[1]),np.array(external_dataset[2])])   
                    print np.shape(val_scores), np.shape(external_dataset[3])
                    #Use the prediction scores and exerimental measurements to calculate AUROC and SRCC
                    test_label = [1 if aff > 1-log(500)/log(50000) else 0 for aff in external_dataset[3]]
                    fpr, tpr, thresholds = roc_curve(test_label,val_scores)
                    roc_auc = auc(fpr, tpr)  
                    foutput(dataset_date+"\t"+dataset+"\t"+allele+"\t"+str(len_pep)+"\t"\
                                   +str(ss.spearmanr(val_scores,external_dataset[3])[0])+"\t"+str(roc_auc),
                                   output_path)

    
