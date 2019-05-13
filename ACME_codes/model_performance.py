from foutput import *
from model_eval import *

def model_performance(models, test_dict, global_args, alleles = None):
    '''
    Evaluate the allele-specific performance of our model on the test data
    of different alleles
    Args:
        1. models: A model or an ensemble of keras models
        2. test_dict: One dict in test_dicts. test_dicts should be one of the
            outputs of preparing_data()
        3. alleles: If you want to test on specific alleles, speficy the list 
                of alleles of interest here.
    Return values:
        1. performance_dict: A dictionary. Keys are alleles. Values are the 
                    performances on the corresponding alleles.
    '''
    [blosum_matrix, aa, main_dir, output_path] = global_args
    performance_dict = {}
    for allele in sorted(test_dict.keys()):
        #See if we are testing on specified alleles
        if alleles is not None and allele not in alleles:
            continue
        #We only test on datasets with >= 10 data
        if len(test_dict[allele]) < 10:
            continue
        foutput(allele+" "+str(len(test_dict[allele])), output_path)
        print allele
        [test_pep, test_mhc, test_target] = [[i[j] for i in test_dict[allele]] for j in range(3)] 
        #Evaluate the performance of our model.
        pcc, roc_auc, max_acc = model_eval(models,[np.array(test_pep),np.array(test_mhc)], np.array(test_target))
        performance_dict[allele] = [pcc, roc_auc, max_acc]
        
    return performance_dict