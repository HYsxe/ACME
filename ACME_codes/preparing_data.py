from sklearn.model_selection import KFold
import random
import re
import numpy as np

def preparing_data(data_dict, n_splits = 5, test_len = 9):
    '''
    Partition the data for n-fold training or cross-validation
    Args:
        1. data_dict: Dictionary containing the binding data
                     should be the output of read_binding_data()
                     or read_external()
        2. n_splits: Data are split into n_splits parts
                    in cross-validation, one part is used for validation
                    and the rest are used for training.
    Return values:
        1.  training_data: Contains n_splits lists. Each one is the training
                            set for 1 cross-validtion round.
        2.  test_dicts: Dictionarys corresponding to each training set.
                        Keys are alleles, values are the test data for this allele.
        3.  test_len: Length of the peptides to be included in the test set
    Note: KFold needs to be included from sklearn.model_selection
    '''  
    training_data = []
    test_dicts= []
    cross_validation = KFold(n_splits = n_splits)  
    
    #For each partition, initialize the container of data and target
    for split in range(n_splits):
        training_data.append([])
        test_dicts.append({})    
    
    #For each allele, sort the corresponding data into training sets and test sets.
    for allele in data_dict.keys():
        allele_data = data_dict[allele]
        random.shuffle(allele_data)
        allele_data = np.array(allele_data)
        split = 0
        #We only include the alleles with >= data
        if len(allele_data)< 100:
            continue
        #Partition of data
        for training_indices, test_indices in cross_validation.split(allele_data):
            training_data[split].extend(allele_data[training_indices])
            allele_test_data = allele_data[test_indices]
            #We only select peptides of a specific length to be included in the test set
            allele_test_data = [i for i in allele_test_data if i[3] == test_len]
            test_dicts[split][allele] = allele_test_data
            split += 1 

    for split in range(n_splits):        
        random.shuffle(training_data[split])

    return training_data, test_dicts