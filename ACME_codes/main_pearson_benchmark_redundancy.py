from allele_seq import *
from read_external_train import *
from pseudo_seq import *
from read_external_test import *
from foutput import *

def main_pearson_benchmark_redundancy(global_args):
    [blosum_matrix, aa, main_dir, output_path] = global_args
    path_train = main_dir+ "binding_data/binding_data.txt"
    path_seq = main_dir+ "HLA_A_B.txt"
    seq_dict = allele_seq(path_seq)
    pseq_dict = pseudo_seq(seq_dict, global_args)
    path_external = main_dir+ "binding_data/external_training_set.txt"#This is the one used in NetMHCpan4    
    Pearson_dict = read_external_train(path_external, pseq_dict,global_args)
    for allele in Pearson_dict.keys():
        #Extract the seuqence in the form of strings
        Pearson_dict[allele] = [dt[4] for dt in Pearson_dict[allele]]
    dataset_dates = ["20170901","20170323","20161209","20160503","20160219","20150807","20150731","20150717","20150619", "20150515"]

    print("dataste_date \t dataset \t allele \t ")    
    for dataset_date in dataset_dates:
        #Read the data of external dataset
        path_external = main_dir + "IEDB_benchmarking_datasets/"+dataset_date+".txt"
        external_dict = read_external_test(path_external,pseq_dict,global_args)
        #Test the model on these datasets
        for dataset in external_dict.keys():
            for allele in external_dict[dataset].keys():
                print(allele)
                for len_pep in external_dict[dataset][allele].keys():
                    all_pep = 0
                    overlap_pep = 0
                    for pep in external_dict[dataset][allele][len_pep][0]:
                        #print(pep)
                        all_pep += 1
                        if allele in Pearson_dict.keys() and pep in Pearson_dict[allele]:
                            overlap_pep += 1
                    print(all_pep, overlap_pep)

