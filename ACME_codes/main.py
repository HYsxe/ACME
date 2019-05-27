import numpy as np
import random
from sklearn.metrics import roc_curve, auc
import re
import copy
import scipy.stats as ss
from math import *
import keras
import pickle
from keras.models import Sequential
from keras.layers import *
from keras.models import Model
import tensorflow as tf
from keras import backend as K
from sklearn.model_selection import KFold
config = tf.ConfigProto()
config.gpu_options.per_process_gpu_memory_fraction = 0.3
config.gpu_options.visible_device_list = "3"
sess = tf.Session(config=config)
K.set_session(sess)
from main_model_training import main_model_training
from read_blosum import read_blosum
from main_cross_validation import *
from main_external_testing import *
from main_test_attention import *
from main_motif import *
from main_leave_one_out import *
from main_MHC_clustering import *
from main_pearson_benchmark_redundancy import *
from main_anti_anchor import *
from main_cross_validation_without_attention import *
from main_cross_validation_without_CNN import *
from main_binding_prediction import *
#Output
from foutput import *

def main(func_ind):
    #Path to the Attention_CNN folder
    main_dir = "/home/huyan/ACME_revised/"
    aa={"A":0,"R":1,"N":2,"D":3,"C":4,"Q":5,"E":6,"G":7,"H":8,"I":9,"L":10,"K":11,"M":12,"F":13,"P":14,"S":15,"T":16,"W":17,"Y":18,"V":19}     
    #Load the blosum matrix for encoding
    path_blosum = main_dir + r"blosum50.txt"
    blosum_matrix = read_blosum(path_blosum)
    if func_ind == 0:
        #Train new models and save them to a folder
        output_path = main_dir + "models/model_training.txt"
        main_model_training([blosum_matrix, aa, main_dir,output_path])
    elif func_ind == 1:
        #5-fold cross valudation training
        output_path = main_dir + "results/cross_validation_9mer.txt"
        main_cross_validation([blosum_matrix, aa, main_dir, output_path])
    elif func_ind == 2:
        #Testing the trained models on external datasets
        output_path = main_dir + "results/testing_on_IEDB_benchmark_datasets.txt"
        main_external_testing([blosum_matrix, aa, main_dir, output_path])
    elif func_ind == 3:
        #Test whether or not the attention module can detect important residues
        output_path = main_dir + "results/attention_test.txt"
        main_test_attention([blosum_matrix, aa, main_dir, output_path])  
    elif func_ind == 4:
        #Generate binding motifs for binders
        output_path = main_dir + "results/binder_motif.txt"
        main_motif([blosum_matrix, aa, main_dir, output_path], mode = 'binder')    
    elif func_ind == 5:
        #Generate binding motifs for non-binders
        output_path = main_dir + "results/non_binder_motif.txt"
        main_motif([blosum_matrix, aa, main_dir, output_path], mode = 'nonbinder')   
    elif func_ind == 6:
        #Making predictions for alleles with no training data
        output_path_1 = main_dir + "results/ACME_leave_one_out.txt"
        output_path_2 = main_dir + "results/ACME_leave_one_out_training.txt"
        main_leave_one_out([blosum_matrix, aa, main_dir, output_path_1, output_path_2]) 
    elif func_ind == 7:
        #Test whether the anti-anchor residues can impair binding
        output_path = main_dir + "results/anti_anchor.txt"
        main_anti_anchor([blosum_matrix, aa, main_dir, output_path])  
    elif func_ind == 8:
        #Cross validation using a model without the attention module
        output_path = main_dir + "results/cross_validation_without_attention_9mer.txt"
        main_cross_validation_without_attention([blosum_matrix, aa, main_dir, output_path])  
    elif func_ind == 9:
        #Cross validation using a model without the convolutional module
        output_path = main_dir + "results/cross_validation_without_CNN_9mer.txt"
        main_cross_validation_without_CNN([blosum_matrix, aa, main_dir, output_path])  
    elif func_ind == 10:
        #Calculating the reference information (RI) of different MHC alleles.
        output_path = main_dir + "results/MHC_RI.txt"
        main_MHC_clustering([blosum_matrix, aa, main_dir, output_path])   

main(0)
        
