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
config.gpu_options.visible_device_list = "0"
sess = tf.Session(config=config)
K.set_session(sess)
from read_blosum import read_blosum
from main_binding_prediction import *
#Output
from foutput import *

#Path to the Attention_CNN folder
main_dir = "/home/user/ACME_revised/"
aa={"A":0,"R":1,"N":2,"D":3,"C":4,"Q":5,"E":6,"G":7,"H":8,"I":9,"L":10,"K":11,"M":12,"F":13,"P":14,"S":15,"T":16,"W":17,"Y":18,"V":19}     
#Load the blosum matrix for encoding
path_blosum = main_dir + r"blosum50.txt"
blosum_matrix = read_blosum(path_blosum) 
#Making binding predictions
output_path = main_dir + "results/binding_prediction.txt"
main_binding_prediction([blosum_matrix, aa, main_dir, output_path])  
        
