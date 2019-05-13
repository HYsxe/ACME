import numpy as np
import keras
import pickle
from keras.models import Sequential
from keras.layers import *
from keras.models import Model
import tensorflow as tf
from keras import backend as K
from foutput import *
from model_eval import *
from math import log
from model_performance import *

def cross_validation_training_without_CNN_fc(training_data, test_dict, validation_data, validation_target, global_args):
    '''
    Cross validation training using only the attention module
    Args:
        1. training_data: Data for training, should be one split of 
                    the output of preaparing_data()
        2. test_dict: Test data for evaluating the models performance in 
                    cross-validation.This should be one dict in test_dicts.
                    The latter is the output of preparing_data()                    
        3. validation_data: Data for validation, should be the output of 
                    read_binding_data_val()            
        4. validation_target: Target for validation, should be the output of 
                    read_binding_data_val()  
    Return values:
        1. models: A list of trained models.
    '''
    [blosum_matrix, aa, main_dir, output_path] = global_args
    [training_pep, training_mhc, training_target] = [[i[j] for i in training_data] for j in range(3)]
    validation_pep, validation_mhc = [i[0] for i in validation_data], [i[1] for i in validation_data] 
    #Specifying the hyperparameters: number of filters and size of network layers.
    filters, fc_size = 512, 256 
    foutput(str(filters), output_path)
    #Now train the new model
    kernel_size = 3
    models = []
    while len(models) < 5:
        inputs_1 = Input(shape = (np.shape(training_pep[0])[0],20))
        inputs_2 = Input(shape = (np.shape(training_mhc[0])[0],20))
        #Initial feature extraction using a convolutional layer
        pep_conv = Conv1D(filters,kernel_size,padding = 'same',activation = 'relu',strides = 1)(inputs_1)
        mhc_conv_1 = Conv1D(filters,kernel_size,padding = 'same',activation = 'relu',strides = 1)(inputs_2)
        #The attention module
        mhc_attention_weights = Flatten()(TimeDistributed(Dense(1))(mhc_conv_1))
        pep_attention_weights = Flatten()(TimeDistributed(Dense(1))(pep_conv))
        mhc_attention_weights = Activation('softmax')(mhc_attention_weights)
        pep_attention_weights = Activation('softmax')(pep_attention_weights)        
        mhc_conv_permute = Permute((2,1))(mhc_conv_1)
        pep_conv_permute = Permute((2,1))(pep_conv)
        mhc_attention = Dot(-1)([mhc_conv_permute, mhc_attention_weights])
        pep_attention = Dot(-1)([pep_conv_permute, pep_attention_weights])
        #Concatenating the output of the two modules
        merge = Concatenate()([mhc_attention,pep_attention])
        #Output of the model
        fc = Dense(fc_size,activation = "relu")(merge)
        out = Dense(1,activation = "sigmoid")(fc)
        model = Model(inputs=[inputs_1, inputs_2],outputs=out)  
        model.summary()
        model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mse'])
        poor_init = False
        for n_epoch in range(39):
            #Start training
            model.fit([np.array(training_pep),np.array(training_mhc)], np.array(training_target), batch_size=64,\
                       epochs = 1)
            #Evaluate the performance of the model (on the validation set) after each training epoch. (note that this is not cross-validation)
            pcc, roc_auc, max_acc = model_eval(model,[np.array(validation_pep),np.array(validation_mhc)],np.array(validation_target))
            foutput(str(n_epoch)+"\t"+str(pcc)+"\t"+str(roc_auc)+"\t"+str(max_acc), output_path)
            #If the pcc after the first epoch is very low, the network is unlikely to be trained well, so start again
            if n_epoch == 0 and not (pcc > 0):
                poor_init = True
                break
        if poor_init:
            continue
        if pcc > 0:
            foutput("Model Adopted", output_path)
            models.append(model)        
    performance_dict = model_performance(models, test_dict, global_args)

    return performance_dict
