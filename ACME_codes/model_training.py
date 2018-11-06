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

def model_training(training_data, validation_data, validation_target, global_args, n_estimators):
    '''
    Training new models
    Args:
        1. training_data: Data for training, should be one split of 
                    the output of preaparing_data()
        2. validation_data: Data for validation, should be the output of 
                    read_binding_data_val()            
        3. validation_target: Target for validation, should be the output of 
                    read_binding_data_val()  
        4. n_estimators: Number of estimators in the ensemble
    Return values:
        1. models: A list of trained models.
    '''
    [blosum_matrix, aa, main_dir, output_path] = global_args
    [training_pep, training_mhc, training_target] = [[i[j] for i in training_data] for j in range(3)]
    validation_pep, validation_mhc = [i[0] for i in validation_data], [i[1] for i in validation_data] 
    #Specifying the hyperparameters: number of filters and size of network layers.
    filters, fc1_size, fc2_size, fc3_size= 128, 256, 64, 2
    foutput(str(filters)+"\t"+str(fc1_size)+"\t"+str(fc2_size)+"\t"+str(fc3_size), output_path)
    #Now train the new model
    kernel_size = 3
    models = []
    while len(models) < n_estimators:
        inputs_1 = Input(shape = (np.shape(training_pep[0])[0],20))
        inputs_2 = Input(shape = (np.shape(training_mhc[0])[0],20))
        #Initial feature extraction using a convolutional layer
        pep_conv = Conv1D(filters,kernel_size,padding = 'same',activation = 'relu',strides = 1)(inputs_1)
        pep_maxpool = MaxPooling1D()(pep_conv)
        mhc_conv_1 = Conv1D(filters,kernel_size,padding = 'same',activation = 'relu',strides = 1)(inputs_2)
        mhc_maxpool_1 = MaxPooling1D()(mhc_conv_1)
        #The convolutional module
        mhc_conv_2 = Conv1D(filters,kernel_size,padding = 'same',activation = 'relu',strides = 1)(mhc_maxpool_1)
        mhc_maxpool_2 = MaxPooling1D()(mhc_conv_2)
        flat_pep_0 = Flatten()(pep_conv)
        flat_pep_1 = Flatten()(pep_conv)
        flat_pep_2 = Flatten()(pep_conv)
        flat_mhc_0 = Flatten()(inputs_2)
        flat_mhc_1 = Flatten()(mhc_maxpool_1)
        flat_mhc_2 = Flatten()(mhc_maxpool_2)
        cat_0 = Concatenate()([flat_pep_0, flat_mhc_0])
        cat_1 = Concatenate()([flat_pep_1, flat_mhc_1])
        cat_2 = Concatenate()([flat_pep_2, flat_mhc_2])        
        fc1_0 = Dense(fc1_size,activation = "relu")(cat_0)
        fc1_1 = Dense(fc1_size,activation = "relu")(cat_1)
        fc1_2 = Dense(fc1_size,activation = "relu")(cat_2)
        merge_1 = Concatenate()([fc1_0, fc1_1, fc1_2])
        fc2 = Dense(fc2_size,activation = "relu")(merge_1)
        fc3 = Dense(fc3_size,activation = "relu")(fc2)
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
        merge_2 = Concatenate()([mhc_attention,pep_attention,fc3])
        #Output of the model
        out = Dense(1,activation = "sigmoid")(merge_2)
        model = Model(inputs=[inputs_1, inputs_2],outputs=out)  
        model.summary()
        model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mse'])
        poor_init = False
        for n_epoch in range(12):
            #Starts training
            model.fit([np.array(training_pep),np.array(training_mhc)], np.array(training_target), batch_size=64,\
                       epochs = 1)  
            #Evaluate the performance of the model (on the validation set) after each training epoch. (note that this is not cross-validation)
            pcc, roc_auc, max_acc = model_eval(model,[np.array(validation_pep),np.array(validation_mhc)],np.array(validation_target))
            #If the pcc after the first epoch is very low, the network is unlikely to be trained well, so start again
            if n_epoch == 0 and not(pcc > 0.75):
                poor_init = True
                break
            foutput(str(n_epoch)+"\t"+str(pcc)+"\t"+str(roc_auc)+"\t"+str(max_acc), output_path)
        if poor_init:
            continue
        #Only select the models with good performances on the validation set (note that this is not cross-validation)
        if pcc > 0.84:
            foutput("Model Adopted", output_path)
            models.append(model)      
    return models