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

def simplified_panANN_training(training_data, validation_data, validation_target, global_args):
    '''
    Training a simplified pan specific CNN model to compare with pan specific ANN model
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
    #Now train the new model
    models = []
    for n_nodes in [60, 70]:
        #A single-layered neural network
        inputs_1 = Input(shape = (np.shape(training_pep[0])[0],20))
        inputs_2 = Input(shape = (np.shape(training_mhc[0])[0],20))
        flat_1 = Flatten()(inputs_1)
        flat_2 = Flatten()(inputs_2)
        merge_1 = Concatenate()([flat_1, flat_2])
        #Output of the model
        fc1 = Dense(n_nodes, activation = "relu")(merge_1)
        out = Dense(1, activation = "sigmoid")(fc1)
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
            print(pcc, roc_auc, max_acc)
            if n_epoch == 0 and not(pcc > 0.60):
                poor_init = True
                break
            foutput(str(n_epoch)+"\t"+str(pcc)+"\t"+str(roc_auc)+"\t"+str(max_acc), output_path)
        if poor_init:
            continue
        #Only select the models with good performances on the validation set (note that this is not cross-validation)
        if pcc > 0.80:
            foutput("Model Adopted", output_path)
            models.append(model)      
    return models