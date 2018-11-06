import scipy.stats as ss
import numpy as np
from scoring import *
from sklearn.metrics import roc_curve, auc
from math import log

def model_eval(model,test_data,test_target):
    '''
    Evaluate the performance of a model (or an ensemble of models)
    on an dataset
    Args:
        1. model: A keras model or a list of models 
        2,3. test_data, test_target: data and labels for testing
    Return values:
        1.PCC: Pearson correlation coefficient between the actual target 
                values and the predicted values
        2.roc_auc: AUROC of prediction
        3.acc: Accuracy of prediction
    '''  
    #Using the model(s) to make predictions
    if type(model) != type([]):
        #If it's not an ensemble
        scores= np.transpose(model.predict(test_data))[0]
    else:
        #Otherwise
        scores = scoring(model, test_data)

    #Compute the Pearson correlation coefficient
    pcc = ss.pearsonr(test_target,scores)    
    
    # Compute ROC curve and area the curve
    test_label = [1 if aff > 1 - log(500) / log(50000) else 0 for aff in test_target]
    fpr, tpr, thresholds = roc_curve(test_label,scores)
    roc_auc = auc(fpr, tpr)
    
    #Compute the accuracy
    threshold = 1 - log(500) / log(50000) 
    predictions = [0 if score < threshold else 1 for score in scores]
    accurate = [1 if predictions[i] == test_label[i] else 0 for i in range(len(predictions))]
    acc = np.sum(accurate)/float(len(accurate))  

    return pcc[0], roc_auc, acc