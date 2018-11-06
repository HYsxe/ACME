import numpy as np

def scoring(models, data):
    '''
    Use an ensemble of models to yield prediction scores for the data
    Args:
        1. models: Model ensemble
        2. data: Data to be predicted
    Return values:
        1. probas_: Prediction scores
    '''
    probas_ = [np.transpose(model.predict(data))[0] for model in models]
    probas_ = [np.mean(scores) for scores in zip(*probas_)]
    return probas_    
