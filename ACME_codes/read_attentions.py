#Here we adopted the code from https://github.com/philipperemy/keras-visualize-activations/blob/master/read_activations.py

import keras.backend as K
import numpy as np
import matplotlib.pyplot as plt
path_img = r"/home/huyan/peptide_splicing/workstation/attention/output_attention/"

def get_attentions(model, model_inputs, print_shape_only=False, layer_name=None):
    activations = []
    inp = model.input

    model_multi_inputs_cond = True
    if not isinstance(inp, list):
        # only one input! let's wrap it in a list.
        inp = [inp]
        model_multi_inputs_cond = False

    outputs = [layer.output for layer in model.layers if
               layer.name == layer_name or layer_name is None]  # all layer outputs

    funcs = [K.function(inp + [K.learning_phase()], [out]) for out in outputs]  # evaluation functions

    if model_multi_inputs_cond:
        list_inputs = []
        list_inputs.extend(model_inputs)
        list_inputs.append(0.)
    else:
        list_inputs = [model_inputs, 0.]

    # Learning phase. 0 = Test mode (no dropout or batch normalization)
    # layer_outputs = [func([model_inputs, 0.])[0] for func in funcs]
    layer_outputs = [func(list_inputs)[0] for func in funcs]
    attentions = []
    layer_activations = layer_outputs[27]
    return layer_activations

