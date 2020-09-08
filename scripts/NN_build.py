# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 12:30:47 2019

@author: meine
"""
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras import optimizers
from tensorflow.keras import backend as K
import numpy as np
import pandas as pd
import os
#FUNCTIONS

"""
Function to to generate noise in input data,
noise_0 (Contamination) is changing 0 -> 1,
noise_1 (False Ommision) is changing 1 -> 0
"""
def noise_data(i, noise_0, noise_1):
    temp = i.copy()
    a=np.arange(len(temp))[temp==0]
    b=np.arange(len(temp))[temp!=0]
    n0 = int(len(a)*noise_0)
    n1 = int(len(b)*noise_1)
    np.random.shuffle(a)
    np.random.shuffle(b)
    s0 = a[0:n0]
    s1 = b[0:n1]
    temp[s0]= 1
    temp[s1]= 0
    o = temp
    return o

"""
Function to calculate weighted loss dI is the model.input
"""
def custom_weighted_loss(dI, bias):
    if not (K.int_shape(dI)[0] == 'None'):
        bf = K.ones((1, 1970))
        w = bf - dI
    def custom_loss(y_true, y_pred):
        loss = w * K.binary_crossentropy(y_true, y_pred)
        return bias*(1-y_true)*loss+(1-bias)*y_true*loss
    return custom_loss
"""
Most important function, creates actual NN,
INPUT:
data       = training input
labels     = training labels (truth)
val_data   = validation input
val_labels = validation labels (truth)
nlayers    = number of hidden layers (layers that are not input or output).
nnodes     = number of nodes per layer,
nepochs    = how often the NN needs to loop over all the data
b_size     = batch_size (number of training examples that are simultaneously evaluated)
dropout    = ''
save       = ?save model (output is path+NN/.)
segment    = which segment of data is currently evaluated (used for saving) default = 42


OUTPUT:
model   = Neural network
history = History of training process (e.g. binary accuracy and loss over the different epochs )
"""
def create_neural_network(data, labels, nlayers=1, nnodes=512, nepochs=7, b_size=32, dropout=0.1, bias_0=0.3, save=False, name='noname'):
    model = Sequential()
    model.add(Dense(nreactions, input_shape=(nreactions,), activation='relu'))
    for _ in range(nlayers):
        model.add(Dense(nnodes, activation='relu'))
        model.add(Dropout(dropout))
    model.add(Dense(nreactions, activation='sigmoid'))

    model.compile(optimizers.Adam(lr=0.005, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.01, amsgrad=True),
                  loss=custom_weighted_loss(model.input, bias_0),
                  metrics=['binary_accuracy'])

    if(save):
        model.save(output_path+"%s.h5"%name)
    return (model)


#PARAMETERS

#Training data parameters
nuplo = 15                    #number of copies of train data
max_con_train = 0.0          #noise in the train data 0->1 fr
max_for_train = 0.3          #noise in the train data 1->0


#NN parameters

nlayers =  1                #number of layers in the model
nnodes = 512                #number of nodes in every layer
nepochs = 7                #number of epochs (training) (needs to be an integer)
b_size  = 50               #batch_size
dropout = 0.1              #dropout
bias_0 = 0.3            #bias towards negative class


path = os.getcwd()
split_path = path.split('/')
data_path = '/'+os.path.join(*split_path[:-1])+'/files/'
output_path = ''


#Load and shuffle dataset
metadata = pd.read_csv(data_path+'new_metadata.csv', index_col=0 )  #Full dataset
metadata = metadata[metadata.columns[metadata.sum()>100]]           #Drop models < 100 reactions
genus_ids = list(metadata.columns)                                  #list of genus ids
np.random.shuffle(genus_ids)                                        #shuffle the ids
metadata = metadata[genus_ids]                                      #shuffle dataset


#Load data into numpy array
reaction_ids = list(data.index)
matrix = np.asarray(data).T
ngenus, nreactions = matrix.shape

#Train
train_data =  matrix
train_ids = genus_ids
test_ids = genus_lists[seg]

ndata = np.zeros((len(train_data)*nuplo, nreactions))
for i in range(len(train_data)):
    for j in range(nuplo):
        con_train = np.random.uniform(0, max_con_train)
        for_train = np.random.uniform(0, max_for_train)

        ndata[nuplo*i+j] = noise_data(train_data[i],con_train, for_train)

train_labels = np.repeat(np.copy(train_data), nuplo, axis = 0)
network = create_neural_network(ndata, train_labels, nlayers, nnodes, nepochs, b_size, dropout, bias_0, True, 'NN_full')
type(network)
