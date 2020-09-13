# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 12:30:47 2019

@author: meine
"""
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout
from tensorflow.keras import optimizers
from tensorflow.keras import backend as K
from tensorflow import compat
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

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
def custom_weighted_loss(dI):
    if not (K.int_shape(dI)[0] == 'None'):
        bf = K.ones((1, 1970))
        w = bf - dI
    def custom_loss(y_true, y_pred):
        return w * K.binary_crossentropy(y_true, y_pred)
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
history = <..> of training process (e.g. binary accuracy and loss over the different epochs )
"""
def create_neural_network(data, labels, val_data, val_labels, nlayers=3, nnodes=32, nepochs=10, b_size=32, dropout=0, save=True, segment=42):
    model = Sequential()
    model.add(Dense(nreactions, input_shape=(nreactions,), activation='relu'))
    for _ in range(nlayers):
        model.add(Dense(nnodes, activation='relu'))
        model.add(Dropout(dropout))
    model.add(Dense(nreactions, activation='sigmoid'))

    model.compile(optimizers.Adam(lr=0.005, beta_1=0.9, beta_2=0.999, epsilon=1e-08, decay=0.01, amsgrad=True),
                  loss=custom_weighted_loss(model.input),
                  metrics=['binary_accuracy'])



    history = model.fit(ndata, train_labels, epochs = nepochs, shuffle=True, batch_size = b_size, validation_data=(val_data, val_labels), verbose=2)
#    if(save):
#        model.save(path+"output/NN/networks/network%i.h5"%segment)
    return (model, history)


#splits list l in n equal parts
def split_list(l, n):
    a = np.asarray(l)
    sa = np.array_split(a, n)
    o = [list(i) for i in sa]

    return list(o)

#flattens list
def flatten_list(l):
    return [i for j in l for i in j]

#creates a list contain all elements of l except i
def exclude(l, i):
    return l[:i] + l[i+1:]

#functions to calculate f1-score
def calc_cm(t, p, m):
    ti = 1-t
    pi = 1-t
    TP = sum(p[t*m==1])
    FP = sum(p[ti*m==1])
    FN = sum(pi[t*m==1])
    TN = sum(pi[ti*m==1])

    return (TP, FP, FN, TN)

def calc_precision(cm):
    return cm[0]/(cm[0]+cm[1])

def calc_sensitivity(cm):
    return cm[0]/(cm[0]+cm[2])

def calc_f1score(t, p, m):
    cm = calc_cm(t, p, m)
    p = calc_precision(cm)
    s = calc_sensitivity(cm)
    return 2.*(p*s)/(p+s)

#plots accuracy over epochs for train and validation data based on network history
def plot_network_accuracy(h, seg):
    plt.plot(h.history['binary_accuracy'])
    plt.plot(h.history['val_binary_accuracy'])
    plt.title('model accuracy segment %i'%seg)
    plt.ylabel('binary accuracy')
    plt.xlabel('epoch')
    plt.legend(['train', 'test'], loc='upper left')
#    plt.savefig(path+'output/NN/figures/model_accuracy_s%i_d%.2f.png'%(seg, dropout))
    plt.show()

#plots loss over epochs for train and validation data based on network history
def plot_network_loss(h, seg):
    plt.plot(h.history['loss'])
    plt.plot(h.history['val_loss'])
    plt.title('model loss segment %i'%seg)
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'test'], loc='upper left')
#    plt.savefig(path+'output/NN/figures/model_loss_s%i_d%.2f.png'%(seg, dropout))
    plt.show()


#PARAMETERS

#Training parameters
nuplo = 5                    #number of copies of train data
nsegments = 10
max_con_train = 0.0          #noise in the train data 0->1 fr
max_for_train = 0.1          #noise in the train data 1->0
max_con_test  = 0.0          #noise in the test data 0->1
max_for_test  = 0.1         #noise in the test data 1->0

#NN parameters

nlayers =  8                #number of layers in the model
nnodes = 512                  #number of nodes in every layer
nepochs = 10                    #number of epochs (training) (needs to be an integer)
b_size  = 100                    #batch_size
dropout = 0.01                     #dropout

compat.v1.disable_eager_execution()
path = os.getcwd()
split_path = path.split('/')
data_path = '/'+os.path.join(*split_path[:-1])+'/files/'


#<Curate> and shuffle dataset
metadata = pd.read_csv(data_path+'new_metadata.csv', index_col=0 ) #Full dataset
metadata = metadata[metadata.columns[metadata.sum()>100]]      #Drop models < 100 reactions
genus_ids = list(metadata.columns)  #list of genus ids
np.random.shuffle(genus_ids)        #shuffle the ids
metadata = metadata[genus_ids]      #shuffle dataset


#Split of validation data
ncolumns = len(metadata.columns)
val_data = metadata.iloc[:, int(0.9*ncolumns):]
data = metadata.iloc[:, :int(0.9*ncolumns)]

#Load data into numpy array
reaction_ids = list(data.index)
matrix = np.asarray(data).T
ngenus, nreactions = matrix.shape

#Split dataset into multiple segments
data_sets = np.array_split(matrix, nsegments)
genus_lists = split_list(genus_ids, nsegments)

#Train and evaluate

for seg in range(nsegments):
    print("seqment: %i"%(seg))
    train_data = np.concatenate(exclude(data_sets, seg))
    test_data = data_sets[seg]
    train_ids = flatten_list(exclude(genus_lists, seg))
    test_ids = genus_lists[seg]

    ndata = np.zeros((len(train_data)*nuplo, nreactions))
    for i in range(len(train_data)):
        for j in range(nuplo):
            con_train = np.random.uniform(0, max_con_train)
            for_train = np.random.uniform(0, max_for_train)

            ndata[nuplo*i+j] = noise_data(train_data[i],con_train, for_train)

    ntest =  np.zeros((len(test_data), nreactions))
    for i in range(len(test_data)):
        con_test = np.random.uniform(0, max_con_test)
        for_test = np.random.uniform(0, max_for_test)

        ntest[i]= noise_data(test_data[i],con_test, for_test)

    train_labels = np.repeat(np.copy(train_data), nuplo, axis = 0)
    test_labels = np.copy(test_data)
    network, history = create_neural_network(ndata, train_labels, ntest, test_labels, nlayers, nnodes, nepochs, b_size, dropout,False, seg)
#    plot_network_accuracy(history, seg)
#    plot_network_loss(history, seg)


    m = 1-ntest
    p = network.predict(test_data)
    bp = np.round(p)

    m2 = np.ones(test_data.shape)

    f1 = calc_f1score(test_labels, bp, m)
    f1_2 = calc_f1score(test_labels, bp, m2)

    print('f1score masked: %.2f| f1-score complete: %.2f'%(f1, f1_2))
