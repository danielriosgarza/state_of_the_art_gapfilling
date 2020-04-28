#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 13:19:13 2019

@author: meine
"""
import numpy as np
import pandas as pd
import os
import tensorflow as tf
from datetime import date

today = date.today()



#splits list l in n equal parts
def split_list(l, n):
    a = np.asarray(l)
    sa = np.array_split(a, n)
    o = [list(i) for i in sa]
        
    return list(o)


path = os.getcwd()
split_path = path.split('/')
data_path = '/'+os.path.join(*split_path[:-1])+'/files/'


reaction_ids = np.load(data_path+'rxn_vector.npy')
genus_ids =  np.load(data_path+'id_vecotor.npy')
df_input = pd.read_csv(data_path+'incomplete_30_0.csv', index_col=0)

df_input = df_input[genus_ids]

genus_lists = split_list(genus_ids, 5)


prediction = np.zeros(df_input.shape)

for seg in range(5):
#   load network
    network = tf.keras.models.load_model(path+'output/NN/model%i_1_layers_b30.h5'%(seg), custom_objects={"custom_loss": 'binary_crossentropy'})
    input_ids = genus_lists[seg]
    input_data = np.asarray(df_input[input_ids]).T
    t_result = network.predict(input_data)        
    prediction[genus_lists[seg]] = t_result.T
prediction.save(data_path+"prediction_"+today+"_1l_30_0.csv")


