#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 13:19:13 2019

@author: meine
"""
import numpy as np
import pandas as pd
import tensorflow as tf
from datetime import date



#splits list l in n equal parts
def split_list(l, n):
    a = np.asarray(l)
    sa = np.array_split(a, n)
    o = [list(i) for i in sa]
        
    return list(o)


path = "/home/meine/Git/"
today = str(date.today())

truth = pd.read_csv(path+"data/dataset_genera.csv", index_col=0)

reaction_ids = list(truth.index)
genus_ids =  list(truth.columns)
df_input = pd.read_csv(path+'data/Incomplete_dfs/incomplete_10_0.csv', index_col=0)

df_input = df_input[genus_ids]

ncolumns, nrows = truth.shape

genus_lists = split_list(genus_ids, 5)


df_prediction = pd.DataFrame(index=reaction_ids, columns=genus_ids)

for seg in range(5):
    network = tf.keras.models.load_model(path+'output/NN/model%i_1_layers_b30.h5'%(seg), custom_objects={"custom_loss": 'binary_crossentropy'})
    input_ids = genus_lists[seg]
    seg_truth = np.asarray(truth[input_ids]).T
    input_data = np.asarray(df_input[input_ids]).T
    t_result = network.predict(input_data)        
    df_prediction[genus_lists[seg]] = t_result.T
df_prediction.to_csv(path+"output/NN/predictions/prediction_1l_"+today+"_10_0.csv")


