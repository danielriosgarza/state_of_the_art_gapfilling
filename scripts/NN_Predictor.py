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
single_input = False

#Load input
path = os.getcwd()
data_path = os.path.join(os.path.split(path)[0],'files')
input_path = os.path.join(data_path, 'example_binary.npy')
df_input = np.load(input_path)

#load reaction ids
reaction_ids = np.load(os.path.join(data_path,'rxn_vector.npy')).astype('str')

#test for single input (trips up NN)
if np.ndim(df_input) == 1:
    input_data = np.tile(df_input,[2,1])
else:
    input_data = df_input
prediction = np.zeros(input_data.shape)

#   load network and make prediction
network = tf.keras.models.load_model(os.path.join(data_path,'NN_full.h5'), custom_objects={"custom_loss": 'binary_crossentropy'})
t_result = network.predict(input_data)
prediction = np.asarray(t_result)

#ugly i know. but do not have a solution at this point
if np.ndim(df_input) == 1:
    prediction  = prediction[0].T
else:
    prediction = prediction.T
    
np.save(os.path.join(data_path,"prediction_example.npy"), prediction)
