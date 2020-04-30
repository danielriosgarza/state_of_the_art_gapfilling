#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 12:19:41 2019

@author: meine
"""

from datetime import date
import gapfill_function as gap
import pandas as pd
import numpy as np
from Reaction_class_08_27 import Reaction
import os

#function to create a vector where every id has the same cost
def create_costs_fixed(ids, cost):
    temp_costs = {}
    for i in ids:
        temp_costs[i] = cost
    return temp_costs

#Set path names
path = os.getcwd()
split_path = path.split('/')
data_path = '/'+os.path.join(*split_path[:-1])+'/files/'
reaction_db = data_path+'modelseed_reactions.tsv'
model_path = '/'+os.path.join(*split_path[:-1])+'/test_model/'
#output_path = path+'/output'


#Get predictions NN
#predictions = np.load(data_path+"/prediction_1l_2020-02-17_30_0.npy")
#model_ids = np.load("id_vector.npy")
nn_reaction_ids = np.load("rxn_vector.npy")

#PLACEHOLDER predictions
predictions = np.random.uniform(len(nn_reaction_ids))
model_ids = ['ungapfilled_model']

#Read reaction database used for gapfilling
biochem = gap.read_biochemistry_table(reaction_db)
db_reactions = Reaction(biochem_input=biochem) #Load all reactions that will be considered for gap filling.

def_cost = 50 #Set a cost that each reaction should get when it is not present in your vector.
def_nn_cost = -np.log(0.5) #Set a standard cost for NN_set reactions (used for latendresse split)

p_linear = 1 - predictions #Linear translation prob -> cost
p_binary = p_linear.apply(np.round) #binarisation of p_linear costs
p_log = (predictions+10**-12).apply(lambda x:-np.log(x)) #log translation prob -> cost

p_binary_2 = p_binary.replace(1, def_cost) #replaces all 1 costs in NN_binary with the high def cost

#Dictionary of NN_method cost vectors
NN_methods_dic = {'NN_linear':p_linear, 'NN_log':p_log, 'NN_binary':p_binary, 'NN_binary_2':p_binary_2}



#default latendresse
costs_0 = {}
costs_1 = create_costs_fixed(nn_reaction_ids, cost = def_nn_cost)

#methods used for gapfilling
methods = ['latendresse','latendresse_split', 'NN_linear', 'NN_log', 'NN_binary', 'NN_binary_2']

#optional shuffling of model_ids in case of sampling
#np.random.shuffle(model_ids)

#output    
gapfill_result = []

#loop for the actual gapfilling
for i, A in enumerate(model_ids):
    A_complete = Reaction(model=model_path+'complete/'+A+'.sbml')
    A_incomplete = Reaction(model = model_path+'incomplete/0.2deleted/repl0/'+A+'.sbml')
    A_reaction_ids = set(A_complete.reactions.keys())
    A_incomplete_reaction_ids = set(A_incomplete.reactions.keys())
    deleted_reactions = A_reaction_ids.difference(A_incomplete_reaction_ids)
    all_reactions = A_complete
    all_reactions.reactions = all_reactions.add_dict(all_reactions.reactions, db_reactions.reactions)
    for m in methods:
        print('Gapfilling model %i: id =  %s, method = %s'%(i, A, m))
        if m == 'latendresse':
            candidate_reactions = costs_0
        elif m == 'latendresse_split':
            candidate_reactions = costs_1
        else:
            candidate_reactions = NN_methods_dic[m][A].to_dict()
        result = gap.gapfill(all_input_reactions=all_reactions, compare_to_added_reactions=deleted_reactions, input_reaction_ids=A_incomplete_reaction_ids, input_candidate_reactions=candidate_reactions, biomass_name='bio1', medium = 'complete', default_cost = def_cost)
        result.append(A)
        result.append(m)
        gapfill_result.append(result)

gapfill_result_df = pd.DataFrame(data=gapfill_result, columns=['delta','TP','FP','TN','FN','reactions','id','method'])
#gapfill_result_df.to_csv(output_path + today + "_20_1l_b30.csv")


