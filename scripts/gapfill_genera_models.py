#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 16:30:38 2019

@author: jan
"""
import gc
import numpy as np
from Reaction_class_08_27 import Reaction
import gapfill_function
from copy import deepcopy
import pandas as pd
from datetime import date
import os

'''
Set directory paths.
Make sure to run this script from git/state_of_the_art_gapfilling/scripts.
'''
current_dir = os.getcwd()
split_dir = current_dir.split('/')
files_folder = '/'+os.path.join(*split_dir[:-1])+'/files/'
data_folder = '/'+os.path.join(*split_dir[:-2])+'/data/'

#Define functions
def read_list(location):
    f = file(location, 'r')
    for line in f:
        id_list = eval(line.strip())
    return id_list

def get_distances(tree, X, target_list):
    """
    Returns a list of distances from X to every node in target list.
    """
    distances = {}
    for node in target_list:
        if node != X:
            d = tree.get_distance(X, node)
            distances[node] = d
        else:
            distances[node] = 0
            
    return distances

def get_cost(distance_to_X, reaction_presence):
    """
    For a given reaction, based on its presence in other models in the tree, returns a cost for that reaction.
    The cost is higher when the reaction is less present in the tree (more abesent), 
    and every node in the tree is scaled inversely to their distance.
    """
    inverse_distance_to_X = {k : 1/v for k,v in distance_to_X.items()}
    sum_absent = 0
    for id_ in distance_to_X:
        if reaction_presence[id_] == 0:
            sum_absent += inverse_distance_to_X[id_]
            
    sum_all = sum(inverse_distance_to_X.values())
    cost = sum_absent/float(sum_all)
    return cost

def cost_function(model, id_array, rxn_vector, dist_matrix, rxn_matrix, all_model_reactions):
    '''
    Function to get candidate reactions using the cost function.
    '''
    candidate_reactions = {}
    model_index = np.where(id_array==model)[0][0]
    distances = dist_matrix[model_index]
    dist_dict = {id_array[i]:distances[i] for i in range(len(distances))}
    del dist_dict[model] #Delete self from dict
    for reaction in all_model_reactions.reactions.keys():
        rxn_index = np.where(rxn_vector==reaction)[0][0]
        presences = rxn_matrix.T[rxn_index]
        pres_dict = {id_array[i]:presences[i] for i in range(len(presences))}
        del pres_dict[model] #Delete self from dict
        candidate_reactions[reaction] = get_cost(dist_dict, pres_dict)
    return candidate_reactions

def func(x, a, b, c):
    return a * np.exp(-b * x) + c

def emp_likelihood_func(x, parameters = [1.28467195,  1.91736455, -0.00487022], maxi=1.27980173):
    y = func(x, *parameters)
    y_squished = y/maxi #Divide by max possible y value (x=0) such that the cost will never be negative.
    cost = 1. - y_squished
    return cost

def min_dist(model, id_array, rxn_vector, min_dist_matrix):
    '''
    Function to get candidate reactions where the cost is taken from empirical cost function.
    '''
    model_index = np.where(id_array==model)[0][0]
    min_dists = min_dist_matrix[model_index]
    min_dist_dict = {rxn_vector[i]:min_dists[i] for i in range(len(rxn_vector))}
    candidate_reactions = {key:emp_likelihood_func(value) for key,value in min_dist_dict.items()}
    return candidate_reactions

id_list = read_list(data_folder+'genera_models/complete/ids_in_tree.txt')
model_path = data_folder+'genera_models/complete/'
biochem = gapfill_function.read_biochemistry_table(files_folder+'modelseed_reactions.tsv')

#Load reactions
db_reactions = Reaction(biochem_input=biochem)
all_model_reactions = Reaction(model_folder=model_path, model_list=id_list)
#Only keep rxn, delete bio and EX reactions from all_model_reactions
for reaction in all_model_reactions.reactions.keys():
    if 'rxn' not in reaction:
        del all_model_reactions.reactions[reaction]


#Import reaction presence absence matrix and distance matrix
id_array = np.load(files_folder+'id_vector.npy')
rxn_vector = np.load(files_folder+'rxn_vector.npy')
rxn_matrix = np.load(files_folder+'rxn_presence_matrix_complete.npy')
dist_matrix = np.load(files_folder+'dist_matrix.npy')
min_dist_matrix = np.load(files_folder+'min_dist_matrix.npy')


models_to_gapfill = [id_ for id_ in id_list]
models_to_gapfill = models_to_gapfill[1830:] #Skip malformed XML file 1829 of 0.2 deleted repl 0.
#Setup gap fill methods, which models to gap fill
#method = ['control', 'min_dist', 'cost']
method = ['flat_control']
deleted = [0.3]
r = 0 #Replicate
cost_db_rxn = [40]
Results = []

for m in method:
    for d in deleted:
        for model in models_to_gapfill:
            for c in cost_db_rxn:
                if d == 0.0:
                    model_map = data_folder+'genera_models/complete/'
                elif d == 0.1:
                    model_map = data_folder+'genera_models/incomplete/0.1deleted/repl' + str(r) + '/'
                elif d == 0.2:
                    model_map = data_folder+'genera_models/incomplete/0.2deleted/repl' + str(r) + '/'
                elif d == 0.3:
                    model_map = data_folder+'genera_models/incomplete/0.3deleted/repl' + str(r) + '/'
                
                Model_reactions = Reaction(model =  model_map + model + '.sbml')
                Model_reaction_ids = set(Model_reactions.reactions.keys())
                all_reactions = deepcopy(Model_reactions)
                all_reactions.reactions = all_reactions.add_dict(all_reactions.reactions, all_model_reactions.reactions)
                all_reactions.reactions = all_reactions.add_dict(all_reactions.reactions, db_reactions.reactions)
                #Method for getting candidate_reactions
                if m == 'control':
                    candidate_reactions = {}
                    c = 1 #All reactions will get a default cost of 1.
                elif m == 'min_dist':
                    candidate_reactions = min_dist(model, id_array, rxn_vector, min_dist_matrix)
                elif m == 'cost':
                    candidate_reactions = cost_function(model, id_array, rxn_vector, dist_matrix, rxn_matrix, all_model_reactions)
                elif m == 'random_control':
                    candidate_reactions = min_dist(model, id_array, rxn_vector, min_dist_matrix)
                    costs = candidate_reactions.values()
                    np.random.shuffle(costs)
                    reactions = candidate_reactions.keys()
                    candidate_reactions = {reactions[i]:costs[i] for i in range(len(costs))} #Randomly shuffled costs
                elif m == 'flat_control':
                    candidate_reactions = {i: 1.0 for i in all_model_reactions.reactions}
                    
                complete_reactions = Reaction(model = model_path + model + '.sbml')
                complete_reaction_ids = set(complete_reactions.reactions.keys())
                deleted_reaction_ids = complete_reaction_ids.difference(Model_reaction_ids)
                Result = gapfill_function.gapfill(all_reactions, Model_reaction_ids, candidate_reactions, 'bio1', compare_to_added_reactions = deleted_reaction_ids, default_cost=c)
                Result.append(m)
                Result.append(d)
                Result.append(model)
                Result.append(c)
                Results.append(Result)
                gc.collect()
    
            today = str(date.today())
            Results_df = pd.DataFrame(np.array(Results, dtype=object),    
                                      columns = ['delta', 'TP', 'FP', 'TN', 'FN', 'rxns_added', 'method', 'deleted', 'model', 'cost_db_rxn'])
            Results_df.to_csv(files_folder + today + 'gapfill_result.csv')

