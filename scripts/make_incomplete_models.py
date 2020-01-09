#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 10:19:24 2019

@author: jan
For all genus models, make 10 incomplete models, for 10, 20 and 30 percent deleted genes.
Save all these models in directories.
"""
import numpy as np
import cobra
import os

'''
Set directory paths.
Make sure to run this script from git/state_of_the_art_gapfilling/scripts.
'''
current_dir = os.getcwd()
split_dir = current_dir.split('/')
files_folder = '/'+os.path.join(*split_dir[:-1])+'/files/'
data_folder = '/'+os.path.join(*split_dir[:-2])+'/data/'

def delete_fraction_of_genes(model, fraction):
    """
    For a given model, delete given fraction of genes, and all reactions that only have deleted genes associated.
    Returns the incomplete model and lists of deleted genes, deleted reactions.
    """
    incomplete_model = model.copy()
    incomplete_model_genes = list(incomplete_model.genes) #Get a list of genes present in the model.
    np.random.shuffle(incomplete_model_genes) #Shuffle the list to ensure genes are deleted randomly.
    index = int(len(incomplete_model_genes)*(fraction)) #Get the index that matches the fraction of genes to delete.
    genes_to_delete = incomplete_model_genes[0:index]
    reactions_to_delete = [] #List all the reactions that will be deleted from the model.
    for reaction in incomplete_model.reactions:
        if reaction.genes <= set(genes_to_delete) and len(reaction.genes) != 0: #If all genes from a reaction are in the genes to delete list, select that reaction.
            reactions_to_delete.append(reaction.id) #Reactions that do not have genes will not be selected.
            reaction.delete() #delete reaction
            
    incomplete_model.repair()
    return incomplete_model, genes_to_delete, reactions_to_delete

#Read id_list
f = file('/home/jan/genera_models/ids_in_tree.txt', 'r')
id_list = eval(f.read())
f.close()

model_folder = data_folder+'genera_models/'
output_folder = data_folder+'genera_models/incomplete/'

fraction = [0.1, 0.2, 0.3]
repl = 10

#Create dirs
for f in fraction:
    if not os.path.isdir(output_folder + str(f) + 'deleted'):
        os.mkdir(output_folder + str(f) + 'deleted')
    for r in range(repl):
        if not os.path.isdir(output_folder + str(f) + 'deleted/repl' + str(r)):
            os.mkdir(output_folder + str(f) + 'deleted/repl' + str(r))

for strain in id_list: 
    co_model = cobra.io.read_sbml_model(model_folder + strain + '.sbml')
    for f in fraction:
        for r in range(repl):
            incomplete_co_model, g_to_del, r_to_del = delete_fraction_of_genes(co_model, f)
            cobra.io.write_sbml_model(incomplete_co_model, output_folder + str(f) + 'deleted/repl' + str(r) + '/' + strain + '.sbml')
