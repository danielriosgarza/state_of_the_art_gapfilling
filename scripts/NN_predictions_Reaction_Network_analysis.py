#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 14:46:44 2019

@author: meine
"""
import numpy as np
import pandas as pd
import os
import networkx as nx
from Reaction_class_08_27 import Reaction, read_biochemistry_table
import itertools as it
from sklearn.metrics import confusion_matrix as cm

from scipy.stats import pearsonr

"""
Functions
"""

#Auxilliary function to normalise a dictionary

def norm_dic(d):
    o = {}
    vmin = min(d.values())
    vmax = max(d.values())
    for k, v in d.iteritems():
        o[k] = (v - vmin)/(vmax - vmin)
    return o

#Function to determine if set a and set b have at least one common member

def common_member(a, b): 
    a_set = set(a) 
    b_set = set(b) 
    if (a_set & b_set): # Returns the set of elements that are in both a and b, if empty no common members
        return True 
    else: 
        return False



path = os.getcwd()
split_path = path.split('/')
data_path = '/'+os.path.join(*split_path[:-1])+'/files/'
#output_path = "" OUTPUT PATH DISABLED

df_truth = np.load(data_path+"dataset_genera.npy") #load truth

df_predictions =np.load(data_path+"prediction_1l_2020-02-17_30_0.npy") #load prediction

genus_ids = np.load(data_path+"id_vector.npy")  #get list of genera
reaction_ids = np.load(data_path+"rxn_vector.npy") #get list of reactions
t = df_truth.T          #transfer truth df to numpy array
p = df_predictions.T #idem for predictions, with same order of genera
nrows, ncolumns = t.shape  #get number of genera and reactions

bp = np.round(p, 0)       #create binarized predictions



dic_TP = {}
dic_TN = {}
dic_FN = {}
dic_FP = {}


df_input = pd.read_csv(data_path+'incomplete_30_0.csv', index_col=0) #get df used as input
m = 1-np.asarray(df_input).T                                 #get mask (all 0s in input)


#Calculate confusion matrix with m as weights

for p,q,k,w in zip(bp, t, reaction_ids, m):
    c = cm(p, q, sample_weight = w)
    dic_TP[k] = c[1,1]
    dic_TN[k] = c[0,0]
    dic_FN[k] = c[0,1]
    dic_FP[k] = c[1,0]
    
    
#Calculate Recall precision and f1

score_df  = pd.DataFrame(data={'TP':dic_TP,'FP':dic_FP, 'TN':dic_TN, 'FN':dic_FN})

score_df['sensitivity'] = score_df['TP']/(score_df['TP']+score_df['FN'])

score_df['precision'] = score_df['TP']/(score_df['TP']+score_df['FP'])

score_df.plot.scatter(x='sensitivity', y='precision')

score_df['f1score'] = 2.*(score_df['precision']*score_df['sensitivity'])/(score_df['precision']+score_df['sensitivity'])


#get modelSEED reaction database
reaction_dict = read_biochemistry_table('/home/meine/Old_Scripts/reactions.tsv')
reaction_db = Reaction(biochem_input=reaction_dict)

#Define common metabolites
common_metabolites = ['cpd00067_c', 'cpd00067_e']
for i in range(22):
    if i < 9:
        common_metabolites.append('cpd0000%i_c'%(i+1))
        common_metabolites.append('cpd0000%i_e'%(i+1))
    else:
        common_metabolites.append('cpd000%i_c'%(i+1))
        common_metabolites.append('cpd000%i_e'%(i+1))
        
    

#create dictionary of reactions and metabolites
dict_metabolites = {}
for r  in reaction_db.reactions:
    me = reaction_db.reactions[r]['metabolites'].keys()
    for met in common_metabolites:
        if met in me: me.remove(met)
    dict_metabolites[r] = me

#count the frequency of metabolites
count_metabolites = {}
for r in reaction_db.reactions:
    mes = reaction_db.reactions[r]['metabolites'].keys()
    for me in mes:
        if me not in count_metabolites.keys():
            count_metabolites[me] = 1
        else:
            count_metabolites[me] += 1




#Create graph with metabolites
G = nx.DiGraph()
for reaction in reaction_ids:
    r = reaction_db.reactions[reaction]    
    G.add_node(reaction,nt='reaction')
    
    for metabolite in dict_metabolites[reaction]:
        if not G.has_node(metabolite):
            G.add_node(metabolite,nt='metabolite')
        
        if r['metabolites'][metabolite]<0:
            G.add_edge(metabolite, reaction)
        else:
            G.add_edge(reaction, metabolite)

#nx.write_gml(G,Output_path + "network_with_metabolites.gml') #Ouput disabled

#connect reactions directly
reaction_network = []
for i,j in it.combinations(reaction_ids,2):
    a = dict_metabolites[i]
    b = dict_metabolites[j]
    if common_member(a, b):
        reaction_network.append([i, j])
          

#Create graph without metabolites            
G2 = nx.Graph()
for i in reaction_network:
    a = i[0]
    b = i[1]
    if not G2.has_node(a):
        G2.add_node(a, f1=score_df['f1score'])
    if not G2.has_node(b):
        G2.add_node(b, f1=score_df['f1score'])
    G2.add_edge(a, b)


#nx.write_gml(G2,output_path+'network_without_metabolites.gml') #Output disabled


#Define charcteristics
degree = dict(nx.degree(G2))
clustering = dict(nx.clustering(G2))
centrality = dict(nx.degree_centrality(G2))


#plot characteristics
network_df = score_df.loc[list(nx.nodes(G2)), ['precision','sensitivity','f1score']]
network_df['degree'] = degree.values()
network_df['betweenness_centrality'] = centrality.values()
network_df['clustering'] = clustering.values()
network_df.plot.scatter(x='degree', y='f1score')
network_df.plot.scatter(x='clustering', y='f1score')
network_df.plot.scatter(x='betweenness_centrality', y='f1score', xlim=(0, 0.07))

temp = network_df.dropna()
pearsonr(temp['betweenness_centrality'], temp['f1score'])
pearsonr(temp['clustering'], temp['f1score'])


