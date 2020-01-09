#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 16:40:25 2019

@author: meine
"""

import ete3
import os
import numpy as np
import pandas as pd
from matplotlib import colors, cm
from datetime import date

import matplotlib.pyplot as plt



path = os.getcwd()
split_path = path.split('/')
data_path = '/'+os.path.join(*split_path[:-1])+'/files/'


df_truth = pd.read_csv(data_path+"new_metadata.csv", index_col=0)
p = np.load(data_path+"prediction_2019-12-10_30_0.npy").T
m = 1-np.asarray(pd.read_csv(data_path+"incomplete_30_0.csv")).T
t = np.asarray(df_truth).T

genus_ids = list(df_truth.columns)
reaction_ids = list(df_truth.index)

nrows, ncolumns = t.shape

bp = np.round(p, 0)
bp_i = 1-bp


dic_TP = {}
dic_TN = {}
dic_FN = {}
dic_FP = {}





for pi,ni,j,k,w in zip(bp, bp_i, t, genus_ids, m):
    dic_TP[k] = sum(pi[j*w==1])
    dic_TN[k] = sum(ni[(1-j)*w==1])
    dic_FN[k] = sum(ni[j*w==1])
    dic_FP[k] = sum(pi[(1-j)*w==1])


score_df  = pd.DataFrame(data={'TP':dic_TP,'FP':dic_FP, 'TN':dic_TN, 'FN':dic_FN})

score_df['nreactions'] = sum(t.T)

score_df.plot.scatter(x='FP', y='TP',c='nreactions', cmap='autumn')

score_df['sensitivity'] = score_df['TP']/(score_df['TP']+score_df['FN'])
score_df.plot.scatter(x='nreactions', y='sensitivity')
#score_df['sensitivity'].hist(bins=50)

score_df['specificity'] = score_df['TN']/(score_df['TN']+score_df['FP'])
score_df.plot.scatter(x='nreactions', y='specificity')
#score_df['specificity'].hist(bins=50)

score_df['precision'] = score_df['TP']/(score_df['TP']+score_df['FP'])
score_df.plot.scatter(x='nreactions', y='precision')
#score_df['precision'].hist(bins=50)

score_df.plot.scatter(x='sensitivity', y='precision', c='nreactions', cmap='cool')

score_df['f1score'] = 2.*(score_df['precision']*score_df['sensitivity'])/(score_df['precision']+score_df['sensitivity'])
score_df.plot.scatter(x='nreactions', y='f1score')
#score_df['f1score'].hist(bins=50)


contree = ete3.Tree(data_path+'concatenated_msa_trimmed.contree')
ts = ete3.TreeStyle()
ts.mode = 'c'

cmap = cm.cool

for n in contree.traverse():
   nstyle = ete3.NodeStyle()
   if n.name in score_df['f1score']:
       rgba = cmap(score_df['f1score'][n.name])
   else:
       rgba = (1, 1, 1, 1.0)
   nstyle["bgcolor"] = colors.to_hex(rgba)
   n.set_style(nstyle)

contree.show(tree_style=ts)   

data = score_df['f1score'].dropna()


# Plot histogram.
n, bins, patches = plt.hist(data, 100, color='green')
bin_centers = 0.5 * (bins[:-1] + bins[1:])

# scale values to interval [0,1]
col = bin_centers - min(bin_centers)
col /= max(col)

for c, p in zip(col, patches):
    plt.setp(p, 'facecolor', cmap(c))
plt.xlabel("f1score")
plt.ylabel("count")
#plt.savefig(path+'output/Trees/Hist_'+str(date.today())+'.svg')
plt.show()


bl_dic = {}
for n in contree.traverse():
    if n.name in score_df.index:
        bl_dic[n.name] = n.dist
    
branch_length = pd.Series(bl_dic)    
score_df['branch_length'] = branch_length
score_df.plot.scatter(x='branch_length', y='f1score')
