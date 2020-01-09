#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 16:52:43 2019

@author: jan
"""
import urllib
import matplotlib.pyplot as plt
import os

'''
Set directory paths.
Make sure to run this script from git/state_of_the_art_gapfilling/scripts.
'''
current_dir = os.getcwd()
split_dir = current_dir.split('/')
files_folder = '/'+os.path.join(*split_dir[:-1])+'/files/'
data_folder = '/'+os.path.join(*split_dir[:-2])+'/data/'

#Map strain ids to strain names
id_to_name = {}
best_per_genus_tsv = file(files_folder+'best_per_genus.20190620.tsv', 'r')
for line in best_per_genus_tsv:
    split_line = line.strip().split('\t')
    strain_id = split_line[0]
    strain_name = split_line[8]
    id_to_name[strain_id] = strain_name

def get_taxonomy(tax_id):
    '''
    For a given NCBI taxonomy id, return the a list with [phylum, class, order, family].
    Using urllib and the ncbi taxonomy browser.
    '''
    page = urllib.urlopen('https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id='+tax_id)
    lines = page.readlines()
    tax = [0,0,0,0,0]
    for line in lines:
        if 'superkingdom' in line:
            split_line = line.split(';')
            for l in split_line:
                if 'TITLE' in l and 'superkingdom' in l or 'ALT' in l and 'superkingdom' in l:
                    tax[0] = l.split('>')[1].replace('<', '').replace('/a','').replace('/A', '')
                elif 'TITLE' in l and 'phylum' in l or 'ALT' in l and 'phylum' in l:
                    tax[1] = l.split('>')[1].replace('<', '').replace('/a','').replace('/A', '')
                elif 'TITLE' in l and 'class' in l or 'ALT' in l and 'class' in l:
                    tax[2] = l.split('>')[1].replace('<', '').replace('/a','').replace('/A', '')
                elif 'TITLE' in l and 'order' in l or 'ALT' in l and 'order' in l:
                    tax[3] = l.split('>')[1].replace('<', '').replace('/a','').replace('/A', '')
                elif 'TITLE' in l and 'family' in l or 'ALT' in l and 'family' in l:
                    tax[4] = l.split('>')[1].replace('<', '').replace('/a','').replace('/A', '')
            return tax
        
id_list = id_to_name.keys()

ncbi_id_list = [id_.split('.')[0] for id_ in id_list]

taxonomy = {}
for id_ in ncbi_id_list:
    taxonomy[id_] = get_taxonomy(id_)
    if len(taxonomy.keys())%100==0:
        print len(taxonomy.keys())
        
#Rerun taxonomy with new function version for all previously unfound taxonomies.
unfound_tax = [id_ for id_ in taxonomy.keys() if taxonomy[id_] == [0,0,0,0,0]]
for id_ in unfound_tax:
    taxonomy[id_] = get_taxonomy(id_)

#For each unique phylum, class, store how often it is found in the tree.
phyla = {}
for taxo in taxonomy.values():
    if taxo[1] not in phyla:
        phyla[taxo[1]] = 1
    else:
        phyla[taxo[1]] +=1
        
classes = {}
for tax in taxonomy.values():
    if tax[2] not in classes:
        classes[tax[2]] = 1
    else:
        classes[tax[2]] +=1
        
plt.hist(phyla.values())
colors = ['darkred', 'darkorange', 'forestgreen', 'royalblue', 'indigo', 'cyan', 'pink']
phyla_colors = {}
for phylum in phyla:
    if phyla[phylum] > 49 and phylum != 0:
        phyla_colors[phylum] = colors[0]
        colors = colors[1:]
    elif phylum == 0:
        phyla_colors[phylum] = 'linen'
    else:
        phyla_colors[phylum] = 'gray'
#Write phyla_colors to a file
file_name = files_folder+'phyla_colors.txt'
f = file(file_name, 'w')
f.write(str(phyla_colors))
f.close()

#Write taxonomy to a file
file_name = files_folder+'taxonomy.txt'
f = file(file_name, 'w')
f.write(str(taxonomy))
f.close()
