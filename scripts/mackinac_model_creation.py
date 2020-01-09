#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:25:22 2019

@author: jan
This script generates and stores sbml models based on an input list of patric genome ids.
"""
#Imports
import mackinac
import cobra

#Define functions
def read_list(location, skip_first_line=False):
    '''
    Used to read a list of all strain_ids, to load all models.
    Expected input: strain_id +\n, no header.
    '''
    f = file(location, 'r')
    templist = []
    
    if skip_first_line:
        f.readline()
        
    for line in f:
        templist.append(line.split(',')[0].strip('"')) #Works for both id_list and PATRIC_genome.csv files.
        
    for i in templist:
        if len(i) < 3: #To remove possible empty lines from the list
            templist.remove(i)
            
    f.close()            
    return templist

#Script
path = '/home/jan/sbml_models/'
input_list = read_list('/home/jan/Downloads/PATRIC_genome.csv', skip_first_line=True)
modelled_strain_ids = []

model_dict = {}
gapfilled_model_dict = {}
mackinac.get_token('JBaijens')
for patric_id in input_list[43:]: #CHANGE THIS IF RESTARTING
    try:
        mackinac.reconstruct_modelseed_model(patric_id)
    except:
        pass
    model_dict[patric_id] = mackinac.create_cobra_model_from_modelseed_model(patric_id)
    tempmodel = model_dict[patric_id]
    if len(tempmodel.genes) > 299:
        tempname = path +patric_id+'.sbml'
        cobra.io.write_sbml_model(tempmodel, tempname)
        modelled_strain_ids.append(patric_id)
        print 'Added model %s to model_dict.' %patric_id
        try:
            mackinac.gapfill_modelseed_model(patric_id)
        except:
            pass
        gapfilled_model_dict[patric_id] = mackinac.create_cobra_model_from_modelseed_model(patric_id)
        tempmodel = gapfilled_model_dict[patric_id]
        tempname = path + 'gapfilled/'+patric_id+'.sbml'
        cobra.io.write_sbml_model(tempmodel, tempname)
        print 'Added model %s to gapfilled_model_dict.' %patric_id

f = file(path+'modelled_strain_ids.txt', 'w')
for strain_id in modelled_strain_ids:
    f.write(strain_id+'\n')
f.close()
print 'Done!'
print 'Number of models created is %i.' %len(modelled_strain_ids)

#for model in model_dict:
#    tempmodel = model_dict[model]
#    tempname = '/home/jan/Documents/sbml_models/'+model+'.sbml'
#    cobra.io.write_sbml_model(tempmodel, tempname)
#    
#for model in gapfilled_model_dict:
#    tempmodel = gapfilled_model_dict[model]
#    tempname = '/home/jan/Documents/sbml_models/gapfilled/'+model+'.sbml'
#    cobra.io.write_sbml_model(tempmodel, tempname)