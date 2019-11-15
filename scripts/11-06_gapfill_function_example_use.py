#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 12:36:37 2019

@author: jan

This file shows some examples of how to use the gapfill function, which can be imported from gapfill_function.py.
The gapfill function is based on the FastGapFilling algorithm, 
see Latendresse, M. (2014). Efficiently gap-filling reaction networks. BMC bioinformatics, 15(1), 225.
"""
from Reaction_class_08_27 import Reaction
import gapfill_function

'''
Example 1. Gapfilling a model

Load the model'''
model_location = '/home/jan/Documents/gapfill_function/ungapfilled_model.sbml'
ungapfilled_model = Reaction(model=model_location) 
#incomplete_model now contains names, bounds and stoichiometric information of all reactions in the model.
#These are all stored in a dictionary, called incomplete_model.reactions.

ungapfilled_reaction_ids = set(ungapfilled_model.reactions) #Set of reaction ids in ungapfilled_model.

'''Load the database with reactions'''
biochem = gapfill_function.read_biochemistry_table('/home/jan/Downloads/modelseed_reactions.tsv')
db_reactions = Reaction(biochem_input=biochem)
#db_reactions now contains names, bounds and stoichiometric information of all reactions in the database.

'''Combine all reactions into one object'''
all_reactions = ungapfilled_model
all_reactions.reactions = all_reactions.add_dict(all_reactions.reactions, db_reactions.reactions) #add_dict returns all items from both input dictionaries

'''Gapfill'''
gapfilled_model_location = '/home/jan/Documents/gapfill_function/gapfilled_model.sbml'
result = gapfill_function.gapfill(all_reactions, ungapfilled_reaction_ids, {}, 'bio1', write_sbml=gapfilled_model_location)
#Certain reactions can be given a custom cost for gap filling, by adding them to the now empty dictionary {}. See example 2.
#The biomass name must be specified, here it is 'bio1'.
#write_sbml is off on default, but here it was given a location to write the gapfilled model to.

#result[0] is the delta of the gapfill solution from the latendresse algorithm.
#result[1] tells you if medium defining reactions are checked. 'checked_EX' means they are checked, and unusable ones are deleted.
#result[3] lists all the reaction ids that are in the gap filled model.

'''
Example 2. Comparing two different methods

Here, an ungapfilled model is made more incomplete, by removing some reactions.
The resulting incomplete_model is then gapfilled with two different methods.
For each method, the reactions that are added are compared with the reactions that were originally deleted, to give TP, FP, TN and FN.
These values may then be used to compare different methods in their ability to gapfill with the correct reactions.

Load the model'''
model_location = '/home/jan/Documents/gapfill_function/ungapfilled_model.sbml'
ungapfilled_model = Reaction(model=model_location) 
#incomplete_model now contains names, bounds and stoichiometric information of all reactions in the model.
#These are all stored in a dictionary, called incomplete_model.reactions.

ungapfilled_reaction_ids = set(ungapfilled_model.reactions) #Set of reaction ids in ungapfilled_model.
deleted_reactions = set(ungapfilled_model.reactions.keys()[0:10]) #Select 10 reactions for deletion
incomplete_reaction_ids = ungapfilled_reaction_ids.difference(deleted_reactions) #Delete reactions from model.

'''Load the database with reactions'''
biochem = gapfill_function.read_biochemistry_table('/home/jan/Downloads/modelseed_reactions.tsv')
db_reactions = Reaction(biochem_input=biochem)
#db_reactions now contains names, bounds and stoichiometric information of all reactions in the database.

'''Combine all reactions into one object'''
all_reactions = ungapfilled_model
all_reactions.reactions = all_reactions.add_dict(all_reactions.reactions, db_reactions.reactions) #add_dict returns all items from both input dictionaries

'''Create custom costs'''
candidate_reactions = {rxn_id: 0.0 for rxn_id in list(deleted_reactions)} #Dict mapping rxn_id : cost

'''Gapfill'''
results = []
#Gapfill with default costs (method 1)
result = gapfill_function.gapfill(all_reactions, incomplete_reaction_ids, {}, 'bio1', compare_to_added_reactions=deleted_reactions, output_reactions=False)
result.append('control')
results.append(result)

#Gapfill with custom costs (method 2)
result = gapfill_function.gapfill(all_reactions, incomplete_reaction_ids, candidate_reactions, 'bio1', compare_to_added_reactions=deleted_reactions, output_reactions=False)
result.append('custom_costs')
results.append(result)

#The number of TP, FP, TN or FN can now be compared between the different methods
#Because we added an extra label to the results, we can keep track of which results belong to which methods, even when gapfilling many models.

'''
Notable extra features:
    medium may be customised, by providing a dictionary mapping 
        {EX_rxn_id : {'lower_bound' : lower_bound, 'upper_bound' : upper_bound, 'metabolites' : metabolite }}
        For more info see gapfill_function.create_EX_reactions()
    default_cost may be customised by setting default cost = some_cost, where some_cost is a float or int.
    
'''
    


