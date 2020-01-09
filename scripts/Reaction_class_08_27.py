#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 09:12:52 2019

@author: jan
"""

#imports
import cobra

import numpy as np

import os

import ast

#import gurobipy as gu

class Reaction:
    
    def __init__(self, reaction_database=None, model_folder=None, model_list=None, model=None, biochem_input=None, dbtype=2):
        
        self.reaction_database = self.__get_reactions_from_database(reaction_database, dbtype)
        
        self.model_folder = model_folder
        
        self.model_list = model_list
        
        self.model = self.__get_model(model)
        
        self.biochem_input = self.__get_reactions_from_biochem_input(biochem_input)
        
        self.reactions = self.__get_reactions()
        
        print 'Reaction class version 19_08_27'
    
    def __get_reactions_from_database(self, reaction_database, dbtype):
        '''
        db type 1: ModelSEED
        db type 2: python_dict_string
        '''
        if reaction_database is None:
            return None
        
        elif dbtype==2:
            f= file(reaction_database)
            database = ast.literal_eval(f.read())
            f.close()
            return database
         
    def __get_model(self, model):
        
        if model is None:
            return None
        
        else:
            self.model=cobra.io.read_sbml_model(model)
            return self.model
    
    def __get_reactions_from_model(self, model):
        '''
        reactions[reaction_id]={upper_bound:, lower_bound:, metabolites:{metabolite:stoichiometry}}
        '''
        reactions={}
        for reaction in model.reactions:
            reactions[reaction.id] = {'lower_bound':reaction.lower_bound/1000.0, 'upper_bound':reaction.upper_bound/1000.0}
            mets = reaction.metabolites
            reactions[reaction.id]['metabolites']={i.id:mets[i] for i in mets}
            
        return reactions
    
    def __get_reactions_from_biochem_input(self, biochem_input):
        '''
        reactions[reaction_id]={upper_bound:, lower_bound:, metabolites:{metabolite:stoichiometry}}
        '''
        if biochem_input is None:
            return None
        
        reactions = {}
        for reaction in biochem_input:
            reaction_id = str(reaction.replace("_c", "")+"_c")
            direction = biochem_input[reaction][1]
            if direction == '=':
                reactions[reaction_id] = {'lower_bound':-1.0, 'upper_bound':1.0}
                
            elif direction == '>':
                reactions[reaction_id] = {'lower_bound':0.0, 'upper_bound':1.0}
                
            elif direction == '<': #Flip the direction of this reaction.
                reactions[reaction_id] = {'lower_bound':0.0, 'upper_bound':1.0}
            
            if reaction_id in reactions: #Check if reaction was added to the dict (could be skipped if direction is not specified correctly)..
                correct_reaction = True
                mets = {}
                metabolites = biochem_input[reaction][0].split(';')
                for i in metabolites:
                    metabolite = i.split(':')
                    if len(metabolite) > 1: #Check if metabolite is specified.
                        stoc = metabolite[0]
                        cpd = metabolite[1]
                        loc = metabolite[2]
                        if int(loc) == 0:
                            name = cpd + '_c'
                            
                        elif int(loc) == 1:
                            name = cpd + '_e'
                        
                        else:
                            correct_reaction = False
                            
                        if direction == '<':
                            mets[name] = -float(stoc) #Flip the sign of the stoichiometry of reversed reactions.
                                
                        else:
                            mets[name] = float(stoc)
                                
                if correct_reaction: #Check if location is specified correctly for all metabolites.   
                    reactions[reaction_id]['metabolites'] = {i:mets[i] for i in mets}
                    
                else:
                    del reactions[reaction_id]
        
        return reactions
    
    def add_dict(self, dict_x, dict_y, suffix = '_a'):
        '''
        Add all items from dict_y to dict_x. Same keys mapping to different values will be added with a suffix plus a counter.
        '''
        x = {k:v for k,v in dict_x.items()}
        y = {k:v for k,v in dict_y.items()}
        for key in y:
            if key not in x:
                x[key] = y[key]
                    
            else:
                if y[key] == x[key]:
                    pass
                    
                else:
                    reaction_added = False
                    counter = -1
                    while reaction_added == False:
                        counter += 1
                        new_name = key+suffix+str(counter)
                        if new_name not in x:
                            x[new_name] = y[key]
                            reaction_added = True
                            
                        else:
                            if y[key] == x[new_name]:
                                reaction_added = True
    
        return x
    
    def __get_reactions(self):
        '''
        reactions[reaction_id]={upper_bound:, lower_bound:, metabolites:{metabolite:stoichiometry}}
        '''
        reaction_dict={}
        if self.model is not None:
            model_reactions = self.__get_reactions_from_model(self.model)
            reaction_dict = self.add_dict(reaction_dict, model_reactions)
         
        if self.model_folder is not None:
            if self.model_list is None:
                list_of_models=os.listdir(self.model_folder)
                
            else:
                list_of_models = [i.replace('.sbml', '') +'.sbml' for i in self.model_list] #Ensure .sbml is behind every model ID.
                
            for i in list_of_models:
                if '.sbml' in i:
                    mod=cobra.io.read_sbml_model(self.model_folder+i)
                    mod_reac=self.__get_reactions_from_model(mod)
                    reaction_dict= self.add_dict(reaction_dict, mod_reac)
        
        if self.reaction_database is not None:
            reaction_dict = self.add_dict(reaction_dict, self.reaction_database)
            
        if self.biochem_input is not None:
            reaction_dict = self.add_dict(reaction_dict, self.biochem_input)
            
        self.reactions = reaction_dict
        return self.reactions
    
    def get_gurobi_reaction_dict(self, list_of_reactions=None):
        '''
        Get a reaction dict as input for gurobi model, for reactions on list_of_reactions, or all reactions in self.reactions.
        
        Reaction_dict = { reaction_id : (lower_bound, upper_bound) }
        '''
        if self.reactions is None:
            return None
        
        if list_of_reactions is None:
            return {i:(self.reactions[i]['lower_bound'],self.reactions[i]['upper_bound'] ) for i in self.reactions}
            
        else:
            return {i:(self.reactions[i]['lower_bound'],self.reactions[i]['upper_bound'] ) for i in list_of_reactions}
    
    def get_gurobi_metabolite_dict(self, list_of_reactions=None):
        '''
        Get a metabolite dict as input for gurobi model, for reactions on list_of_reactions, or all reactions in self.reactions.
        
        Metab_dict = { metabolite_id : { reaction_id : stoichiometry } }
        '''
        if self.reactions is None:
            return None
        
        metab_dict = {}
        if list_of_reactions is None:
            for reaction in self.reactions:
                for met in self.reactions[reaction]['metabolites']:
                    if met not in metab_dict:
                        metab_dict[met] = {reaction : self.reactions[reaction]['metabolites'][met]}
                        
                    else:
                        metab_dict[met][reaction] = self.reactions[reaction]['metabolites'][met]
        
        else:
            for reaction in list_of_reactions:
                for met in self.reactions[reaction]['metabolites']:
                    if met not in metab_dict:
                        metab_dict[met] = {reaction : self.reactions[reaction]['metabolites'][met]}
                        
                    else:
                        metab_dict[met][reaction] = self.reactions[reaction]['metabolites'][met]
        
        return metab_dict
    
    def split_bidirectional_reaction(self, reaction):
        
        if float(reaction['lower_bound']) == -1 and float(reaction['upper_bound']) == 1: #Check if bidirectional.
            split_reaction_f = {k:v for k,v in reaction.items()} #Deep copy of reaction.
            split_reaction_f['lower_bound'] = 0.0 #Make reaction unidirectional.
            split_reaction_r = {k:v for k,v in split_reaction_f.items()} #Deep copy of forward reaction.
            split_reaction_r['metabolites'] = {i:-reaction['metabolites'][i] for i in reaction['metabolites']}
            return split_reaction_f, split_reaction_r
        
        else:
            return reaction, False
        
    def split_all_bidirectional_reactions(self, reaction_dictionary=None):
        '''
        Split all bidirectional reactions in a dictionary (of similar structure to self.reactions), or all reactions in self.reactions.
        Splitted reactions become two unidirectional reactions, reaction_id and reaction_id_r.
        When using a dictionary, this does not change the input dictionary, only the output dictionary has splitted reactions.
        When using self.reactions, these are changed.
        '''
        if reaction_dictionary is None:
            self_reaction_list = self.reactions.keys()
            for reaction in self_reaction_list:
                r1, r2 = self.split_bidirectional_reaction(self.reactions[reaction])
                if r2:
                    r2_name = reaction + '_r'
                    self.reactions[reaction] = r1
                    self.reactions[r2_name] = r2
                    
            print 'Done splitting all bidirectional reactions in self.reactions.'
                    
        else:
            split_reaction_dictionary = {k:v for k,v in reaction_dictionary.items()} #Deep copy of reaction_dictionary.
            for reaction in reaction_dictionary:
                r1, r2 = self.split_bidirectional_reaction(split_reaction_dictionary[reaction])
                if r2:
                    r2_name = reaction + '_r'
                    split_reaction_dictionary[reaction] = r1
                    split_reaction_dictionary[r2_name] = r2
                    
            return split_reaction_dictionary
        
    def split_reaction(self, reaction):
        '''
        Splits a reaction into a forward and reverse version (regardeless of original directionality), 
        where both have bounds (0., 1.) but the sign of involved metabolites is flipped for the _r reaction.
        '''
        split_reaction_f = {k:v for k,v in reaction.items()}
        split_reaction_f['lower_bound'] = 0.0
        split_reaction_f['upper_bound'] = 1.0
        split_reaction_r = {k:v for k,v in split_reaction_f.items()}
        split_reaction_r['metabolites'] = {i:-reaction['metabolites'][i] for i in reaction['metabolites']}
        return split_reaction_f, split_reaction_r
        
    def split_ALL_reactions(self, reaction_dictionary=None):
        '''
        Split all reactions in a dictionary into a forward and reverse reaction, regardless of original directionality.
        This means that all reactions will become bidirectional.
        When using a dictionary, the input dictionary does not change.
        When using self.reactions, these are changed.
        '''
        if reaction_dictionary is None:
            self_reaction_list = self.reactions.keys()
            for reaction in self_reaction_list:
                r1, r2 = self.split_reaction(self.reactions[reaction])
                r2_name = reaction + '_r'
                self.reactions[reaction] = r1
                self.reactions[r2_name] = r2
                
            print 'Done splitting ALL reactions in self.reactions'
            
        else:
            split_reaction_dict = {k:v for k,v in reaction_dictionary.items()} #Deep copy of reaction_dictionary, to leave inpun unchanged.
            for reaction in reaction_dictionary:
                r1, r2 = self.split_reaction(self.split_reaction_dict[reaction])
                r2_name = reaction + '_r'
                self.reactions[reaction] = r1
                self.reactions[r2_name] = r2
                
            return split_reaction_dict
            
    def write_database(self, file_name):
        
        if self.reactions is None:
            return None
            
        f = file(file_name, 'w')
        f.write(str(self.reactions))
        f.close()
                    
def read_biochemistry_table(location):
    biochem = {}
    f = file(location, 'r')
    f.readline() #Skip the header
    for line in f:
        splitted_line = line.strip().split("\t")
        biochem[splitted_line[0]] = [splitted_line[4], splitted_line[9], splitted_line[6]]
        
    return biochem
