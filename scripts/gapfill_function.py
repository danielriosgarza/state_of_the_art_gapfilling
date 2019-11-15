#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 14:51:22 2019

@author: jan
"""
import gurobipy as gu
from Reaction_class_08_27 import Reaction
import cobra
from copy import deepcopy

def read_biochemistry_table(location):
    """
    Reads a ModelSEED/PATRIC reaction database table, using the id, stoichiometry, direction and equation columns.
    Returns a dictionary where id maps to a list of the other columns.
    """
    biochem = {}
    f = file(location, 'r')
    f.readline() #Skip the header
    for line in f:
        splitted_line = line.strip().split("\t")
        biochem[splitted_line[0]] = [splitted_line[4], splitted_line[9], splitted_line[6]]
        
    return biochem

def create_EX_reactions(metab_dict):
    '''
    For each metabolite in a metab_dict, create a reaction that exclusively exports that metabolite.
    '''
    EX_reactions = {}
    for metab in metab_dict:
        EX_reaction_name = 'EX_' + metab
        if '_e' in metab: #Only create EX reactions for extracellular reactions.
            EX_reactions[EX_reaction_name] = {'lower_bound': -1., 'upper_bound': 1., 'metabolites': {metab:-1.}} #Reversible reaction (Complete medium).
        
    return EX_reactions

def build_gurobi_model(reaction_dict, metabolite_dict, objective_name, model_reactions, candidate_reactions, export_reactions, delta):
    """
    Builds a gurobi model with a variable for each reaction in reaction dict,
    constrains variables using the metabolite dict, and sets the objective of the model.
    The objective is based on the candidate dict, method, objective name and delta.
    Candidate dict includes all reactions that are considered for gap filling, mapping to their cost.
    The objective of the model is defined as:
        MAXIMIZE: delta * objective_name (biomass reaction) - sum(cost * variable in candidate reactions)
    """
    #define the model
    m= gu.Model('mipl')
    #define variables
    for i in reaction_dict:
        if i==objective_name: #Why the if statements, they all do the same
            b=m.addVar(vtype= gu.GRB.CONTINUOUS, name = i, lb = reaction_dict[i][0], ub = reaction_dict[i][1])
            
        elif i in candidate_reactions:
            b=m.addVar(vtype= gu.GRB.CONTINUOUS, name = i, lb = reaction_dict[i][0], ub = reaction_dict[i][1])
            
        elif i in model_reactions:
            b=m.addVar(vtype= gu.GRB.CONTINUOUS, name = i, lb = reaction_dict[i][0], ub = reaction_dict[i][1])
        
        else:
            b=m.addVar(vtype= gu.GRB.CONTINUOUS, name = i, lb = reaction_dict[i][0], ub = reaction_dict[i][1])
            
    m.update()
    #set the objective function
    var = m.getVars()
    coef=[]
    for i in var:
        if i.VarName==objective_name:
            coef.append(delta)
            
        elif i.VarName in candidate_reactions:
            cost = candidate_reactions[i.VarName]
            coef.append(-cost) #Costs are assumed to be positive numbers in input.
            
        elif i.VarName in model_reactions:
            coef.append(0)
            
        elif i.VarName in export_reactions:
            coef.append(0)
            
        else:
            print 'Cannot set objective for %s, not objective_name, model, export or candidate reaction!' %i.VarName

    m.setObjective(gu.LinExpr(coef, var), gu.GRB.MAXIMIZE)
    m.update()    
    #set the stoichiometric constraints for metabolites    
    for i in metabolite_dict:
        var = metabolite_dict[i].keys()
        var = [m.getVarByName(z) for z in var]
        coef = metabolite_dict[i].values()
        m.addConstr(gu.LinExpr(coef, var), 'E', 0, i)  
        
    m.update()
    m.setParam('OutputFlag', False) #Keep gurobi silent
    m.optimize()
    return m

def get_minimum(R, R_flux, R_cost, D, output):
    
    if output == 'min_cost':
        sum_cost = [sum(cost) for cost in R_cost]
        if sum_cost != []:
            min_cost = min(sum_cost)
            min_cost_index = sum_cost.index(min_cost)
            return R[min_cost_index], D[min_cost_index]
    
    if output == 'min_reactions':
        len_reactions = [len(e) for e in R]
        if len_reactions != []:
            min_reactions = min(len_reactions)
            min_reactions_index = len_reactions.index(min_reactions)
            return R[min_reactions_index], D[min_reactions_index]
    
    if output == 'min_flux':
        sum_flux = [sum(flux) for flux in R_flux]
        if sum_flux != []:
            min_flux = min(sum_flux)
            min_flux_index = sum_flux.index(min_flux)
            return R[min_flux_index], D[min_flux_index]
    
    else:
        return [], 0
    
def latendresse_gapfill(all_split_reactions, N, M, split_EX_reactions, B, output):
    """
    Function based on Latendresse BMC Bioinformatics 2014.
    Input:
        all_split_reactions is Reaction class object containing data for all reactions.
        N is reactions in input model.
        M is dictionary of all candidate reactions mapping to their cost.
        split_EX_reactions are all environment defining reactions that the model can use, 
            these are not associated to a cost, hence not found in M.
        B is the name of the biomass reaction.
    Output:
        List of the minimum set of candidate reactions to add to gap fill the model.
    """
    R = []
    R_flux = []
    R_cost = []
    D = []
    alpha = 0
    beta = 2 * len(M)
    reaction_dict = all_split_reactions.get_gurobi_reaction_dict(all_split_reactions.reactions.keys())
    metabolite_dict = all_split_reactions.get_gurobi_metabolite_dict(all_split_reactions.reactions.keys())
    while abs(alpha - beta) > 1:
        delta = int((alpha + beta) / 2.0) #Take the floor function of the mean of alpha and beta.
        print 'Delta is %i.' %delta
        gu_model = build_gurobi_model(reaction_dict, metabolite_dict, B, N, M, split_EX_reactions, delta)
        if gu_model.getVarByName(B).X > 0:
            print 'Flux through biomass reaction is %f.' %gu_model.getVarByName(B).X
            beta = delta
            
            #Store results      
            R.append([var.VarName for var in gu_model.getVars() if var.VarName in N or gu_model.getVarByName(var.VarName).X > 0])
            R_flux.append([gu_model.getVarByName(e).X for e in M if gu_model.getVarByName(e).X > 0])
            R_cost.append([M[e] for e in M if gu_model.getVarByName(e).X > 0])
            D.append(delta)
           
        else:
            print 'Flux through biomass reaction is %f.' %gu_model.getVarByName(B).X
            alpha = delta
            
    minimum_set = get_minimum(R, R_flux, R_cost, D, output) #List that has minimum nr of reactions, sum of cost or sum of flux, dependend on output.
    return  minimum_set

def get_classification(all_reaction_ids, deleted_reactions, gapfill_result, do_not_count):
    
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    
    for rxn in all_reaction_ids:
        if rxn not in do_not_count: #Only count reactions that are not in this set.
            if rxn in deleted_reactions: #rxn IS deleted when making the model incomplete, should be predicted positive.
                if rxn in gapfill_result:
                    tp += 1
                    
                else:
                    fn += 1
                    
            else: #rxn IS negative, should be predicted negative.
                if rxn in gapfill_result:
                    fp += 1
                    
                else:
                    tn += 1
            
    return tp, fp, tn, fn

def make_cobra_metabolites(metab_dict):
    '''
    Make a cobra metabolite object for each metabolite in metab_dict.
    Return a dict mapping the metab_id (str) : object.
    '''
    cobra_metabs = {}
    for metab in metab_dict:
        cobra_metabs[metab] = cobra.Metabolite(metab)
        if '_c' in metab:
            cobra_metabs[metab].compartment = 'c'
        
        elif '_e' in metab:
            cobra_metabs[metab].compartment = 'e'
            
    return cobra_metabs
        
def make_cobra_reaction(reaction_dict, cobra_metabs, rxn):
    '''
    Make a cobra reaction object from a rxn, setting the bounds, and adding metabolites.
    Return the cobra reaction object.
    '''
    reaction = cobra.Reaction(rxn)
    reaction.name = rxn
    reaction.lower_bound = reaction_dict[rxn]['lower_bound']
    reaction.upper_bound = reaction_dict[rxn]['upper_bound']
    metabolites = {cobra_metabs[k]:v for k,v in reaction_dict[rxn]['metabolites'].items()}
    reaction.add_metabolites(metabolites)
    return reaction

def make_cobra_model(reaction_dict, metab_dict, reactions_in_model, objective_name):
    '''
    Make a cobra model of all reactions in reactions_in_model, 
    using a reaction class instance and a metab dict containing all metabs used in the reactions.
    '''
    #Make cobra metabolites and reaction objects.
    cobra_metabs = make_cobra_metabolites(metab_dict)
    cobra_reactions = [make_cobra_reaction(reaction_dict, cobra_metabs, e) for e in reactions_in_model]
    #Make cobra model.
    cobra_model = cobra.Model('tempmodel')
    cobra_model.add_reactions(cobra_reactions)
    cobra_model.objective = objective_name
    return cobra_model

def EX_reaction_fva_check(cobra_model, reactions_to_check):
    '''
    For all reactions in reactions_to_check, look at the ranges of flux variability.
    If both minimum and maximum are 0, it cannot be used without decreasing biomass, so it should be deleted.
    If minimum is -1 and max is 1, then the flux can be any value without decreasing biomass, so it might also be deleted.
    Return the longest possible list without decreasing biomass.
    '''
    fva = cobra.flux_analysis.flux_variability_analysis(cobra_model)
    min0_max0 = []
    minneg1_max1 = []
    for i in reactions_to_check:
        if fva['minimum'][i] == 0 and fva['maximum'][i] == 0:
            min0_max0.append(i)
        if fva['minimum'][i] == -1 and fva['maximum'][i] == 1:
            minneg1_max1.append(i)
    cobra_model.remove_reactions(min0_max0)
    cobra_model.repair()
    obj = cobra_model.optimize().objective_value
    cobra_model.remove_reactions(minneg1_max1)
    cobra_model.repair()
    new_obj = cobra_model.optimize().objective_value
    if new_obj == obj:
        return min0_max0 + minneg1_max1
    else:
        if obj > 0:
            return min0_max0
        else:
            return []

def gapfill(all_input_reactions, input_reaction_ids, input_candidate_reactions, biomass_name, medium = 'complete', default_cost = 1,
            compare_to_added_reactions = False, write_sbml = False, output_reactions = True, latendresse_result_selection = 'min_cost'):
    '''
    Gapfill an incomplete model using the following input:
        - all_input_reactions 
            This is a Reaction class object containing the chemical information of all reactions 
            (bounds, stoichiometry and metabolites).
            
        - input_reaction_ids 
            This is a set containing all the reaction ids in the input model for gap filling.
            
        - input_candidate_reactions 
            This is a dictionary mapping reaction_ids to their cost during gap filling. When reactions are not present in 
            candidate_reactions, their cost will be default_cost.
            
        - biomass_name 
            This is the id correspoding to the biomass reaction in all_split_reactions.
            
        - medium 
            This determines the EX_reactions that the model can use. As default a complete medium is assumed, where the model
            can use all extracellular metabolites that are defined in all_reactions. (This means they are available extra cellularly, 
            not that a transport reaction is added by default!) A reaction class can be given to medium for custom media.
            
        - default_cost
            Every reaction in all_reactions that is not in candidate_reactions or A_incomplete_reaction_ids will get this cost for
            gap filling.
    
    Output:
            
        - compare_to_added_reactions 
            This evalutates all reactions by comparing reactions added during gap filling to the reactions that are given as a set to 
            compare_to_added_reactions. Reactions in compare_added_reactions that are found back in the reactions added during gapfilling
            are counted as TP, added reactions that were not in compare_to_added_reactions are counted as FP, reactions that were not
            in compare_to_added_reactions and not added are counted as TN and reactions that were in compare_to_added_reactions but not added are counted as FN.
            The TP, FP, TN, FN are then returned.
            
        - write_sbml 
            This can be given a file location to write an the gap filled model to a sbml file.
            
        - output_reactions 
            This returns the set of reaction ids that were added during gap filling.
            
        - latendresse_result_selection 
            This is the parameter that specifies which of the gapfill solutions should be returned.
            It can either be the set that has the lowest cost ("min_cost") [default], lowest number of reactions ("min_reactions"), 
            or the lowest flux ("min_flux").
            
    Reactions from all_reactions are added to candidate_reactions if they are not in A_incomplete_reaction_ids or candidate_reactions. 
    This ensures that all reactions are considered as candidates during gap filling, while only reactions given as input with 
    candidate_reactions have a custom cost.
    During gap filling reactions are made unidirectional, so bidirectional reactions are split into a forward and reverse version.
    EX_reactions are added to the model. These specifiy which metabolites are avaiable in the medium. When no custom medium is used,
    the model will be gap filled on a complete medium. The model is gap filled, and the result is changed back into bidirectional reactions.
    Reactions added during gap filling are stored in the set added_reactions. This can be used to calculate TP, FP, TN and FN, where
    added EX reactions are not counted.
    '''
    all_reactions = deepcopy(all_input_reactions)
    candidate_reactions = deepcopy(input_candidate_reactions)
    #Add reactions from all_reactions to candidate_reactions, with cost = default_cost.
    for reaction in all_reactions.reactions:
        if reaction not in input_reaction_ids:
            if reaction not in candidate_reactions:
                candidate_reactions[reaction] = default_cost
        else: 
            if reaction in candidate_reactions:
                del candidate_reactions[reaction] #Delete reaction from candidate_reactions if it is present in the starting model.
        
    #Split bidirectional reactions into a forward and reverse reaction.
    all_split_reactions = Reaction()
    all_split_reactions.reactions = all_reactions.split_all_bidirectional_reactions(reaction_dictionary=all_reactions.reactions)
    #Add reverse reactions to A_split_incomplete_reaction_ids and candidate_reactions.
    A_split_incomplete_reaction_ids = set()
    for reaction in all_split_reactions.reactions:
        forward_version = reaction.replace('_r', '')
        if forward_version in input_reaction_ids:
            A_split_incomplete_reaction_ids.add(reaction)
        if forward_version in candidate_reactions and '_r' in reaction: #If forward version of a reverse reaction is in candidate_reactions.
            candidate_reactions[reaction] = candidate_reactions[forward_version] #Give reverse reaction same cost as forward version.
    
    #Create medium-defining EX reactions or use EX_reactions from custom medium.
    if medium == 'complete':
        metabolites = all_split_reactions.get_gurobi_metabolite_dict()
        EX_reactions = Reaction()
        EX_reactions.reactions = create_EX_reactions(metabolites)
        split_EX_reactions = EX_reactions.split_all_bidirectional_reactions(reaction_dictionary=EX_reactions.reactions)
        all_split_reactions.reactions = all_split_reactions.add_dict(all_split_reactions.reactions, split_EX_reactions)
        all_reactions.reactions = all_reactions.add_dict(all_reactions.reactions, EX_reactions.reactions)
    else:
        EX_reactions = Reaction()
        EX_reactions.reactions = medium
        split_EX_reactions = EX_reactions.split_all_bidirectional_reactions(reaction_dictionary=EX_reactions.reactions)
        
    
    #Run gapfilling algorithm
    split_gapfill_result, delta = latendresse_gapfill(all_split_reactions, A_split_incomplete_reaction_ids, candidate_reactions, 
                                               split_EX_reactions, biomass_name, latendresse_result_selection)
    gapfill_result = set([r.replace('_r', '') for r in split_gapfill_result])
    added_reactions = gapfill_result.difference(input_reaction_ids) #All reactions that are added to the model during gapfilling.
    #Create cobra model
    metab_dict = all_reactions.get_gurobi_metabolite_dict(all_reactions.reactions.keys())
    cobra_model = make_cobra_model(all_reactions.reactions, metab_dict, gapfill_result, biomass_name)
    objective_value = cobra_model.optimize().objective_value
    print 'Objective value is %f.' %objective_value

    #Check EX reaction essentiality
    EX_added = [e for e in list(added_reactions) if 'EX' in e]
    EX_to_delete = EX_reaction_fva_check(cobra_model, EX_added)
    gapfill_result_checked_EX = [reaction for reaction in list(gapfill_result) if reaction not in EX_to_delete]
    #Create cobra model without non-essential EX reactions
    cobra_model = make_cobra_model(all_reactions.reactions, metab_dict, gapfill_result_checked_EX, biomass_name)
    objective_value_checked_EX = cobra_model.optimize().objective_value
    print 'Objective value checked EX is %f.' %objective_value_checked_EX
    #Return output
    Output = [delta]
    if compare_to_added_reactions is not False:
        TP, FP, TN, FN = get_classification(all_reactions.reactions, compare_to_added_reactions, added_reactions, EX_reactions.reactions)
        Output += [TP, FP, TN, FN]
        
    if write_sbml: #Needs to be added
        if objective_value_checked_EX > 0:
            output_cobra_model = make_cobra_model(all_reactions.reactions, metab_dict, gapfill_result_checked_EX, biomass_name)
            cobra.io.write_sbml_model(output_cobra_model, write_sbml)
            Output.append('checked_EX')
        else:
            output_cobra_model = make_cobra_model(all_reactions.reactions, metab_dict, gapfill_result, biomass_name)
            cobra.io.write_sbml_model(output_cobra_model, write_sbml)
            Output.append('no_check_EX')
    
    if output_reactions:
        if objective_value_checked_EX > 0:
            Output.append(gapfill_result_checked_EX)
        else:
            Output.append(gapfill_result)
        
    return Output


    
        
        
        
    
    
    
