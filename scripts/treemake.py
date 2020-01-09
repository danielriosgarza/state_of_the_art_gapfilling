#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 13:49:35 2019

@author: jan
"""
import os
import matplotlib.pyplot as plt
import numpy as np

'''
Get a list of Patric IDs that should be included in the tree.
'''
model_folder = '/home/jan/genera_models/'
id_list = [e.replace('.sbml', '') for e in os.listdir(model_folder) if '.sbml' in e]

'''
Make a bash script to get .faa files for each id in id_list.
'''

output_path = '/home/jan/trees/genera_tree/'
os.mkdir(output_path+'faa')
ftp_bash_script = file(output_path + 'faa/get_faa.sh', 'w')
for strain in id_list:
    ftp_bash_script.write('wget ftp://ftp.patricbrc.org/genomes/' + strain + '/' + strain + '.PATRIC.faa\n')

ftp_bash_script.close()   

found_faa = [e.replace('.PATRIC.faa', '') for e in os.listdir(output_path+ 'faa') if 'PATRIC' in e]
missing_faa = [e for e in id_list if e not in found_faa]
re_ftp_bash_script = file(output_path + 'faa/re_get_faa.sh', 'w')
for strain in missing_faa:
    re_ftp_bash_script.write('wget ftp://ftp.patricbrc.org/genomes/' + strain + '/' + strain + '.PATRIC.faa\n')
    
re_ftp_bash_script.close()
'''
Change >Name of each protein in .faa files
'''
for strain in found_faa:
    fin = file(output_path + 'faa/' + strain + '.PATRIC.faa', 'r')
    fout = file(output_path + 'faa/' + strain + '.faa', 'w')
    for line in fin:
        if '>' in line:
            new_line = line.split(' ')[0]
            fout.write(new_line + '\n')
            
        else:
            fout.write(line)
    
    fin.close()
    fout.close()
            
'''
Do hmmsearch for each proteome.
Use --tblout to save a parseble table which is used to determine the best hit.
Use -o to output a file from which the alignment from the protein sequence to the profile can be parsed.
'''
profiles = '/home/jan/Downloads/genes.hmm'
#os.mkdir(output_path + 'hmm')
hmm_bash = file(output_path + 'hmm/get_hmm.sh', 'w')
for strain in found_faa:
    hmm_bash.write('hmmsearch --tblout ' + output_path + 'hmm/' + strain + '.txt -o '+output_path+'hmmalign/'+strain+' --cut_tc --cpu 8 ' + profiles + ' ' + output_path + 'faa/' + strain + '.faa\n')
hmm_bash.close()

'''
Find best hit for each gene, save target name.
Parse both output files to get a dict that maps gene name to protein sequence for every strain.
'''

def find_best_hits(hmm_file_location):
    hits = {}
    evalues = {}
    hmm_file = file(hmm_file_location, 'r')
    for line in hmm_file:
        if line[0] == '#':
            pass
        else:
            split_line = line.split()
            target = split_line[0]
            query_name = split_line[2]
#            print target
#            print query_name
            if query_name not in hits:
                hits[query_name] = target
                evalues[query_name] = split_line[4]
    return hits, evalues

#nr_of_genes = 71
hits = {}
e_values = {}
for strain in found_faa:
    hits[strain], e_values[strain] = find_best_hits(output_path + 'hmm/' + strain + '.txt')

profile_seqs = {strain:{} for strain in hits}
for strain in profile_seqs:
    prot_seq = ''
    protein_name = 'tempname'
    gene_map = {v:k for k,v in hits[strain].items()}
    hmm_output = file(output_path+'hmmalign/'+strain,'r')
    for line in hmm_output:
        if line[0] != '#':
            if line[0] == '>':
                if protein_name in gene_map:
                    gene_name = gene_map[protein_name]
                    profile_seqs[strain][gene_name] = prot_seq #Store protein sequence to profile_seqs dict under gene key.
                protein_name = line.strip().split()[1]
                prot_seq = ''
            elif protein_name in line:
                prot_seq += line.strip().split()[2].replace('-', '').upper() #Leave out gaps in sequence
                
#Throw out bad hits
bad_hits = []
for strain in profile_seqs.keys():
    for gene in profile_seqs[strain].keys():
        if len(profile_seqs[strain][gene]) < 10:
            del profile_seqs[strain][gene]
            bad_hits.append((strain, gene))
#all_genes_found = [strain for strain in hits.keys() if len(hits[strain]) == nr_of_genes]
nr_genes_found = {strain:len(hits[strain]) for strain in hits.keys()}

genes = {}
for strain in hits:
    for gene in hits[strain]:
        if gene not in genes:
            genes[gene] = {strain}
        else:
            genes[gene].add(strain)

nr_times_found = {gene:len(genes[gene]) for gene in genes.keys()}

#Check the hits
nr_genome = np.asarray(nr_genes_found.values())
plt.hist(nr_genome)
nr_gene = np.asarray(nr_times_found.values())
plt.hist(nr_gene)

'''
Make multifasta file for each gene
'''

def insert_newlines(string, every=60):
    return '\n'.join(string[i:i+every] for i in xrange(0, len(string), every))

for gene in genes:
    gene_mfa = file(output_path+'mfa/'+gene, 'w')
    for strain in profile_seqs:
        if gene in profile_seqs[strain]:
            gene_mfa.write('>'+strain + '\n')
            sequence_w_newlines = insert_newlines(profile_seqs[strain][gene])
            gene_mfa.write(sequence_w_newlines.strip()+'\n') #.strip() to aviod double \n\n.
    gene_mfa.close()
    
'''
Run a multiple sequence aligner over each mfa file.
'''
if not os.path.isdir(output_path+'msa'):
    os.mkdir(output_path+'msa')
    
run_msa = file(output_path+'mfa/run_msa.sh', 'w')
for gene in genes:
    command = './clustalo-1.2.4-Ubuntu-x86_64 -i ' + output_path + 'mfa/' + gene + ' -o ' + output_path + 'msa_08-06/' + gene + ' --dealign --threads $(nproc)\n'
    run_msa.write(command)
run_msa.close()

def check_alignment_len(path_to_file):
    f = file(path_to_file, 'r')
    lens = []
    counter = 0
    for line in f:
        if '>' in line:
            if counter != 0:
                lens.append(counter)
            counter = 0
        else:
            counter += len(line.strip())
    return lens

msa_lens = {}
for gene in genes:
    msa_lens[gene] = set(check_alignment_len(output_path+'msa_trimmed/'+gene))

'''
Write a bash script that trims every gene msa file.
'''

trimbash = file(output_path+'msa_08-06/run_trimal.sh', 'w')
for gene in genes:
    trimbash.write('trimal -in '+gene+' -out '+output_path+'msa_trimmed/'+gene+' -automated1\n')
trimbash.close()

'''
Write one msa file where the seqs of all genes are concatenated per strain.
'''
#Test if each strain has a sequence for enough genes to be considered.
usable_strains = [strain for strain in hits if len(hits[strain].values()) > 10]

for strain in hits:
    if len(hits[strain].values()) < 10:
        print strain
        
strains_found_per_gene = {gene:[] for gene in genes}
for gene in genes:
    f = file(output_path+'msa_trimmed/'+gene, 'r')
    for line in f:
        if '>' in line:
            strains_found_per_gene[gene].append(line.strip().split()[0][1:])

'''
Add gaps for genes that were not found.
'''

def get_alignment_len(path_to_file):               
    f = file(path_to_file, 'r')
    f.readline() #Skip first >line
    alignment_len = 0
    for line in f:
        if '>' not in line:
            alignment_len += len(line.strip())
        else:
            break
    f.close()
    return alignment_len

for gene in genes:
    gene_len = get_alignment_len(output_path+'msa_trimmed/'+gene)
    gap = insert_newlines('-'*gene_len)
    for strain in usable_strains:
        if strain not in strains_found_per_gene[gene]:
            with open(output_path+'msa_trimmed/'+gene, 'a') as f:
                f.write('>'+strain+'\n')
                f.write(gap+'\n')
                f.close()
                
concatenated_seqs = {strain:'' for strain in usable_strains}
concatenation_len = {strain:0 for strain in usable_strains}
for gene in genes:
    f = file(output_path+'msa_trimmed/'+gene, 'r')
    for line in f:
        if '>' in line:
            strain = line.strip().split()[0][1:]
        else:
            if strain in concatenated_seqs:
                concatenated_seqs[strain] += line.strip()
                concatenation_len[strain] += len(line.strip())

museal = file(output_path+'msa_trimmed/concatenated_msa_trimmed', 'w')
for strain in concatenated_seqs:
    museal.write('>'+strain+'\n')
    seq = insert_newlines(concatenated_seqs[strain]).strip()
    museal.write(seq+'\n')
museal.close()


'''
Run iqtree on concatenated_msa_trimmed.
With iqtree -s concatenated_msa_trimmed -bb 1000 -alrt 1000 -nt AUTO
'''
