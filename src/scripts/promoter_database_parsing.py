#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 12:00:45 2020
@author: maalcantar
"""

import os
import argparse
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from promoter_parsing_util import *

"""
EPDnew
Contains sequences for 15 different eukaryotic species
1. C. elegans (worm)
2. D. melanogaster (fly)
3. D. rerio (zebrafish)
4. R. norbegicus (rat)
5. Z. mays (maize)
6. H. sapien (humans)
7. G. gallus (chicken)
8. S. pombe (fission yeast)
9. A. thaliana (arabidopsis)
10. M. mulatta (monkey)
11. P. falciparum (malaria-causing parasite)
12. M. musculus (mouse)
13. S. cerevisiae (yeast)
14. A. mellifera (bee)
15. C. familiaris (dog)
"""

def parse_EPDnew(save_csv=False):
    """
    parse through EPDnew database and apply QC / filtering

    inputs
    save_csv: boolean argument indicating whether csv files should be saved

    outputs
    EPDnew_promoters_df: dataframe containing EPDnew promoters and metadata,
            if save_csv, then csv file is saved to data directory

    """

    # hard-coding directory containing data of interest
    directory_str = '../../data/EPDnew'
    directory = os.fsencode(directory_str)

    # initializing lists that will contain information of interest
    promoter_organism = []
    promoter_database = []
    promoter_sequences = []
    promoter_gene = []
    promoter_range = []
    promoter_notes = []

    # loop through all files
    # each file contains promoters for a particular organism and motif (e.g., h. sapien TATA promoters)
    for file in os.listdir(directory):

        # establishing file path and name
        filename = os.fsdecode(file)
        filepath = directory_str + '/' + filename

        # initializing gene number and sequence (the latter helps with appending sequences to list)
        counter_gene = 0
        sequence = ''

        # ignore hidden files
        if filename[0] != '.':

            with open(filepath) as fp:

                line = fp.readline()

                while line:

                    # The '>' signifies the start of a new gene
                    if line[0] == '>':

                        # append last completed sequence if we've already been through 1 iteration
                        if counter_gene >= 1:

                            promoter_sequences.append(sequence)

                        sequence = ''
                        counter_gene += 1

                        line_split = line.split(' ')

                        # extracting gene regulated by the current promoter
                        promoter_gene.append(line_split[1])

                        # edge case where gene does not have a specified range
                        try:
                            range_index = line.split(' ').index('range')
                            range_temp = line.split(' ')[range_index+2] + ' to ' + line.split(' ')[range_index+6].strip('.\n')
                            promoter_range.append(range_temp)
                        except ValueError:
                            promoter_range.append('')

                        # append organism name, database where promoter was obtained, and notes (e.g., TATA or nonTATA promoter)
                        promoter_organism.append(filename.split('_')[2] + '_' + filename.split('_')[3])
                        promoter_database.append(filename.split('_')[1])
                        promoter_notes.append(filename.split('_')[4].replace(' ',''))

                    else: # keep appending nucleotides until we reach a new sequence
                        sequence = sequence + line
                    line = fp.readline()
            promoter_sequences.append(sequence)

        # organizing data into a dictionary
        EPDnew_data_dict = {
            'organism': promoter_organism,
            'database/source': promoter_database,
            'DNA sequence':promoter_sequences,
            'regulated gene': promoter_gene,
            'range (with respect to TSS)': promoter_range,
            'sigma factor/motif': promoter_notes
                            }

    # organizing data into dataframe
    EPDnew_promoters_df = pd.DataFrame.from_dict(EPDnew_data_dict,orient='index').transpose()
    EPDnew_promoters_df = QC_DNA_sequences(EPDnew_promoters_df)

    # save as csv
    if save_csv:
        EPDnew_promoters_df.to_csv('../../data/parsed_promoter_data/20191125_EPDnew_promoters.csv',index=False)

    return(EPDnew_promoters_df)

# creating stacked barplot showing number of promoters from each organism and their specific motifs
def plot_EPDnew_data(motifs_occurances_df, save_plot=False):

    """
    create bar plot showing EPDnew promoters stratified by organism and motif

    # inputs
    motifs_occurances_df: dataframe containing EPDnew organisms x motifs (the output of EPDnew_motifs function)
    save_plot: boolean indicating whether plot should be saved

    # outputs: if save_plot, then a bar plot is saved
    """

    dark_blue = '#4D98FA'
    dark_yellow = '#FAAE65'
    light_blue = '#60B5E6'
    purple = '#CC7DAA'
    plot_colors = [dark_blue, dark_yellow, light_blue, purple]

    # create horizontal barplot that specifies promoter by organism and motif
    motifs_occurances_df.iloc[:,:-1].plot.barh(stacked=True, color = plot_colors ,figsize=(10,7))
    plt.xlabel('Count', size=24,fontname="Arial"),
    plt.ylabel('Organism', size=24,fontname="Arial")
    plt.title('Promoter counts by organism and motif',fontname="Arial", size=24)
    plt.xticks(fontsize = 16,fontname="Arial")
    plt.yticks(fontsize = 16,fontname="Arial")
    plt.legend(loc='lower right', prop={'size': 16})
    plt.tight_layout()

    if save_plot:
        plt.savefig('../../figs/promoter_database_summaries/20191127_EPDnew_promoters_by_organisms_and_motifs.pdf')
        plt.savefig('../../figs/promoter_database_summaries/20191127_EPDnew_promoters_by_organisms_and_motifs.png', dpi=400)
        plt.savefig('../../figs/promoter_database_summaries/20191127_EPDnew_promoters_by_organisms_and_motifs.svg')

"""
RegulonDB
Contains sequences for E. coli K12, stratified by cognate sigma factor
"""

def parse_RegulonDB(save_csv=False):

    """
    parse through RegulonDB database and apply QC / filtering

    inputs
    save_csv: boolean argument indicating whether csv files should be saved

    outputs
    regulonDB_promoters_df: dataframe containing RegulonDB promoters and metadata,
            if save_csv, then csv file is saved to data directory

    """

    # Download regulon DB promoter database and organize dataframe
    regulonDB_promoters_df = pd.read_csv('../../data/RegulonDB/20191127_PromoterSet.csv', header=36)
    regulonDB_promoters_df = rename_promoter_df_columns(regulonDB_promoters_df,promoter_df_or_type='RegulonDB')
    regulonDB_promoters_df = QC_DNA_sequences(regulonDB_promoters_df)

    if save_csv:
        regulonDB_promoters_df.to_csv('../../data/parsed_promoter_data/20191202_RegulonDB_promoters.csv',index=False)

    return(regulonDB_promoters_df)

def plot_RegulonDB_data(Ecoli_sigma_factor_occurances, save_plot=False):

    """
    create bar plot showing RegulonDB promoters stratified by cognate sigma factor

    # inputs
    Ecoli_sigma_factor_occurances: dataframe containing RegulonDB organisms x motifs (the output of RegulonDB_motifs function)
    save_plot: boolean indicating whether plots should be saved

    # outputs: if save_plot, then a bar plot is saved
    """

    # plot promoter counts by sigma factor
    Ecoli_sigma_factor_occurances.astype(int).plot.barh(figsize=(10,7), color = '#C4C4C4')
    plt.xlabel('Count', size=24,fontname="Arial")
    plt.ylabel('E. coli sigma factor', size=24,fontname="Arial")
    plt.title('E. coli promoter counts by sigma factor', size=24,fontname="Arial")
    plt.xticks(fontsize = 16,fontname="Arial") # work on current fig
    plt.yticks(fontsize = 16,fontname="Arial") # work on current fig
    plt.legend(loc='lower right', prop={'size': 16})
    plt.tight_layout()

    if save_plot:
        plt.savefig('../../figs/promoter_database_summaries/20191127_RegulonDB_promoters_by_ecoli_sigmaFactor.pdf')
        plt.savefig('../../figs/promoter_database_summaries/20191127_RegulonDB_promoters_by_ecoli_sigmaFactor.png', dpi=400)
        plt.savefig('../../figs/promoter_database_summaries/20191127_RegulonDB_promoters_by_ecoli_sigmaFactor.svg')

"""
DBTBS
Contains sequences for B. subtilis, stratified by cognate sigma factor
"""

def parse_DBTBS(save_csv=False):

    """
    parse through DBTBS database and apply QC / filtering

    inputs
    save_csv: boolean argument indicating whether csv files should be saved

    outputs
    DBTBS_promoters_df: dataframe containing DBTBS promoters and metadata,
            if save_csv, then csv file is saved to data directory

    """

    # merge Bacillus subtilus promoter excel sheets into one comprehensive dataframe
    DBTBS_all_promoter_sheets_df = pd.read_excel('../../data/DBTBS/20191129_DBTBS_promoter_sequences.xlsx', sheet_name=None, ignore_index=True)
    DBTBS_promoters_df = pd.concat(DBTBS_all_promoter_sheets_df.values()).fillna('')# #dropna(thresh=1)
    DBTBS_promoters_df.reset_index(drop=True, inplace=True)
    DBTBS_promoters_df = rename_promoter_df_columns(DBTBS_promoters_df, promoter_df_or_type='DBTBS')

    DBTBS_promoters_df = QC_DNA_sequences(DBTBS_promoters_df)

    # save to csv
    if save_csv:
        DBTBS_promoters_df.to_csv("../../data/parsed_promoter_data/20191129_DBTBS_promoter_sequences.csv",index=False)

    return(DBTBS_promoters_df)

def plot_DBTBS_data(Bsubtilis_sigma_factor_occurances, save_plot=False):

    """
    create bar plot showing DBTBS promoters stratified by cognate sigma factor

    # inputs
    Bsubtilis_sigma_factor_occurances: dataframe containing DBTBS organisms x motifs (the output of DBTBS_motifs function)
    save_plot: boolean indicating whether plot should be saved
    # outputs
     if save_plot, then a bar plot is saved
    """

    Bsubtilis_sigma_factor_occurances.astype(int).plot.barh(figsize=(10,7), color = '#6E8DAB')
    plt.xlabel('Count', size=24,fontname="Arial"),
    plt.ylabel('B. subtilis sigma factor', size=24,fontname="Arial")
    plt.title('B. subtilis promoter counts by sigma factor', size=24,fontname="Arial")
    plt.xticks(fontsize = 16) # work on current fig
    plt.yticks(fontsize = 16) # work on current fig
    plt.legend(loc='lower right', prop={'size': 16})
    plt.tight_layout()

    if save_plot:
        plt.savefig('../../figs/promoter_database_summaries/20191127_DBTBS_promoters_by_Bsubtilis_sigmaFactor.pdf')
        plt.savefig('../../figs/promoter_database_summaries/20191127_DBTBS_promoters_by_Bsubtilis_sigmaFactor.png', dpi=400)
        plt.savefig('../../figs/promoter_database_summaries/20191127_DBTBS_promoters_by_Bsubtilis_sigmaFactor.svg')

"""
inducible promoters from bacterial and fungal species
"""

def parse_bacterial_inducible_promoters(save_csv=True):

    """
    parse through bacterial inducible promoters obtained from Meyer et al 2019 and apply QC / filtering

    inputs
    save_csv: boolean argument indicating whether csv files should be saved

    outputs
    inducible_promoters_bact_df: dataframe containing bacterial inducible promoters and metadata,
            if save_csv, then csv file is saved to data directory

    """

    # load raw inducible promoter dataframe
    inducible_promoters_bact_df = pd.read_excel('../../data/inducible_promoters/20191202_inducible_promoters_bacteria.xlsx')

    inducible_promoters_bact_df = QC_DNA_sequences(inducible_promoters_bact_df)

    if save_csv:
        inducible_promoters_bact_df.to_csv("../../data/parsed_promoter_data/20191129_bacterial_inducible_promoters.csv",index=False)


    return(inducible_promoters_bact_df)

def parse_fungal_inducible_promoters(save_csv=False):

    """
    parse through fungal inducible promoters obtained from various sources and apply QC / filtering

    inputs
    save_csv: boolean argument indicating whether csv files should be saved

    outputs
    inducible_promoters_fungus_df: dataframe containing DBTBS promoters and metadata,
            if save_csv, then csv file is saved to data directory

    """

    # load raw inducible promoter dataframe
    inducible_promoters_fungus_df = pd.read_excel('../../data/inducible_promoters/20191210_inducible_promoters_fungus.xlsx')

    inducible_promoters_fungus_df = QC_DNA_sequences(inducible_promoters_fungus_df)

    if save_csv:
        inducible_promoters_fungus_df.to_csv("../../data/parsed_promoter_data/20200128_fungal_inducible_promoter_sequences.csv",index=False)

    return(inducible_promoters_fungus_df)

def inducible_promoters_count(inducible_promoters_bact_df_to_concat, inducible_promoters_fungus_df_to_concat, save_csv=False):

    """
    stratify inducible promoters by organism (e.g., fungal, bacterial, mammalian)

    inputs
    inducible_promoters_bact_df_to_concat: dataframe for bacterial inducible promoters
    inducible_promoters_fungus_df_to_concat: dataframe for fungal inducible promoters
    save_csv: boolean indicating whether csv should be saved

    outputs: dataframe containing inducible promoter counts
    """

    bact_inducible_promoters = list(inducible_promoters_bact_df_to_concat['promoter'])
    bact_inducible_promoter_varients_count = 0

    # counting promoter varients
    for induc_prom in bact_inducible_promoters:
        if 'varient' in induc_prom:
            bact_inducible_promoter_varients_count += 1

    fungus_inducible_promoter_varients_count = 0

    mammalian_inducible_promoter_varients_count = 0


    inducible_promoters_counts_df = pd.DataFrame(columns = ["Promoters total", "Original promoters", "Varients"])
    inducible_promoters_counts_df.loc['Mammalian', 'Promoters total'] = 0
    inducible_promoters_counts_df.loc['Bacteria', 'Promoters total'] = inducible_promoters_bact_df_to_concat.shape[0]
    inducible_promoters_counts_df.loc['Fungi', 'Promoters total'] = inducible_promoters_fungus_df_to_concat.shape[0]

    inducible_promoters_counts_df.loc['Mammalian', 'Original promoters'] = 0 - mammalian_inducible_promoter_varients_count
    inducible_promoters_counts_df.loc['Bacteria', 'Original promoters'] = inducible_promoters_bact_df_to_concat.shape[0] - bact_inducible_promoter_varients_count
    inducible_promoters_counts_df.loc['Fungi', 'Original promoters'] = inducible_promoters_fungus_df_to_concat.shape[0] - fungus_inducible_promoter_varients_count

    inducible_promoters_counts_df.loc['Mammalian', 'Varients'] = mammalian_inducible_promoter_varients_count
    inducible_promoters_counts_df.loc['Bacteria', 'Varients'] = bact_inducible_promoter_varients_count
    inducible_promoters_counts_df.loc['Fungi', 'Varients'] = fungus_inducible_promoter_varients_count

    inducible_promoters_counts_df = inducible_promoters_counts_df.sort_values(by="Promoters total")
    inducible_promoters_counts_df = inducible_promoters_counts_df.astype('int')

    if save_csv:
        inducible_promoters_counts_df.to_csv('../../data/parsed_promoter_data/20200128_inducible_promoter_counts.csv',index=False)
    return(inducible_promoters_counts_df)

def plot_inducible_promoters(inducible_promoters_counts_df, save_plot=True):

    """
    create bar plot showing DBTBS promoters stratified by cognate sigma factor

    # inputs
    inducible_promoters_counts_df: dataframe containing DBTBS organisms x motifs (the output of DBTBS_motifs function)
    save_plot: boolean indicating whether plot should be saved

    # outputs
     if save_plot, then a bar plot is saved
    """

    black = '#000000'
    dark_slate_gray = '#2F4F4F'
    plot_colors = [black, dark_slate_gray]


    inducible_promoters_counts_df.iloc[:,1:].plot.barh(stacked=True, color = plot_colors, figsize=(10,7))
    plt.xlabel('Count', size=24,fontname="Arial")
    plt.ylabel('Organism', size=24,fontname="Arial")
    plt.title('Inducible promoter counts',fontname="Arial", size=24)
    plt.xticks(fontsize = 16,fontname="Arial")
    plt.yticks(fontsize = 16,fontname="Arial")
    plt.legend(loc='lower right', prop={'size': 16})
    plt.tight_layout()

    if save_plot:
        plt.savefig('../../figs/promoter_database_summaries/20191210_inducible_promoter_counts.pdf')
        plt.savefig('../../figs/promoter_database_summaries/20191210_inducible_promoter_counts.png', dpi=400)
        plt.savefig('../../figs/promoter_database_summaries/20191210_inducible_promoter_counts.svg')


def main():
    parser = argparse.ArgumentParser()

    # run EPDnew functions
    EPDnew_promoters_df = parse_EPDnew(save_csv=True)
    motifs_occurances_df = EPDnew_motifs(EPDnew_promoters_df)
    plot_EPDnew_data(motifs_occurances_df, save_plot=True) # plot EPDnew data
    # edit EPDnew dataframe such that it is compatible for concatenation
    EPDnew_promoters_df_to_concat = reorganize_promoter_df_columns(EPDnew_promoters_df)

    # run RegulonDB functions
    regulonDB_promoters_df = parse_RegulonDB(save_csv=True)
    Ecoli_sigma_factor_occurances = RegulonDB_motifs(regulonDB_promoters_df, save_csv=True)
    plot_RegulonDB_data(Ecoli_sigma_factor_occurances, save_plot=True)
    # edit RegulongDB dataframe such that it is compatible for concatenation
    regulonDB_promoters_df_to_concat = reorganize_promoter_df_columns(regulonDB_promoters_df)

    # run DBTBS functions
    DBTBS_promoters_df = parse_DBTBS(save_csv=True) #Bsubtilis_sigma_factor_occurances
    Bsubtilis_sigma_factor_occurances = DBTBS_motifs(DBTBS_promoters_df, save_csv=True)
    plot_DBTBS_data(Bsubtilis_sigma_factor_occurances, save_plot=True)

    # edit DBTBS dataframe such that it is compatible for concatenation
    DBTBS_promoters_df_to_concat = reorganize_promoter_df_columns(DBTBS_promoters_df)

    # run bacterial inducible promoters functions
    inducible_promoters_bact_df = parse_bacterial_inducible_promoters(save_csv=True)
    inducible_promoters_bact_df_to_concat = reorganize_promoter_df_columns(inducible_promoters_bact_df)

    # run fungal inducible promoters functions
    inducible_promoters_fungus_df = parse_fungal_inducible_promoters(save_csv=True)
    inducible_promoters_fungus_df_to_concat = reorganize_promoter_df_columns(inducible_promoters_fungus_df)

    # run functions for inducible promoter counts and plot these counts
    inducible_promoters_counts_df = inducible_promoters_count(inducible_promoters_bact_df_to_concat, inducible_promoters_fungus_df_to_concat)
    plot_inducible_promoters(inducible_promoters_counts_df, save_plot=True)

    # concatenate promoters
    promoters_df = pd.concat([EPDnew_promoters_df_to_concat,
               regulonDB_promoters_df_to_concat,
               DBTBS_promoters_df_to_concat,
               inducible_promoters_bact_df_to_concat,
                inducible_promoters_fungus_df_to_concat],
                sort=False).reset_index().drop('index', axis=1)

    promoters_df.to_csv("../../data/parsed_promoter_data/20191203_promoters.csv",index=False)

main()
