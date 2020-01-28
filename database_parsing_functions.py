#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 13:21:30 2020
@author: malcantar
"""

# importing packages required to run notebook
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# parse through EPDnew database
def parse_EPDnew(save_csv=True):

    #####------parsing through EPDnew database and acquring promoter sequences------#####

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
        counter_gene = 0 # counts number of genes included
        sequence = '' # initialize empty sequence || need to initialize here too so we don't get an error

        if filename[0] != '.': # ignore hidden files

        # loop through lines in a file
            with open(filepath) as fp:

                line = fp.readline()

                while line:

                    # The '>' signifies the start of a new gene
                    if line[0] == '>':

                        # append last completed sequence if we've already been through 1 iteration
                        if counter_gene >= 1:

                            promoter_sequences.append(sequence)
                        # initilizing empty sequence
                        sequence = ''
                        # increasing gene count
                        counter_gene += 1
                        # splitting line to extact important values
                        line_split = line.split(' ')

                        # extracting gene regulated by the current promoter
                        promoter_gene.append(line_split[1])

                        # edge case where gene does not have a specified range
                        try:
                            range_index = line.split(' ').index('range')
                            range_temp = line.split(' ')[range_index+2] + ' to ' + line.split(' ')[range_index+6].strip('.\n')
                            promoter_range.append(range_temp)
                        except ValueError:
                            promoter_range.append('NA')

                        # append organism name, database where promoter was obtained, and notes (e.g., TATA or nonTATA promoter)
                        promoter_organism.append(filename.split('_')[2] + '_' + filename.split('_')[3])
                        promoter_database.append(filename.split('_')[1])
                        promoter_notes.append(filename.split('_')[4].replace(' ',''))

                    else: # keep appending nucleotides until we reach a new sequence
                        sequence = sequence + line
                    line = fp.readline()

            # removing linebreaks from sequences
            promoter_sequences.append(sequence)
            promoter_sequences = [seq.strip('\n') for seq in promoter_sequences] # remove edge '\n'
            promoter_sequences = [seq.replace('\n','') for seq in promoter_sequences] # remove other '\n'

        # organizing data into a dictionary
        EPDnew_data_dict = {
            "organism": promoter_organism,
            "promoter database": promoter_database,
            "sequence":promoter_sequences,
            "regulated gene": promoter_gene,
            "range (with respect to TSS)": promoter_range,
            "motifs": promoter_notes
                            }

        # organizing data into dataframe
    EPDnew_promoters_df = pd.DataFrame.from_dict(EPDnew_data_dict,orient='index').transpose()

    #####------organize data by specific motifs and organisms------#####


    # initialize list the will eventually contain the information to be plotted
    organism_motifs_dict = dict()

    # organisms and motifs of interest
    organisms_list = list(set(EPDnew_promoters_df['organism']))
    motifs_list = list(set(EPDnew_promoters_df['motifs']))

    # calculate number of promoters that were obtained for each organism (and further stratify by motif)
    for motif in motifs_list:
        motif_occurances_list = []
        for organism in organisms_list:
            organism_motifs_count = EPDnew_promoters_df[(EPDnew_promoters_df['motifs'] == motif) & (EPDnew_promoters_df['organism'] == organism)].shape[0]
            motif_occurances_list.append(organism_motifs_count)
        organism_motifs_dict.update({motif: motif_occurances_list})

    # storing motifs in a form useful for plotting
    # creating dataframe of organism x motif
    motifs_occurances_df = pd.DataFrame.from_dict(organism_motifs_dict,orient='index').transpose()

    # convert h_sapien --> h. sapien (for all organisms)
    organisms_replace_underscore = [org.replace("_",". ") for org in organisms_list]
    reindex_motifs_occuranges_dict = dict(zip(list(motifs_occurances_df.index), organisms_replace_underscore))
    motifs_occurances_df = motifs_occurances_df.rename(index = reindex_motifs_occuranges_dict)

    # find total occurance occurance of promoters for each organism and then organize by motif
    motifs_occurances_df['total'] = motifs_occurances_df[list(motifs_occurances_df.columns)].sum(axis=1)
    motifs_occurances_df = motifs_occurances_df.sort_values(by="total")
    motifs_occurances_df = motifs_occurances_df[['nonTATA', 'TATA','Inr-pfa', 'nonInr-pfa', 'total']]

    # find total number of each kind of motif
    motifs_occurances_df.loc['total'] = np.nan
    motifs_occurances_df.loc['total', 'nonTATA'] = motifs_occurances_df.loc[:, 'nonTATA'].sum()
    motifs_occurances_df.loc['total', 'TATA'] = motifs_occurances_df.loc[:, 'TATA'].sum()
    motifs_occurances_df.loc['total', 'Inr-pfa'] = motifs_occurances_df.loc[:, 'Inr-pfa'].sum()
    motifs_occurances_df.loc['total', 'nonInr-pfa'] = motifs_occurances_df.loc[:, 'nonInr-pfa'].sum()
    motifs_occurances_df.loc['total', 'total'] = motifs_occurances_df.loc[:, 'total'].sum()

    # convert dataframe datatype to integer
    motifs_occurances_df = motifs_occurances_df.astype('int')

    # save as csv and feather
    if save_csv:
        EPDnew_promoters_df.to_csv('../../data/parsed_promoter_data/20191125_EPDnew_promoters.csv')
        motifs_occurances_df.to_csv('../../data/parsed_promoter_data/20191125_EPDnew_motifs_occurances.csv')
#     if save_feather:
#         EPDnew_promoters_df.to_feather('../../data/parsed_promoter_data/20191125_EPDnew_promoters.feather')
#         motifs_occurances_df.reset_index().to_feather('../../data/parsed_promoter_data/20191125_EPDnew_motifs_occurances.feather')


    return(EPDnew_promoters_df, motifs_occurances_df)
