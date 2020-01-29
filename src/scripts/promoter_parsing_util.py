#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 12:01:42 2020
@author: maalcantar
"""

import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

### HELPER FUNCTIONS ###

## HELPER FUNCTIONS FOR QC AND FILTERING ##

def calc_range(promoter_sequence_df):

    """
    calculate range of provided promoter sequences, with respect to the TSS.
    this is done based on the delimiter '/', which surrounds the TSS (+1).

    inputs
    promoter_sequence_df: dataframe containing promoter sequences. sequences must be under column
    named "range (with respect to TSS)".

    outputs
    promoter_sequence_df: dataframe containing promoter sequences with updated "range (with respect to TSS)".
    """

    # loop through dataframe and update range if there is a '/' marking the TSS
    for idx, sequence in promoter_sequence_df['DNA sequence'].iteritems():

        if '/' in sequence:
                seq_split = sequence.split('/')
                TSS_upstream_len = str(len(seq_split[0]) - 1)
                TSS_downstream_len = str(len(seq_split[2]) + 1)
                range_string = '-' + TSS_upstream_len + ' to ' + TSS_downstream_len
                promoter_sequence_df.loc[idx, 'range (with respect to TSS)'] = range_string

    return(promoter_sequence_df)

def QC_DNA_sequences(promoter_sequence_df):

    """
    perform quality control (QC) and filter dna promoter sequences. QC involves the following:
    i) remove unwanted char acters: line-break character '\n', spaces, and brackets
    ii) use the "calc_range" function to update missing ranges with a known TSS
        '/' characters denoting TSS are removed after this step
    iii) replace NaN elements with empty string
    iv) removing sequences that don't contain canonical nucleotides: a, t, c, g
    this function also prints out the number of sequences filtered

    inputs
    promoter_sequence_df: dataframe containing promoter sequences. sequences must be under column
    named "range (with respect to TSS)".

    outputs
    promoter_sequence_df: dataframe that has undergone QC and filtering
    """

    original_sequence_number = promoter_sequence_df.shape[0]

    if list(promoter_sequence_df['database/source'])[0] == 'EPDnew':
        db = 'EPDnew'
    elif list(promoter_sequence_df['database/source'])[0] == 'RegulonDB':
        db = 'RegulonDB'
    elif list(promoter_sequence_df['database/source'])[0] == 'DBTBS':
        db = 'DBTBS'
    elif list(promoter_sequence_df['database/source'])[0] == 'Meyer et al 2019 (PMID: 30478458)':
        db = 'bacterial inducible promoters'
    else:
        db = 'fungal inducible promoters'

    # remove unwanted characters from DNA sequences
    sequences = list(promoter_sequence_df['DNA sequence'])
    sequences = [str(seq).strip('\n') for seq in sequences] # remove edge '\n'
    sequences = [seq.replace('\n','') for seq in sequences] # remove other '\n'
    sequences = [seq.replace(' ','') for seq in sequences] # remove empty spaces
    sequences = [seq.lower() for seq in sequences] # convert sequences to lowercase
    sequences = [seq.replace('{', '').replace('}', '') for seq in sequences]

    promoter_sequence_df['DNA sequence'] = sequences
    # update range (with respect to TSS)
    promoter_sequence_df = calc_range(promoter_sequence_df)
    sequences = [seq.replace('/', '') for seq in sequences]

    # replace old sequences with QC sequences
    promoter_sequence_df['DNA sequence'] = sequences

    # filter rows that don't only contain nucleotide sequence
    promoter_sequence_df = promoter_sequence_df.fillna('').replace('NaN','').replace('ND','').replace('None','')
    promoter_sequence_df = promoter_sequence_df[(~promoter_sequence_df['DNA sequence'].str.contains('[^atcg]')) & (promoter_sequence_df['DNA sequence'].str.len() > 0)]
    promoter_sequence_df.reset_index()
    new_sequence_number = promoter_sequence_df.shape[0]

    # display number of filtered sequences
    print(str(original_sequence_number-new_sequence_number) + ' sequences filtered from ' + db + ' dataset')

    return(promoter_sequence_df)

def rename_promoter_df_columns(promoter_sequence_df, promoter_df_or_type):

    """
    standardize column names for dataframe containing promoter sequences and metadata.
    new columns are also added, according the selected database

    inputs
    promoter_sequence_df: dataframe containing promoter sequences and metedata
    promoter_df_or_type: string indicating promoter database or type.
        valid string inputs are: "RegulonDB", "DBTBS"
                string input is NOT case-sensitive

    outputs
    promoter_sequence_df: dataframe containing promoter sequences with updated column names
    """

    # standardized column names for RegulonDB dataframe
    newregulonDB_columns_list = ["RegulonDB ID",
                                "promoter",
                                "DNA strand",
                                "location of TSS on genome",
                                "sigma factor/motif",
                                "DNA sequence",
                                "evidence",
                                "evidence confidence"]

    # standardized column names for DBTBS dataframe
    newDBTBS_columns_list = ['operon',
                             'regulated gene',
                             'absolute position',
                             'range (with respect to TSS)',
                             'DNA sequence',
                             'experimental evidence',
                             'sigma factor/motif']

    # standardize column names and add new columns, as needed
    if promoter_df_or_type.lower() == 'regulondb':
        promoter_sequence_df_columns_dict = dict(zip(list(promoter_sequence_df.columns), newregulonDB_columns_list))
        promoter_sequence_df = promoter_sequence_df.rename(columns=promoter_sequence_df_columns_dict)
        promoter_sequence_df['database/source'] = "RegulonDB"
        promoter_sequence_df['organism'] = 'e_coli'
        promoter_sequence_df['range (with respect to TSS)'] = '-60 to +20'

    elif promoter_df_or_type.lower() == 'dbtbs':

        promoter_sequence_df_columns_dict = dict(zip(list(promoter_sequence_df.columns), newDBTBS_columns_list))
        promoter_sequence_df = promoter_sequence_df.rename(columns=promoter_sequence_df_columns_dict)
        promoter_sequence_df['organism'] = 'b. subtilis'
        promoter_sequence_df['database/source'] = 'DBTBS'

        # standardize range syntax for DBTBS (position before +1 is now 0, as opposed to -1)
        DBTBS_range_list = list(promoter_sequence_df['range (with respect to TSS)'])
        range_new_list = []
        for range_orig in DBTBS_range_list:
            range_orig_splt = str(range_orig).split(':')
            if len(range_orig_splt) > 1:
                if int(range_orig_splt[0]) < 0:
                    range_orig_lower_range = str(int(range_orig_splt[0]) + 1)
                else:
                    range_orig_lower_range = range_orig_splt[0]

                range_orig_upper_range = range_orig_splt[1]
                range_new = range_orig_lower_range + ' to ' + range_orig_upper_range
            else:
                range_new = range_orig
            range_new_list.append(range_new)
        # DBTBS_range_list_renamed = [str(seq_range).replace(':', ' to ') for seq_range in DBTBS_range_list]
        # promoter_sequence_df['range (with respect to TSS)'] = DBTBS_range_list_renamed

        promoter_sequence_df['range (with respect to TSS)'] = range_new_list
    else:
        print('Error! Please enter a valid promoter dataframe or type: \'RegulonDB\' or \'DBTBS\'')

    return(promoter_sequence_df)

def reorganize_promoter_df_columns(promoter_sequence_df):

    """
    reorganize promoter dataframe columns such that they are consistent between databases.
    this will involve dropping some unnecessary columns

    inputs
    promoter_sequence_df: dataframe containing promoter sequences
        this dataframe will ideally have been passed through the rename_promoter_df_columns function if needed

    outputs
    promoter_sequence_df: dataframe with reordered columns
    """

    columns_to_conserve = ['organism', 'database/source', 'DNA sequence', 'regulated gene',
                           'range (with respect to TSS)', 'motifs', 'sigma factor/motif', 'inducer/repressor', 'promoter']

    columns_to_add = list(set(columns_to_conserve) - set(promoter_sequence_df.columns))
    columns_to_drop = list(set(promoter_sequence_df.columns) - set(columns_to_conserve))

    # if new columns are added, populate with empty strings
    for col_to_add in columns_to_add:
        promoter_sequence_df[col_to_add] = ''

    promoter_sequence_df = promoter_sequence_df.drop(columns_to_drop, axis=1)
    promoter_sequence_df = promoter_sequence_df[columns_to_conserve]

    return(promoter_sequence_df)

## functions for finding number of motifs in dataset ##

def EPDnew_motifs(EPDnew_promoters_df, save_csv=False):

    """
    stratify EPDnew promoters by motif and organism

    inputs
    EPDnew_promoters_df: EPDnew promoter dataframe,
    save_csv: boolean indicating whether to save csv

    outputs
    motifs_occurances_df: dataframe containing matrix of organism x motif (nonTATA, TATA, Inr-pfa, nonInr-pfa, total)

    """

    organism_motifs_dict = dict()
    organisms_list = list(set(EPDnew_promoters_df['organism']))
    motifs_list = list(set(EPDnew_promoters_df['sigma factor/motif']))

    # calculate number of promoters that were obtained for each organism (and further stratify by motif)
    for motif in motifs_list:
        motif_occurances_list = []
        for organism in organisms_list:
            organism_motifs_count = EPDnew_promoters_df[(EPDnew_promoters_df['sigma factor/motif'] == motif) & (EPDnew_promoters_df['organism'] == organism)].shape[0]
            motif_occurances_list.append(organism_motifs_count)
        organism_motifs_dict.update({motif: motif_occurances_list})

    # creating dataframe of organism x motif (a form that is useful for plotting)
    motifs_occurances_df = pd.DataFrame.from_dict(organism_motifs_dict,orient='index').transpose()

    # replace underscors with periods for organism name
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

    if save_csv:
        motifs_occurances_df.to_csv('../../data/parsed_promoter_data/20191125_EPDnew_motifs_occurances.csv',index=False)

    return(motifs_occurances_df)

def RegulonDB_motifs(regulonDB_promoters, save_csv=False):

    """
    stratify RegulonDB promoters by sigma factor that binds promoter

    inputs
    regulonDB_promoters: RegulonDB promoter dataframe,
    save_csv: boolean indicating whether to save csv

    outputs
    Ecoli_sigma_factor_occurances: dataframe containing matrix of sigma factor x number of promoters
    """

    # define E. coli sigma factors
    sigma_factors_Ecoli = ['unknown',
                           'Sigma70',
                           'Sigma54',
                           'Sigma38',
                           'Sigma32',
                           'Sigma28',
                           'Sigma24',
                           'Sigma19']

    # dictionary that will contain promoter sigma factor frequencies
    Ecoli_sigma_factor_count_dict = dict()

    # determine counts for each sigma factor
    for sigma_factor in sigma_factors_Ecoli:

        if sigma_factor == 'unknown':
            SF_count = int(regulonDB_promoters[regulonDB_promoters['sigma factor/motif'].str.contains(sigma_factor)].shape[0]) + \
                        int(regulonDB_promoters[regulonDB_promoters['sigma factor/motif'] == ''].shape[0])
        else:
            SF_count = regulonDB_promoters[regulonDB_promoters['sigma factor/motif'].str.contains(sigma_factor)].shape[0]

        Ecoli_sigma_factor_count_dict.update({sigma_factor: SF_count})

    # build dataframe containing sigma factor count information
    Ecoli_sigma_factor_occurances = pd.DataFrame.from_dict(Ecoli_sigma_factor_count_dict,orient='index')
    Ecoli_sigma_factor_occurances.loc['total',:] = regulonDB_promoters.shape[0]
    Ecoli_sigma_factor_occurances = Ecoli_sigma_factor_occurances.rename(columns = {0:"E. coli"}).sort_values(by="E. coli")
    Ecoli_sigma_factor_occurances = Ecoli_sigma_factor_occurances.astype('int')

    if save_csv:
        Ecoli_sigma_factor_occurances.to_csv('../../data/parsed_promoter_data/20191202_RegulonDB_sigma_factor_occurances.csv', index=False)

    return(Ecoli_sigma_factor_occurances)

def DBTBS_motifs(DBTBS_promoters_df, save_csv=False):

    """
    stratify DBTBS promoters by sigma factor that binds promoter


    inputs
    DBTBS_promoters_df: DBTBS promoter dataframe,
    save_csv: boolean indicating whether to save csv

    outputs
    Bsubtilis_sigma_factor_occurances: dataframe containing matrix of sigma factor x number of promoters
    """

    sigma_factors_Bsubtilis = list(set(DBTBS_promoters_df['sigma factor/motif']))
    Bsubtilis_sigma_factor_count_dict = dict()

    # loop through sigma factors to find counts of each sigma factor type
    for sigma_factor in sigma_factors_Bsubtilis:
        SF_count = DBTBS_promoters_df[DBTBS_promoters_df['sigma factor/motif'].str.contains(sigma_factor[0:5])].shape[0]
        Bsubtilis_sigma_factor_count_dict.update({sigma_factor: SF_count})

    Bsubtilis_sigma_factor_occurances = pd.DataFrame.from_dict(Bsubtilis_sigma_factor_count_dict,orient='index')
    Bsubtilis_sigma_factor_occurances.loc["total",:] = DBTBS_promoters_df.shape[0]
    Bsubtilis_sigma_factor_occurances = Bsubtilis_sigma_factor_occurances.rename(columns = {0:"B. subtilis"}).sort_values(by="B. subtilis")
    Bsubtilis_sigma_factor_occurances = Bsubtilis_sigma_factor_occurances.astype('int')

    if save_csv:
        Bsubtilis_sigma_factor_occurances.to_csv("../../data/parsed_promoter_data/20191129_DBTBS_sigma_factor_occurances.csv",index=False)

    return(Bsubtilis_sigma_factor_occurances)
