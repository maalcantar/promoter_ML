#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Jan 30 11:59:41 2020
@author: maalcantar
"""

import os
import random
import argparse
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from promoter_df_adjustment_util import *
from promoter_parsing_util import *

random.seed(777)

def create_promoter_df_for_ml(positive_promoter_df,
                              negative_permuted_promoter_df,
                              random_df,
                              #percentage_permuted,
                              random_idx_start,organism, seed=777):

    """
    create dataframe cotaining promoter sequences needed for machine learning
        that is, this dataframe with contain: i) positive promoters ii) negative
        promomters (that were permuted from positive ones)
        iii) random negative promoters. For each organism, we will have an equal split of organism and motif

    inputs
    positive_promoter_df: dataframe containing positive promoter sequences
    negative_permuted_promoter_df: dataframe containing negative promoters (obtained by permuting positive promoters)
    random_df: random negative promoter sequences
    random_idx_start: where to start indexing negative random num_sequences_to_fill

    taken out: percentage_permuted: percent of permuted promoters that should constitute the negative subset

    outputs
    promoters_for_ML_df: dataframe containing promoters ready for ML (see function description for more detail)
    random_idx_start: integer indicating the index for the next available random sequence

    """
    random.seed(777)
    # specify the number of permuted promoter sequences we want
    # we will randomly sample this many permuted sequences from the given pool
    promoters_for_ML = pd.DataFrame()
    num_sequences_to_fill = list(range(0,positive_promoter_df.shape[0]))
    num_sequences = len(num_sequences_to_fill)

    # create a dataframe with the selected permuted promoters
    negative_set_permuted_df = negative_permuted_promoter_df.copy()#.iloc[permuted_promoters_to_pull_idx,:]
    random_set_df = random_df.copy().iloc[0:num_sequences,:]
    random_set_df['organism'] = [organism] * num_sequences
    # since we continually sample from the random promoters, we need to keep track
    # of how many we have used
    random_idx_start = random_idx_start + num_sequences#num_random_promoter

    # put together real promoters, permuted sequences, and random sequences
    promoters_for_ML_df = pd.concat([promoters_for_ML,
           positive_promoter_df,
           negative_set_permuted_df,
           random_set_df],
            sort=False).reset_index().drop('index', axis=1)
    return(promoters_for_ML_df, random_idx_start)

def concat_intragenic_region(prom_df, prom_for_ML, organism_name):

    """
    concatenate intragenic regions with positive and negative promoter sets for specified organism_motifs_dict

    inputs
    prom_df: original containing only positive promoter sequences
    prom_for_ML: dataset with positive and negative promoter num_sequences
    organism_name: character string indicating organism whose intragenic regions
    are being concatenated. Currently supported inputs are: human, athaliana, ecoli_intragenic_seq_df

    output
    prom_for_ML: dataframe with new intragenic regions concatenated to positive
    and negative promoter sequences
    """

    genome_transcrips_file = '../../data/parsed_genome_transcripts/' + organism_name + '_transcriptome_trimmed.csv'
    intragenic_seq_df = pd.read_csv(genome_transcrips_file)
    organism_name_in_df = list(intragenic_seq_df['organism'])[0]

    prom_seq_num = prom_df.copy()[(prom_df['organism'] == organism_name_in_df)].shape[0]
    num_intra_seqs = intragenic_seq_df.shape[0]
    random_indices = random.sample(range(0, num_intra_seqs), prom_seq_num)
    intragenic_seq_df = intragenic_seq_df.copy().iloc[random_indices,:]
    prom_for_ML = pd.concat([prom_for_ML,
       intragenic_seq_df],
        sort=False).reset_index().drop('index', axis=1)

    return(prom_for_ML)

def main():

    promoters_all_df = pd.read_csv('../../data/promoters_for_ML/promoters_all.csv', low_memory = False).fillna('')

    # separate EPDnew dataframe by real promoter and permuted promoter (for negative set)
    EPDnew_df = promoters_all_df.copy()[(promoters_all_df['database/source'] == 'EPDnew') |
                                      (promoters_all_df['database/source'] == 'EPDnew (permuted)')
                                      ]

    RegulonDB_df = promoters_all_df.copy()[(promoters_all_df['database/source'] == 'RegulonDB') |
                                      (promoters_all_df['database/source'] == 'RegulonDB (permuted)')
                                      ]

    # load in random sequences
    random_df = promoters_all_df.copy()[promoters_all_df['database/source'] == 'random']

    eukaryotic_organisms = set(EPDnew_df['organism'])
    eukaryotic_motifs = ['Inr-pfa', 'nonInr-pfa','TATA', 'nonTATA'] # set(EPDnew_df['sigma factor/motif'])

    EPDnew_promoters_for_ML_df = pd.DataFrame()
    random_set_counter = 0

    # for each organism / motif, create dataframe with positive and negative promoter set
    # we do this for each organism / motif to ensure sequences are equally split up
    for organism in eukaryotic_organisms:
        for motif in eukaryotic_motifs:
            EPDnew_promoter_positive_df = EPDnew_df.copy()[(EPDnew_df['organism'] == organism) &
                                                            (EPDnew_df['sigma factor/motif'] == motif) &
                                                            ((EPDnew_df['database/source'] == 'EPDnew'))]

            EPDnew_promoter_negative_permuted_df = EPDnew_df.copy()[(EPDnew_df['organism'] == organism) &
                                                                    (EPDnew_df['sigma factor/motif'] == motif) &
                                                            (       (EPDnew_df['database/source'] == 'EPDnew (permuted)'))]
            EPDnew_promoters_for_ML_temp_df, random_set_counter = create_promoter_df_for_ml(EPDnew_promoter_positive_df,
                                          EPDnew_promoter_negative_permuted_df,
                                          random_df,
                                          # percentage_permuted,
                                          random_set_counter,organism)


            EPDnew_promoters_for_ML_df = pd.concat([EPDnew_promoters_for_ML_df,
               EPDnew_promoters_for_ML_temp_df],
                sort=False).reset_index().drop('index', axis=1)
  # intragenic sequences
    # human_prom_seq_num = EPDnew_df.copy()[(EPDnew_df['organism'] == 'h_sapien')].shape[0]
    # human_intragenic_seq_df = pd.read_csv('../../data/parsed_genome_transcripts/human_transcriptome_trimmed.csv')
    # num_intra_seqs = human_intragenic_seq_df.shape[0]
    # random_indices = random.sample(range(0, num_intra_seqs), human_prom_seq_num)
    # human_intragenic_seq_df = human_intragenic_seq_df.copy().iloc[random_indices,:]
    # EPDnew_promoters_for_ML_df = pd.concat([EPDnew_promoters_for_ML_df,
    #    human_intragenic_seq_df],
    #     sort=False).reset_index().drop('index', axis=1)
    EPDnew_promoters_for_ML_df = concat_intragenic_region(EPDnew_df,
                                                          EPDnew_promoters_for_ML_df,
                                                          organism_name='human')
    EPDnew_promoters_for_ML_df = concat_intragenic_region(EPDnew_df,
                                                          EPDnew_promoters_for_ML_df,
                                                          organism_name='athaliana')

    # create positive / negative set promoter sataframe for RegulonDB sequences
    RegulonDB_promoters_for_ML_df = pd.DataFrame()

    RegulonDB_promoter_positive_df = RegulonDB_df.copy()[RegulonDB_df['database/source'] == 'RegulonDB']
    RegulonDB_promoter_negative_permuted_df = RegulonDB_df.copy()[RegulonDB_df['database/source'] == 'RegulonDB (permuted)']

    # trim random sequences to 80bp such that they are comptabile with RegulonDB
    num_RegulonSB_sequences = RegulonDB_promoter_positive_df.shape[0]
    random_df_RegulonDB = random_df.copy().iloc[random_set_counter:random_set_counter+num_RegulonSB_sequences,:]
    random_seq_RegulonDB_list = list(random_df_RegulonDB['DNA sequence'])
    random_seq_RegulonDB_trimmed = [random_seq[0:80] for random_seq in random_seq_RegulonDB_list]
    random_df_RegulonDB['DNA sequence'] = random_seq_RegulonDB_trimmed

    # intragenic regions
    # ecoli_intragenic_seq_df = pd.read_csv('../../data/parsed_genome_transcripts/ecoli_transcriptome_trimmed.csv')
    # num_intra_seqs = ecoli_intragenic_seq_df.shape[0]
    # random_indices = random.sample(range(0, num_intra_seqs), num_RegulonSB_sequences)
    # ecoli_intragenic_seq_df = ecoli_intragenic_seq_df.copy().iloc[random_indices,:]
    #
    # RegulonDB_promoters_for_ML_df, random_set_counter = create_promoter_df_for_ml(RegulonDB_promoter_positive_df,
    #                               RegulonDB_promoter_negative_permuted_df,
    #                               random_df_RegulonDB,
    #                               # percentage_permuted,
    #                               random_set_counter,'e_coli')

    RegulonDB_promoters_for_ML_df = concat_intragenic_region(RegulonDB_promoter_positive_df,
                                                              RegulonDB_promoters_for_ML_df,
                                                              organism_name='ecoli')

    # RegulonDB_promoters_for_ML_df = pd.concat([RegulonDB_promoters_for_ML_df,
    #     ecoli_intragenic_seq_df],
    #      sort=False).reset_index().drop('index', axis=1)

    EPDnew_promoters_for_ML_df.to_csv('../../data/promoters_for_ML/EPDnew_promoters_for_ML.csv', index=False)
    RegulonDB_promoters_for_ML_df.to_csv('../../data/promoters_for_ML/RegulonDB_promoters_for_ML.csv', index=False)

main()
