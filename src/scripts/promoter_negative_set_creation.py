#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Jan 29 16:45:22 2020
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

def mutate_sequence(sequence, mutation_rate, seed = 777):

    """
    mutate nucleotides in a given sequence

    inputs
    sequence: string containing nucleotide sequence -- it is assumed that the sequence only contains 'atcg'
    mutation_rate: float indicating desired mutation rate (should be values between 0 and 1)
    """
    random.seed(seed)
    bases = 'actg'
    new_sequence = []
    for base in sequence:
        if random.random() < mutation_rate:
            new_base = bases.strip(base)[random.randint(0, 2)]
            new_sequence.append(new_base)
        else:
            new_sequence.append(base)
    return(''.join(new_sequence))


def create_permuted_set(promoter_df, number_of_splits, percentage_to_conserve, mutation_rate, seed=777):

    """
    create negative promoter set from existing promoter sequences
    future iterations may want to pull genomic sequences

    inputs
    promoter_df: dataframe containing promoter sequences
    number_of_splits: integer indicating how many fragments should be created from the main sequence
    percentage_to_conserve: float indicating what percentage of sequences should not be permuted

    outputs
    negative_promoter_set: dataframe containing negative promoter set

    """
    random.seed(seed)
    promoter_sequences = promoter_df['DNA sequence']
    negative_sequences = []

    for sequence in promoter_sequences:

        # define length of nucleotide segments
        seq_length = int(len(sequence))
        segements_length = int(seq_length / number_of_splits)
        segment_indices = list(range(0,number_of_splits,1))

        # segment sequence
        sequence_segments = [sequence[start_idx:start_idx+segements_length] for start_idx in range(0, seq_length, segements_length)]

        # define which positions should be conserved and which should be mutated
        num_positions_permute = round(number_of_splits*(1-percentage_to_conserve))
        num_positions_conserve = number_of_splits - num_positions_permute
        permuted_indices = random.sample(segment_indices,num_positions_permute)
        conserved_indices = list(set(segment_indices) - set(permuted_indices))
        conserved_indices.sort()
        permuted_indices.sort()
        conserved_segments = [sequence_segments[conserved_idx] for conserved_idx in conserved_indices]
        permuted_segments = [sequence_segments[permuted_idx] for permuted_idx in permuted_indices]

        # permute segments
        random.shuffle(permuted_segments)

        # mutate permuted segmenes
        permuted_segments = [mutate_sequence(segment, mutation_rate, seed) for segment in permuted_segments]

        # place segments back into sequence
        mutated_permuted_sequence = [None] * len(sequence_segments)
        # populate permuted sequences
        for idx_to_insert, permuted_segment in zip(permuted_indices,permuted_segments):
            mutated_permuted_sequence[idx_to_insert] = permuted_segment

        for idx_to_insert, conserved_segment in zip(conserved_indices,conserved_segments):
            mutated_permuted_sequence[idx_to_insert] = conserved_segment

        negative_sequences.append(''.join(mutated_permuted_sequence))

    # create dataframe with negative promoters
    negative_sequence_df = promoter_df.copy()
    negative_sequence_df['DNA sequence'] = negative_sequences
    negative_sequence_df['class'] = 0
    negative_sequence_df['database/source'] = list(negative_sequence_df['database/source'])[0] +' (permuted)'
    return(negative_sequence_df)


def generate_random_sequence(num_sequences, sequence_length, seed=777):

    """
    generates a specified number of random sequences

    inputs
    num_sequences: integer indicating the number of random sequences to be generated
    sequence_length: integer indicating length of sequences to be generated

    outputs
    sequences: list of random sequences
    """
    random.seed(seed)
    bases = ['a', 'c', 't', 'g']
    random_sequences = []

    # create extra sequences than desired in case there is any redundancy
    # (though that is extremely unlikely)
    for idx in range(0,round(num_sequences*1.2)):
        random_sequence = [random.choice(bases) for _ in range(sequence_length)]
        random_sequences.append(''.join(random_sequence))
    random_sequences_final = list(set(random_sequences[0:num_sequences]))
    random_sequences_df = pd.DataFrame()
    random_sequences_df['DNA sequence'] = random_sequences_final
    random_sequences_df['database/source'] = 'random'
    random_sequences_df = reorganize_promoter_df_columns(random_sequences_df)
    random_sequences_df['class'] = 0
    return(random_sequences_df)



def main():
    parser = argparse.ArgumentParser()

    promoter_df = pd.read_csv('../../data/parsed_promoter_data/20191203_promoters.csv',low_memory=False).fillna('')
    # promoter_df = trim_promter_seq(promoter_df, [-249,50], 'EPDnew')
    promoter_df = trim_promter_seq(promoter_df, [-59,20], 'RegulonDB')
    num_promoter_sequences = int(promoter_df.shape[0])
    EPDnew_promoters_df = promoter_df.copy()
    EPDnew_promoters_df = EPDnew_promoters_df[EPDnew_promoters_df['database/source'] == 'EPDnew']
    EPDnew_negative_promoters_df = create_permuted_set(EPDnew_promoters_df, 20, 0.4, 0.2, seed=777)

    RegulonDB_promoters_df = promoter_df.copy()
    RegulonDB_promoters_df = RegulonDB_promoters_df[RegulonDB_promoters_df['database/source'] == 'RegulonDB']
    RegulonDB_negative_promoters_df = create_permuted_set(RegulonDB_promoters_df, 8, 0.4, 0.2, seed=777)

    random_sequences_df = generate_random_sequence(num_promoter_sequences,600)

    promoters_all_df = pd.concat([promoter_df,
               EPDnew_negative_promoters_df,
               RegulonDB_negative_promoters_df,
               random_sequences_df],
                sort=False).reset_index().drop('index', axis=1)

    promoters_all_df.to_csv('../../data/promoters_for_ML/promoters_all.csv')

main()
