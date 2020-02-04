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

from promoter_parsing_util import *

random.seed(777)
def parse_human_transcriptome(save_csv = True):

    """
    create a dataframe containing nucleotide sequences from human human_transcriptome

    inputs:
    save_csv: boolean indicating whether csv files should be saved

    outputs:
    None
    """

    # define location of human transcriptome
    DIRECTORY_STR = '../../data/genome_transcripts/'
    filename = 'human_GRCh38_latest_rna.txt'
    filepath = DIRECTORY_STR + filename

    # information to parse
    sequences = []
    gene = []
    seq_length = []

    sequence = ''
    counter_gene = 0

    with open(filepath) as fp:
        line = fp.readline()

        while line:

            # The '>' signifies the start of a new gene
            if line[0] == '>':

                # append last completed sequence if we've already been through 1 iteration
                if counter_gene >= 1:
                    sequences.append(sequence)
                    seq_length.append(str(len(sequence)) + ' bp')

                sequence = ''
                counter_gene += 1
                gene_name = ' '.join(line.split(',')[0].split(' ')[3:])

                # extracting gene being encoded
                gene.append(gene_name)

            else: # keep appending nucleotides until we reach a new sequence
                sequence = sequence + line
            line = fp.readline()

    sequences.append(sequence)
    seq_length.append(str(len(sequence)) + ' bp')

    organism = ['h_sapien'] * len(sequences)
    database = ['human genes (NCBI)']* len(sequences)
    SF_motif = [''] * len(sequences)
    seq_class = [0] * len(sequences)

    # organizing data into a dictionary
    human_intragenic_sequences_data_dict = {
        'organism': organism,
        'database/source': database,
        'DNA sequence': sequences,
        'regulated gene': gene,
        'range (with respect to TSS)': seq_length,
        'sigma factor/motif': SF_motif
                        }

    # organizing data into dataframe
    human_transcriptome_df = pd.DataFrame.from_dict(human_intragenic_sequences_data_dict,orient='index').transpose()
    human_transcriptome_df = human_transcriptome_df[~human_transcriptome_df['regulated gene'].str.contains('uncharacterized')]
    human_transcriptome_df = QC_DNA_sequences(human_transcriptome_df)
    human_transcriptome_df = reorganize_promoter_df_columns(human_transcriptome_df)
    human_transcriptome_df = human_transcriptome_df[human_transcriptome_df['DNA sequence'].str.len() > 600]
    human_transcriptome_df['class'] = 0
    human_transcriptome_df = human_transcriptome_df.reset_index().drop('index', axis=1)

    trimmed_sequences = []
    for sequence in human_transcriptome_df['DNA sequence']:
        highest_start = len(sequence) - 600
        seq_start = random.randint(0,highest_start)
        trimmed_seq = sequence[seq_start:seq_start+600]
        trimmed_sequences.append(trimmed_seq)
    human_transcriptome_trimmed_df = human_transcriptome_df.copy().reset_index().drop('index', axis=1)
    human_transcriptome_trimmed_df['DNA sequence'] = trimmed_sequences
    human_transcriptome_trimmed_df['range (with respect to TSS)'] = ['600 bp'] * len(trimmed_sequences)

    if save_csv:
        human_transcriptome_df.to_csv('../../data/parsed_genome_transcripts/human_transcriptome.csv', index=False)
        human_transcriptome_trimmed_df.to_csv('../../data/parsed_genome_transcripts/human_transcriptome_trimmed.csv', index=False)

def main():
    parse_human_transcriptome(save_csv = True)
main()
