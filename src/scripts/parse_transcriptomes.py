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

def parse_ecoli_transcriptome(save_csv = True):

    """
    create a dataframe containing nucleotide sequences from e coli gene sequences

    inputs:
    save_csv: boolean indicating whether csv files should be saved

    outputs:
    None
    """
    # Download regulon DB promoter database and organize dataframe
    ecoli_gene_sequences_df = pd.read_csv('../../data/genome_transcripts/e_coli_Gene_sequence.csv', header=40)
    new_columns_list = ["RegulonDB ID",
                    "regulated gene",
                    "gene left end position in the genome",
                    "gene right end position in the genome",
                    "DNA strand where the gene is coded",
                    "product type",
                    "product name",
                    "start codon sequence",
                    "stop codon sequence",
                    "DNA sequence",
                    "all bnumber related to gene",
                    "other databases id  related to gene",
                    "other databases id  related to gene (2)"]
    old_columns_list = list(ecoli_gene_sequences_df.columns)
    reindex_dict = dict(zip(old_columns_list, new_columns_list))
    ecoli_gene_sequences_df = ecoli_gene_sequences_df.rename(columns=reindex_dict)
    ecoli_gene_sequences_df = ecoli_gene_sequences_df.fillna('')
    ecoli_gene_sequences_df = ecoli_gene_sequences_df[~ecoli_gene_sequences_df['product name'].str.contains('putative')]
    ecoli_gene_sequences_df = ecoli_gene_sequences_df[["DNA sequence", 'regulated gene']]
    ecoli_gene_sequences_df['database/source'] = ['e_coli genes (RegulonDB)'] * ecoli_gene_sequences_df.shape[0]
    ecoli_gene_sequences_df['organism'] = 'e_coli'
    ecoli_gene_sequences_df = reorganize_promoter_df_columns(ecoli_gene_sequences_df)

    ecoli_gene_sequences_df = QC_DNA_sequences(ecoli_gene_sequences_df)
    ecoli_gene_sequences_df = ecoli_gene_sequences_df[ecoli_gene_sequences_df['DNA sequence'].str.len() > 160]
    ecoli_gene_sequences_df = ecoli_gene_sequences_df.reset_index()
    ecoli_gene_sequences_df['class'] = 0

    final_seq_len = 80
    new_trimmed_seqs = []
    origin_gene = []
    for seq, gene in zip(ecoli_gene_sequences_df['DNA sequence'],ecoli_gene_sequences_df['regulated gene']) :
        num_iterations = int(np.floor(len(seq)/final_seq_len))
        for num_iter in range(0,num_iterations):
            new_trimmed_seqs.append(seq[final_seq_len*num_iter:(num_iter+1)*final_seq_len])
            origin_gene.append(gene)
    trimmed_ecoli_dict = {
                    'DNA sequence': new_trimmed_seqs,
                    'regulated gene': origin_gene,
                    'database/source': ['e_coli genes (RegulonDB)']*len(origin_gene)
                    }
    ecoli_gene_sequences_trimmed_df = pd.DataFrame(trimmed_ecoli_dict)
    ecoli_gene_sequences_trimmed_df['organism'] = "e_coli"
    ecoli_gene_sequences_trimmed_df = reorganize_promoter_df_columns(ecoli_gene_sequences_trimmed_df)
    ecoli_gene_sequences_trimmed_df = QC_DNA_sequences(ecoli_gene_sequences_trimmed_df)
    ecoli_gene_sequences_trimmed_df['class'] = 0

    if save_csv:
        ecoli_gene_sequences_df.to_csv('../../data/parsed_genome_transcripts/ecoli_transcriptome.csv', index=False)
        ecoli_gene_sequences_trimmed_df.to_csv('../../data/parsed_genome_transcripts/ecoli_transcriptome_trimmed.csv', index=False)
def parse_athaliana_transcriptome(save_csv=True):

    """
    create a dataframe containing nucleotide sequences from plant cDNA

    inputs:
    save_csv: boolean indicating whether csv files should be saved

    outputs:
    None
    """

    # define location of human transcriptome
    DIRECTORY_STR = '../../data/genome_transcripts/'
    filename = 'a_thaliana.txt'
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
                gene_name = line.split('|')[-1].split(',')[0][0:-5]

                # extracting gene being encoded
                gene.append(gene_name)

            else: # keep appending nucleotides until we reach a new sequence
                sequence = sequence + line
            line = fp.readline()

    sequences.append(sequence)
    seq_length.append(str(len(sequence)) + ' bp')

    organism = ['a_thaliana'] * len(sequences)
    database = ['a thaliana genes (plantgDB)']* len(sequences)
    SF_motif = [''] * len(sequences)
    seq_class = [0] * len(sequences)

    # organizing data into a dictionary
    plant_intragenic_sequences_data_dict = {
        'organism': organism,
        'database/source': database,
        'DNA sequence': sequences,
        'regulated gene': gene,
        'range (with respect to TSS)': seq_length,
        'sigma factor/motif': SF_motif
                        }

    # organizing data into dataframe
    plant_transcriptome_df = pd.DataFrame.from_dict(plant_intragenic_sequences_data_dict,orient='index').transpose()
    plant_transcriptome_df = plant_transcriptome_df[~plant_transcriptome_df['regulated gene'].str.contains('uncharacterized')]
    plant_transcriptome_df = QC_DNA_sequences(plant_transcriptome_df)
    plant_transcriptome_df = reorganize_promoter_df_columns(plant_transcriptome_df)
    plant_transcriptome_df = plant_transcriptome_df[plant_transcriptome_df['DNA sequence'].str.len() > 600]
    plant_transcriptome_df['class'] = 0
    plant_transcriptome_df = plant_transcriptome_df.reset_index().drop('index', axis=1)

    # trim sequences to match required length
    trimmed_sequences = []
    for sequence in plant_transcriptome_df['DNA sequence']:
        highest_start = len(sequence) - 600
        seq_start = random.randint(0,highest_start)
        trimmed_seq = sequence[seq_start:seq_start+600]
        trimmed_sequences.append(trimmed_seq)
    plant_transcriptome_trimmed_df = plant_transcriptome_df.copy().reset_index().drop('index', axis=1)
    plant_transcriptome_trimmed_df['DNA sequence'] = trimmed_sequences
    plant_transcriptome_trimmed_df['range (with respect to TSS)'] = ['600 bp'] * len(trimmed_sequences)

    if save_csv:
        plant_transcriptome_df.to_csv('../../data/parsed_genome_transcripts/athaliana_transcriptome.csv', index=False)
        plant_transcriptome_trimmed_df.to_csv('../../data/parsed_genome_transcripts/athaliana_transcriptome_trimmed.csv', index=False)


def main():
    parse_human_transcriptome(save_csv = True)
    parse_ecoli_transcriptome(save_csv = True)
    parse_athaliana_transcriptome(save_csv=True)
main()
