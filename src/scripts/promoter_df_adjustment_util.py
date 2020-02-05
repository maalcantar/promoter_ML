#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 17:21:34 2020
@author: maalcantar
"""

import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def trim_promter_seq(promoter_df, new_range, dataset_to_trim):

    """
    reduce size of the promoter sequences according to newly specified range

    inputs
    promoter_df: pandas dataframe
        dataframe containing promoter sequences
    new_range: 2-element integer list
        list containing new range for promoter sequences (e.g., [-249, 50])
    dataset_to_trim: character string
        string indicating which dataset to trim
        valid inputs are: 'EPDnew', 'RegulonDB', 'DBTBS', 'bacterial inducible promoters',
        bacterial inducible promoters

    returns
        promoter_df: pandas dataframe
            dataframe with trimmed promoter sequences
    """

    # split dataframe based on which dataset we need to trim
    if dataset_to_trim == 'EPDnew':
        promoter_df_to_trim = promoter_df.copy()[promoter_df['database/source'] == 'EPDnew']
        promoter_df_not_to_trim = promoter_df.copy()[promoter_df['database/source'] != 'EPDnew']

    elif dataset_to_trim == 'RegulonDB':
        promoter_df_to_trim = promoter_df.copy()[promoter_df['database/source'] == 'RegulonDB']
        promoter_df_not_to_trim = promoter_df.copy()[promoter_df['database/source'] != 'RegulonDB']

    elif dataset_to_trim == 'DBTBS':
        promoter_df_to_trim = promoter_df.copy()[promoter_df['database/source'] == 'DBTBS']
        promoter_df_not_to_trim = promoter_df.copy()[promoter_df['database/source'] != 'DBTBS']

    elif dataset_to_trim == 'bacterial inducible promoters':
        promoter_df_to_trim = promoter_df.copy()[promoter_df['database/source'] == 'Meyer et al 2019 (PMID: 30478458)']
        promoter_df_not_to_trim = promoter_df.copy()[promoter_df['database/source'] != 'Meyer et al 2019 (PMID: 30478458)']

    elif dataset_to_trim == 'fungal inducible promoters':
        promoter_df_to_trim = promoter_df.copy()[((promoter_df['database/source'] != 'EPDnew') &
                                                 (promoter_df['database/source'] != 'RegulonDB') &
                                                 (promoter_df['database/source'] != 'DBTBS') &
                                                (promoter_df['database/source'] != 'Meyer et al 2019 (PMID: 30478458)'))]


        promoter_df_not_to_trim = promoter_df.copy()[~((promoter_df['database/source'] != 'EPDnew') &
                                                 (promoter_df['database/source'] != 'RegulonDB') &
                                                 (promoter_df['database/source'] != 'DBTBS') &
                                                (promoter_df['database/source'] != 'Meyer et al 2019 (PMID: 30478458)'))]

    new_sequence_list = []
    new_range_list = []

    # change sequence length
    # this is done by setting the lower index at: (new lower bound) - (old lower bound)
    # and upper index at: (sequence length) -( (old upper bound) - (new upper bound) )
    for seq, seq_range in zip(list(promoter_df_to_trim['DNA sequence']), list(promoter_df_to_trim['range (with respect to TSS)'])):
        range_orig_splt = seq_range.split(' to ')
        lower_range_orig = int(range_orig_splt[0])
        upper_range_orig = int(range_orig_splt[1])

        if new_range[0] < lower_range_orig or new_range[1] > upper_range_orig:
            print('Requested new range is larger than original sequence length.')
            new_sequence_list.append(seq)
            new_range_list.append(seq_range)
        else:
            new_lower_idx = int(new_range[0]) - lower_range_orig
            new_upper_idx = upper_range_orig - int(new_range[1])
            seq_new = seq[new_lower_idx:len(seq) - new_upper_idx]
            new_sequence_list.append(seq_new)
            new_range_list.append(str(new_range[0]) + ' to ' + '+'+str(new_range[1]))

    promoter_df_to_trim['DNA sequence'] = new_sequence_list
    promoter_df_to_trim['range (with respect to TSS)'] = new_range_list

    promoter_df = pd.concat([promoter_df_to_trim,
           promoter_df_not_to_trim],
            sort=False).reset_index().drop('index', axis=1)

    return(promoter_df)
