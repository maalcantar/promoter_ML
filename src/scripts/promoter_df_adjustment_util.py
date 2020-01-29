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
    promoter_df: dataframe containing promoter sequences
    new_range: list containing new range for promoter sequences (e.g., [-249, 50])

    outputs:
    promoter_df: dataframe with updated promoter sequences
    """
