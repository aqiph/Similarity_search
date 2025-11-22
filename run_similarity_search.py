#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 11:14:00 2023

@author: guohan
"""

import sys, time, os
path_list = sys.path
module_path = '/Users/guohan/Documents/Codes/Similarity_search'
if module_path not in sys.path:
    sys.path.append(module_path)
    print('Add module path')

from similarity_search import similarity_search_multiple_query, select_analogs, plot_distribution


def run_similarity_search_multiple_query(method = 'fingerprint'):
    """
    Perform similarity search for multiple query SMILES in input_file_query
    """
    input_file_library = 'tests/library.csv'
    # input_file_library = '/Users/guohan/Documents/Projects/Datasets/HTS/Combination/forGNN/HTS_forGNN_446663.csv'
    input_file_query = 'tests/query_cmps.csv'
    kwargs = {'id_column_name_query': 'ID',
              'smiles_column_name_query': 'SMILES',
              'output_folder': 'similarity_search_results_FP_MolGen',
              'output_option': 'satisfied'}

    ### Method: 'fingerprint' ###
    if method == 'fingerprint':
        kwargs['similarity_cutoff'] = 0.3
        similarity_search_multiple_query(input_file_library, input_file_query, method='fingerprint', **kwargs)

    ### Method: 'mcs' ###
    elif method == 'mcs':
        kwargs['mcs_match'] = 'exact'
        kwargs['similarity_cutoff'] = 0.3
        similarity_search_multiple_query(input_file_library, input_file_query, method='mcs', **kwargs)

    ### Method: 'substructure' ###
    elif method == 'substructure':
        kwargs['substructure_method'] = 'SMARTS'
        similarity_search_multiple_query(input_file_library, input_file_query, method='substructure', **kwargs)


def run_select_analogs():
    """
    Select analogs
    """
    input_file_query = 'tests/query_cmps.csv'
    analogs_dir = 'tests/similarity_search_results'
    ### Select analogs above similarity cutoff
    # analog_method = 'cutoff'
    # similarity_cutoff = 0.4
    # select_analogs(input_file_query, analogs_dir, analog_method, similarity_cutoff=similarity_cutoff, deduplication=False)
    ### Select the rank-n most similar analog
    # analog_method = 'rank'
    # similarity_rank = 1
    # select_analogs(input_file_query, analogs_dir, analog_method, similarity_rank=similarity_rank, deduplication=False)
    ### Select top-n most similar analog
    analog_method = 'topN'
    similarity_topN = 3
    select_analogs(input_file_query, analogs_dir, analog_method, similarity_topN=similarity_topN, deduplication=False)


def run_plot_distribution():
    input_file = 'tests/query_cmps_Top1_12.csv'
    plot_distribution(input_file)



if __name__ == '__main__':
    # run_similarity_search_multiple_query()
    run_select_analogs()
    # run_plot_distribution()