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
    id_column_name_query = 'ID'
    smiles_column_name_query = 'Cleaned_SMILES'

    output_folder = 'tests/similarity_search_results'
    output_option = 'satisfied'

    ### Method: 'fingerprint' ###
    if method == 'fingerprint':
        similarity_search_multiple_query(input_file_library, input_file_query, method='fingerprint',
                                         similarity_cutoff=0.3, output_folder=output_folder, output_option=output_option)
    ### Method: 'mcs' ###
    elif method == 'mcs':
        similarity_search_multiple_query(input_file_library, input_file_query, method='mcs',
                                         mcs_match='exact', similarity_cutoff=0.3, output_folder=output_folder, output_option=output_option)
    ### Method: 'substructure' ###
    elif method == 'substructure':
        similarity_search_multiple_query(input_file_library, input_file_query, method='substructure',
                                         substructure_method='SMARTS',output_folder=output_folder, output_option=output_option)


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