#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 16:10:00 2023

@author: guohan

1. Similarity search
2. Select analogs based on similarity cutoff or ranking
3. Plot similarity score distribution

"""


import os, time
from tqdm import tqdm
import pandas as pd
import numpy as np
from rdkit import Chem

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_size(10)

from utils.molecular_description import get_fingerprint, cal_fingerprint_distance, cal_MCS
from utils.tools import remove_unnamed_columns


### Similarity search ###
def is_similar_by_fingerprint(fp_1, fp_2, similarity_cutoff=0.7):
    """
    Determine whether fp_1 is similar to fp_2
    :param fp_1: fingerprint for molecule 1
    :param fp_2: fingerprint Mol object for molecule 2
    :param similarity_cutoff: float, if similarity score >= similarity_cutoff, smiles_1 is considered to be similar to smiles_2
    :return: bool, whether smiles_1 is considered to be similar to smiles_2
    """
    dist = cal_fingerprint_distance(fp_1, fp_2, fp_method='ecfp4')
    similarity_score = 1.0 - dist
    return similarity_score >= similarity_cutoff, similarity_score


def is_similar_by_mcs(mol_1, mol_2, match='exact', similarity_cutoff=0.7):
    """
    Determine whether smiles_1 is similar to smiles_2 based on MCS
    :param mol_1: RDKit Mol object for molecule 1
    :param mol_2: RDKit Mol object for molecule 2
    :param match: str, specify 'match' to use different comparison functions, allowed values include 'exact', 'anyAtom'
    :param similarity_cutoff: float, if similarity score >= similarity_cutoff, smiles_1 is considered to be similar to smiles_2
    :return: bool, whether smiles_1 is considered to be similar to smiles_2
    """
    similarity_score = cal_MCS(mol_1, mol_2, match)
    return similarity_score >= similarity_cutoff, similarity_score


def is_similar_by_substructure(mol, mol_sub, substructure_method='SMILES'):
    """
    Determine whether smiles_sub is a substructure of smiles
    :param mol: RDKit Mol object
    :param mol_sub: RDKit Mol object for substructure
    :param substructure_method: str, method used to convert smiles_sub to mol_sub, allowed values include 'SMILES', 'SMARTS'
    :return: bool, whether smiles_sub is a substructure of smiles
    """
    # # check SMILES
    # if not smiles:
    #     print(f"Error: Invalid SMILES {smiles}")
    #     return False
    # try:
    #     mol = Chem.MolFromSmiles(smiles)
    # except Exception:
    #     print(f"Error: Invalid SMILES {smiles}")
    #     return False
    # if not mol:
    #     print(f"Error: Invalid SMILES {smiles}")
    #     return False

    # determine whether smiles_sub is a substructure of smiles
    if substructure_method == 'SMILES':
        pass
    elif substructure_method == 'SMARTS':
        smiles_sub = Chem.MolToSmiles(mol_sub)   # Change molecule before get mol as a pattern, not using the Kekule form
        mol_sub = Chem.MolFromSmarts(smiles_sub)
    else:
        raise Exception('Error: Invalid substructure method.')
    hasSubstructure = mol.HasSubstructMatch(mol_sub)

    # compute MCS similarity score for compounds containing the given substructure
    if hasSubstructure:
        similarity_score = cal_MCS(mol, mol_sub, match='exact')
    else:
        similarity_score = 0.0

    return hasSubstructure, similarity_score


def similarity_search_multiple_query(input_file_library, input_file_query, method, **kwargs):
    """
    Find SMILES in input_file_library that are similar to the SMILES in input_file_query.
    Allowed methods include fingerprint, Maximum Common Structure (mcs), and substructure.
    :param input_file_library: str, path of the input library.
    :param input_file_query: str, path of the input query SMILES.
    :param method: str, method for similarity search, allowed values include 'fingerprint', 'mcs' and 'substructure'.
    The method 'substructure' will select compounds that contain the given substructure.
    :param smiles_column_name_library: str, name of the SMILES column in input_file_library.
    :param id_column_name_query: str, name of the ID column in input_file_query.
    :param smiles_column_name_query: str, name of the SMILES column in input_file_query.
    :param similarity_cutoff: float, similarity cutoff for 'fingerprint' and 'mcs' method.
    :param mcs_match: str, specify 'mcs_match' to use different comparison functions for MCS, allowed values include 'exact', 'anyAtom'.
    :param substructure_method: str, method used to convert smiles_query to mol_query, allowed values include 'SMILES', 'SMARTS'.
    :param output_folder: str, path of the output file.
    :param output_file: str, name of the output file.
    :param output_option: str, options to output data, allowed values include 'satisfied', 'not_satisfied' and 'all'.
    """
    folder = kwargs.get('output_folder', os.getcwd())
    if not os.path.exists(folder):
        os.makedirs(folder)
        print(f'Folder {folder} created.')

    # read library SMILES
    smiles_column_name_lib = kwargs.get('smiles_column_name_library', 'Cleaned_SMILES')
    df_lib = pd.read_csv(input_file_library)
    COLUMNS = df_lib.columns.tolist()
    SMILES_lib = df_lib[smiles_column_name_lib].values.tolist()
    Mols_lib = [Chem.MolFromSmiles(smiles) for smiles in SMILES_lib]

    # read query SMILES
    id_column_name_query = kwargs.get('id_column_name_query', 'ID')
    smiles_column_name_query = kwargs.get('smiles_column_name_query', 'Cleaned_SMILES')
    df_query = pd.read_csv(input_file_query)
    IDs_query = [str(id) for id in df_query[id_column_name_query].values.tolist()]
    SMILES_query = df_query[smiles_column_name_query].values.tolist()
    Mols_query = [Chem.MolFromSmiles(smiles) for smiles in SMILES_query]

    # fingerprint for similarity
    FP_lib, FP_query = [], []
    if method == 'fingerprint':
        FP_lib = [get_fingerprint(mol, fp_method='ecfp4') for mol in Mols_lib]
        FP_query = [get_fingerprint(mol, fp_method = 'ecfp4') for mol in Mols_query]

    # run similarity search for each query compounds
    total_time = []
    for i in tqdm(range(len(IDs_query)), desc='Processing'):
        time.sleep(0.0001)
        start_time = time.time()
        id_query = IDs_query[i]
        print(f'Start: {id_query}')

        is_similar, Similarity_Score = [], []
        if method == 'fingerprint':
            similarity_cutoff = kwargs.get('similarity_cutoff', 0.7)
            fp_query = FP_query[i]
            for fp_lib in FP_lib:
                is_sim, sim_score = is_similar_by_fingerprint(fp_lib, fp_query, similarity_cutoff)
                is_similar.append(is_sim)
                Similarity_Score.append(sim_score)
        elif method == 'mcs':
            match = kwargs.get('mcs_match', 'exact')
            similarity_cutoff = kwargs.get('similarity_cutoff', 0.7)
            mol_query = Mols_query[i]
            for mol_lib in Mols_lib:
                is_sim, sim_score = is_similar_by_mcs(mol_lib, mol_query, match, similarity_cutoff)
                is_similar.append(is_sim)
                Similarity_Score.append(sim_score)
        elif method == 'substructure':
            substructure_method = kwargs.get('substructure_method', 'SMILES')
            mol_query = Mols_query[i]
            for mol_lib in Mols_lib:
                is_sim, sim_score = is_similar_by_substructure(mol_lib, mol_query, substructure_method)
                is_similar.append(is_sim)
                Similarity_Score.append(sim_score)
        else:
            raise Exception('Error: Invalid similarity search method.')
        df_sim = pd.DataFrame(df_lib)
        df_sim['is_similar'] = is_similar
        df_sim['Similarity_Score'] = Similarity_Score

        # prepare output file
        output_option = kwargs.get('output_option', 'satisfied')
        if output_option == 'satisfied':
            df_sim = pd.DataFrame(df_sim[df_sim['is_similar']], columns=COLUMNS+['Similarity_Score'])
        elif output_option == 'not_satisfied':
            df_sim = pd.DataFrame(df_sim[~df_sim['is_similar']], columns=COLUMNS+['Similarity_Score'])
        elif output_option == 'all':
            df_sim = pd.DataFrame(df_sim, columns=COLUMNS + ['Similarity_Score', 'is_similar'])
        else:
            raise Exception('Error: Invalid output option.')

        # write to file
        if method in {'fingerprint', 'mcs', 'substructure'}:
            df_sim.sort_values(by=['Similarity_Score'], ascending=False, inplace=True)
        df_sim = df_sim.reset_index(drop=True)
        print('Number of rows:', df_sim.shape[0])
        df_sim = remove_unnamed_columns(df_sim)
        output_file = os.path.join(folder, id_query)
        df_sim.to_csv(output_file + f'_{df_sim.shape[0]}.csv')

        end_time = time.time()
        elapsed_time = end_time - start_time
        total_time.append(elapsed_time)
        print(f'Done, time = {elapsed_time:.2f}')

    print(f'Mean time is {np.mean(np.array(total_time)):.4f}, std is {np.std(np.array(total_time)):.4f}')


### Select analogs ###
def select_analogs(input_file_query, analogs_dir, analog_method, **kwargs):
    """
    Select analogs based on three criteria: 1. similarity scores are above the given similarity score cutoff;
    2. the rank n most similar analogs; or 3. top n most similar analogs.
    :param input_file_query: str, path of the input query SMILES.
    :param analogs_dir: str, path of the directory containing analogs.
    :param analog_method: str, method for selecting analogs, allowed values include 'cutoff', 'rank' and 'topN'.
    :param similarity_cutoff: float, similarity cutoff for selecting analogs.
    :param similarity_rank: int, rank of the selected analog.
    :param similarity_topN: int, number of top selected analogs.
    :param deduplication: bool, whether to deduplicate the selected analogs
    """
    # output file
    output_file = os.path.splitext(os.path.abspath(input_file_query))[0]

    # get query IDs and SMILES
    df_query = pd.read_csv(input_file_query)
    IDs_query = df_query['ID'].tolist()
    SMILES_query = df_query['SMILES'].tolist()

    # get analogs
    if analog_method == 'cutoff':
        similarity_cutoff = kwargs.get('similarity_cutoff', 0.5)
        print(f'Get analogs whose similarity scores are above {similarity_cutoff}.')
        output_file = f'{output_file}_cutoff{similarity_cutoff}'

        dfs = []
        files = os.listdir(analogs_dir)
        for i, id in enumerate(IDs_query):
            analogs_file = [file for file in files if file.startswith(f'{id}_')][0]
            df_analogs = pd.read_csv(f'{analogs_dir}/{analogs_file}')
            df_analogs['ID'] = df_analogs['ID'].apply(lambda analog_id: str(analog_id))

            df_selected_analogs = df_analogs[df_analogs['Similarity_Score'] >= similarity_cutoff]

            if df_selected_analogs.empty:
                continue
            df_selected_analogs = pd.DataFrame(df_selected_analogs, columns=['ID', 'SMILES', 'Similarity_Score'])
            df_selected_analogs.rename(columns={'ID': 'Analog_ID', 'SMILES': 'Analog_SMILES'}, inplace=True)
            df_selected_analogs['Query_ID'] = str(id)
            df_selected_analogs['Query_SMILES'] = SMILES_query[i]
            dfs.append(df_selected_analogs)

    elif analog_method == 'rank':
        similarity_rank = kwargs.get('similarity_rank', 1)
        assert similarity_rank >= 1, 'Error: Invalid rank.'
        print(f'Get the rank {similarity_rank} most similar analog.')
        output_file = f'{output_file}_rank{similarity_rank}'

        dfs = []
        files = os.listdir(analogs_dir)
        for i, id in enumerate(IDs_query):
            analogs_file = [file for file in files if file.startswith(f'{id}_')][0]
            df_analogs = pd.read_csv(f'{analogs_dir}/{analogs_file}')
            df_analogs['ID'] = df_analogs['ID'].apply(lambda analog_id: str(analog_id))

            df_selected_analogs = pd.DataFrame()
            if df_analogs.shape[0] >= similarity_rank:
                df_analogs = df_analogs.sort_values(by=['Similarity_Score'], ascending=[False], ignore_index=True)
                df_selected_analogs = df_analogs.loc[[similarity_rank - 1]]

            if df_selected_analogs.empty:
                continue
            df_selected_analogs = pd.DataFrame(df_selected_analogs, columns=['ID', 'SMILES', 'Similarity_Score'])
            df_selected_analogs.rename(columns={'ID': 'Analog_ID', 'SMILES': 'Analog_SMILES'}, inplace=True)
            df_selected_analogs['Query_ID'] = str(id)
            df_selected_analogs['Query_SMILES'] = SMILES_query[i]
            dfs.append(df_selected_analogs)

    elif analog_method == 'topN':
        similarity_topN = kwargs.get('similarity_topN', 1)
        assert similarity_topN >= 1, 'Error: Invalid N value.'
        print(f'Get top {similarity_topN} most similar analogs.')
        output_file = f'{output_file}_top{similarity_topN}'

        dfs = []
        files = os.listdir(analogs_dir)
        for i, id in enumerate(IDs_query):
            analogs_file = [file for file in files if file.startswith(f'{id}_')][0]
            df_analogs = pd.read_csv(f'{analogs_dir}/{analogs_file}')
            df_analogs['ID'] = df_analogs['ID'].apply(lambda analog_id: str(analog_id))

            num_selected_analogs = min(similarity_topN, df_analogs.shape[0])
            df_analogs = df_analogs.sort_values(by=['Similarity_Score'], ascending=[False], ignore_index=True)
            df_selected_analogs = df_analogs.loc[0:(num_selected_analogs - 1)]

            if df_selected_analogs.empty:
                continue
            df_selected_analogs = pd.DataFrame(df_selected_analogs, columns=['ID', 'SMILES', 'Similarity_Score'])
            df_selected_analogs.rename(columns={'ID': 'Analog_ID', 'SMILES': 'Analog_SMILES'}, inplace=True)
            df_selected_analogs['Query_ID'] = str(id)
            df_selected_analogs['Query_SMILES'] = SMILES_query[i]
            dfs.append(df_selected_analogs)

    else:
        raise Exception('Error: Invalid method for selecting analogs.')

    # concat results
    if len(dfs) >= 1:
        df_concat = pd.concat(dfs, ignore_index=True, sort=False)
    else:
        print(f'No analog found!')
        return

    # remove duplicates
    deduplication = kwargs.get('deduplication', False)
    if deduplication:
        df_concat = df_concat.sort_values(by=['Similarity_Score'], ascending=[False], ignore_index=True)
        df_concat = df_concat.drop_duplicates('Analog_ID', keep='first', ignore_index=True)

    # write output file
    df_concat = df_concat.reset_index(drop=True)
    print('Number of selected analogs:', df_concat.shape[0])
    df_concat = remove_unnamed_columns(df_concat)
    df_concat.to_csv(f'{output_file}_{df_concat.shape[0]}.csv')


### Plot similarity score distribution ###
def plot_distribution(input_file):
    """
    Plot similarity score distribution
    :param input_file: str, path of the input file
    """
    # files
    output_file = os.path.splitext(os.path.abspath(input_file))[0]
    df = pd.read_csv(input_file)

    # plot distribution
    values = df['Similarity_Score'].tolist()
    output_file = f'{output_file}_similarity_distribution.pdf'

    plt.figure(1)
    plt.hist(values, 10, range=(0.0, 1.0))
    plt.xlabel('Similarity Score', fontproperties=font)
    plt.ylabel('Counts', fontproperties=font)
    plt.xticks(fontproperties=font)
    plt.yticks(fontproperties=font)

    plt.savefig(output_file, dpi=300)
    plt.show()
    plt.close()



if __name__ == '__main__':
    ### Similarity search for multiple SMILES ###
    input_file_library = 'tests/library.csv'
    # input_file_library = '/Users/guohan/Documents/Projects/Datasets/HTS/Combination/forGNN/HTS_forGNN_446663.csv'
    input_file_query = 'tests/query_cmps.csv'
    output_folder = 'tests/similarity_search_results'
    output_option = 'satisfied'

    ### Method: 'fingerprint' ###
    # similarity_search_multiple_query(input_file_library, input_file_query, method='fingerprint',
    #                                  similarity_cutoff=0.3, output_folder=output_folder, output_option=output_option)
    ### Method: 'mcs' ###
    # similarity_search_multiple_query(input_file_library, input_file_query,method='mcs',
    #                                  mcs_match='exact', similarity_cutoff=0.3, output_folder=output_folder, output_option=output_option)
    ### Method: 'substructure' ###
    similarity_search_multiple_query(input_file_library, input_file_query, method='substructure',
                                     substructure_method='SMARTS',output_folder=output_folder, output_option=output_option)


    ### Select analogs ###
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
    # analog_method = 'topN'
    # similarity_topN = 3
    # select_analogs(input_file_query, analogs_dir, analog_method, similarity_topN=similarity_topN, deduplication=False)


    ### Plot similarity score distribution ###
    # input_file = 'tests/query_cmps_rank1_5.csv'
    # plot_distribution(input_file)


