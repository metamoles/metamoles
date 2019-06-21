# rdkit imports
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs


# housekeeping imports
import pandas as pd


def fingerprint_products(input_df):  # fingerprints all products in a given df
    """
    Args:
        input_df : pandas dataframe with smiles in column named 'SMILES'

    Returns:
        input_df : with added mol and fingerprint columns
    """

    mol_list = []
    fp_list = []

    for index, row in input_df.iterrows():
        # get mols from SMILES and add mols to list
        mol_list.append(Chem.rdmolfiles.MolFromSmiles(row['SMILES']))
        fp_list.append(FingerprintMols.FingerprintMol(Chem.rdmolfiles.MolFromSmiles(
            row['SMILES'])))  # get fingerprints from mols and and fingerprints to list

    input_df.insert(1, column='Mol', value=mol_list)
    input_df.insert(2, column='Fingerprint', value=fp_list)

    return input_df


def input_data(input_df):  # cleans input df and returns neccessary elements
    '''DocString'''

    for index, row in input_df.iterrows():

        smiles = row['SMILES']
        if len(smiles) <= 2:

            input_df.drop(index, inplace=True)

    return input_df


def sim_i_j(row_i, row_j):
    """DocString"""
    return DataStructs.FingerprintSimilarity(
        row_i['Fingerprint'],
        row_j['Fingerprint'],
        metric=DataStructs.TanimotoSimilarity)


def sim_i_all(input_df, index_i, row_i, metric):
    """DocString"""
    for index_j, row_j in input_df.iterrows():
        if index_j < index_i:  # skip redundant rows
            continue
        elif index_i == index_j:  # autocorrelate rows
            metric.loc[index_i, index_j] = 1
        else:
            # fill matrix with calculated similarity at two positions at once
            metric.loc[index_i, index_j] = sim_i_j(row_i, row_j)
            metric.loc[index_j, index_i] = metric.loc[index_i, index_j]
    return


def sim_metric(input_df):
    """DocString"""
    metric = pd.DataFrame()
    for index_i, row_i in input_df.iterrows():
        sim_i_all(input_df, index_i, row_i, metric)
    return metric


def add_prom_tag(df, cut_off):
    """Cut-off is inclusive"""
    tag = []
    for index, row in df.iterrows():
        if row['Dist'] >= cut_off:
            tag.append(0.0)
        elif row['Dist'] < cut_off:
            tag.append(1.0)

    df['Promiscuous'] = tag

    return df


def split_by_enzyme(input_df):
    '''DocString'''

    unique_enzymes = set(input_df['enzyme'].unique())

    enzyme_df_list = []

    for entry in unique_enzymes:  # for each unique enzyme in the input dataframe...

        # ...initialize a new dataframe with the same columns as the input dataframe...
        enzyme_df = pd.DataFrame(columns=input_df.columns)

        # ...iterate through the input dataframe...
        for index, row in input_df.iterrows():

            # ... and add product rows that correspond to the unique enzyme entry...
            if row['enzyme'] == entry:
                enzyme_df.loc[index] = row

        # ...then add the completed dataframe of unique enzyme products to a list
        enzyme_df_list.append(enzyme_df)

    return enzyme_df_list  # return list of dataframes


def main(input_df):
    '''DocString'''

    # expand input df: generate mols from SMILES then generate fingerprints
    # from mols, adding columns for each
    input_df = fingerprint_products(input_df)

    # input_df.drop(columns=['Mol'])

    # split expanded df by rows, grouped by enzyme entry (1.1.1.110 etc), into
    # a list of dataframes
    enzyme_df_list = split_by_enzyme(input_df)

    # outer = []

    for enzyme_df in enzyme_df_list:  # loop through list of enzyme dataframes

        enzyme_df['Dist'] = ''  # initialize distance column

        metric = sim_metric(enzyme_df)  # get similarity matrix dataframe

        vals = metric.values  # use np array of similarity matrix

        start_at = 1  # skip autocorrelation

        dist_list = []  # initialize list

        if len(vals) == 1:

            dist_list.append(vals)  # add distance value to list

        elif len(vals) > 1:
            for i in range(len(vals) - 1):  # row of matrix except for last row

                for j in range(
                        start_at,
                        len(vals)):  # col of matrix skipping first column

                    dist_list.append(vals[i][j])  # add distance value to list

                start_at += 1  # start at higher index to skip redundancy

        # outer.append(dist_list)
        avg_dist = sum(dist_list) / len(dist_list)  # compute average distance
        # outer.append(avg_dist)
        for index, row in enzyme_df.iterrows(
        ):  # loop through enzyme dataframe
            # add averaged distance to each product row of enzyme dataframe
            enzyme_df['Dist'].loc[index] = avg_dist

    # concatenate enzyme dataframes into master_df
    master_df = pd.concat(enzyme_df_list)

    return master_df
