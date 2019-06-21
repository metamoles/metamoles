import pandas as pd
import pubchempy as pc
import gzip
import Bio
from Bio.KEGG import Compound
from Bio.KEGG import REST
from Bio.KEGG import Enzyme

import re
import scipy as sp
import pandas as pd
import numpy as np

#temp start here

def create_kegg_df(file_path: str, kegg_db: str):
    """
    create_kegg_df() parses a gzipped text file of KEGG records using the
        Biopython.Bio.KEGG package, and returns the data as a pandas dataframe.

        NOTE: only 'enzyme' and 'compound' KEGG records are supported

    Args:
        file_path (str): filepath string pointing to gzipped text file
            of KEGG records
        kegg_db (str): either 'enzyme' or 'compound' keyword argument

    Returns:
        pandas.DataFrame: containing all Enzyme.Record() fields in columns
    """
    supported_dbs = ['enzyme', 'compound']

    if kegg_db == 'enzyme':
        parser = Enzyme
    elif kegg_db == 'compound':
        parser = Compound
    else:
        raise ValueError('supported kegg_db values include: {}'.format(supported_dbs))

    field_list = [method for method in dir(parser.Record()) if not method.startswith('_')]
    data_matrix = []

    with gzip.open(file_path, 'rt') as file:
        for record in parser.parse(file):
            data_matrix.append([getattr(record, field) for field in field_list])

    kegg_df = pd.DataFrame(data_matrix, columns=field_list)
    return kegg_df

def parse_compound_ids(field: str):
    """
    parse_compound_ids() uses regular expressions to extract the KEGG compound
        IDs from a product or substrate field in a KEGG record field

    Args:
        field (str): name of field that contains KEGG compound IDs in a string

    Returns:
        list: contains parsed KEGG compound IDs
    """

    cpd_list = []
    regex = 'CPD:(C\d+)'
    # matches 'CPD:' chars exactly and captures 'C' + any following digits (\d+)
    for entry in field:
        ids = re.findall(regex, str(entry), re.IGNORECASE)
        for i in ids:
            cpd_list.append(i)

    return cpd_list

#temp stop here

def explode_dataframe(dataframe: pd.DataFrame, explosion_function,
                        explosion_target_field: str, fields_to_include: list):
    """
    explode_dataframe() applies the input explosion_function to the target
        field in each row of a dataframe. Each item in the output of the
        explosion_function is an anchor for a new row in the new dataframe. All
        of the supplied fields_to_include are added to the explosion item,
        and appended to the new dataframe row.

    Args:
        dataframe (pandas.DataFrame): input dataset
        explosion_function (function): function to be applied to target
            column in dataframe
        explosion_target_field (str): name of field in dataframe to which the
            explosion funciton will be applied
        fields_to_include (list): a list of strings that denote the columns of
            the input dataframe to be included in the output

    Returns:
        pandas.DataFrame: new exploded dataframe
    """
    new_rows = []
    for _, row in dataframe.iterrows():
        explosion_list = explosion_function(row[explosion_target_field])
        for item in explosion_list:
            row_data = [row[field] for field in fields_to_include]
            row_data.append(item)
            new_rows.append(row_data)

    fields_to_include.append(explosion_target_field)
    new_df = pd.DataFrame(new_rows, columns=fields_to_include)

    return new_df

def remove_cofactors(master_df: pd.DataFrame, master_cpd_field: str,
                     cofactor_df: pd.DataFrame, cofactor_field: str,
                     drop_na=True):
    """
    remove_cofactors() should be used to clean the dataset of cofactors. These
        will be included in the KEGG records as substrates and products, but
        are not actually products in the reaction

    Args:
        master_df (pandas.DataFrame): input dataset
        master_cpd_field (str): field that contains products
        cofactor_df (pandas.DataFrame): contains cofactors to be removed
        cofactor_field (str): field that contains cofactors
        drop_na (bool): default True

    Returns:
        pandas.DataFrame: cleaned data without cofactor entries
    """
    cofactor_list = parse_compound_ids(cofactor_df[cofactor_field])
    bool_mask = [False if cpd in cofactor_list else True for cpd in master_df[master_cpd_field]]
    clean_df = master_df[bool_mask]
    clean_df = clean_df.drop_duplicates()

    if drop_na:
        clean_df = clean_df[clean_df[master_cpd_field] != 'NA']
    else:
        pass

    return clean_df

def binarize_enzyme_class(dataframe, column):
    """
    binarize_enzyme_class() converts the enzyme class into binary dummy variables
        that are appended onto the input dataframe

    Args:
        dataframe (pandas.DataFrame): input dataset
        column (str): column name containing kegg enzyme id

    Returns:
        pandas.DataFrame: with seven columns appended for the seven enzyme classes
    """
    dataframe['enzyme_class'] = [row[column][0] for _, row in dataframe.iterrows()]
    dataframe = pd.get_dummies(dataframe, columns=['enzyme_class'])
    return dataframe

def create_negative_matches(dataframe: pd.DataFrame,
                            enzyme_field: str, compound_field: str):
    """
    create_negative_matches() returns two dataframes.
        One dataframe is positive data that contains all the enzyme-compound
        pairs that exist in the input dataset.
        The second data frame is negative data made from matching all
        enzyme-compound pairs that do not exist in the dataset.

    Args:
        dataframe (pandas.DataFrame): input dataset
        enzyme_field (str): column in dataframe that contains enzyme ids
        compound_field (str): column in dataframe that contains compound ids

    Returns:
        pandas.DataFrame: positive data
            (contains fields ['enzyme', 'product', 'reacts'])
        pandas.DataFrame: negative data
            (contains fields ['enzyme', 'product', 'reacts'])
    """
    unique_enzymes = set(dataframe[enzyme_field].unique())
    # set of all unique enzymes in provided dataframes
    unique_cpds = set(dataframe[compound_field].unique())
    # set of all unique compounds in provided dataframe

    positive_data = []
    negative_data = []
    # initialize empty lists

    for enzyme in unique_enzymes:
    # iterate through unique enzyme set
        working_prods = set(dataframe[dataframe[enzyme_field] == enzyme][compound_field].unique())
        # unique set of all products reported to reaction with this enzyme in provided dataset
        non_working_prods = (unique_cpds - working_prods)
        # set math of all remaining products in the dataset minus those reported to react

        reactions = [{'reacts':1.0, 'enzyme':enzyme, 'product':product} for product in working_prods]
        # create new entry for each positive reaction
        non_reactions = [{'reacts':0.0, 'enzyme':enzyme, 'product':product} for product in non_working_prods]
        # create new entry for each negative reaction

        positive_data.extend(reactions)
        # add positive reactions to master list
        negative_data.extend(non_reactions)
        # add negative reactions to master list

    positive_df = pd.DataFrame(positive_data)
    negative_df = pd.DataFrame(negative_data)

    return positive_df, negative_df

def remove_single_cpd_rows(dataframe, enzyme_col, smiles_col):
    """
    remove_single_cpd_rows() is meant to be a pre-processing function prior to passing a dataframe to the
        calculate_dist() function

    Args:
        dataframe (pandas.Dataframe): input dataset
        enzyme_col (str): name for column that contains kegg enzyme ids
        smiles_col (str): name for column that contains smiles string

    Returns:
        pandas.Dataframe: output dataframe with rows removed in which there was only one product paired with
            the enzyme entry, enzyme_col renamed 'entry', and smiles_col renamed 'SMILES'
    """
    dataframe = dataframe.rename(columns={enzyme_col:'entry', smiles_col:'SMILES'})
    counts_df = dataframe.groupby('entry').count()
    singles_df = counts_df[counts_df['SMILES'] == 1]
    singles = singles_df.index.tolist()
    bool_mask = [False if row['entry'] in singles else True for _, row in dataframe.iterrows()]
    clean_df = dataframe[bool_mask]
    return clean_df

def join_pubchem_ids(master_df, pubchem_df, master_join_key, pubchem_join_key,
                        pubchem_id_field):
    """
    join_pubchem_ids() takes an input dataframe containing a column of KEGG
        compound ids, and a second dataframe containing KEGG compound ids and
        their corresponding PubChem ids. The function parses the PubChem ids
        from the correct column, and joins these onto the input dataframe

    Args:
        master_df (pandas.DataFrame): input dataset
        pubchem_df (pandas.DataFrame): dataframe containing PubChem ids
        master_join_key (str): field in master_df with KEGG compound ids
        pubchem_join_key (str): field in pubchem_df with KEGG compound ids
        pubchem_id_field (str): field in pubchem_df with PubChem ids

    Returns:
         pandas.DataFrame:
    """
    pubchem_id_data = []

    for _, row in pubchem_df.iterrows():
        pubchem_id = parse_pubchem_ids(row[pubchem_id_field])
        join_key = row[pubchem_join_key]
        entry = {'pubchem_id': pubchem_id, pubchem_join_key: join_key}
        pubchem_id_data.append(entry)

    join_df = pd.DataFrame(pubchem_id_data)
    master_df = master_df.merge(join_df, left_on=master_join_key,
                                right_on=pubchem_join_key)
    master_df = master_df.drop(columns=pubchem_join_key)

    return master_df

def sid_to_smiles(sid):
    """Takes a PubChem SID. Returns the associated isomeric SMILES string and PubChem CID.

    Args:
        sid : The PubChem SID number.

    Returns:
        str: isomeric smiles.
        int: Pubchem CID number.

    """

    substance = pc.Substance.from_sid(sid)
    cid = substance.standardized_cid
    compound = pc.get_compounds(cid)[0]

    return compound.isomeric_smiles, cid

def kegg_df_to_smiles(kegg_df, column_name):
    """
    Args:
        kegg_df : pandas dataframe with SID numbers in the third column

    Returns:
        kegg_df : modified with a fourth column containing CID and fifth column containing SMILES
        unsuccessful_list : list of SIDs for which no CID or SMILES were found

    """

    res = []
    cid_list = []
    unsuccessful_list = []

    for i in range(len(kegg_df)):
        # cell index of desired SID
        sid = kegg_df.loc[i, column_name]
        try:
            smile_result = sid_to_smiles(sid)[0]
            res.append(smile_result)
            cid_result = sid_to_smiles(sid)[1]
            cid_list.append(cid_result)
        except BaseException:
            res.append('none')
            cid_list.append('none')
            unsuccessful_list.append(sid)
            pass

    kegg_df.insert(0, column='CID', value=cid_list)
    # Change this 2 to the number where the smiles column should be
    kegg_df.insert(1, column='SMILES', value=res)
    # kegg_df.to_csv(r'../datasets/df_cleaned_kegg_with_smiles.csv')

    return kegg_df, unsuccessful_list
