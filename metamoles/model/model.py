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

import sklearn
from sklearn import linear_model
from sklearn.model_selection import train_test_split

def pair_query_compound(master_df, enzyme_col, pubchem_col, smiles_col, pubchem_sid):
    """
    pair_query_compound_with_enzymes() queries pubchem to get a SMILES string from an input pubchem_sid,
        then pairs that query compound with each unique enzyme id in the master dataframe

    Args:
        master_df (pandas.DataFrame): master dataframe containing enzyme ids
        enzyme_col (str): column containing enzyme id
        pubchem_col (str): column containing pubchem sid
        smiles_col (str): column containing SMILES string
        pubchem_sid (str): query PubChem sid

    Returns:
        pandas.DataFrame: with rows added to include query compound
    """
    master_df = master_df[[enzyme_col, pubchem_col, smiles_col]]
    new_pairs = []
    smiles, _ = sid_to_smiles(pubchem_sid)
    if len(smiles) == 0:
        raise 'query compound SMILES string could not be retrieved'
    else:
        pass
    unique_enzymes = master_df[enzyme_col].unique().tolist()
    for enzyme in unique_enzymes:
        pair = {enzyme_col:enzyme, pubchem_col:pubchem_sid, smiles_col:smiles}
        new_pairs.append(pair)
    new_paris_df = pd.DataFrame(new_pairs)
    output_df = pd.concat((master_df, new_paris_df), axis=0, sort=False)
    return output_df

def query_model(master_df, query_sid):
    """
    NOTE: Fields containing enzyme, compound PubChem sid, and SMILES string must be named
        ['entry', 'PubChem', 'SMILES'] respectively
    """
    # get query SMILES string & pair query compound with each unique enzyme in the master DataFrame
    updated_df = pair_query_compound(master_df, 'entry', 'PubChem', 'SMILES', query_sid)
    # calculate molecular distances between products of the same enzyme
    distance_df = calculate_dist(updated_df)
    # remove any rows that are not the query compound
    reduced_df = distance_df[distance_df['PubChem'] == query_sid]
    # get dummy variables to represent enzyme class
    query_df = binarize_enzyme_class(reduced_df, 'entry')
    query_df = query_df.reset_index(drop=True)
    # add in compound features with RDKit
    cpd_query_df = create_cpd_info(query_df)

    # re-instantiate log reg model
    feature_df = master_df[['dist', 'enzyme_class_1', 'enzyme_class_2', 'enzyme_class_3',
           'enzyme_class_4', 'enzyme_class_5', 'enzyme_class_6', 'enzyme_class_7',
           'n_O', 'n_N', 'n_S', 'n_X', 'DoU']]
    features = np.array(feature_df) #shape balance array for regression
    reactions = list(master_df['reacts'])
    feature_train, feature_test, reaction_train, reaction_test = train_test_split(features, reactions,
                                                      test_size=0.20, random_state=42)
    model_1 = linear_model.LogisticRegression(solver='liblinear', penalty='l1', random_state=1, class_weight='balanced')
    model_1.fit(feature_train, np.ravel(reaction_train))

    # select query features
    query_feat_df = query_df[['dist', 'enzyme_class_1', 'enzyme_class_2', 'enzyme_class_3',
           'enzyme_class_4', 'enzyme_class_5', 'enzyme_class_6', 'enzyme_class_7',
           'n_O', 'n_N', 'n_S', 'n_X', 'DoU']]
    # query reactive enzymes
    predictions = model_1.predict(query_feat_df)
    pred = model_1.predict_proba(query_feat_df)

    # write results to a DataFrame
    prediction_values = pd.DataFrame(pred)
    model_descriptive_df = pd.DataFrame()
#     model_descriptive_df['0']=prediction_values[0]
    model_descriptive_df['p_reacts']=prediction_values[1]
    prediction_df = pd.merge(model_descriptive_df, query_df, left_index=True, right_index=True)
    # sort DataFrame
    prediction_df = prediction_df.sort_values(by=['p_reacts'], ascending=False)
    # reset index in output dataframe
    prediction_df = prediction_df.reset_index(drop=True)
    # add rank to dataframe
    prediction_df['rank'] = prediction_df.index + 1
    # return DataFrame
    return prediction_df
