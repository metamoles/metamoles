import numpy as np
import pandas as pd
import math
import ast

df_cpd = pd.read_csv('df_cpd.csv', index_col = 0)
df_rxn = pd.read_csv('parsed_rxns.csv', index_col = 0)
df_enz = pd.read_csv('df_enzrxns.csv', index_col = 0)
df_cpd = df_cpd.set_index(keys ='UNIQUE-ID')
df_rxn = df_rxn.set_index(keys = 'UNIQUE-ID')
df_enz = df_enz.set_index(keys = 'UNIQUE-ID')

def recover_list(df, column):
    """This function will recover a list formatted string read from .csv into a list"""
    assert type(df[column][0]) != type([]), "TypeError: The data type is already a list, it should not be converted again"
    replacement = []
    for index, row in df.iterrows():

        data = []

        if type(row[column]) == type('string'):
            data = ast.literal_eval(row[column])
        else:
            pass
        replacement.append(data)
    df[column] = replacement
    return

##################################################
# Change PubChemID into int type in df_cpd
PubChemID_int = df_cpd['PubChemID'].fillna(0).astype(int)
df_cpd['PubChemID'] = PubChemID_int

# Recover list format of df_rxn
rxn_list_fix = ['EC-NUMBER', 'ERXN-NUMBER', 'SUBSTRATES', 'PRODUCTS']
for col in rxn_list_fix:
    recover_list(df_rxn, col)

# Recover list format of df_enz
enz_list_fix = ['REACTION', 'ALTERNATIVE-SUBSTRATES', '^SUBSTRATE', 'KM', 'KCAT', 'VMAX']
for col in enz_list_fix:
    recover_list(df_enz, col)

#################################################

def get_inchi(ID):

    """This function accept UNIQUE-ID and return InChI string of a certain compound"""

    inchi = df_cpd['INCHI'][ID]

    return inchi

def get_smiles(ID):

    """This function accept UNIQUE-ID and return SMILES string of a certain compound"""

    smiles = df_cpd['SMILES'][ID]

    return smiles

def get_pubchem(ID):

    """This function accept UNIQUE-ID and return InChI string of a certain compound"""
    if ID in df_cpd['PubChemID']:
        pubchem = df_cpd['PubChemID'][ID]
    else:
        pubchem = ID

    return pubchem

def add_pubchem_column(df_rxn):

    # Start from df_rxn and rerun the master dataframe again
    subs_id = []
    pdts_id = []

    for index, row in df_rxn.iterrows():

        subs = []
        for item in row['SUBSTRATES']:
            subs.append(get_pubchem(item))
        subs_id.append(subs)

        pdts = []
        for item in row['PRODUCTS']:
            pdts.append(get_pubchem(item))
        pdts_id.append(pdts)

    df_rxn['SUBSTRATES_PubChemID'] = subs_id
    df_rxn['PRODUCTS_PubChemID'] = pdts_id

    return df_rxn
    

def df_ec_number(df_rxn):

    EC = []
    rxn = []

    for index, row in df_rxn.iterrows():

        if len(row['EC-NUMBER']) > 1:
            for i in range(len(row['EC-NUMBER'])):
                EC.append(row['EC-NUMBER'][i])
                rxn.append(index)
        elif len(row['EC-NUMBER']) == 1:
            EC.append(row['EC-NUMBER'][0])
            rxn.append(index)
        else:
            EC.append('No_Data')
            rxn.append(index)

        df_master = pd.DataFrame({'EC-NUMBER' : EC,
                              'UNIQUE-ID' : rxn})

        rxn_num = []
        subs = []
        pdts = []
        gibbs = []

    for index, row in df_master.iterrows():
        ID = row['UNIQUE-ID']
        rxn_num.append(df_rxn['ERXN-NUMBER'][ID])
        subs.append(df_rxn['SUBSTRATES'][ID])
        pdts.append(df_rxn['PRODUCTS'][ID])
        gibbs.append(df_rxn['GIBBS'][ID])

    df_master['ERXN-NUMBER'] = rxn_num
    df_master['SUBSTRATES'] = subs
    df_master['PRODUCTS'] = pdts
    df_master['GIBBS'] = gibbs
    df_master.head()

    df_sorted = df_master.sort_values(by=['EC-NUMBER'])
    df_sorted.reset_index(inplace=True, drop=True)

    EC_a = 'EC-1'

    EC = []
    ID = []
    erxn = []
    subs = []
    pdts = []
    gibbs = []
    counter = 0

    ID_temp = []
    erxn_temp = []
    subs_temp = []
    pdts_temp = []
    gibbs_temp = []

    for index, row in df_sorted.iterrows():

        if row['EC-NUMBER'] == EC_a:
            ID_temp.append(row['UNIQUE-ID'])
            erxn_temp.append(row['ERXN-NUMBER'])
            subs_temp.append(row['SUBSTRATES'])
            pdts_temp.append(row['PRODUCTS'])
            gibbs_temp.append(row['GIBBS'])
            counter += 1

        elif counter == 0:
            ID.append(row['UNIQUE-ID'])
            erxn.append(row['ERXN-NUMBER'])
            subs.append(row['SUBSTRATES'])
            pdts.append(row['PRODUCTS'])
            gibbs.append(row['GIBBS'])

            EC.append(EC_a)
            EC_a = row['EC-NUMBER']
        else:
            ID.append(ID_temp)
            erxn.append(erxn_temp)
            subs.append(subs_temp)
            pdts.append(pdts_temp)
            gibbs.append(gibbs_temp)

            EC.append(EC_a)
            counter = 0
            EC_a = row['EC-NUMBER']

    df_sorted_master = pd.DataFrame({'EC-NUMBER' : EC,
                                'UNIQUE-ID' : ID,
                                'ERXN-NUMBER' : erxn,
                                'SUBSTRATES' : subs,
                                'PRODUCTS' : pdts,
                                'GIBBS' : gibbs})

    return df_sorted_master
