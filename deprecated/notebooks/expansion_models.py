import numpy as np
import pandas as pd
import pubchempy as pc

#rdkit imports
import rdkit
from rdkit import Chem
from rdkit.Chem.EState import Fingerprinter
from rdkit.Chem.rdmolops import RDKFingerprint
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from rdkit.Avalon.pyAvalonTools import GetAvalonFP

def check_for_known(all_enz_df, mol_fingerprint, threshold):
    
    bad_cols = ['Known', 'Mol', 'Fingerprint']
    
    for col in bad_cols:
        if col in all_enz_df.columns:
            all_enz_df.drop(columns=col, inplace=True)
        else:
            continue
            
    fingerprint = mol_fingerprint
    
    input_df = fingerprint_products(all_enz_df) # fingerprint the input dataframe and return it
    
    input_df['Known'] = '' # initialize similarity column
    
    for index, row in input_df.iterrows():
        similarity = DataStructs.FingerprintSimilarity(fingerprint, row['Fingerprint'], 
                                                       metric=DataStructs.TanimotoSimilarity)
        input_df['Known'].loc[index] = similarity
    
    
    known_df = input_df[input_df['Known'] >= threshold]
    
    if len(known_df) > 0:
        known_df.sort_values(by='Known', ascending=False, inplace=True)
        result = known_df
    else:
        #call to promiscuous search code here
        result = print('No known enzymes. Beginning promiscuous search.')
        
    return result


def cid_to_smiles(cid):

    try:
        compound = pc.get_compounds(cid)[0]
        smiles = compound.canonical_smiles
    except BaseException:
        pass
    
    return smiles, cid


def dist_for_expansion(prom_df, pubchem_cid, num_similar):
    
    #####this probably isn't necessary because in the final version, this will all be run in one sitting and not saved and loaded etc
    #if 'Fingerprint' and 'Mol' in prom_df.columns:
     #   prom_df = prom_df.drop(columns=['Fingerprint','Mol'])
    #else:
     #   continue
    
    bad_cols = ['ExpanDist']
    
    for col in bad_cols:
        if col in prom_df.columns:
            prom_df.drop(columns=col, inplace=True)
        else:
            continue
    
        
    smiles, _ = cid_to_smiles(pubchem_cid)
    
    if len(smiles) == 0:
        raise 'query compound SMILES string could not be retrieved'
    else:
        pass
    
    mol = Chem.rdmolfiles.MolFromSmiles(smiles)
    fingerprint = FingerprintMols.FingerprintMol(mol)
    
    
    prom_df['ExpanDist'] = ''
    
    for index, row in prom_df.iterrows():
        comp = DataStructs.FingerprintSimilarity(row['Fingerprint'], fingerprint, metric=DataStructs.TanimotoSimilarity)
        prom_df['ExpanDist'].loc[index] = comp

    prom_df.sort_values(by='ExpanDist', ascending=False, inplace=True)
    
    if num_similar == 'None':
        n = 20
    else:
        n = num_similar
    
    selected = prom_df.iloc[:n,:].copy()
    
    return selected

def fingerprint_products(input_df): #fingerprints all products in a given df
    '''DocString'''
    
    mol_list = []
    fp_list = []
    
    for index, row in input_df.iterrows():
        mol_list.append(Chem.rdmolfiles.MolFromSmiles(row['SMILES'])) #get mols from SMILES and add mols to list
        fp_list.append(FingerprintMols.FingerprintMol(Chem.rdmolfiles.MolFromSmiles(row['SMILES']))) #get fingerprints from mols and and fingerprints to list
        
    input_df.insert(1, column='Mol', value=mol_list)
    input_df.insert(2, column='Fingerprint', value= fp_list)
            
    return input_df
