import pubchempy as pc

# rdkit imports
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs


def check_for_known(all_enz_df, mol_fingerprint, threshold):
    """
    Args:
        all_enz_df : pandas dataframe with all known enzymatic reactions
        mol_fingerprint : RDKit fingerprint of query molecule
        threshold (float, 0-1.0) : cut-off value for similarity.
                                    Recommend 0.95 or above for truly similar compounds.

    Returns:
        known_df : df with known enzymatic transformations of this molecule
                    ordered from most similar molecule to least similar
        or
        printed message : indication to begin searching through promiscuous enzymes

    """
    # drop any columns leftover from the last iteration
    bad_cols = ['Known', 'Mol', 'Fingerprint']

    for col in bad_cols:
        if col in all_enz_df.columns:
            all_enz_df.drop(columns=col, inplace=True)
        else:
            continue

    fingerprint = mol_fingerprint

    # fingerprint the input dataframe and return it
    input_df = fingerprint_products(all_enz_df)

    input_df['Known'] = ''

    for index, row in input_df.iterrows():
        # get similarity between query compound and each molecule 
        # in df
        similarity = DataStructs.FingerprintSimilarity(
            fingerprint, row['Fingerprint'], metric=DataStructs.TanimotoSimilarity)
        input_df['Known'].loc[index] = similarity

    # select only reactions that have product similarity to query
    # greater than the threshold
    known_df = input_df[input_df['Known'] >= threshold]

    if len(known_df) > 0:
        known_df.sort_values(by='Known', ascending=False, inplace=True)
        result = known_df
    else:
        # call to promiscuous search code here
        result = print('No known enzymes. Begin promiscuous search.')

    return result


def cid_to_smiles(cid):
    """
    Args:
        cid : cid to query for smiles
        
    Returns:
       smiles (str) : smiles string retrieved from pubchem
        cid : as input
    
    """
    
    try:
        # query pubchem for information about this cid
        compound = pc.get_compounds(cid)[0]
        smiles = compound.canonical_smiles
    except BaseException:
        pass

    return smiles, cid


def dist_for_expansion(prom_df, pubchem_cid, num_similar):
    """
    Args:
        prom_df : pandas dataframe with all promiscuous enzymatic reactions
        pubchem_cid : cid of query molecule
        num_similar (int) : number of similar compounds to retrieve for tree search.
                            If input is 'None', default is 20 molecules.
                            A higher number will make search time longer.

    Returns:
        selected (df) : df with selected promiscuous enzyme pairs to use in tree search
    """
    
    # drop column leftover from previous iteration
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
    
    # get molecule and fingerprint for query compound
    # using rdkit
    mol = Chem.rdmolfiles.MolFromSmiles(smiles)
    fingerprint = FingerprintMols.FingerprintMol(mol)

    prom_df['ExpanDist'] = ''

    for index, row in prom_df.iterrows():
        # get similarity between query compound and each molecule 
        # in df
        comp = DataStructs.FingerprintSimilarity(
            row['Fingerprint'], fingerprint, metric=DataStructs.TanimotoSimilarity)
        prom_df['ExpanDist'].loc[index] = comp

    prom_df.sort_values(by='ExpanDist', ascending=False, inplace=True)

    if num_similar == 'None':
        n = 20
    else:
        n = num_similar

    selected = prom_df.iloc[:n, :].copy()

    return selected


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
