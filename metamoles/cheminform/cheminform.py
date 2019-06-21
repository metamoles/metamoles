import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.Descriptors import MolWt
from rdkit import DataStructs

def input_data(input_df): #cleans input df and returns neccessary elements
    """From the input dataframe, removes rows that do not contain product
    SMILES strings. Returns the cleaned dataframe"""
    for index, row in input_df.iterrows():

        if row['SMILES'] == 'none':

            input_df.drop(index, inplace=True)

    return input_df


def fingerprint_products(input_df): #fingerprints all products in a given df
    """From the input dataframe, makes a list of rdkit Mol objects and makes a
    list of rdkit fingerprints generated from those Mol objects. Inserts both
    lists as new columns and returns the expanded dataframe."""
    mol_list = []
    fp_list = []

    for index, row in input_df.iterrows():
        mol_list.append(Chem.rdmolfiles.MolFromSmiles(row['SMILES'])) #get mols from SMILES and add mols to list
        fp_list.append(FingerprintMols.FingerprintMol(Chem.rdmolfiles.MolFromSmiles(row['SMILES']))) #get fingerprints from mols and and fingerprints to list

    input_df['Mol'] = mol_list
    input_df['Fingerprint'] = fp_list

    return input_df

def sim_i_j(row_i, row_j):
    """For two given rows of a dataframe, use the rdkit fingerprints to compute
    TanimotoSimilarity and return the resulting float"""
    return DataStructs.FingerprintSimilarity(row_i['Fingerprint'], row_j['Fingerprint'], metric=DataStructs.TanimotoSimilarity)

def sim_i_all(input_df, index_i, row_i, metric):
    """From the input dataframe, check the passed indexes against the DataFrame,
    and construct a new dataframe which is the similarity matrix of all of the
    products contained in the dataframe."""
    for index_j, row_j in input_df.iterrows():
        if index_j < index_i: #skip redundant rows
            continue
        elif index_i == index_j: #autocorrelate rows
            metric.loc[index_i, index_j] = 1
        else:
            metric.loc[index_i, index_j] = sim_i_j(row_i, row_j) #fill matrix with calculated similarity at two positions at once
            metric.loc[index_j, index_i] = metric.loc[index_i, index_j]
    return

def sim_metric(input_df):
    """From an input_df, use sim_i_j and sim_i_all to build and return a
    similarity matrix dataframe."""
    metric = pd.DataFrame()
    for index_i, row_i in input_df.iterrows():
        sim_i_all(input_df, index_i, row_i, metric)
    return metric

def calculate_dist(input_df):
    """Main method, takes an input dataframe and builds and returns a master
    dataframe which is the original dataframe, with three additional columns,
    an rdkit Mol column, an rdkit Fingerprint column, and a column which
    describes the average distance of a product row to all the products of the
    associated enzyme entry. Requires the KEGG enzyme entry column to be named 'entry'
	and the SMILES string column to be named 'SMILES' """

    master_df = fingerprint_products(input_data(input_df))    #expand input df: generate mols from SMILES then generate fingerprints from mols, adding columns for each
    # enzyme_df_list = split_by_enzyme(input_df)    #split expanded df by rows, grouped by enzyme entry (1.1.1.110 etc), into a list of dataframes
    unique_enzymes = set(master_df['entry'].unique()) # create set of unique enzymes
    dist_lookup = {} # initialize master dist list
    for enzyme in unique_enzymes:    #loop through list of enzyme dataframes
        # enzyme_df['Dist'] = '' #initialize distance column
        enzyme_df = master_df[master_df['entry'] == enzyme]
        metric = sim_metric(enzyme_df) #get similarity matrix dataframe
        vals = metric.values #use np array of similarity matrix
        start_at = 1 #skip autocorrelation
        dist_list =[] #initialize list
        for i in range(len(vals)-1): #row of matrix except for last row
            for j in range(start_at, len(vals)): #col of matrix skipping first column
                dist_list.append(vals[i][j]) #add distance value to list
            start_at += 1 #start at higher index to skip redundancy
        avg_dist = sum(dist_list)/len(dist_list) #compute average distance
        dist_lookup[enzyme] = avg_dist
        # for _, row in enzyme_df.iterrows():    #loop through enzyme dataframe
        #     # enzyme_df['Dist'].loc[index] = avg_dist #add averaged distance to each product row of enzyme dataframe
    master_df['dist'] = [dist_lookup[row['entry']] for _, row in master_df.iterrows()]
    return master_df

def count_C(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
def count_O(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
def count_N(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
def count_P(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
def count_S(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
def count_X(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9 or atom.GetAtomicNum() == 17 or
        atom.GetAtomicNum() == 35 or atom.GetAtomicNum() == 53)
def count_H(mol):
    H = 0
    for i in range(mol.GetNumAtoms()):
        H += mol.GetAtomWithIdx(i).GetTotalNumHs(includeNeighbors=True)
    return H

def cpd_inform(SMILES):
    """
    A function for getting compound information from SMILES string it received
        a SMILES string and return a dictionary of information consisted of
        number of C, H, O , N, P, S, X, Degree of Unsaturation
        and Molecular Weight
    """
    info = []
    mol = Chem.rdmolfiles.MolFromSmiles(SMILES)
    info.append(float(count_C(mol)))
    info.append(float(count_H(mol)))
    info.append(float(count_O(mol)))
    info.append(float(count_N(mol)))
    info.append(float(count_P(mol)))
    info.append(float(count_S(mol)))
    info.append(float(count_X(mol)))
    info.append((2*info[0] + 2 + info[3] + info[4] - info[6] - info[1])/2) # it is (2*C + 2 + N + P - X - H)/2
    info.append(MolWt(mol))
    return info

# Create a function that create a new column of chemical information

def create_cpd_info(input_df, col_name='SMILES'):

    """
    Receive a DataFrame and return a dataframe with additional columns
        named n_C, n_H, ..., DoU, and MW
    """
    n_C = []
    n_H = []
    n_O = []
    n_N = []
    n_P = []
    n_S = []
    n_X = []
    DoU = []
    MW = []

    for _, row in input_df.iterrows():
        mol = row[col_name]
        info = cpd_inform(mol)
        n_C.append(info[0])
        n_H.append(info[1])
        n_O.append(info[2])
        n_N.append(info[3])
        n_P.append(info[4])
        n_S.append(info[5])
        n_X.append(info[6])
        DoU.append(info[7])
        MW.append(info[8])

    input_df['n_C'] = n_C
    input_df['n_H'] = n_H
    input_df['n_O'] = n_O
    input_df['n_N'] = n_N
    input_df['n_P'] = n_P
    input_df['n_S'] = n_S
    input_df['n_X'] = n_X
    input_df['DoU'] = DoU
    input_df['MW'] = MW

    return input_df
