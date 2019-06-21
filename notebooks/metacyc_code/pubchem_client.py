import pubchempy as pc


def cid_df_to_smiles(df, cid_colname):
    """
    Args:
        df : pandas dataframe with CID numbers
        cid_colname (str) : name of column that contains PubChem SID numbers

    Returns:
        df : modified with column SMILES
        unsuccessful_list : list of CIDs for which no SMILES were found

    """

    res = []
    unsuccessful_list = []

    for index, row in df.iterrows():
        cid = row[cid_colname]
        try:
            # pubchempy calls to get compound info
            compound = pc.get_compounds(cid)[0]
            smiles = compound.canonical_smiles
            res.append(smiles)

        except BaseException:
            res.append('none')
            unsuccessful_list.append(cid)
            pass

    df['SMILES'] = res
    # df.to_csv(r'../datasets/df_cleaned_kegg_with_smiles.csv')

    return df, unsuccessful_list
