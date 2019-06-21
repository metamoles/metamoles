import os

import pandas as pd
from pandas.util.testing import assert_frame_equal
from rdkit import Chem
from rdkit.Chem import AllChem

import metamoles
from metamoles import curation

data_path = os.path.join(metamoles.__path__[0], 'data')

def test_explode_dataframe():
    """Unit tests for explode_dataframe()
    """
    df = curation.create_kegg_df(data_path + "/test_kegg_enzyme_records.txt.gz","enzyme")
    exploded_df = curation.explode_dataframe(df, curation.parse_compound_ids,
                                    'product', ['entry'])
    assert exploded_df.shape == (11, 2)
    assert exploded_df['product'].tolist() == ['C00071','C00004','C00080','C01450',
        'C00071','C00005','C00080','C00441','C00004','C00005', 'C00080']
    return

# test_remove_cofactors --- is not available

def test_binarize_enzyme_class():
     """Unit tests for binarize_enzyme_class()
     """
     df = curation.create_kegg_df(data_path + "/test_kegg_enzyme_records.txt.gz","enzyme")
     actual = curation.binarize_enzyme_class(df,column="entry").enzyme_class_1.tolist()
     assert actual == [1, 1, 1]
     return

def test_create_negative_matches():
    """Unit tests for create_negative_matches()
    """
    df = curation.create_kegg_df(data_path + "/test_kegg_enzyme_records.txt.gz","enzyme")
    exploded_df = curation.explode_dataframe(df, curation.parse_compound_ids,
                                    'product', ['entry'])
    pos_df, neg_df = curation.create_negative_matches(exploded_df, 'entry', 'product')
    assert pos_df.shape == (11, 3)
    assert neg_df.shape == (7, 3)
#    assert pos_df['enzyme'].tolist() == ['1.1.1.2', '1.1.1.2', '1.1.1.2', '1.1.1.1', '1.1.1.1', '1.1.1.1', '1.1.1.1', '1.1.1.3', '1.1.1.3', '1.1.1.3', '1.1.1.3'], "ValueError enzyme listing is wrong"
#    assert neg_df['product'].tolist() == ['C00004', 'C01450', 'C00441', 'C00441', 'C00005', 'C01450', 'C00071'], "ValueError compound listing is wrong"
    return

def test_remove_single_cpd_rows():
     """Unit tests for remove_single_cpd_rows()
     """
     smilesdf = pd.read_csv(data_path + "/playground_df_cleaned_kegg_with_smiles.csv")
     assert curation.remove_single_cpd_rows(smilesdf,"entry","SMILES").shape == (50, 6)
     return

def test_join_pubchem_ids():
     """Unit tests for join_pubchem_ids()
     """
     pubchemdf = curation.create_kegg_df(data_path + "/test_kegg_compound_records.txt.gz","compound")
     masterdf =  curation.create_kegg_df(data_path + "/test_kegg_enzyme_records.txt.gz","enzyme")
     return

def test_sid_to_smiles():
    """Unit test for pubchem_client.py sid_to_smiles."""
    sids = ['3489', '3990']
    expected = ['C(CO)N', 'C1CSSC1CCCCC(=O)O']
    actual = []
    for sid in sids:
        result_smile = curation.sid_to_smiles(sid)
        assert len(
            result_smile) >= 1, 'SMILES string is very short. Check SMILES.'
        isinstance(result_smile, str), 'SMILES not returned as string.'
        actual.append(result_smile[0])
    assert expected == actual, 'Actual SMILES are not the expected SMILES.'
    return

def test_kegg_df_to_smiles():
    """Unit test for pubchem_client.py kegg_df_to_smiles."""
    test_frame = pd.DataFrame([['space fill', 'ethanolamine', '3489'], [
                              'space fill', 'pyruvate', '3324']], columns=['Filler', 'Compound Name', 'SID'])
    expected_frame = pd.DataFrame([[int(700),
                                    'C(CO)N',
                                    'space fill',
                                    'ethanolamine',
                                    '3489'],
                                   [int(1060),
                                    'CC(=O)C(=O)O',
                                    'space fill',
                                    'pyruvate',
                                    '3324']],
                                  columns=['CID',
                                           'SMILES',
                                           'Filler',
                                           'Compound Name',
                                           'SID'])
    result_frame = curation.kegg_df_to_smiles(test_frame, 'SID')
    assert_frame_equal(
        result_frame[0], expected_frame), 'Did not generate expected df.'
    return
