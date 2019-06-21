
import pandas as pd

from pandas.util.testing import assert_frame_equal

import pubchem_client


def test_cid_df_to_smiles():
    """Unit test for pubchem_client.py kegg_df_to_smiles."""

    test_frame = pd.DataFrame([['EC-1.1.1.321', 'CPD-685', 1, 5363397],
 ['EC-1.1.1.111', '1-INDANOL', 1, 22819],
 ['EC-1.21.99.M2', '4-HYDROXYPHENYLACETATE', 1, 4693933],
 ['EC-1.21.99.M2', 'Cl-', 1, 312]], columns=['EC', 'Compound', 'Reacts', 'CID'])

    expected_frame = pd.DataFrame([['EC-1.1.1.321', 'CPD-685', 1, 5363397, 'CC(=CCO)CCC=C(C)CO'],
 ['EC-1.1.1.111', '1-INDANOL', 1, 22819, 'C1CC2=CC=CC=C2C1O'],
 ['EC-1.21.99.M2',
  '4-HYDROXYPHENYLACETATE',
  1,
  4693933,
  'C1=CC(=CC=C1CC(=O)[O-])O'],
 ['EC-1.21.99.M2', 'Cl-', 1, 312, '[Cl-]']],
                                  columns=['EC',
                                           'Compound',
                                           'Reacts',
                                           'CID',
                                           'SMILES',
                                           ])
    cid_colname = 'CID'
    result_frame = pubchem_client.cid_df_to_smiles(test_frame, cid_colname)

    assert_frame_equal(
        result_frame[0], expected_frame), 'Did not generate expected df.'

    return

