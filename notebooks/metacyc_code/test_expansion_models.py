import pandas as pd

# rdkit imports
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols

from pandas.util.testing import assert_frame_equal

import expansion_models


def test_check_for_known():
    """Testing check_for_known function with CIDs = [243, 985]"""

    expected_df = pd.DataFrame([['EC-1.14.14.80',
                                 'PALMITATE',
                                 504166,
                                 'CCCCCCCCCCCCCCCC(=O)[O-]',
                                 0.714717543728323,
                                 1.0],
                                ['EC-1.14.14.80',
                                 'STEARIC_ACID',
                                 3033836,
                                 'CCCCCCCCCCCCCCCCCC(=O)[O-]',
                                 0.714717543728323,
                                 1.0]],
                               columns=['enzyme',
                                        'product',
                                        'PubChemID',
                                        'SMILES',
                                        'Dist',
                                        'Known'])

    test_df = pd.DataFrame([['EC-1.14.14.77',
                             '2-METHYL-3-PHYTYL-14-NAPHTHOQUINONE',
                             5280483,
                             'CC1=C(C(=O)C2=CC=CC=C2C1=O)CC=C(C)CCCC(C)CCCC(C)CCCC(C)C',
                             1.0],
                            ['EC-1.14.14.80',
                             'CPD-10515',
                             25201835,
                             'CCCCCCCCC(C(CCCCCCCC(=O)[O-])O)O',
                             0.714717543728323],
                            ['EC-1.14.14.80',
                             'PALMITATE',
                             504166,
                             'CCCCCCCCCCCCCCCC(=O)[O-]',
                             0.714717543728323],
                            ['EC-1.14.14.80',
                             'OLEATE-CPD',
                             5460221,
                             'CCCCCCCCC=CCCCCCCCC(=O)[O-]',
                             0.714717543728323],
                            ['EC-1.14.14.80',
                             'STEARIC_ACID',
                             3033836,
                             'CCCCCCCCCCCCCCCCCC(=O)[O-]',
                             0.714717543728323]],
                           columns=['enzyme',
                                    'product',
                                    'PubChemID',
                                    'SMILES',
                                    'Dist'])

    test_smiles = ['C1=CC=C(C=C1)C(=O)O', 'CCCCCCCCCCCCCCCC(=O)O']

    fingerprint_list = []
    for smile in test_smiles:
        mol = Chem.rdmolfiles.MolFromSmiles(smile)
        fingerprint = FingerprintMols.FingerprintMol(mol)
        fingerprint_list.append(fingerprint)

    for mol_fingerprint in fingerprint_list:
        test_result = expansion_models.check_for_known(
            test_df, mol_fingerprint, 0.95)

        if test_result is not type(pd.core.frame.DataFrame):
            continue
        else:
            test_result.drop(columns=['Fingerprint', 'Mol'], inplace=True)
            assert test_result == expected_df

    return


def test_dist_for_expansion():
    test_df = pd.DataFrame([['EC-1.14.14.77',
                             '2-METHYL-3-PHYTYL-14-NAPHTHOQUINONE',
                             5280483,
                             'CC1=C(C(=O)C2=CC=CC=C2C1=O)CC=C(C)CCCC(C)CCCC(C)CCCC(C)C',
                             1.0],
                            ['EC-1.14.14.80',
                             'CPD-10515',
                             25201835,
                             'CCCCCCCCC(C(CCCCCCCC(=O)[O-])O)O',
                             0.714717543728323],
                            ['EC-1.14.14.80',
                             'PALMITATE',
                             504166,
                             'CCCCCCCCCCCCCCCC(=O)[O-]',
                             0.714717543728323],
                            ['EC-1.14.14.80',
                             'OLEATE-CPD',
                             5460221,
                             'CCCCCCCCC=CCCCCCCCC(=O)[O-]',
                             0.714717543728323]],
                           columns=['enzyme',
                                    'product',
                                    'PubChemID',
                                    'SMILES',
                                    'Dist'])

    expected = pd.DataFrame([['EC-1.14.14.80',
                              'PALMITATE',
                              504166,
                              'CCCCCCCCCCCCCCCC(=O)[O-]',
                              0.714717543728323,
                              1.0],
                             ['EC-1.14.14.80',
                              'OLEATE-CPD',
                              5460221,
                              'CCCCCCCCC=CCCCCCCCC(=O)[O-]',
                              0.714717543728323,
                              0.75],
                             ['EC-1.14.14.80',
                              'CPD-10515',
                              25201835,
                              'CCCCCCCCC(C(CCCCCCCC(=O)[O-])O)O',
                              0.714717543728323,
                              0.7377049180327869],
                             ['EC-1.14.14.77',
                              '2-METHYL-3-PHYTYL-14-NAPHTHOQUINONE',
                              5280483,
                              'CC1=C(C(=O)C2=CC=CC=C2C1=O)CC=C(C)CCCC(C)CCCC(C)CCCC(C)C',
                              1.0,
                              0.3515625]],
                            columns=['enzyme',
                                     'product',
                                     'PubChemID',
                                     'SMILES',
                                     'Dist',
                                     'ExpanDist'])

    fingers = expansion_models.fingerprint_products(test_df)
    actual = expansion_models.dist_for_expansion(fingers, 985, 4)
    actual.drop(columns=['Mol', 'Fingerprint'], inplace=True)
    assert_frame_equal(actual.reset_index(drop=True),
                       expected.reset_index(drop=True), check_dtype=False)

    return
