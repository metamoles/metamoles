import pandas as pd

import mol_sim


def test_input_data():
    '''Tests input_data function in mol_sim.py'''

    input_df = pd.read_csv('mol_sim_test_df.csv')
    test_df = mol_sim.input_data(input_df)

    assert isinstance(test_df, pd.DataFrame), """TypeError,
    function should return a pandas dataframe"""

    return '1/1 tests successful'


def test_fingerprint_products():
    '''Tests fingerprint_products function in mol_sim.py'''

    input_df = pd.read_csv('mol_sim_test_df.csv')
    test_df = mol_sim.input_data(input_df)

    assert isinstance(mol_sim.fingerprint_products(
        test_df), pd.DataFrame), """TypeError,
    function should return a pandas dataframe"""

    return '1/1 tests successful'


def test_split_by_enzyme():
    '''Tests split_by_enzyme function in mol_sim.py'''

    input_df = pd.read_csv('mol_sim_test_df.csv')
    test_df = mol_sim.fingerprint_products(mol_sim.input_data(input_df))

    assert isinstance(mol_sim.split_by_enzyme(test_df), list), """TypeError,
    function should return a pandas dataframe"""

    return '1/1 tests successful'


def test_sim_i_j():
    '''Tests sim_i_j function in mol_sim.py'''

    input_df = pd.read_csv('mol_sim_test_df.csv')
    test_df = mol_sim.fingerprint_products(mol_sim.input_data(input_df))

    A = test_df.iloc[0]

    assert mol_sim.sim_i_j(A, A) == 1, "Self correlation is broken"

    return '1/1 tests successful'


def test_sim_i_all():
    '''Test sim_i_all functionin mol_sim.py'''

    input_df = pd.read_csv('mol_sim_test_df.csv')
    test_df = mol_sim.fingerprint_products(mol_sim.input_data(input_df))
    metric = pd.DataFrame()

    assert metric.empty, """ShapeError, input metric dataframe
    should be initialized as empty"""

    for index, row in test_df.iterrows():
        assert mol_sim.sim_i_all(
            test_df, index, row, metric) is None, """OutputError, function
        shouldn't return anything"""
        assert metric[index].all(
        ) >= 0 and metric[index].all() <= 1.0, """ValueError,
        metric should be between 0 and 1"""

    return "3/3 Tests successful"


def test_sim_metric():
    '''Test sim_i_all functionin mol_sim.py'''

    input_df = pd.read_csv('mol_sim_test_df.csv')
    test_df = mol_sim.fingerprint_products(mol_sim.input_data(input_df))
    assert isinstance(mol_sim.sim_metric(test_df), pd.DataFrame), """TypeError,
    function should return a dataframe"""
    assert mol_sim.sim_metric(
        test_df).isnull().values.any() == False, """ValueError,
    function-generated dataframe should not contain null values"""

    return "2/2 Tests successful"


def test_main():

    df = pd.read_csv('mol_sim_test_df.csv')
    test_df = mol_sim.main(df)

    assert isinstance(test_df, pd.DataFrame), """TypeError,
    function should return a dataframe"""

    return "1/1 Tests successful"
