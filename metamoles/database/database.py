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

def create_kegg_df(file_path: str, kegg_db: str):
    """
    create_kegg_df() parses a gzipped text file of KEGG records using the
        Biopython.Bio.KEGG package, and returns the data as a pandas dataframe.

        NOTE: only 'enzyme' and 'compound' KEGG records are supported

    Args:
        file_path (str): filepath string pointing to gzipped text file
            of KEGG records
        kegg_db (str): either 'enzyme' or 'compound' keyword argument

    Returns:
        pandas.DataFrame: containing all Enzyme.Record() fields in columns
    """
    supported_dbs = ['enzyme', 'compound']

    if kegg_db == 'enzyme':
        parser = Enzyme
    elif kegg_db == 'compound':
        parser = Compound
    else:
        raise ValueError('supported kegg_db values include: {}'.format(supported_dbs))

    field_list = [method for method in dir(parser.Record()) if not method.startswith('_')]
    data_matrix = []

    with gzip.open(file_path, 'rt') as file:
        for record in parser.parse(file):
            data_matrix.append([getattr(record, field) for field in field_list])

    kegg_df = pd.DataFrame(data_matrix, columns=field_list)
    return kegg_df

def select_promiscuous_enzymes(enzyme_df: pd.DataFrame):
    """
    select_promiscuous_enzymes() down selects the enzymes from an input
        dataframe to include only those in which 2 or more reactions are
        associated with their KEGG record.

    Args:
        enzyme_df (pandas.DataFrame): must contain at least fields
            ['reaction', 'entry', 'product', 'substrate']

    Returns:
        pandas.DataFrame: pandas dataframe containing only the fields
            ['reaction', 'entry', 'product', 'substrate'] and containing only
            rows of enzymes with two or more reactions in their KEGG record
    """

    promiscuous_df = enzyme_df[[True if len(rxn) > 1 else False for rxn in enzyme_df['reaction']]]
    compact_promiscuous_df = promiscuous_df[['entry','reaction','product','substrate']]

    return compact_promiscuous_df

def parse_compound_ids(field: str):
    """
    parse_compound_ids() uses regular expressions to extract the KEGG compound
        IDs from a product or substrate field in a KEGG record field

    Args:
        field (str): name of field that contains KEGG compound IDs in a string

    Returns:
        list: contains parsed KEGG compound IDs
    """

    cpd_list = []
    regex = 'CPD:(C\d+)'
    # matches 'CPD:' chars exactly and captures 'C' + any following digits (\d+)
    for entry in field:
        ids = re.findall(regex, str(entry), re.IGNORECASE)
        for i in ids:
            cpd_list.append(i)

    return cpd_list

def parse_pubchem_ids(field: str):
    """
    parse_pubchem_ids() uses regular expressions to extract the PubChem
        compound IDs from a field in a record

    Args:
        field (str): name of a pandas.DataFrame field containing PubChem
            compound IDs in a string

    Returns:
        str: extracted pubchem_id
    """

    regex = "'PubChem', \[\'(\d+)\'\]\)"
    # matches "'PubChem', ['" characters exactly, then captures any following
    # digits (\d+), before another literal "']" character match

    ids = re.findall(regex, str(field), re.IGNORECASE)
    if len(ids) > 0:
        pubchem_id = ids[0]
    else:
        pubchem_id = ''

    return pubchem_id

def parse_reaction_ids(df: pd.DataFrame):
    """
    parse_reaction_ids() parses the list of reaction numbers from a
        dataframe containing a column labeled 'reaction'

    Args:
        df (pandas.DataFrame): must contain column 'reaction'

    Returns:
        list: contains parsed KEGG reaction IDs
    """
    reaction_list = []
    for index,row in df.iterrows():
        for reaction in row['reaction']:
            reaction_split = reaction.split("[RN:")[-1]
            if reaction_split.startswith("R") and not reaction_split.startswith("RN"):
                for i in reaction_split[:-1].split(" "):
                    reaction_list.append(i)
    # note - this approach is case sensitive and will miss reaction numbers
    # labeled with a lowercase r, or with incorrect spacing. A regex extraction
    # may be one option for improved robustness, though likely sacrifices speed.
    return reaction_list

def parse_reversible_reactions(reaction_list: list):
    """
    parse_reversible_reactions() queries the KEGG database with the input
        reaction list, and parses the results for all reactions that have been
        annotated with "<=>" in the reaction equation, which suggests that the
        catalyzed reaction is reversible

    Args:
        reaction_list (list): contains KEGG reaction IDs (e.g. 'R00709')

    Returns:
        list: contains KEGG IDs of reversible reactions
    """

    reversible_reaction = []
    for reaction in reaction_list:
        reaction_file = REST.kegg_get(reaction).read()
        for i in reaction_file.rstrip().split("\n"):
            if i.startswith("EQUATION") and "<=>" in i:
                reversible_reaction.append(reaction)
    return reversible_reaction

def combine_substrates_products(df: pd.DataFrame):
    """
    combine_substrates_products() is for use with a collection of enzymes
        in which it is understood that they are capable of catalyzing both the
        forward and reverser reactions. In this case, both the substrates and
        the products should be considered as bioreachable products.
        This function parses the list of substrates and products from their
        respective fields in the input dataframe, and returns a new dataframe
        with the combined substrates and products in a column labeled 'product'

    WARNING: combine_substrates_products() should not be run multiple times on
        the same dataframe becuase it will will append duplicate substrates

    Args:
        df (pandas.DataFrame): must contain fields
            ['entry', 'substrate', 'product']

    Returns:
        pandas.DataFrame: contains only fields ['entry', 'product']
    """

    rowindex = np.arange(0,len(df))
    df_with_ordered_index = df.set_index(rowindex)

    newdf = df_with_ordered_index
    # should this be a .copy()?

    for index,row in df_with_ordered_index.iterrows():
        productlist = row['product']
        substratelist = row['substrate']
        newdf.iloc[index,2] = productlist + substratelist

    return newdf[["entry","product"]]
