import gzip
import sklearn
import Bio
import rdkit
import re
import scipy as sp
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pubchempy as pc
from Bio.KEGG import Compound
from Bio.KEGG import REST
from Bio.KEGG import Enzyme
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.EState import Fingerprinter
from rdkit.Chem import Descriptors
from rdkit.Chem import rdFMCS
from rdkit.Chem.rdmolops import RDKFingerprint
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
from rdkit.Avalon.pyAvalonTools import GetAvalonFP
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import MolWt
from sklearn import linear_model
from sklearn.model_selection import train_test_split
