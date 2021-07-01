import cobra
import numpy
import pandas

from multiprocessing import Process

from sympy import Min, Max, Add, Mul, Symbol
from sympy.parsing.ast_parser import parse_expr

import driven
from driven.data_sets import ExpressionProfile
from optlang.symbolics import Zero
from six import iteritems
from gimme import gimme_mod

import os

# Get the current working directory
cwd = os.getcwd()

# Print the current working directory
print("Current working directory: {0}".format(cwd))


# import human genome model
print('importing human genome model...')
#model = cobra.io.read_sbml_model('data/HUMAN-GEM.xml')

num_genes = len(model.genes)
num_samples = 226

# import transcript table
print('importing transcript table...')
transcript_df = pandas.read_csv('../id_tpm.csv', index_col=0, header=0)
# transcript_id = pandas.read_csv('../tpm_table.csv', index_col=0, header=0)

# make sure the index is string, not int, for integration with GIMME
# print('converting transcript indices...')
# new_index = ['gene'+str(i) for i in transcript_df.index.tolist()]
# transcript_df.index = new_index

# read in the metadata which contains sample IDs
print('importing patient metadata...')
metadata = pandas.read_csv('data/SraRunTable.txt', index_col=0, header=0)