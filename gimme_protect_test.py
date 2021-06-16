import cobra
from cobra.test import create_test_model
import numpy
import pandas

import driven
from driven.data_sets import ExpressionProfile
from optlang.symbolics import Zero
from six import iteritems
from gimme import gimme_mod

# import human genome model
print('importing human genome model...')
model = cobra.io.read_sbml_model('data/HUMAN-GEM.xml')

print(model.genes)

num_genes = len(model.genes)
num_samples = 226

# import transcript table
print('importing transcript table...')
transcript_df = pandas.read_csv('../tpm_table.csv', index_col=0, header=0)

# make sure the index is string, not int, for integration with GIMME
print('converting transcript indices...')
new_index = ['gene'+str(i) for i in transcript_df.index.tolist()]
transcript_df.index = new_index

# read in the metadata which contains sample IDs
print('importing patient metadata...')
metadata = pandas.read_csv('../protect/SraRunTable.txt', index_col=0, header=0)

# transform the simulated counts using the rank-based approach
# print('transforming simulated counts...')
# for sample,subdata in transcript_df.groupby(by=transcript_df.columns,axis=1):
#     # within the sample, sort transcripts by counts
#     transcripts_ranked = subdata[sample].sort_values()
#     # create a vector of ranks from 0 to the length of transcripts
#     rank_vec = [x for x in range(0,len(transcripts_ranked))]
#     rank_vec = pandas.Series(rank_vec)
#     rank_vec.index = transcripts_ranked.index
#     # convert the ranks to a percentile
#     percentiles = (rank_vec + 1)/len(rank_vec)
#     # reorder the percentiles to match the original dataframe
#     percentiles = percentiles.reindex(transcript_df.index)
#     # multiply the original counts by the percentiles
#     transcript_df[sample] = transcript_df[sample].multiply(percentiles)
#     if sample == 'V1':
#     	print(transcripts_ranked)
#     	print(rank_vec)
#     	print(percentiles)
#     	print(transcript_df)
#     else:
#     	pass

# Now normalize each transcript by dividing by the max percentile-normalized
# value for that transcript across all samples
# print('normalizing transcripts...')
# transcript_df = transcript_df.div(transcript_df.max(axis=1), axis=0)
# transcript_df.fillna(0, inplace=True)
# print(transcript_df)

# print(transcript_df.index.values)
# print(transcript_df.columns.values)
# print(transcript_df.values)