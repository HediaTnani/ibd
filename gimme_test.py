import cobra
from cobra.test import create_test_model
import numpy
import pandas

import driven
from driven.data_sets import ExpressionProfile
from optlang.symbolics import Zero
from six import iteritems
from gimme import gimme_mod

# create a test model. By default, this is the E. coli core model.
model = create_test_model()

# create a dataframe containing simulated trancript abundances
num_genes = len(model.genes)
num_samples = 10
# driven expects a table of genes (rows) x samples (columns):
transcript_df = pandas.DataFrame(numpy.random.randint(low=0,high=100,size=(num_genes,num_samples)))
transcript_df.index = [gene.id for gene in model.genes]

# transform the simulated counts using the rank-based approach
for sample,subdata in transcript_df.groupby(by=transcript_df.columns,axis=1):
    # within the sample, sort transcripts by counts
    transcripts_ranked = subdata[sample].sort_values()
    # create a vector of ranks from 0 to the length of transcripts
    rank_vec = [x for x in range(0,len(transcripts_ranked))]
    rank_vec = pandas.Series(rank_vec)
    rank_vec.index = transcripts_ranked.index
    # convert the ranks to a percentile
    percentiles = (rank_vec + 1)/len(rank_vec)
    # reorder the percentiles to match the original dataframe
    percentiles = percentiles.reindex(transcript_df.index)
    # multiply the original counts by the percentiles
    transcript_df[sample] = transcript_df[sample].multiply(percentiles)

# Now normalize each transcript by dividing by the max percentile-normalized
# value for that transcript across all samples
transcript_df = transcript_df.div(transcript_df.max(axis=1), axis=0)

# create the expression profile object using Driven
exp_prof = ExpressionProfile(identifiers=transcript_df.index.values,
                             conditions=transcript_df.columns.values,
                             expression=transcript_df.values)

# rename the genes to have a non-numeric character so they don't get evaluated prematurely by driven.
# This is a problem with the AST parser that cobrapy uses to represent gene identifiers.
rename_dict = {gene.id:'gene'+gene.id for gene in model.genes}
cobra.manipulation.modify.rename_genes(model,rename_dict)
model.repair()
# update the GPRs with the 'gene' prefix that we added to the gene IDs
for reaction in model.reactions:
    old_rule = reaction.gene_reaction_rule
    split_by_white = old_rule.split(' ')
    new_rule = ''
    for entry in split_by_white:
        if entry.isdigit():
            new_rule += ' gene'+entry
        else:
            new_rule += ' ' + entry
    reaction.gene_reaction_rule = new_rule


gimme_solutions = {}
count = 0
with model:
    for sample in transcript_df.columns:
        constrained_model,gimme_solution = gimme_mod(model,
                                                      exp_prof,
                                                      condition=sample,
                                                      cutoff = 1.0,
                                                      fraction_of_optimum = 0.1,
                                                      max_penalty=1.0)
        #print(gimme_solution.fluxes)
        print(str(gimme_solution.objective_value))
        gimme_solutions[count] = gimme_solution
        count = count+1

# print the fluxes in the solution for the first sample
#print(gimme_solutions[transcript_df.columns[0]].fluxes)

#print(gimme_solutions[0].fluxes)
#print(gimme_solutions[1].fluxes)

#print('writing fluxes to csv...')
#gimme_solutions.fluxes.to_csv('data/fluxes/fluxes.csv', sep = '\t')



