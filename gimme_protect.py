import cobra
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

num_genes = len(model.genes)
num_samples = 226

# import transcript table
print('importing transcript table...')
transcript_df = pandas.read_csv('../tpm_table_old.csv', index_col=0, header=0)
transcript_id = pandas.read_csv('../tpm_table.csv', index_col=0, header=0)

# make sure the index is string, not int, for integration with GIMME
print('converting transcript indices...')
new_index = ['gene'+str(i) for i in transcript_df.index.tolist()]
transcript_df.index = new_index

# read in the metadata which contains sample IDs
print('importing patient metadata...')
metadata = pandas.read_csv('../protect/SraRunTable.txt', index_col=0, header=0)

# transform the simulated counts using the rank-based approach
print('transforming simulated counts...')
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
print('normalizing transcripts...')
transcript_df = transcript_df.div(transcript_df.max(axis=1), axis=0)
transcript_df.fillna(0, inplace=True)
transcript_df.reindex(transcript_id.index)


# create the expression profile object using Driven
print('generating expression profile...')
exp_prof = ExpressionProfile(identifiers=transcript_df.index.values,
                             conditions=transcript_df.columns.values,
                             expression=transcript_df.values)

# rename the genes to have a non-numeric character so they don't get evaluated prematurely by driven.
# This is a problem with the AST parser that cobrapy uses to represent gene identifiers.
print('renaming genes...')
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

print('calling gimme function...')
gimme_solutions = {}
with model:
    for sample in transcript_df.columns:
        constrained_model,gimme_solution = gimme_mod(model,
                                                      exp_prof,
                                                      condition=sample,
                                                      cutoff = 1.0,
                                                      fraction_of_optimum = 0.1,
                                                      max_penalty=1.0)
        print(str(gimme_solution.objective_value))
        print(gimme_solution.fluxes)
        print(constrained_model)
        print(sample)
        gimme_solutions[sample] = gimme_solution

# write output to csv
print('writing fluxes to csv...')

count = 0
while count < num_samples:
    gimme_solutions[transcript_df.columns[count]].fluxes.to_csv('data/fluxes/patient_'+str(count+1)+'.tsv', sep='\t')

    with open('data/consistency.txt', 'a+') as consistency_scores:
        consistency_scores.write('patient_'+str(count+1)+':\t'+str(gimme_solution.objective_value)+'\n')

    count = count + 1

