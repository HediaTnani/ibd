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


def gimme_thread(columns):
    # import human genome model
    # print('importing human genome model...')
    model = cobra.io.read_sbml_model('data/Human-GEM.xml')

    num_genes = len(model.genes)
    num_samples = 226

    # import transcript table
    # print('importing transcript table...')
    transcript_df = pandas.read_csv('../id_tpm.csv', index_col=0, header=0)
    # transcript_id = pandas.read_csv('../tpm_table.csv', index_col=0, header=0)

    # make sure the index is string, not int, for integration with GIMME
    # print('converting transcript indices...')
    # new_index = ['gene'+str(i) for i in transcript_df.index.tolist()]
    # transcript_df.index = new_index

    # read in the metadata which contains sample IDs
    # print('importing patient metadata...')
    metadata = pandas.read_csv('data/SraRunTable.txt', index_col=0, header=0)

    # transform the simulated counts using the rank-based approach
    # print('transforming simulated counts...')
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
    # print('normalizing transcripts...')
    transcript_df = transcript_df.div(transcript_df.max(axis=1), axis=0)
    transcript_df.fillna(0, inplace=True)
    # transcript_df.reindex(transcript_id.index)

    # print(transcript_df.columns)


    # create the expression profile object using Driven
    # print('generating expression profile...')
    exp_prof = ExpressionProfile(identifiers=transcript_df.index.values,
                                 conditions=transcript_df.columns.values,
                                 expression=transcript_df.values)

    # print(Symbol(model.genes[1].id))
    # print(exp_prof.to_reaction_dict("SRR6467548", model))

    # rename the genes to have a non-numeric character so they don't get evaluated prematurely by driven.
    # This is a problem with the AST parser that cobrapy uses to represent gene identifiers.
    # print('renaming genes...')
    # rename_dict = {gene.id:'gene'+gene.id for gene in model.genes}
    # cobra.manipulation.modify.rename_genes(model,rename_dict)
    # model.repair()
    # # update the GPRs with the 'gene' prefix that we added to the gene IDs
    # for reaction in model.reactions:
    #     old_rule = reaction.gene_reaction_rule
    #     split_by_white = old_rule.split(' ')
    #     new_rule = ''
    #     for entry in split_by_white:
    #         if entry.isdigit():
    #             new_rule += ' gene'+entry
    #         else:
    #             new_rule += ' ' + entry
    #     reaction.gene_reaction_rule = new_rule

    gimme_solutions = {}



    with model:
        for sample in transcript_df.columns[columns]:
            print('calling gimme on '+str(sample)+'...')
            constrained_model,gimme_solution,condition = gimme_mod(model,
                                                          exp_prof,
                                                          condition=sample,
                                                          cutoff = 1.0,
                                                          fraction_of_optimum = 0.1,
                                                          max_penalty=1.0)

            gimme_solution.fluxes.to_csv('data/fluxes/'+str(sample)+'.tsv', sep='\t')
            with open('data/consistency.txt', 'a+') as consistency_scores:
                consistency_scores.write(str(sample)+':\t'+str(gimme_solution.objective_value)+'\n')
            
            # print(coefficients)
            gimme_solutions[sample] = gimme_solution


if __name__ == "__main__":

    p1 = Process(target=gimme_thread, args=(numpy.linspace(0,13,14)))
    p2 = Process(target=gimme_thread, args=(numpy.linspace(14,27,14)))
    p3 = Process(target=gimme_thread, args=(numpy.linspace(28,41,14)))
    p4 = Process(target=gimme_thread, args=(numpy.linspace(42,55,14)))
    p5 = Process(target=gimme_thread, args=(numpy.linspace(56,69,14)))
    p6 = Process(target=gimme_thread, args=(numpy.linspace(70,83,14)))
    p7 = Process(target=gimme_thread, args=(numpy.linspace(84,97,14)))
    p8 = Process(target=gimme_thread, args=(numpy.linspace(98,111,14)))
    p9 = Process(target=gimme_thread, args=(numpy.linspace(112,125,14)))
    p10 = Process(target=gimme_thread, args=(numpy.linspace(126,139,14)))
    p11 = Process(target=gimme_thread, args=(numpy.linspace(140,153,14)))
    p12 = Process(target=gimme_thread, args=(numpy.linspace(154,167,14)))
    p13 = Process(target=gimme_thread, args=(numpy.linspace(168,181,14)))
    p14 = Process(target=gimme_thread, args=(numpy.linspace(182,195,14)))
    p15 = Process(target=gimme_thread, args=(numpy.linspace(196,209,14)))
    p16 = Process(target=gimme_thread, args=(numpy.linspace(210,225,14)))

    p1.start()
    p2.start()
    p3.start()
    p4.start()
    p5.start()
    p6.start()
    p7.start()
    p8.start()
    p9.start()
    p10.start()
    p11.start()
    p12.start()
    p13.start()
    p14.start()
    p15.start()
    p16.start()

    p1.join()
    p2.join()
    p3.join()
    p4.join()
    p5.join()
    p6.join()
    p7.join()
    p8.join()
    p9.join()
    p10.join()
    p11.join()
    p12.join()
    p13.join()
    p14.join()
    p15.join()
    p16.join()

    print("Done!")

# write output to csv
# print('writing fluxes to csv...')

# count = 0
# while count < num_samples:
#     gimme_solutions[transcript_df.columns[count]].fluxes.to_csv('data/fluxes/patient_'+str(count+1)+'.tsv', sep='\t')

#     with open('data/consistency.txt', 'a+') as consistency_scores:
#         consistency_scores.write('patient_'+str(count+1)+':\t'+str(gimme_solution.objective_value)+'\n')

#     count = count + 1

