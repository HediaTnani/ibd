import cobra
from cobra.test import create_test_model
import numpy
import pandas
import sympy

import driven
from driven.data_sets import ExpressionProfile
from optlang.symbolics import Zero
from six import iteritems

# installation instructions:
# in a virtual environment (e.g., conda), install cobrapy via pip:
# pip install cobra
#
# Then install the devel branch of my fork of Driven:
# pip install -e git+https://github.com/gregmedlock/driven@devel#egg=driven


# This is a modified version of the standard GIMME algorithm from driven.
# In this version, genes that don't have transcripts detected in the data
# are assigned the max flux penalty. When the original algorithm was
# published, most transcriptomics data was being generated from microarrays
# that use a pre-constructed set of nucleotide probes. If a gene did not have
# a probe on the microarray, there was no way to tell if transcripts were
# present or not, so the authors decided to have zero penalty for those genes.
def gimme_mod(model, expression_profile, cutoff, fraction_of_optimum=0.9,
          condition=0, max_penalty=1):
    r"""
    Modified version of the GIMME algorithm which applies the maximum
    flux penalty to reactions that have no associated transcripts in the
    dataset.
    
    Parameters
    ----------
    model: cobra.Model
        A constraint-based model to perform GIMME on.
    expression_profile: ExpressionProfile
        An expression profile to integrate in the model.
    cutoff: float
        The cutoff value to be defined by the user.
    fraction_of_optimum: float
        The fraction of the Required Metabolic Functionalities.
    condition: str or int, optional (default 0)
        The condition from the expression profile.
        If None, the first condition is used.
    max_penalty: float
        The maximum penalty possible given the users preprocessing
        of transcriptomics data. This penalty will be applied to all
        reactions without any associated transcripts.
    Returns
    -------
    context-specific model: cobra.Model
    solution: cobra.Solution
    Notes
    -----
    The formulation for obtaining the Inconsistency Score is given below:
    minimize: \sum c_i * |v_i|
    s.t.    : Sv = 0
              a_i <= v_i <= b_i
    where   : c_i = {x_cutoff - x_i where x_cutoff > x_i
                     0 otherwise} for all i
    References
    ----------
    .. [1] Becker, S. and Palsson, B. O. (2008).
           Context-specific metabolic networks are consistent with experiments.
           PLoS Computational Biology, 4(5), e1000082.
           doi:10.1371/journal.pcbi.1000082
    """
    with model:

        # print('performing first optimization...')

        solution = model.slim_optimize() # returns the flux through the objective
        # print("1")
        prob = model.problem # extracts the optimization problem
        # print("2")
        rxn_profile = expression_profile.to_reaction_dict(condition, model)
        # print("3")

        if model.objective_direction == 'max':
            fix_obj_const = prob.Constraint(model.objective.expression,
                                            lb=fraction_of_optimum * solution,
                                            name="RMF")
        else:
            fix_obj_const = prob.Constraint(model.objective.expression,
                                            ub=fraction_of_optimum * solution,
                                            name="RMF")
        model.add_cons_vars(fix_obj_const)
        # print("4")

        # print('generating coefficients...')

        coefficients = {rxn_id: cutoff - expression
                        for rxn_id, expression in iteritems(rxn_profile)
                        if expression is not None and cutoff > expression}
        
        obj_vars = []
        for rxn_id, coefficient in iteritems(coefficients):
            rxn = model.reactions.get_by_id(rxn_id)
            obj_vars.append((rxn.forward_variable, coefficient))
            obj_vars.append((rxn.reverse_variable, coefficient))

        # print('adding penalties to reactions...')

        # Add the max penalty to all reactions in the model
        # that are not already in the rxn_profile (e.g., no expression)
        for reaction in model.reactions:
            if reaction.id not in rxn_profile.keys():
                obj_vars.append((reaction.forward_variable,max_penalty))
                obj_vars.append((reaction.reverse_variable,max_penalty))
            
        model.objective = prob.Objective(Zero, sloppy=True, direction="min")

        # print(model.objective)
        model.objective.set_linear_coefficients({v: c for v, c in obj_vars})

        # print(model.objective)

        # print('performing second optimization...')

        sol = model.optimize()

        return model, sol, coefficients

