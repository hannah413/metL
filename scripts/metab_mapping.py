"""Functions for mapping human-readable metabolite names to ChEBI IDs
as well as to different conformations/ions/etc. of same metabolite

NOTE: MUST USE PYTHON2 UNTIL NEIL FIGURES OUT COMPATIBILITY :(
USE CONDA ENV METABP27
For interactive work, type: python2.7 -m IPython

"""

import libchebipy as lc
import pandas as pd
import os

# define locations of useful files
dir = os.path.dirname(__file__)
inputdir = os.path.join(dir,'../data/input/')
compounds_fl = inputdir + 'compounds.tsv'

# todo: follow edges in ChEBI to retrieve approximate matches
# todo: with option to write out if user is curious about which ChEBIs were used
# todo: NOTE: metabs received here should be derived like this: metabs = assay_results.inde.dropna().tolist()
def get_chebi_matches(metabs, exact=True):
    """Queries ChEBI API for exact or approximate ChEBI ID matches to provided list
    of metabolites. Note: capitalization is irrelevant, even if exact==True.

    Set exact to False if more inclusive search is desired.
        Uses ChEBI API's non-exact search. Not perfect.
        i.e. includes residues and derivatives related to compound
            i.e. methionine sulfone (CHEBI:132188) and methionine sulfoximine (CHEBI:47833)
            will be returned if queried with methionine (CHEBI:16811)

    Returns
        metab_map (dict):
            key: user-provided metab name
            value: list (of any length) of derived ChEBI IDs
        failed_to_map (list):
            user provided metabolite names for which there were no exact or approximate
            matches in ChEBI
    """
    metab_map = {m: set() for m in metabs}
    failed_to_map = set()

    print("Communicating with ChEBI API\n"
          "This will take a few moments")
    for name in metabs:
        # todo: follow up with Neil re: incompatibility with python3 :(
        # todo: try chebi_obj.getId() and getName() when up and running
        chebi_ID_obj = lc.search(name, exact=exact)
        if len(chebi_ID_obj) > 0:
            [metab_map[name].add('CHEBI:' +
                                    str(chebi_ID_obj[i]._ChebiEntity__chebi_id))
                                    for i in range(len(chebi_ID_obj))]
        elif len(chebi_ID_obj) == 0:
            del metab_map[name]
            failed_to_map.add(name)

    # scavenge for any left-behind metabs for which we have ChEBI matches in the
    # PathwayCommons compounds.tsv
    addl_chebis, failed_to_map = search_exact_matches_compoundstsv(failed_to_map)
    # this kind of dictionary unpacking to combine the 2 only works in python3 :(
    metab_map = {**metab_map, **addl_chebis}

    return metab_map, failed_to_map


def search_compounds_tsv(metabs):
    """Find exact ChEBI ID matches for metabolite names in locally stored
    Pathway Commons compounds.tsv

    Generally used for metabolites that failed to map using the ChEBI API

    Returns
        metab_map (dict):
            key: user-provided metab name
            value: set (of any length) of derived ChEBI IDs
        failed_to_map (set):
            user provided metabolite names for which there were no exact or approximate
            ChEBI matches in compounds.tsv
    """
    compounds_df = pd.read_csv(compounds_fl,
                               sep='\t',
                               usecols=[2,5],
                               header=0,
                               index_col=0)

    # build a dict with metab names as keys and list of corresponding chebis as values
    # note: value doesn't have to be list as it's only ever length of 1, but it is a list here
    # to be consistent with metab_map in other functions
    metab_map = {metab: [compounds_df[compounds_df['NAME'] == metab].index[0]] for metab in
                set(metabs) & set(compounds_df['NAME'])}

    failed_to_map = set(metabs) - set(metab_map.keys())

    return metab_map, failed_to_map