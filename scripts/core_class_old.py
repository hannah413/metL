"""Defines classes for organizing and analyzing the user's assay results"""

from scripts.metab_mapping import *
import pandas as pd
import numpy as np
import os

# define locations of useful files
dir = os.path.dirname(__file__)
inputdir = os.path.join(dir,'../data/input/')
compounds_fl = inputdir + 'compounds.tsv'

# todo: make an option where users provide ChEBI IDs instead of whatever they want

class UserAssayData(object):
    # todo: do attributes refer only to those in the __init__(X)?
    """Base class for metabolomic data

    Attributes:
        metabs (frozenset)
            user provided (human-readable) metabolites
        metabs_map (dict)
            key (str): i for i in metabs
            value (list of str): ChEBI IDs that map to key
        expl_groups (dict)
            key (str): i for i in experimental groups
            value (list of str): member IDs
                values should match colnames in user-provided assay data.
        input_data_path (str)
            Path to input data file(s)
        assay_results (pd.DataFrame)
            Results at location = input_data_path
        map_exact (bool)
            Boolean determining whether exact or approximate ChEBI ID matches will be used
    """

    def __init__(self, expl_groups, input_data_path, map_exact):

        self.input_data_path = input_data_path
        self.assay_results = pd.read_csv(input_data_path,
                                         sep='\t',
                                         header=0,
                                         index_col=0,
                                         na_values='nd')

        # todo: how to get this dict of info from the user?
        self.expl_groups = pass

        self.metabs = frozenset(self.assay_results.index)
        if not metabs:
            raise UserAssayDataError("The user must provide names of measured metabolites\n"
                                     "as first column of assay_results file!")

        metabs_chebi, failed_to_map = get_chebi_matches(self.metabs, map_exact)
        if not metabs_chebi:
            raise UserAssayDataError("Failed to map metabolite names to ChEBI IDs")

        # todo: calculate spearman correlations for every pair?
        # self.filtered_assay_results = filtering.filter_spearman(assay_results)


class UserAssayDataError(Exception):
    """Defines exceptions caused by incorrectly-formatted input file"""
    pass

class SimpleAssayData(UserAssayData):
    pass

class ComplexAssayData(UserAssayData):
    pass

class TimeCourseData(object):
    pass

class PerturbationData(object):
    pass

class SimpleTimeCourseData(SimpleAssayData, TimeCourseData):
    pass

class ComplexTimeCourseData(ComplexAssayData, TimeCourseData):
    pass

class SimplePerturbationData(SimpleAssayData, Perturbation):
    pass

class ComplexPerturbationData(ComplexAssayData, Perturbation):
    pass

class SimpleTimeCoursePerturbation

