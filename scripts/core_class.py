"""Defines classes for organizing and analyzing the metabolomic abundance data.

Author: Hannah Manning
Date: February, 2018
"""

from scripts import file_conversion, filtering, metab_mapping, formatting
from scipy import stats as st
import pandas as pd
import itertools
import os

# define locations of useful files
dir = os.path.dirname(__file__)
inputdir = os.path.join(dir,'../data/input/')
outputdir = os.path.join(dir,'../data/output/')
compounds_fl = inputdir + 'compounds.tsv'
full_sif = inputdir + 'used-to-produce.sif'


# todo: expand and use exceptions
class VisualizerError(Exception):
    """Defines generic errors occuring in ComparativeVisualizer
    """
    pass


class UserInputDataError(Exception):
    """Defines errors in input data file type.
    """
    pass


class ComparativeVisualizer(object):
    """
    Object which visualizes differences in metabolomic abundance data
    between two experimental groups. The terms 'analytes' and 'metabolites'
    are used interchangeably here.

    Attributes:
        abundance_filepath (str)
            Full path to metabolomic abundance data.
        abundance_df (pandas.DataFrame)
            Dataframe composed of abundance_filepath data.
            shape = (n_analytes, n_samples)
        group1 & group2 (frozenset of str)
            Experimental groups used to stratify columns (by sample name)
            in abundance_df (i.e. frozenset("cell_line1", "cell_line2", ...))
        metabs (frozenset of str)
            metabolic analytes in abundance_df.
        mwu_p (float)
            Maximum Mann Whitney U p-value. Used as a threshold to filter out
            metabolites who are not significantly different between the two
            experimental groups.
        c1 (float)
            Spearman correlation threshold used to filter out edges between
            metabolites that directly interact.
            (i.e. between A and B in "A used-to-produce B").
        c2 (float)
            Spearman correlation threshold used to filter out edges between
            metabolites that have an intermediate actor.
            (i.e. between A and C in "A used-to-produce B used-to-produce C").
        identifier (str)
            User-provided string to be used in file titles.
        mwu_df (pandas.DataFrame)
            Dataframe of shape (n_metabs, ['stat', 'pval']]. Contains Mann
            Whitney U test statistics and p-values.
        corr_df (pandas.DataFrame)
            Dataframe of shape (n_metabs, n_metabs) containing Spearman rank
            order correlation coefficients and p-values among all pairs of
            metabolites. Each cell contains a tuple (coef, pval).

    """

    def __init__(self, abundance_filepath, group1, group2,
                 mwu_p, c1, c2, identifier):

        # todo: allow csv, excel, etc. (bring in file_conversion.py functions)
        if not abundance_filepath.endswith('.tsv'):
            raise UserInputDataError("abundance_filepath must be .tsv format")
        else:
            self.abundance_filepath = abundance_filepath
            self.abundance_df = pd.read_csv(abundance_filepath,
                                           sep='\t', index_col=0,
                                           header=0, na_values='nd')

        self.group1 = group1
        self.group2 = group2
        self.metabs = frozenset(self.abundance_df.index)

        self.mwu_p = mwu_p
        self.c1 = c1
        self.c2 = c2

        self.run_name = "{}_mwu_{}__c1_{}__c2_{}".format(identifier,
                                                         str(mwu_p), str(c1), str(c2))

    def _map_names(self):
        """Maps human-readable ID names to ChEBI IDs."""
        # todo: check in on Neil and ChEBI API;
        # if not updated by March, convert code to python2
        # for now only working with dataframes containing ChEBI IDs.
        pass


    # todo: add optional normalizing functions


    def _calc_signif(self):
        """Applies Mann Whitney U test to each metabolite (between
        two experimental groups). Tests likelihood that the two independent
        samples came from different populations and returns a dataframe of
        shape [n_metabs, [stat,pval]].
        """
        mwu_df = pd.DataFrame(index=self.metabs,columns=['stat', 'pval'])
        for metab in self.metabs:
            grp1_data = self.abundance_df.loc[metab][list(self.group1)]
            grp2_data = self.abundance_df.loc[metab][list(self.group2)]
            mwu = st.mannwhitneyu(grp1_data,
                                  grp2_data)
            mwu_df.loc[metab] = [mwu.statistic, mwu.pvalue]
        self.mwu_df = mwu_df


    def _calc_foldchange(self):
        """"""
        # initialize empty df
        fc_df = pd.DataFrame(index=list(self.metabs),
                             columns=["grp1_mean", "grp2_mean", "fold_change"])

        # first calc means for each metabolite in each group
        fc_df['grp1_mean'] = self.abundance_df.loc[:][list(self.group1)].mean(axis=1)
        fc_df['grp2_mean'] = self.abundance_df.loc[:][list(self.group2)].mean(axis=1)

        # then calc fold_change between them
        fc_df['fold_change'] = fc_df['grp1_mean'] / fc_df['grp2_mean'] - 1
        self.fc_df = fc_df

    # todo: consider using pearson instead
    def _correlate(self):
        """Produces Spearman rank-order correlations among all pairs of metabolites.
        """

        # derive all possible pairs of metabolites
        pairs = itertools.combinations(self.metabs, 2)

        # initialize an empty dataframe of metabs x metabs
        corr_df = pd.DataFrame(index=list(self.metabs),
                               columns=list(self.metabs))

        for pair in pairs:
            metab1 = self.abundance_df.loc[pair[0]]
            metab2 = self.abundance_df.loc[pair[1]]

            # calculate spearman's correlation
            spearmanr = st.spearmanr(metab1, metab2, nan_policy='omit')
            corr_df[pair[0], pair[1]] = (spearmanr.correlation, spearmanr.pvalue)
        self.corr_df = corr_df


    def _build_ggm(self):
        """Produces partial correlations via a Gaussian Graphical Model (GGM).
        """

        pass


    def _filter(self):
        """Uses Mann Whitney U and Spearman values calculated in _calc_signif()
        and _correlate() to remove unwanted edges. Writes SIF to file.
        """

        # read the used-to-produce sif in as a series
        # format: {producerA: [produced1, produced2, ...]}
        ref_sif = file_conversion.sif_to_series(full_sif)

        # build both_sif
        self.both_sif = filtering.keep_if_both(sif_series=ref_sif,
                                               mwu_df=self.mwu_df, mwu_p=self.mwu_p,
                                               corr_df=self.corr_df, c1=self.c1,
                                               run_name=self.run_name)

        # no correlations possible in "either" sif (analytes missing)
        # builds on both_sif
        self.either_sif = filtering.keep_if_either(sif_series=ref_sif,
                                                   mwu_df=self.mwu_df, mwu_p=self.mwu_p,
                                                   both_sif_series=self.both_sif)

        # correlations possible again in dist2
        # builds on either_sif
        self.dist2_sif = filtering.keep_dist2(sif_series=ref_sif,
                                                  mwu_df=self.mwu_df, mwu_p=self.mwu_p,
                                                  corr_df=self.corr_df, c1=self.c1, c2=self.c2,
                                                  either_sif_series=self.either_sif)


    def _format(self):
        """Builds a format file for the output SIF."""

        formatting.specify_chebi_formatting(sif_path=outputdir + run_name + "_BOTH.sif",
                                 mwu_df=self.mwu_df, mwu_p = self.mwu_p,
                                 corr_df=self.corr_df, fc_series=self.fc_df['fold_change'])

        formatting.specify_chebi_formatting(sif_path=outputdir + run_name + "_EITHER.sif",
                                 mwu_df=self.mwu_df, mwu_p=self.mwu_p,
                                 corr_df=self.corr_df, fc_series=self.fc_df['fold_change'])

        formatting.specify_chebi_formatting(sif_path=outputdir + run_name + "_DIST2.sif",
                                 mwu_df=self.mwu_df, mwu_p=self.mwu_p,
                                 corr_df=self.corr_df, fc_series=self.fc_df['fold_change'])


    def run(self):
        self._map_names()
        self._calc_signif()
        self._calc_foldchange()
        self._correlate()
        self._filter()
        self._format()


    def visualize(self):
        """Pass SIF and format file to cytoscape here"""
        pass


class OneToOneVisualizer(ComparativeVisualizer):
    """Visualizes metabolomic abundance differences between two individual samples.
    """
    pass


class OneToManyVisualizer(ComparativeVisualizer):
    """Visualizes metabolomic abundance differences between one sample of interest and a
    given background.
    """
    pass


class TimeCourseVisualizer(ComparativeVisualizer):
    """Visualizes metabolomic abundance time course data with one SIF and N format files.
    """
    pass


class OneToOneFCVisualizer(OneToOneVisualizer):
    """Visualizes metabolic abundance differences between two groups whose fold change is
    given in each cell of the input dataframe
    """
    pass



