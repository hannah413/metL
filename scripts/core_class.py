"""Defines classes for organizing and analyzing the metabolomic abundance data.

Author: Hannah Manning
Date: February, 2018
"""

from scripts.metab_mapping import *
from scripts.file_conversion import *
from scripts.filtering immport *
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


class ComparativeVisualizer(object):
    """
    Object which visualizes differences in metabolomic abundance data
    between two experimental groups. The terms 'analytes' and 'metabolites
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


    # todo: add optional normalizing function(s)


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
        fc_df['grp1_mean'] = self.abundance_df.loc[:][list(group1)].mean(axis=1)
        fc_df['grp2_mean'] = self.abundance_df.loc[:][list(group2)].mean(axis=1)

        # then calc fold_change between them
        fc_df['fold_change'] = fc_df['grp1_mean'] / fc_df['grp2_mean'] - 1
        self.fc_df = fc_df


    def _correlate(self):
        """Produces Spearman rank-order correlations among all pairs of metabolites.
        """

        # derive all possible pairs of metabolites
        pairs = itertools.combinations(metabs, 2)

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

        specify_chebi_formatting(sif_path=outputdir + run_name + "_BOTH.sif",
                                 mwu_df=self.mwu_df, self.fc_df['fold_change'])

        specify_chebi_formatting(sif_path=outputdir + run_name + "_EITHER.sif",
                                 mwu_df=self.mwu_df, self.fc_df['fold_change'])

        specify_chebi_formatting(sif_path=outputdir + run_name + "_DIST2.sif",
                                 mwu_df=self.mwu_df, self.fc_df['fold_change'])


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


# todo: expand and use exceptions
class UserAssayDataError(Exception):
    """Defines exceptions caused by incorrectly-formatted input file"""
    pass

if __name__ == "__main__":

    # take user input
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--assay_file",
                        type=str,
                        default= inputdir + 'chebi_dummy_results.tsv',
                        help="Name of assay results file in input/ dir\n"
                             "Shape should be [n_analytes, n_samples]\n"
                             "(default = 'chebi_dummy_results.tsv')")

    parser.add_argument("-g1", "--group1",
                        type=list,
                        default=['COLO 320 DM', 'MDA-MB-157', 'NOMO-1',
                                 'SNU-1', 'RPMI8226', 'COLO-704'],
                        help="Sample IDs for one experimental group in input file\n"
                             "(default is associated with chebi_dummy_results.tsv)")

    parser.add_argument("-g2", "--group2",
                        type=list,
                        default=['MV-4-11', 'SUDHL8', 'P12',
                                 'A549', 'HPB-ALL', 'ALL-SIL'],
                        help="Sample IDs for one experimental group in input file\n"
                             "(default is associated with chebi_dummy_results.tsv)")

    parser.add_argument("-c1", "--minimum_corr_dist1",
                        type=float,
                        default=0.2,
                        help="Min absolute value of spearman correlation\n"
                             "between two nodes in order for them\n"
                             "to be kept in 'both' relationships\n"
                             "default = 0.5")

    parser.add_argument("-c2", "--minimum_corr_dist2",
                        type=float,
                        default=0.5,
                        help="Min absolute value of spearman correlation\n"
                             "between two nodes in order for them\n"
                             "to be kept in 'distance of 2' relationships\n"
                             "default = 0.5")

    parser.add_argument("-p", "--maximum_mwu_p",
                        type=float,
                        default=0.1,
                        help="The maximum p value allowed for each metabolite\n"
                             "(mann whitney u test between sensitive and resistant)\n"
                             "(default = 0.10)")

    parser.add_argument("-i", "--identifier",
                        type=str,
                        default="UNNAMED",
                        help="")

    # todo: implement
    parser.add_argument("-ll", "--linker_lenience",
                        type=int,
                        default=4,
                        help="Max number of .sif relationships in which an unprofiled\n"
                             "'linker' node may be present in order to be included\n"
                             "in the distance-of-2 network\n"
                             "default = 4")

    # todo: connect to excel & csv conversions in file_conversion.py
    parser.add_argument("-ft", "--filetype",
                        type=str,
                        default=".tsv")

    # parse user input data
    args = parser.parse_args()
    assay_file = args.assay_file
    group1 = args.group1
    group2 = args.group2

    # parse user visualization specifications
    min_corr_both = args.minimum_corr_dist1
    min_corr_dist2 = args.minimum_corr_dist2
    max_p = args.maximum_mwu_p

    # todo: incorporate linker_lenience into this version
    ll = args.linker_lenience

    # instantiate, run, and visualize
    V = ComparativeVisualizer(abundance_filepath = assay_file,
                              group1 = group1, group2 = group2,
                              c1 = min_corr_both, c2 = min_corr_dist2, mwu_p=max_p)
    V.run()
    V.visualize()

