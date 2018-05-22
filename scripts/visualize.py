from scripts import core_class
import argparse


def main():
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
    V = core_class.ComparativeVisualizer(abundance_filepath = assay_file,
                              group1 = group1, group2 = group2,
                              c1 = min_corr_both, c2 = min_corr_dist2, mwu_p=max_p)
    V.run()
    V.visualize()

if __name__ == "__main__":
    main()