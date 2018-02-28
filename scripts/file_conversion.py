"""
Contains functions for file conversions
"""

import xlrd
import csv
import os

# define locations of useful files
dir = os.path.dirname(__file__)
outputdir = os.path.join(dir, '../data/output/')

class FileConversionError(Exception):
    pass


def excel_to_tsv(excel_filepath, sheet_num=0):
    """Creates a tsv file containing the same information as the input
    xlsx file. Deposits tsv in same directory as xlsx file.
    Note: only converts first sheet of Excel file.

    Arguments:
        excel_filepath (str):
            full path to excel file to be converted
        sheet_num (int):
            which sheet of the Excel file to convert
            default = 0 (first)

    Returns:
        tsv_filepath (str)
    """

    # assert proper file format and location
    if not excel_filepath.endswith(".xlsx"):
        raise FileConversionError("excel_filepath not in .xlsx format!")

    if not os.path.isfile(excel_filepath):
        raise FileConversionError("excel_filepath does not exist!")

    tsv_filepath = excel_filepath.replace(".xlsx", "_sheet_{}.tsv".format(str(sheet_num)))

    # open the output tsv
    with open(tsv_filepath, 'w') as tsv_file:
        # define a writer
        wr = csv.writer(tsv_file, delimiter="\t")

        # open the xlsx file
        excel_file = xlrd.open_workbook(excel_filepath)

        # get chosen sheet
        sheet = excel_file.sheet_by_index(sheet_num)

        # write the rows
        for rownum in range(sheet.nrows):
            wr.writerow(sheet.row_values(rownum))

    return tsv_filepath


# todo: add csv to tsv?

# todo: consider moving the following functions to a more appropriate file
def safelist(listable):
    """Returns a list object derived from the list-like object input.
    Note: list() or pd.object.tolist() both fail when the object read in is of length 1.
    """
    if type(listable) == str:
        return [listable]
    else:
        return listable.tolist()


def sif_to_series(sif_fl):
    """Converts SIF file to pd.Series

    Arguments:
        sif_fl (str)
            string indicating full path to SIF file.
    Returns:
        sif_series (pandas Series of lists)
            A reference sif in the format of a pandas Series.
            Each key in the series is a producer in the reference sif.
            Each value is a list of metabolites produced by that producer."""
    sif_df = pd.read_csv(sif_fl, sep='\t', header=None, index_col=0)
    uniq_producers = np.unique(sif_df.index)
    sif_series = pd.Series({producer: safelist(sif_df.ix[producer][2])
                            for producer in uniq_producers})
    return sif_series


def write_out_sif_series(sif_series, output_name, location=outputdir):
    """Writes a sif_series out to .sif file.

    Arguments:
        sif_series (pandas Series of lists)
            A reference sif in the format of a pandas Series.
            Each key in the series is a producer in the reference sif.
            Each value is a list of metabolites produced by that producer.

    """
    with open(location + output_name + ".sif", "w+") as out_fh:
        for producer in sif_series.index:
            for produced in sif_series[producer]:
                out_fh.write("{}\tused-to-produce\t{}".format(producer, produced))


