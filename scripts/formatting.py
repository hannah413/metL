"""Functions for formatting SIF visualizations
"""


def specify_chebi_formatting(sif_path, mwu_df, fc_series):
    """Produces a .format file with the same pre-extension name as the .sif
    file to which it refers.

    Arguments
        sif_path (str)
            full path to the .sif file that will be formatted.
        mwu_df (pandas DataFrame)
            Dataframe of shape (n_metabs, ['stat', 'pval']]. Contains Mann
            Whitney U test statistics and p-values.
        fc_series (pandas Series of floats)
            Fold change values for each metabolite between experimental groups
    """

    format_path = sif_path.replace('sif', 'format')
    format_fh = open(format_path, "w")

    #todo: clean up old format_sif methods
    pass


def calculate_node_color_by_fc(fold_change):
    """
    Calculates node color based on fold change.
    Returns node color in RGB format (i.e. np.array([100, 150, 100]))
    Assumes color is RGB between [0, 0, 0] and [255, 255, 255]
    """
    # todo: fix--FC of 0-1.0 indicates reduction, should be diff color than 1.0+
    # TODO: Deal with (normalize?) fold changes not being out of 1 or -1...
    gray = [220, 220, 220]
    red = [255, 89, 0]
    darkblue = [0, 111, 255]

    white = np.array([255, 255, 255])
    blue = np.array([71, 151, 255])
    orange = np.array([255, 150, 94])


    blue_vector = white - blue
    orange_vector = white - orange

    if fold_change < -1.0:
        rgb = darkblue

    # if the spearman correlation is negative, set the node to a gradation of orange
    if -1.0 <= fold_change < -0.1:
        rgb = blue + blue_vector * (1 + fold_change)
        rgb = [int(i) for i in rgb.tolist()]

    if -0.1 <= fold_change <= 0.1:
        rgb = gray

    # if the spearman correlation is positive, set the node to a gradation of blue
    if 0.1 < fold_change <= 1.0:
        rgb = orange + orange_vector * (1 - fold_change)
        rgb = [int(i) for i in rgb.tolist()]

    if fold_change > 1.0:
        rgb = red

    return rgb