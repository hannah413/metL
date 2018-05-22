"""Functions for formatting SIF visualizations
"""


def specify_chebi_formatting(sif_path, mwu_df, mwu_p, corr_df, fc_series):
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

    color: The background color of the node. Colors are given in RGB.
    textcolor: Color of the node text (name)
    bordercolor: Border color
    borderwidth: Border width
    tooltip: Sets the text that will show up when the mouse is over the node.

    The .format file is tab delimited text file. Example:
    node CXCR4 color 23 65 13
    node RTN4 textcolor 0 0 0
    node IL6ST bordercolor 180 23 14
    node CARD9 borderwidth 2
    node IRS1 tooltip
    """

    format_path = sif_path.replace('sif', 'format')
    format_fh = open(format_path, "w")

    sif_fh = open(sif_path, "r")
    sif = sif_fh.readlines()

    # get a list of all chebis appearing in the SIF to be visualized
    sif_chebis = []
    for relationship in sif:
        relationship = relationship.strip('\n')
        relationship_members = relationship.split('\t')
        if relationship_members[0] not in sif_chebis:
            sif_chebis.append(relationship_members[0])
        if relationship_members[2] not in sif_chebis:
            sif_chebis.append(relationship_members[2])

    # get significantly altered metabs
    sig_chebi_ids = mwu_df[mwu_df['pval'] <= mwu_p].index.tolist()
    present_sig_chebi_ids = []
    for ch_id in sig_chebi_ids:
        if ch_id in sif_chebis:
            present_sig_chebi_ids.append(ch_id)

    # color the borders of significantly altered metabs red
    for ch_id in present_sig_chebi_ids:
        format_fh.write("node\t" + ch_id + "\tbordercolor\t255 0 0\n")
        format_fh.write("node\t" + ch_id + "\tborderwidth\t3\n")

    # color remaining chebis
    for metab in sif_chebis:
        fc = fc_series[metab]
        rgb = calculate_node_color_by_fc(fc)
        rgb_str = "{} {} {}\n".format(rgb[0], rgb[1], rgb[2])
        format_fh.write("node\t" + metab + "\tcolor\t" + rgb_str)

    format_fh.close()


def calculate_node_color_by_fc(fold_change):
    """
    Calculates node color based on fold change.
    Returns node color in RGB format (i.e. np.array([100, 150, 100]))
    Assumes color is RGB between [0, 0, 0] and [255, 255, 255]
    """

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

    # TODO: Deal with (normalize?) fold changes not being out of 1 or -1...
    # todo: fix these ranges and calculations; a fold change of 1 should indicate *no change*
    # if the spearman correlation is negative, set the node to a gradation of blue
    if -1.0 <= fold_change < 0.9:
        rgb = blue + blue_vector * (1 + fold_change)
        rgb = [int(i) for i in rgb.tolist()]

    if 0.9 <= fold_change <= 1.1:
        rgb = gray

    # if the spearman correlation is positive, set the node to a gradation of orange
    if 1.1 < fold_change <= 2.0:
        rgb = orange + orange_vector * (1 - fold_change)
        rgb = [int(i) for i in rgb.tolist()]

    if fold_change > 2.0:
        rgb = red

    return rgb