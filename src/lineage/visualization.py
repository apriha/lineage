""" Chromosome plotting functions.

Notes
-----
Adapted from Ryan Dale's GitHub Gist for plotting chromosome features. [#Dale]_

References
----------
.. [#Dale] Ryan Dale, GitHub Gist,
   https://gist.github.com/daler/c98fc410282d7570efc3#file-ideograms-py

"""

"""
The MIT License (MIT)

Copyright (c) 2016 Ryan Dale

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

"""
MIT License

Copyright (c) 2017 Andrew Riha

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
import logging
import os

from atomicwrites import atomic_write
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import patches

logger = logging.getLogger(__name__)


def plot_chromosomes(one_chrom_match, two_chrom_match, cytobands, path, title, build):
    """Plots chromosomes with designated markers.

    Parameters
    ----------
    one_chrom_match : pandas.DataFrame
        segments to highlight on the chromosomes representing one shared chromosome
    two_chrom_match : pandas.DataFrame
        segments to highlight on the chromosomes representing two shared chromosomes
    cytobands : pandas.DataFrame
        cytobands table loaded with Resources
    path : str
        path to destination `.png` file
    title : str
        title for plot
    build : {37}
        human genome build
    """
    # Height of each chromosome
    chrom_height = 1.25

    # Spacing between consecutive chromosomes
    chrom_spacing = 1

    # Decide which chromosomes to use
    chromosome_list = ["chr%s" % i for i in range(1, 23)]
    chromosome_list.append("chrY")
    chromosome_list.append("chrX")

    # Keep track of the y positions for chromosomes, and the center of each chromosome
    # (which is where we'll put the ytick labels)
    ybase = 0
    chrom_ybase = {}
    chrom_centers = {}

    # Iterate in reverse so that items in the beginning of `chromosome_list` will
    # appear at the top of the plot
    for chrom in chromosome_list[::-1]:
        chrom_ybase[chrom] = ybase
        chrom_centers[chrom] = ybase + chrom_height / 2.0
        ybase += chrom_height + chrom_spacing

    # Colors for different chromosome stains
    color_lookup = {
        "gneg": (202 / 255, 202 / 255, 202 / 255),  # background
        "one_chrom": (0 / 255, 176 / 255, 240 / 255),
        "two_chrom": (66 / 255, 69 / 255, 121 / 255),
        "centromere": (1, 1, 1, 0.6),
    }

    df = _patch_chromosomal_features(cytobands, one_chrom_match, two_chrom_match)

    # Add a new column for colors
    df["colors"] = df["gie_stain"].apply(lambda x: color_lookup[x])

    # Width, height (in inches)
    figsize = (6.5, 9)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    # Now all we have to do is call our function for the chromosome data...
    for collection in _chromosome_collections(df, chrom_ybase, chrom_height):
        xranges, yrange, colors = collection
        ax.broken_barh(xranges, yrange, facecolors=colors)

    # Axes tweaking
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list)
    ax.margins(0.01)
    ax.axis("tight")

    handles = []

    # setup legend
    if len(one_chrom_match) > 0:
        one_chrom_patch = patches.Patch(
            color=color_lookup["one_chrom"], label="One chromosome shared"
        )
        handles.append(one_chrom_patch)

    if len(two_chrom_match) > 0:
        two_chrom_patch = patches.Patch(
            color=color_lookup["two_chrom"], label="Two chromosomes shared"
        )
        handles.append(two_chrom_patch)

    no_match_patch = patches.Patch(color=color_lookup["gneg"], label="No shared DNA")
    handles.append(no_match_patch)

    centromere_patch = patches.Patch(
        color=(234 / 255, 234 / 255, 234 / 255), label="Centromere"
    )
    handles.append(centromere_patch)

    plt.legend(handles=handles, loc="lower right", bbox_to_anchor=(0.95, 0.05))

    ax.set_title(title, fontsize=14, fontweight="bold")
    plt.xlabel("Build " + str(build) + " Chromosome Position", fontsize=10)
    logger.info("Saving {}".format(os.path.relpath(path)))
    plt.tight_layout()

    with atomic_write(path, mode="wb", overwrite=True) as f:
        plt.savefig(f)


def _chromosome_collections(df, y_positions, height, **kwargs):
    """

    Yields data for features that can be added to an Axes
    object using ax.broken_barh.

    Parameters
    ----------

    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.

    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the bars

    height : float
        Height of each bar

    Additional kwargs are passed to ax.broken_barh
    """
    del_width = False
    if "width" not in df.columns:
        del_width = True
        df["width"] = df["end"] - df["start"]
    for chrom, group in df.groupby("chrom"):
        yrange = (y_positions["chr" + chrom], height)
        xranges = group[["start", "width"]].values
        colors = group["colors"].values
        yield xranges, yrange, colors
    if del_width:
        del df["width"]


def _patch_chromosomal_features(cytobands, one_chrom_match, two_chrom_match):
    """Highlight positions for each chromosome segment / feature.

    Parameters
    ----------
    cytobands : pandas.DataFrame
        cytoband table from UCSC
    one_chrom_match : pandas.DataFrame
        segments to highlight on the chromosomes representing one shared chromosome
    two_chrom_match : pandas.DataFrame
        segments to highlight on the chromosomes representing two shared chromosomes

    Returns
    -------
    df : pandas.DataFrame
        the start and stop positions of particular features on each
        chromosome
    """

    def concat(df, chrom, start, end, gie_stain):
        return pd.concat(
            [
                df,
                pd.DataFrame(
                    {
                        "chrom": [chrom],
                        "start": [start],
                        "end": [end],
                        "gie_stain": [gie_stain],
                    }
                ),
            ],
            ignore_index=True,
        )

    chromosomes = cytobands["chrom"].unique()

    df = pd.DataFrame()

    for chromosome in chromosomes:
        chromosome_length = np.max(
            cytobands[cytobands["chrom"] == chromosome]["end"].values
        )

        # get all markers for this chromosome
        one_chrom_match_markers = one_chrom_match.loc[
            one_chrom_match["chrom"] == chromosome
        ]
        two_chrom_match_markers = two_chrom_match.loc[
            two_chrom_match["chrom"] == chromosome
        ]

        # background of chromosome
        df = concat(df, chromosome, 0, chromosome_length, "gneg")

        # add markers for shared DNA on one chromosome
        for marker in one_chrom_match_markers.itertuples():
            df = concat(df, chromosome, marker.start, marker.end, "one_chrom")

        # add markers for shared DNA on both chromosomes
        for marker in two_chrom_match_markers.itertuples():
            df = concat(df, chromosome, marker.start, marker.end, "two_chrom")

        # add centromeres
        for item in cytobands.loc[
            (cytobands["chrom"] == chromosome) & (cytobands["gie_stain"] == "acen")
        ].itertuples():
            df = concat(df, chromosome, item.start, item.end, "centromere")

    return df
