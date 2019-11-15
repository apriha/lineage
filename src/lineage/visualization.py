""" Chromosome plotting functions.

Notes
-----
Adapted from Ryan Dale's GitHub Gist for plotting chromosome features. [1]_

References
----------
.. [1] Ryan Dale, GitHub Gist,
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
Copyright (C) 2017 Andrew Riha

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
import logging
import os

from atomicwrites import atomic_write
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from matplotlib import patches

logger = logging.getLogger(__name__)


def plot_chromosomes(one_chrom_match, two_chrom_match, cytobands, path, title, build):
    """ Plots chromosomes with designated markers.

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
        ax.add_collection(collection)

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

    Yields BrokenBarHCollection of features that can be added to an Axes
    object.

    Parameters
    ----------

    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.

    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection

    height : float
        Height of each BrokenBarHCollection

    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if "width" not in df.columns:
        del_width = True
        df["width"] = df["end"] - df["start"]
    for chrom, group in df.groupby("chrom"):
        yrange = (y_positions["chr" + chrom], height)
        xranges = group[["start", "width"]].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group["colors"], **kwargs
        )
    if del_width:
        del df["width"]


def _patch_chromosomal_features(cytobands, one_chrom_match, two_chrom_match):
    """ Highlight positions for each chromosome segment / feature.

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
        df = df.append(
            {
                "chrom": chromosome,
                "start": 0,
                "end": chromosome_length,
                "gie_stain": "gneg",
            },
            ignore_index=True,
        )

        # add markers for shared DNA on one chromosome
        for marker in one_chrom_match_markers.itertuples():
            df = df.append(
                {
                    "chrom": chromosome,
                    "start": marker.start,
                    "end": marker.end,
                    "gie_stain": "one_chrom",
                },
                ignore_index=True,
            )

        # add markers for shared DNA on both chromosomes
        for marker in two_chrom_match_markers.itertuples():
            df = df.append(
                {
                    "chrom": chromosome,
                    "start": marker.start,
                    "end": marker.end,
                    "gie_stain": "two_chrom",
                },
                ignore_index=True,
            )

        # add centromeres
        for item in cytobands.loc[
            (cytobands["chrom"] == chromosome) & (cytobands["gie_stain"] == "acen")
        ].itertuples():
            df = df.append(
                {
                    "chrom": chromosome,
                    "start": item.start,
                    "end": item.end,
                    "gie_stain": "centromere",
                },
                ignore_index=True,
            )

    return df
