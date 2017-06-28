""" Chromosome plotting functions.

Notes
-----
Adapted from Ryan Dale's GitHub Gist for plotting chromosome features
(see [1]_).

References
----------
..[1] Ryan Dale, GitHub Gist,
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

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection


def plot_chromosomes(markers, cytobands, path, title, build):
    """ Plots chromosomes with designated markers.

    Parameters
    ----------
    markers : list of dicts
        segments to highlight on the chromosomes
    cytobands : pandas.DataFrame
        cytobands table loaded with Resources
    path : str
        path to destination `.png` file
    title : str
        title for plot
    build : {36, 37}
        human genome assembly

    """

    # Height of each chromosome
    chrom_height = 1.25

    # Spacing between consecutive chromosomes
    chrom_spacing = 1

    # Decide which chromosomes to use
    chromosome_list = ['chr%s' % i for i in range(1, 23)]
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
        chrom_centers[chrom] = ybase + chrom_height / 2.
        ybase += chrom_height + chrom_spacing

    # Colors for different chromosome stains
    color_lookup = {
        'gneg': (66/255, 69/255, 121/255),  # background
        'match': (0/255, 176/255, 240/255),
        'centromere': (0.95, 0.95, 0.95, 0.5)
    }

    df = _patch_chromosomal_features(cytobands, markers)

    # Add a new column for colors
    df['colors'] = df['gie_stain'].apply(lambda x: color_lookup[x])

    # Width, height (in inches)
    figsize = (6.5, 9)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    # Now all we have to do is call our function for the chromosome data...
    print("adding chromosomes...")
    for collection in _chromosome_collections(df, chrom_ybase, chrom_height):
        ax.add_collection(collection)

    # Axes tweaking
    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list)
    ax.margins(0.01)
    ax.axis('tight')

    fig.text(0.8, 0.2, "Shared DNA", size=10, va="center", ha="center",
             bbox=dict(boxstyle="square", fc=(0/255, 176/255, 240/255), ec='none'))

    ax.set_title(title, fontsize=14, fontweight='bold')
    plt.xlabel("Build " + str(build) + " Chromosome Position", fontsize=10)
    plt.savefig(path)


def _chromosome_collections(df, y_positions, height,  **kwargs):
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
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        yrange = (y_positions['chr' + chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']


def _patch_chromosomal_features(cytobands, markers):
    """ Highlight positions for each chromosome segment / feature.

    Parameters
    ----------
    cytobands : pandas.DataFrame
        cytoband table from UCSC
    markers : list of dicts
        segments to highlight on the chromosomes

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
            cytobands[cytobands["chrom"] == chromosome]["end"].values)

        # get all markers for this chromosome
        chromosome_markers = [marker for marker in markers if marker["chrom"] == chromosome]

        df = df.append({"chrom": chromosome, "start": 0, "end": chromosome_length,
                        "gie_stain": "gneg"}, ignore_index=True)

        for marker in chromosome_markers:
            df = df.append({"chrom": chromosome, "start": marker["start"], "end": marker[
                "end"],
                            "gie_stain": marker["gie_stain"]}, ignore_index=True)

        for item in cytobands.loc[(cytobands["chrom"] == chromosome) &
                                  (cytobands["gie_stain"] == 'acen')].itertuples():
            df = df.append({"chrom": chromosome, "start": item.start, "end": item.end,
                            "gie_stain": "centromere"}, ignore_index=True)

    return df
