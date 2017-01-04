# -*- coding: utf-8 -*-

"""
@package
@brief
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2015
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Standard library imports
# Third party imports
# Local imports

def fill_between_graph (X, Y, basename="out", img_type="png", title=None, xlabel=None, ylabel=None,
    baseline=0, xsize=15, ysize=10, dpi=100, fill_color='green'):
    """
    Trace a generic fill between graph with matplotlib pyplot
    @param X List of values for x axis
    @param Y List of values for y axis
    @param title Title of graph (facultative)
    @param xlabel Label for x axis (facultative)
    @param ylabel Label for y axis (facultative)
    @param basename Output basename of the image file (Default "out")
    @param img_type Type of the image file (Default "png")
    @param baseline lower value of the colorated area (Default 0)
    @param xsize Width of the graphics (Default 15)
    @param ysize Heigth of the graphics (Default 10)
    @param dpi Resolution of the graphics (Default 100)
    @param fill_color Color of the filled area (Default 'green')
    """

    # Require the Third party package matplotlib
    from matplotlib import pyplot as plt

    # Create a figure object and adding details
    fig = plt.figure(figsize=(xsize, ysize), dpi=dpi)
    if title:
        plt.title(title)
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)

    # Plot an area representing the coverage depth
    plt.fill_between(X, Y, baseline, facecolor=fill_color, alpha=0.5)

    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)

    # Export figure to file
    try:
        fig.savefig(basename+"."+img_type, format = img_type)
    except ValueError as E:
        print (E)
        print ("Saving file as png")
        fig.savefig(basename+".png", format = "png")
