#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.ticker as mtick
from matplotlib.ticker import ScalarFormatter
import seaborn as sns
from numpy import nan_to_num

default_cmap = 'viridis'

__all__ = ["heatmap", "histplot"]


def heatmap(X, Y, Z, ax=None, logy=True, cbar=True, hide_low=True, 
            cmap=default_cmap, fig_kws={}, cbar_kws={}, plot_kws={}, 
            **kwargs):
    """Plot the heatmap of the particle size distribution. 
    All NaN'd values will be converted to zeroes.
    
    :param X: The boundaries of cells in x-axis dimension.
    :type X: array-like
    :param Y: The boundaries of cells in y-axis dimension.
    :type Y: array-like
    :param Z: The color value in the cells.
    :type Z: 2D array-like
    :param ax: A `matplotlib.pyplot` Axes in which to build the 
        plot. If `None`, a new Figure is instantiated. Defaults to 
        `None`.
    :type ax: `matplotlib.pyplot` Axes
    :param logy: Set whether the y axis is on a log scale, 
        defaults to True.
    :type logy: bool
    :param cbar: Set whether to include a color bar in the plot key, 
        defaults to True.
    :type cbar: bool
    :param hide_low: Set whether `Z` values below `cbar_min` should 
        be set to `cbar_min`, defaults to True.
    :type hide_low: bool
    :param cmap: A color map for `matplotlib.pyplot` `pcolormesh`, 
        defaults to `default_cmap`. See `here 
        <https://matplotlib.org/stable/api/_as_gen
        /matplotlib.colors.Colormap.html
        #matplotlib.colors.Colormap>`__.
    :type cmap: string or Colormap
    :param fig_kws: Options for the `matplotlib.pyplot` figure.
        See `here <https://matplotlib.org/stable/api/_as_gen
        /matplotlib.pyplot.figure.html>`__,
        efaults to `dict(figsize=(16,8)`.
    :type fig_kws: dict
    :param cbar_kws: Options for the `matplotlib.pyplot` `colorbar`.
        See `here <https://matplotlib.org/stable
        /api/figure_api.html#matplotlib.figure.Figure.colorbar>`__.
    :type cbar_kws: dict
    :param plot_kws: Options for the `matplotlib.pyplot` Axes.
        See `here <https://matplotlib.org/stable/api/_as_gen
        /matplotlib.pyplot.gca.html>'__,
        Defaults to `dict(alpha=1, edgecolor=None, linewidth=0)`.
    :type plot_kws: dict
    :param cbar_min: The minimum value for the color bar, 
        defaults to the smallest `Z` value
        (unless it is not greater than zero, 
        and then it will default to one).
    :type cbar_min: float
    :param cbar_max: The maximum value for the color bar, 
        defaults to the largest `Z` value.
    :type cbar_max: float
    :return: A plot of the histogram.
    :rtype: `matplotlib.pyplot` Axes
    """
    # Copy to avoid modifying original data
    Z_plot = Z.copy()

    # get rid of NaNs
    Z_plot = nan_to_num(Z_plot)

    # Set the colorbar min and max based on the min and 
    # max of the values
    cbar_min = kwargs.pop(
        'cbar_min', 
        Z_plot.min() if Z_plot.min() > 0.0 else 1.
    )
    cbar_max = kwargs.pop('cbar_max', Z_plot.max())

    if hide_low:
        # Increase values below cbar_min to cbar_min
        below_min = Z_plot < cbar_min
        Z_plot[below_min] = cbar_min

    # Set the plot_kws
    plot_kws = dict(
        dict(norm=LogNorm(vmin=cbar_min, vmax=cbar_max), cmap=cmap), 
        **plot_kws
    )

    # Set the figure keywords
    fig_kws = dict(dict(figsize=(10,5)), **fig_kws)

    if ax is None:
        plt.figure(**fig_kws)
        ax = plt.gca()

    # Plot the data as a pcolormesh
    im = ax.pcolormesh(X, Y, Z_plot, shading='auto', **plot_kws)

    # Set the ylim to match the data
    ax.set_ylim([Y.min(), Y.max()])

    # Set the axis to be log in the y-axis
    if logy:
        ax.semilogy()

        ax.yaxis.set_major_formatter(ScalarFormatter())

    ax.set_ylabel(r"$D_p\;[\mu m]$")

    if cbar:
        # Set the figure keywords
        cbar_kws = dict(
            dict(label=r'$dN/dlogD_p\;[cm^{-3}]$'), 
            **cbar_kws
        )

        clb = plt.colorbar(im, **cbar_kws)

    return ax


def histplot(histogram, bins, ax=None, plot_kws=None, fig_kws=None, 
             **kwargs):
    """
    Plot the histogram in the form of a bar chart.
    
    :param histogram: A histogram indexed by bins.
        If a DataFrame is given the data will be averaged over time.
    :type histogram: Pandas DataFrame, Pandas Series
    :param bins: For each bin a row with the lower boundary, the
        midpoint, and the upper bounary for that bin.
    :type bins: numpy array
    :param ax: `matplotlib.pyplot` Axes in which to build the plot.
        If `None`, a new Figure is instantiated.
    :type ax: `matplotlib.pyplot` Axes
    :param plot_kws: Options for the `matplotlib.pyplot` Axes.
        See `here <https://matplotlib.org/stable/api/_as_gen
        /matplotlib.pyplot.gca.html>'__.
        Defaults to `dict(alpha=1, edgecolor=None, linewidth=0)`.
    :type plot_kws: dict
    :param fig_kws: Options for the `matplotlib.pyplot` figure.
        See `here <https://matplotlib.org/stable/api/_as_gen
        /matplotlib.pyplot.figure.html>`__,
        defaults to `dict(figsize=(16,8)`.
    :type fig_kws: dict
    :return: A plot of the histogram.
    :rtype: `matplotlib.pyplot` Axes
    """
    if isinstance(histogram, pd.DataFrame):
        histogram = histogram.mean().values

    if fig_kws is None:
        fig_kws = dict(figsize=(16,8))

    if plot_kws is None:
        plot_kws = dict(alpha=1, edgecolor=None, linewidth=0)

    if ax is None:
        plt.figure(**fig_kws)
        ax = plt.gca()

    ax.bar(x=bins[:, 0], height=histogram, width=bins[:, -1] - bins[:, 0],
            align='edge', **plot_kws)

    ax.semilogx()

    ax.set_xlabel(r"$D_p\;[\mu m]$")

    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter("%.3g"))

    # create a large number formatter for the y axis
    fmt = mtick.ScalarFormatter()
    fmt.set_powerlimits((-3, 2))

    ax.yaxis.set_major_formatter(fmt)

    return ax
