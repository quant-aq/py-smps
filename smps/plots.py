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
            cmap=default_cmap, fig_kws={}, cbar_kws={}, plot_kws={}, **kwargs):
    """Plot the heatmap of the particle size distribution. All NaN'd values
    will be converted to zeroes.
    """
    # Copy to avoid modifying original data
    Z_plot = Z.copy()

    # get rid of NaNs
    Z_plot = nan_to_num(Z_plot)

    # Set the colorbar min and max based on the min and max of the values
    cbar_min = kwargs.pop('cbar_min', Z_plot.min() if Z_plot.min() > 0.0 else 1.)
    cbar_max = kwargs.pop('cbar_max', Z_plot.max())

    if hide_low:
        # Increase values below cbar_min to cbar_min
        below_min = Z_plot < cbar_min
        Z_plot[below_min] = cbar_min

    # Set the plot_kws
    plot_kws = dict(dict(norm=LogNorm(vmin=cbar_min, vmax=cbar_max), cmap=cmap),
                        **plot_kws)

    # Set the figure keywords
    fig_kws = dict(dict(figsize=(10,5)), **fig_kws)

    if ax is None:
        plt.figure(**fig_kws)
        ax = plt.gca()

    # Plot the data as a pcolormesh
    im = ax.pcolormesh(X, Y, Z_plot, **plot_kws)

    # Set the ylim to match the data
    ax.set_ylim([Y.min(), Y.max()])

    # Set the axis to be log in the y-axis
    if logy:
        ax.semilogy()

        ax.yaxis.set_major_formatter(ScalarFormatter())

    ax.set_ylabel(r"$D_p\;[\mu m]$")

    if cbar:
        # Set the figure keywords
        cbar_kws = dict(dict(label=r'$dN/dlogD_p\;[cm^{-3}]$'), **cbar_kws)

        clb = plt.colorbar(im, **cbar_kws)

    return ax


def histplot(histogram, bins, ax=None, plot_kws=None, fig_kws=None, **kwargs):
    """Plot the histogram in the form of a bar chart."""
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
