"""
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import ScalarFormatter
import seaborn as sns

default_cmap = sns.cubehelix_palette(8, as_cmap=True)

rc_log = {
    'xtick.major.size': 8.0,
    'xtick.minor.size': 5.0,
    'ytick.major.size': 8.0,
    'ytick.minor.size': 5.0
}

def heatmap(X, Y, Z, ax=None, kind='log', cbar=True, cmap=default_cmap,
            fig_kws=None, cbar_kws=None, **kwargs):
    """
    """
    cbar_min = kwargs.pop('cbar_min', Z.min() if Z.min() > 0.0 else 1.)
    cbar_max = kwargs.pop('cbar_max', Z.max())

    if fig_kws is None:
        fig_kws = dict(figsize=(16,8))

    if cbar_kws is None:
        cbar_kws = dict(label='$dN/dlogD_p \; [cm^{-3}]$')

    with sns.axes_style('ticks', rc_log):
        if ax is None:
            plt.figure(**fig_kws)
            ax = plt.gca()

        im = ax.pcolormesh(X, Y, Z, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), cmap=cmap)

        ax.set_ylim([Y.min(), Y.max()])

        if kind == 'log':
            ax.semilogy()

            ax.yaxis.set_major_formatter(ScalarFormatter())

        ax.set_ylabel("$D_p \; [nm]$", fontsize=28)

        if cbar:
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

    with sns.axes_style('ticks', rc_log):
        if ax is None:
            plt.figure(**fig_kws)
            ax = plt.gca()

        ax.bar(bins[:, 0], histogram, bins[:, -1] - bins[:, 0], **plot_kws)

        ax.semilogx()

        ax.set_ylabel("$dN/dlogD_p \; [cm^{-3}]$", fontsize=28)
        ax.set_xlabel("$D_p \; [nm]$", fontsize=28)

        ax.xaxis.set_major_formatter(ScalarFormatter())

    return ax
