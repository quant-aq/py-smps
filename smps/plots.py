"""
"""
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

def heatmap(X, Y, Z, ax=None, log=True, cbar=True, cmap=default_cmap, fig_kws=None, cbar_kws=None, **kwargs):
    """
    """
    if fig_kws is None:
        fig_kws = dict(figsize=(16,8))

    if cbar_kws is None:
        cbar_kws = dict(label='$dN/dlogD_p \; [cm^{-3}]$')

    with sns.axes_style('ticks', rc_log):
        if ax is None:
            plt.figure(**fig_kws)
            ax = plt.gca()

        im = ax.pcolormesh(X, Y, Z, norm=LogNorm(vmin=1., vmax=Z.max()), cmap=cmap)

        ax.set_ylim([Y.min(), Y.max()])

        if log:
            ax.semilogy()

            ax.yaxis.set_major_formatter(ScalarFormatter())

        ax.set_ylabel("$D_p \; (nm)$", fontsize=28)

        if cbar:
            clb = plt.colorbar(im, **cbar_kws)

    return ax
