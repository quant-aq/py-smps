"""Control plot style and scaling using seaborn and the matplotlib rcParams interface"""

import matplotlib as mpl

__all__ = ["set"]

def set(tick_scale=1, rc=dict()):
    """
    """
    rc_log_defaults = {
        'xtick.major.size': 10. * tick_scale,
        'xtick.minor.size': 6. * tick_scale,
        'ytick.major.size': 10. * tick_scale,
        'ytick.minor.size': 6. * tick_scale,
        'xtick.color': '0.0',
        'ytick.color': '0.0',
        'axes.linewidth': 1.75,
        'mathtext.default': 'regular'
    }

    mpl.rcParams.update(dict(rc_log_defaults, **rc))
