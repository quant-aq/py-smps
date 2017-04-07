"""Functions to support fitting multi-modal distributions.
"""
import numpy as np
from scipy.optimize import curve_fit

def dndlogdp(dp, n, gm, gsd):
    """PDF for a lognormal particle size distribution"""
    return  (n / (np.sqrt(2*np.pi)*np.log10(gsd)))* \
            np.exp(-((np.log10(dp) - np.log10(gm))**2) / (2*np.log10(gsd)**2))

def dsdlogdp(dp, n, gm, gsd):
    """PDF for a surface-weighted lognormal particle size distribution"""
    return np.pi * dp**2 * dndlogdp(dp, n, gm, gsd)

def dvdlogdp(dp, n, gm, gsd):
    """PDF for a volume-weighted lognormal particle size distribution"""
    return (np.pi/6.) * (dp**3) * dndlogdp(dp, n, gm, gsd)

def number_weighted_single_mode(dp, n, gm, gsd):
    """"""
    return dndlogdp(dp, n, gm, gsd)

def number_weighted_two_modes(dp, n1, gm1, gsd1, n2, gm2, gsd2):
    """Lognormal Distribution for a 2-mode distribution"""
    N = [n1, n2]
    GM = [gm1, gm2]
    GSD = [gsd1, gsd2]

    s = 0
    for i in range(len(N)):
        s += dndlogdp(dp, N[i], GM[i], GSD[i])

    return s

def number_weighted_three_modes(dp, n1, gm1, gsd1, n2, gm2, gsd2, n3, gm3, gsd3):
    """Lognormal Distribution for a 2-mode distribution"""
    N = [n1, n2, n3]
    GM = [gm1, gm2, gm3]
    GSD = [gsd1, gsd2, gsd3]

    s = 0
    for i in range(len(N)):
        s += dndlogdp(dp, N[i], GM[i], GSD[i])

    return s

def surface_weighted_single_mode(dp, n, gm, gsd):
    """"""
    return dsdlogdp(dp, n, gm, gsd)

def surface_weighted_two_modes(dp, n1, gm1, gsd1, n2, gm2, gsd2):
    """Lognormal Distribution for a 2-mode distribution"""
    N = [n1, n2]
    GM = [gm1, gm2]
    GSD = [gsd1, gsd2]

    s = 0
    for i in range(len(N)):
        s += dsdlogdp(dp, N[i], GM[i], GSD[i])

    return s

def surface_weighted_three_modes(dp, n1, gm1, gsd1, n2, gm2, gsd2, n3, gm3, gsd3):
    """Lognormal Distribution for a 2-mode distribution"""
    N = [n1, n2, n3]
    GM = [gm1, gm2, gm3]
    GSD = [gsd1, gsd2, gsd3]

    s = 0
    for i in range(len(N)):
        s += dsdlogdp(dp, N[i], GM[i], GSD[i])

    return s

def volume_weighted_single_mode(dp, n, gm, gsd):
    """"""
    return dvdlogdp(dp, n, gm, gsd)

def volume_weighted_two_modes(dp, n1, gm1, gsd1, n2, gm2, gsd2):
    """Lognormal Distribution for a 2-mode distribution"""
    N = [n1, n2]
    GM = [gm1, gm2]
    GSD = [gsd1, gsd2]

    s = 0
    for i in range(len(N)):
        s += dvdlogdp(dp, N[i], GM[i], GSD[i])

    return s

def volume_weighted_three_modes(dp, n1, gm1, gsd1, n2, gm2, gsd2, n3, gm3, gsd3):
    """Lognormal Distribution for a 2-mode distribution"""
    N = [n1, n2, n3]
    GM = [gm1, gm2, gm3]
    GSD = [gsd1, gsd2, gsd3]

    s = 0
    for i in range(len(N)):
        s += dvdlogdp(dp, N[i], GM[i], GSD[i])

    return s

class LogNormal(object):
    def __init__(self):
        pass

    def fit(self, X, Y, modes=1, xmin=None, xmax=None, weight='number', fit_kwargs=None, **kwargs):
        """
        """
        if fit_kwargs is None:
            fit_kwargs = dict()

        models = {
            'number': [
                number_weighted_single_mode,
                number_weighted_two_modes,
                number_weighted_three_modes
                ],
            'surface': [
                surface_weighted_single_mode,
                surface_weighted_two_modes,
                surface_weighted_three_modes
                ],
            'volume': [
                volume_weighted_single_mode,
                volume_weighted_two_modes,
                volume_weighted_three_modes
                ]
            }

        p0 = kwargs.pop('p0', [1e5, 1.5, 2, 1e5, 5, 2.5, 1e3, 50, 2.5])

        bounds = kwargs.pop('bounds', (0, [1e9, 1e5, 5]*modes))

        self.model = models[weight][modes-1]

        # Set the initial guesses
        p0 = p0[0:modes*3]

        # Subset the data if xmin or xmax is set
        if xmin is not None:
            X = X[np.where(X >= xmin)]
            Y = Y[np.where(X >= xmin)]

        if xmax is not None:
            X = X[np.where(X <= xmax)]
            Y = Y[np.where(X <= xmax)]

        self.fit_params, cov = curve_fit(self.model, X, Y, p0,
                                bounds=bounds, **fit_kwargs)

        perr = np.sqrt(np.diag(cov))

        summary = "Mode\tN (#/cc)\tGM (nm)\t\tGSD\n"
        for m in range(modes):
            summary += "{}\t{:.2e}".format(m, self.fit_params[m*3])
            summary += "\t{:.2f}".format(self.fit_params[m*3 + 1]*1000.)
            summary += "\t\t{:.2f}".format(self.fit_params[m*3 + 2])
            summary += "\n"

        fittedvalues = self.model(X, *self.fit_params)

        results = dict(
            params=self.fit_params,
            error=perr,
            summary=summary,
            fittedvalues=fittedvalues)

        return results

    def predict(self, X):
        return self.model(X, *self.fit_params)
