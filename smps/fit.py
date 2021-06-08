"""Functions to support fitting multi-modal distributions.
"""
import numpy as np
from scipy.optimize import curve_fit
from statsmodels.iolib.table import SimpleTable
from .utils import Table

def dndlogdp(dp, n, gm, gsd):
    """
    Discrete PDF for a lognormal particle size distribution.
    
    :param dp: Particle diameters at which to evaluate the 
        PDF in building the discrete model.
    :type dp: numpy array
    :param n: Area under the curve of the continuous PDF 
        (the discrete model will be slightly different).
    :type n: float
    :param gm: The geometric mean.
    :type gm: float
    :param gsd: The geometric standard deviation.
    :type gsd: float
    :return: A list of concentrations respective to the 
        list of particle sizes, `dp`.
    :rtype: numpy array
    """
    return  ((n/(np.sqrt(2*np.pi) * np.log10(gsd))) * 
             np.exp(-((np.log10(dp) - np.log10(gm))**2) / 
                    (2*np.log10(gsd)**2)))

def dsdlogdp(dp, n, gm, gsd):
    """
    Discrete PDF for a surface-weighted lognormal 
    particle size distribution.
    
    :param dp: Particle diameters at which to evaluate 
        the PDF in building the discrete model.
    :type dp: numpy array
    :param n: The area under the curve of the continuous PDF 
        (the discrete model will be slightly different).
    :type n: float
    :param gm: The geometric mean.
    :type gm: float
    :param gsd: The geometric standard deviation.
    :type gsd: float
    :return: A list of concentrations respective to the 
        list of particle sizes, `dp`.
    :rtype: numpy array
    """
    return np.pi * dp**2 * dndlogdp(dp, n, gm, gsd)

def dvdlogdp(dp, n, gm, gsd):
    """
    Discrete PDF for a volume-weighted lognormal particle 
    size distribution.
    
    :param dp: Particle diameters at which to evaluate 
        the PDF in building the discrete model.
    :type dp: numpy array
    :param n: The area under the curve of the continuous PDF 
        (the discrete model will be slightly different).
    :type n: float
    :param gm: The geometric mean.
    :type gm: float
    :param gsd: The geometric standard deviation.
    :type gsd: float
    :return: A list of concentrations respective to 
        the list of particle sizes, `dp`.
    :rtype: numpy array
    """
    return (np.pi/6.) * (dp**3) * dndlogdp(dp, n, gm, gsd)

def number_weighted_single_mode(dp, n, gm, gsd):
    """
    Discrete PDF for a lognormal distribution of particle 
    size for a 1-mode distribution.
    
    :param dp: Particle diameters at which to evaluate 
        the PDF in building the discrete model.
    :type dp: numpy array
    :param n: The  area under the curve of the continuous PDF 
        (the discrete model will be slightly different).
    :type n: float
    :param gm: The geometric mean.
    :type gm: float
    :param gsd: The geometric standard deviation.
    :type gsd: float
    :return: A list of concentrations respective to the 
        list of particle sizes, `dp`.
    :rtype: numpy array
    """
    return dndlogdp(dp, n, gm, gsd)

def number_weighted_two_modes(dp, n1, gm1, gsd1, n2, gm2, gsd2):
    """
    Discrete PDF for a lognormal distribution of particle 
    size for a 2-mode distribution.
    
    :param dp: Particle diameters at which to evaluate 
        the PDF in building the discrete model.
    :type dp: numpy array
    :param n1: The area under the curve of the continuous 
        PDF for the first component (the discrete 
        model will be slightly different).
    :type n1: float
    :param gm1: The geometric mean for the first component.
    :type gm1: float
    :param gsd1: The geometric standard deviation for 
        the first component.
    :type gsd1: float
    :param n2: The area under the curve of the continuous 
        PDF for the second component (the discrete 
        model will be slightly different).
    :type n2: float
    :param gm2: The geometric mean for the second component.
    :type gm2: float
    :param gsd2: The geometric standard deviation for 
        the second component.
    :type gsd2: float
    :return: A list of concentrations respective to 
        the list of particle sizes, `dp`.
    :rtype: numpy array
    """
    N = [n1, n2]
    GM = [gm1, gm2]
    GSD = [gsd1, gsd2]

    s = 0
    for i in range(len(N)):
        s += dndlogdp(dp, N[i], GM[i], GSD[i])

    return s

def number_weighted_three_modes(dp, n1, gm1, gsd1, n2, gm2, 
                                gsd2, n3, gm3, gsd3):
    """
    Discrete PDF for a lognormal distribution of particle 
        size for a 3-mode distribution.
    
    :param dp: The particle diameters at which to evaluate 
        the PDF in building the discrete model.
    :type dp: numpy array
    :param n1: The area under the curve of the continuous 
        PDF for the first component (the discrete 
        model will be slightly different).
    :type n1: float
    :param gm1: The geometric mean for the first component.
    :type gm1: float
    :param gsd1: The geometric standard deviation for .
        the first component
    :type gsd1: float
    :param n2: The area under the curve of the continuous 
        PDF for the second component (the discrete 
        model will be slightly different).
    :type n2: float
    :param gm2: The geometric mean for the second component.
    :type gm2: float
    :param gsd2: The geometric standard deviation for 
        the second component.
    :type gsd2: float
    :param n3: The area under the curve of the continuous 
        PDF for the third component (the discrete 
        model will be slightly different).
    :type n3: float
    :param gm3: The geometric mean for the third component.
    :type gm3: float
    :param gsd3: The geometric standard deviation for 
        the third component.
    :type gsd3: float
    :return: A list of concentrations respective to 
        the list of particle sizes, `dp`.
    :rtype: numpy array
    """
    N = [n1, n2, n3]
    GM = [gm1, gm2, gm3]
    GSD = [gsd1, gsd2, gsd3]

    s = 0
    for i in range(len(N)):
        s += dndlogdp(dp, N[i], GM[i], GSD[i])

    return s

def surface_weighted_single_mode(dp, n, gm, gsd):
    """
    Discrete PDF for a lognormal distribution of particle 
    size for a 1-mode distribution weighted by surface.
    
    :param dp: The particle diameters at which to evaluate 
        the PDF in building the discrete model.
    :type dp: numpy array
    :param n: The area under the curve of the continuous 
        PDF (the discrete model will be slightly different).
    :type n: float
    :param gm: The geometric mean.
    :type gm: float
    :param gsd: The geometric standard deviation.
    :type gsd: float
    :return: A list of concentrations respective to 
        the list of particle sizes, `dp`.
    :rtype: numpy array
    """
    return dsdlogdp(dp, n, gm, gsd)

def surface_weighted_two_modes(dp, n1, gm1, gsd1, n2, gm2, gsd2):
    """
    Discrete PDF for a lognormal distribution of particle 
    size for a 2-mode distribution weighted by surface.
    
    :param dp: The particle diameters at which to evaluate 
        the PDF in building the discrete model.
    :type dp: numpy array
    :param n1: The area under the curve of the continuous PDF 
        for the first component (the discrete model will 
        be slightly different).
    :type n1: float
    :param gm1: The geometric mean for the first component.
    :type gm1: float
    :param gsd1: The geometric standard deviation for 
        the first component.
    :type gsd1: float
    :param n2: The area under the curve of the continuous 
        PDF for the second component (the discrete model 
        will be slightly different).
    :type n2: float
    :param gm2: The geometric mean for the second component.
    :type gm2: float
    :param gsd2: The geometric standard deviation for the 
        second component.
    :type gsd2: float
    :return: A list of concentrations respective to the 
        list of particle sizes, `dp`.
    :rtype: numpy array
    """
    N = [n1, n2]
    GM = [gm1, gm2]
    GSD = [gsd1, gsd2]

    s = 0
    for i in range(len(N)):
        s += dsdlogdp(dp, N[i], GM[i], GSD[i])

    return s

def surface_weighted_three_modes(dp, n1, gm1, gsd1, n2, gm2, 
                                 gsd2, n3, gm3, gsd3):
    """
    Discrete PDF for a lognormal distribution of particle size 
    for a 3-mode distribution weighted by surface.
    
    :param dp: The particle diameters at which to evaluate 
        the PDF in building the discrete model.
    :type dp: numpy array
    :param n1: The area under the curve of the continuous PDF 
        for the first component (the discrete model will 
        be slightly different).
    :type n1: float
    :param gm1: The geometric mean for the first component.
    :type gm1: float
    :param gsd1: The geometric standard deviation for 
        the first component.
    :type gsd1: float
    :param n2: The area under the curve of the continuous PDF 
        for the second component (the discrete model 
        will be slightly different).
    :type n2: float
    :param gm2: The geometric mean for the second component.
    :type gm2: float
    :param gsd2: The geometric standard deviation for 
        the second component.
    :type gsd2: float
    :param n3: The area under the curve of the continuous PDF 
        for the third component (the discrete model will be 
        slightly different).
    :type n3: float
    :param gm3: The geometric mean for the third component.
    :type gm3: float
    :param gsd3: The geometric standard deviation for the 
        third component.
    :type gsd3: float
    :return: A list of concentrations respective to the 
        list of particle sizes, `dp`.
    :rtype: numpy array
    """
    N = [n1, n2, n3]
    GM = [gm1, gm2, gm3]
    GSD = [gsd1, gsd2, gsd3]

    s = 0
    for i in range(len(N)):
        s += dsdlogdp(dp, N[i], GM[i], GSD[i])

    return s

def volume_weighted_single_mode(dp, n, gm, gsd):
    """
    Discrete PDF for a lognormal distribution of particle 
    size for a 1-mode distribution weighted by volume.
    
    :param dp: The particle diameters at which to evaluate the 
        PDF in building the discrete model.
    :type dp: numpy array
    :param n: The area under the curve of the continuous PDF 
        (the discrete model will be slightly different).
    :type n: float
    :param gm: The geometric mean.
    :type gm: float
    :param gsd: The geometric standard deviation.
    :type gsd: float
    :return: A list of concentrations respective to the 
        list of particle sizes, `dp`.
    :rtype: numpy array
    """
    return dvdlogdp(dp, n, gm, gsd)

def volume_weighted_two_modes(dp, n1, gm1, gsd1, n2, gm2, gsd2):
    """
    Discrete PDF for a lognormal distribution of particle 
    size for a 2-mode distribution weighted by volume.
    
    :param dp: The particle diameters at which to evaluate the 
        PDF in building the discrete model.
    :type dp: numpy array
    :param n1: The area under the curve of the continuous PDF 
        for the first component (the discrete model will 
        be slightly different).
    :type n1: float
    :param gm1: The geometric mean for the first component.
    :type gm1: float
    :param gsd1: The geometric standard deviation for the 
        first component.
    :type gsd1: float
    :param n2: The area under the curve of the continuous PDF 
        for the second component (the discrete model will be 
        slightly different).
    :type n2: float
    :param gm2: The geometric mean for the second component.
    :type gm2: float
    :param gsd2: The geometric standard deviation for the 
        second component.
    :type gsd2: float
    :return: A list of concentrations respective to the 
        list of particle sizes, `dp`.
    :rtype: numpy array
    """
    N = [n1, n2]
    GM = [gm1, gm2]
    GSD = [gsd1, gsd2]

    s = 0
    for i in range(len(N)):
        s += dvdlogdp(dp, N[i], GM[i], GSD[i])

    return s

def volume_weighted_three_modes(dp, n1, gm1, gsd1, n2, gm2, gsd2,
                                n3, gm3, gsd3):
    """
    Discrete PDF for a log-normal distribution of particle 
    size for a 3-mode distribution weighted by volume.
    
    :param dp: The particle diameters at which to evaluate the 
        PDF in building the discrete model.
    :type dp: numpy array
    :param n1: The area under the curve of the continuous PDF 
        for the first component (the discrete model will be 
        slightly different).
    :type n1: float
    :param gm1: The geometric mean for the first component.
    :type gm1: float
    :param gsd1: The geometric standard deviation for the 
        first component.
    :type gsd1: flaot
    :param n2: The area under the curve of the continuous PDF 
        for the second component (the discrete model will 
        be slightly different).
    :type n2: float
    :param gm2: The geometric mean for the second component.
    :type gm2: float
    :param gsd2: The geometric standard deviation for the 
        second component.
    :type gsd2: float
    :param n3: The area under the curve of the continuous PDF
        for the third component (the discrete model will be
        slightly different).
    :type n3: float
    :param gm3: The geometric mean for the third component.
    :type gm3: float
    :param gsd3: The geometric standard deviation for the
        third component.
    :type gsd3: float
    :return: A list of concentrations respective to the list
        of particle sizes, `dp`.
    :rtype: numpy array
    """
    N = [n1, n2, n3]
    GM = [gm1, gm2, gm3]
    GSD = [gsd1, gsd2, gsd3]

    s = 0
    for i in range(len(N)):
        s += dvdlogdp(dp, N[i], GM[i], GSD[i])

    return s

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

class LogNormal(object):
    """
    `LogNormal` uses non-linear least squares to fit a
    log-normal model, with up to three modes.
    """
    def __init__(self):
        pass

    def fit(self, X, Y, modes=1, xmin=None,
            xmax=None, weight='number', fit_kwargs=None, **kwargs):
        """
        `fit` uses non-linear least squares to fit a
        log-normal model, with up to three modes.
        
        :param X: The independant variable.
        :type X: array-like
        :param Y: The dependent variable.
        :type Y: array-like
        :param modes: Number of modes for the model to 
            have, defaults to 1.
        :type modes: {1,2,3}
        :param xmin: Lower bound by which to crop data.
        :type xmin: float, optional
        :param xmax: Upperbound by which to crop data.
        :type xmax: float, optional
        :param weight: Dimension of interest.
        :type weight: {"number", "surface", "volume"}
        :param fit_kwargs: Options to be passed 
            to ``scipy.optimize.curve_fit``. Read 
            more `here <https://docs.scipy.org/doc/
            scipy/reference/generated/scipy.optimize.
            curve_fit.html>`__.
        :type fit_kwargs: dict, optional
        :param p0: Initial guesses passed to ``scipy.
            optimize.curve_fit`` for the parameters. 
            Must have length 3*`modes`.
        :type p0: list, optional
        :param bounds: Bounds to be passed to ``scipy.
            optimize.curve_fit``. Defaults to `(0, [1e9, 
            1e5, 5]*self.modes)`.
        :type bounds: 2-tuple
        :return: A fitted log-normal model.
        :rtype: `LogNormalFitResults` object
        """
        
        self.modes = modes

        if fit_kwargs is None:
            fit_kwargs = dict()

        p0 = kwargs.pop('p0', [1e5, 1.5, 2, 1e5, 5, 2.5, 1e3, 50, 2.5])

        bounds = kwargs.pop('bounds', (0, [1e9, 1e5, 5]*self.modes))

        self.model = models[weight][self.modes-1]

        # Set the initial guesses
        p0 = p0[0:self.modes*3]

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

        fittedvalues = self.model(X, *self.fit_params)

        return LogNormalFitResults(params=self.fit_params,
                                   error_matrix=perr,
                                   fittedvalues=fittedvalues,
                                   modes=self.modes)


class LogNormalFitResults(object):
    """
    `LogNormalFitResults` holds the results produced by `LogNormal`.
    
    :param params: The parameters of the fitted model.
    :type params: array
    :param error_matrix: The one standard deviation errors associaed 
        with the paramters, given in the same order as `params`. 
    :type error_matrix: array
    :param fittedvalues: The results from evaluating the model at 
        each of the independent variables in the training data.
    :type fittedvalues: array
    :param modes: The number of modes the model has.
    :type modes: {1,2,3}
    """
    def __init__(self, params, error_matrix, fittedvalues,
                 modes, **kwargs):
        self.modes = modes

        self.params = params.reshape(self.modes, 3)
        self.errors = error_matrix.reshape(self.modes, 3)
        self.fittedvalues = fittedvalues

    def summary(self):
        """
        A `Table` object with the a summary of the model.
        """
        # Convert GM from microns to nm
        params = self.params.copy()
        errors = self.errors.copy()

        params[:, 1] *= 1000.
        errors[:, 1] *= 1000.

        _tbl = Table()

        _tbl.add_title(self.__class__.__name__)
        _tbl.add_border("=")
        _tbl.add_header()
        _tbl.add_border("-")

        # Add a row for each mode
        for i in range(self.modes):
            _tbl.add_row(label="Mode {}".format(i),
                         fields=params[i], errors=errors[i])

        _tbl.add_border("-")

        return _tbl

    def predict(self, X, weight='number'):
        """
        Uses the model to predict the outcome for a given input.
        
        :param X: Particle diameters at which to evaluate 
            the model PDF.
        :type X: numpy array
        :param weight: Dimension of concentration to model.
        :type weight: {"number", "surface", "volume"}
        :return: A list of concentrations respective to the 
            list of particle sizes, `X`.
        :rtype: numpy array
        """
        model = models[weight][self.modes-1]

        return model(X, *self.params.flatten())
