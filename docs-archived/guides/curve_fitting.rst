Curve Fitting
=============

There are many reasons one may want to fit their SMPS/OPS/OPC data.
Here, we briefly show how one can use ``py-smps`` to do so. Currently,
there is support for fitting between 1 and 3 modes.

.. code:: ipython3

    import smps
    import seaborn as sns
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick
    import random
    import warnings
    
    warnings.simplefilter(action='ignore', category=FutureWarning)
    
    sns.set("notebook", style='ticks', font_scale=1.25, palette='colorblind')
    smps.set()
    
    %matplotlib inline

Fit a Single-Mode Particle Size Distribution
--------------------------------------------

Here, we will use the provided ``boston`` SMPS dataset as an example.
Let’s begin by loading the data and plotting the particle size
distribution:

.. code:: ipython3

    obj = smps.io.load_sample("boston")
    
    # Plot the histogram
    ax = smps.plots.histplot(
        obj.dndlogdp,
        obj.bins,
        plot_kws=dict(linewidth=.1),
        fig_kws=dict(figsize=(12, 6))
    )
    
    ax.set_title("Wintertime in Cambridge, MA", y=1.02)
    ax.set_ylabel("$dN/dlogD_p \; [cm^{-3}]$")
    
    # remove the spines of the plot
    sns.despine()



.. image:: curve_fitting_files/curve_fitting_3_0.png


Above, we plotted the particle size distribution and can see that during
this time, there was a single dominant mode with a mode particle
diameter of around 50-60 nm. We can then go ahead and fit a simple
1-mode distribution to our data.

.. code:: ipython3

    from smps.fit import LogNormal
    
    # Initiate the class
    model = LogNormal()
    
    # Gather our X and Y values
    X = obj.midpoints
    Y = obj.dndlogdp.mean()
    
    # Fit the data
    results = model.fit(X, Y, modes=1)
    
    # Print out the results
    results.summary()




.. parsed-literal::

                                  LogNormalFitResults                               
    ================================================================================
                              N (cm-3)            GM (nm)               GSD         
    --------------------------------------------------------------------------------
           Mode 0        2.73e+03 (7.3e+00)   56.07 (1.2e-01)      2.10 (5.2e-03)   
    --------------------------------------------------------------------------------



Above, we see the results printed out as a table with the three fit
parameters:

-  the number concentration in particles per cubic centimeter
-  the geometric mean diameter
-  the geometric standard deviation

All three parameters have error estimates as well (standard deviation)
as shown in parens. Now that we’ve successfully fit our data, let’s go
ahead and plot it to make sure it’s correct!

.. code:: ipython3

    ax = smps.plots.histplot(
        obj.dndlogdp,
        obj.bins,
        plot_kws=dict(linewidth=0, alpha=.6, edgecolor=None),
        fig_kws=dict(figsize=(12, 6))
    )
    
    # Plot the fit values
    ax.plot(obj.midpoints, results.fittedvalues, lw=6, label="Fit Data")
    
    ax.set_ylabel("$dN/dlogD_p \; [cm^{-3}]$")
    ax.set_title("Wintertime in Cambridge, MA with Fit Data")
    
    # remove the spines of the plot
    sns.despine()



.. image:: curve_fitting_files/curve_fitting_7_0.png


So, what else is stored alongside the ``fittedvalues`` in the fit
results? Glad you asked! For beginners, you can go ahead and pull the
fit parameters using ``results['params']``. They are stored in format
[``N``, ``GM``, ``GSD``].

.. code:: ipython3

    results.params




.. parsed-literal::

    array([[2.73066908e+03, 5.60740983e-02, 2.10235549e+00]])



You can also go ahead and pull the error associated with those values:

.. code:: ipython3

    results.errors




.. parsed-literal::

    array([[7.30970583e+00, 1.23588578e-04, 5.17954875e-03]])



Upon fitting, an instance of the ``LogNormalFitResults`` class is
returned and has available a couple of useful methods. The first is the
``.summary()`` method we showed above. There is also a ``.predict()``
method so that you can predict values given a fit. It takes two
arguments:

-  ``X`` - an array of values (particle diameters)
-  ``weight`` - one of [``number``, ``surface``, or ``volume``]

.. code:: ipython3

    results.predict(1.)




.. parsed-literal::

    1.8359079994604226



Plot Missing Data
~~~~~~~~~~~~~~~~~

Let’s use the ``predict`` method to fill in the lower portion of the
curve we the SMPS was not scanning. Is this a great idea? Probably not,
but we can still do it anyways!

.. code:: ipython3

    newX = np.logspace(np.log10(.01), np.log10(1), 1000)
    
    # plot the histogram
    ax = smps.plots.histplot(obj.dndlogdp, obj.bins, plot_kws={'linewidth': 0., 'alpha': .5},
                            fig_kws={'figsize': (12, 6)})
    
    # Plot the fit values
    ax.plot(newX, results.predict(newX), lw=6, label="Fit Data")
    
    ax.set_title("Wintertime in Cambridge, MA with Fit Data")
    ax.set_ylabel("$dN/dlogD_p \; [cm^{-3}]$")
    
    # remove the spines of the plot
    sns.despine()



.. image:: curve_fitting_files/curve_fitting_15_0.png


Fit a Multi-Mode Particle Size Distribution
-------------------------------------------

While the existing sample data doesn’t have a strong multi-mode period,
we can mock the data to show the utility of ``py-smps``. **NOTE: If you
are in posession of such a data set and feel like donating its use for
this project, please reach out!**.

First, let’s build a noisy dataset

.. code:: ipython3

    dp = np.logspace(np.log10(1e-4), np.log10(1), 500)
    
    # Sample data pulled from S+P pg371
    N = np.array([9.93e4, 3.64e4])
    GM = np.array([1.3e-3, 20e-3])
    GSD = np.array([10**.245, 10**0.336])
    
    total = 0
    
    for j in range(len(N)):
        total += smps.fit.dndlogdp(dp, N[j], GM[j], GSD[j])
        
    # Let's confuzzle our data
    twisted = total* [random.uniform(0.9, 1.1) for i in range(len(dp))]
    
    with sns.axes_style('ticks'):
        fig, ax = plt.subplots(1, figsize=(12, 6))
    
        ax.plot(dp, twisted, 'o', label="Twisted Data")
    
        ax.set_xlabel("$D_p \; [\mu m]$")
        ax.set_ylabel("$dN/dlogD_p$")
        ax.semilogx()
        
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter("%.4g"))
        ax.legend()
    
    sns.despine()
    plt.show()



.. image:: curve_fitting_files/curve_fitting_17_0.png


Now that we have some mocked data, let’s go ahead and fit it! We’re also
going to need to go ahead and throw some initial guesses in - there need
to be 3xn guesses where n is the number of modes you are fitting. They
should be in format [:math:`N_i`, :math:`GM_i`, :math:`GSD_i`] for i=1
to i=n.

.. code:: ipython3

    model = smps.fit.LogNormal()
    
    X = dp
    Y = twisted
    
    # Let's state some initial guesses
    p0 = [1e5, 1e-3, 2, 3e4, 20e-3, 2]
    
    results = model.fit(X, Y, modes=2, p0=p0)
    
    results.summary()




.. parsed-literal::

                                  LogNormalFitResults                               
    ================================================================================
                              N (cm-3)            GM (nm)               GSD         
    --------------------------------------------------------------------------------
           Mode 0        9.95e+04 (3.5e+02)    1.30 (3.0e-03)      1.76 (4.2e-03)   
           Mode 1        3.64e+04 (4.1e+02)   20.07 (2.0e-01)      2.16 (2.2e-02)   
    --------------------------------------------------------------------------------



Now that we have the results, let’s go ahead and plot them!

.. code:: ipython3

    with sns.axes_style('ticks'):
        fig, ax = plt.subplots(1, figsize=(12, 6))
    
        ax.plot(dp, twisted, 'o', label="Twisted Data")
        ax.plot(dp, results.fittedvalues, lw=6, label="Fitted Values")
    
        ax.set_xlabel("$D_p \; [\mu m]$")
        ax.set_ylabel("$dN/dlogD_p$")
        ax.semilogx()
        
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter("%.4g"))
        ax.legend()
    
    sns.despine()
    plt.show()



.. image:: curve_fitting_files/curve_fitting_21_0.png


Nice! Now, what if we have a dataset like above, but we only want to fit
a portion of it? No worries, just fit 1 mode under specified diameters:

.. code:: ipython3

    model = smps.fit.LogNormal()
    
    X = dp
    Y = twisted
    
    results = model.fit(X, Y, modes=1, xmax=8.5, xmin=0)
    
    print (results.summary())
    
    
    with sns.axes_style('ticks'):
        fig, ax = plt.subplots(1, figsize=(12, 6))
    
        ax.plot(dp, twisted, 'o', label="Twisted Data")
        ax.plot(X[X <= 8.5], results.fittedvalues, lw=6, label="Fitted Values")
    
        ax.set_xlabel("$D_p \; [\mu m]$")
        ax.set_ylabel("$dN/dlogD_p$")
        ax.semilogx()
        
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter("%.4g"))
        ax.legend()
    
    sns.despine()
    plt.show()


.. parsed-literal::

                                  LogNormalFitResults                               
    ================================================================================
                              N (cm-3)            GM (nm)               GSD         
    --------------------------------------------------------------------------------
           Mode 0        1.02e+05 (1.8e+03)    1.32 (1.6e-02)      1.80 (2.1e-02)   
    --------------------------------------------------------------------------------



.. image:: curve_fitting_files/curve_fitting_23_1.png


