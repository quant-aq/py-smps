try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

__version__ = '1.0.1'

setup(
    name='py-smps',
    version=__version__,
    packages=['smps'],
    description='A simple python library to import and visualize data from particle sizing instruments.',
    author='David H Hagan',
    author_email='dhagan@mit.edu',
    license='MIT',
    url='https://github.com/dhhagan/py-smps',
    keywords=['atmospheric chemistry'],
    test_suite='tests',
    install_requires=[
        'pandas',
        'matplotlib',
        'numpy',
        'seaborn',
        'scipy',
        'requests',
        'joblib',
        'statsmodels'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
		'Intended Audience :: Education',
		'Programming Language :: Python :: 2.7',
		'Programming Language :: Python :: 3.3',
		'Programming Language :: Python :: 3.4',
		'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
		'Topic :: Scientific/Engineering :: Atmospheric Science',
		'Topic :: Software Development',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ]
)
