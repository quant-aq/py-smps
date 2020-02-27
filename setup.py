from setuptools import setup
import versioneer

setup(
    name='py-smps',
    version=versioneer.get_version(),
    packages=['smps'],
    description='A simple python library to import and visualize data from particle sizing instruments.',
    author='David H Hagan',
    author_email='david.hagan@quant-aq.com',
    license='MIT',
    url='https://github.com/dhhagan/py-smps',
    keywords=['atmospheric chemistry'],
    test_suite='tests',
    python_requires=">=3.5",
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
		'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
		'Topic :: Scientific/Engineering :: Atmospheric Science',
		'Topic :: Software Development',
        'Topic :: Software Development :: Libraries :: Python Modules'
    ],
    cmdclass=versioneer.get_cmdclass()
)
