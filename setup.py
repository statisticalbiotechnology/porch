from distutils.core import setup

setup(
    name='porch',
    version='0.0.1',
    author='Statistical Biotechnology',
    author_email='lukas.kall@scilifelab.se',
    packages=['porch', 'qvalue', 'example'],
#    scripts=['bin/stowe-towels.py','bin/wash-towels.py'],
#    url='http://pypi.python.org/pypi/TowelStuff/',
    license='LICENSE',
    description='Personal Pathway Analysis',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy",
        "pandas",
        "statsmodels",
        "scikit-learn",
        "requests",
        "biothings_client",
        "bioservices",
        "seaborn",
        "scipy",
        "wpca",
        "xlrd",
        "networkx",
        "lifelines>=0.24",
    ],
)
