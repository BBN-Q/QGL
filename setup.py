from setuptools import setup, find_packages

setup(name='QGL',
      version='2.1',
      url='https://github.com/BBN-Q/QGL',
      packages=find_packages(exclude=["tests"]),
      install_requires=[
        "bbndb >= 0.1",
        "numpy >= 1.11.1",
        "scipy >= 0.17.1",
        "networkx >= 1.11",
        "bqplot >= 0.11.5",
        "sqlalchemy >= 1.2.15"
    ])
