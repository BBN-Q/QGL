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
        "bokeh >= 0.12.13",
        # "pony >= 0.7.4", # This needs to be 0.7.4-dev
    ])
