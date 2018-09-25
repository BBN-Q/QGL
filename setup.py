from setuptools import setup, find_packages

setup(name='QGL',
      version='2.1',
      url='https://github.com/BBN-Q/QGL',
      packages=find_packages(exclude=["tests"]),
      install_requires=[
            "numpy >= 1.11.1",
            "scipy >= 0.17.1",
            "jupyter >= 1.0.0",
            "atom >= 0.4.1",
            "h5py >= 2.6.0",
            "networkx >= 1.11",
            "future >= 0.16",
            "watchdog >= 0.8.3",
            "bokeh >= 0.11",
            "pyyaml >= 3.12"
      ],
      extras_require={"gst": "pygsti>0.9.4"}
)
