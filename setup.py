from setuptools import setup, find_packages

setup(name = 'QGL',
      version = '2.1',
      url = 'https://github.com/BBN-Q/QGL',
      packages = find_packages(exclude=["tests"]),
      package_data = {'QGL': ['config*.json']}
     )
