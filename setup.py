from setuptools import setup, find_packages

setup(name='QGL',
      version='2019.1.2',
      packages=find_packages(exclude=["tests"]),
      url='https://github.com/BBN-Q/QGL',
      download_url='https://github.com/BBN-Q/QGL',
      license="Apache 2.0 License",
      install_requires=[
        "bbndb >= 2019.1.1",
        "numpy >= 1.11.1",
        "scipy >= 0.17.1",
        "networkx >= 1.11",
        "bqplot >= 0.11.5",
        "sqlalchemy >= 1.2.15"
      ],
      description="Quantum Gate Language (QGL) is a domain specific language embedded in python for specifying pulse sequences.",
      long_description_content_type='text/markdown',
      long_description=open('README.md').read(),
      python_requires='>=3.6',
      keywords="quantum qubit experiment configuration gate language"
)

# python setup.py sdist
# python setup.py bdist_wheel
# For testing:
# twine upload --repository-url https://test.pypi.org/legacy/ dist/*
# For distribution:
# twine upload dist/*
# Test with:
# pip install --extra-index-url https://test.pypi.org/simple/ qgl