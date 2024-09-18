from setuptools import setup, find_packages
import glob
import os
import pkg_resources

from kraken_utils import __version__, _program

setup(name='kraken_utils',
      version=__version__,
      packages=find_packages(),
      scripts=[
            'kraken_utils/taxonomy.py',
            'kraken_utils/assignment.py',
            'kraken_utils/report.py'
                ],
      package_data={},
      install_requires=[],
      description='Utility functions to interact with kraken reports and assignment files',
      url='https://github.com/rmcolq/kraken_utils',
      author='Rachel Colquhoun',
      author_email='rachel.colquhoun@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = kraken_utils.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
