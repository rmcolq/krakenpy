from setuptools import setup, find_packages
import glob
import os
import pkg_resources

from pykraken import __version__, _program

setup(name='pykraken',
      version=__version__,
      packages=find_packages(),
      scripts=[
            'pykraken/taxonomy.py',
            'pykraken/assignment.py',
            'pykraken/report.py'
                ],
      package_data={},
      install_requires=[],
      description='Utility functions to interact with kraken reports and assignment files',
      url='https://github.com/rmcolq/pykraken',
      author='Rachel Colquhoun',
      author_email='rachel.colquhoun@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = pangolin_data.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
