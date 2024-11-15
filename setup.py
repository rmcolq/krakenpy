from setuptools import setup, find_packages
import glob
import os
import pkg_resources

setup(name='krakenpy',
      version="1.0.0",
      packages=find_packages(),
      scripts=[
            'krakenpy/taxonomy.py',
            'krakenpy/assignment.py',
            'krakenpy/report.py',
            'krakenpy/merge.py'
                ],
      package_data={},
      install_requires=[],
      description='Utility functions to interact with kraken reports and assignment files',
      url='https://github.com/rmcolq/krakenpy',
      author='Rachel Colquhoun',
      author_email='rachel.colquhoun@ed.ac.uk',
      entry_points={"console_scripts": ["krakenpy = krakenpy.__main__:main"]},
      include_package_data=True,
      keywords=[],
      zip_safe=False)
