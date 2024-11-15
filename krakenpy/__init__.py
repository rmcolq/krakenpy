"""
This file is part of KrakenPy (https://github.com/rmcolq/pykraken).
Copyright 2024 Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""
from pkg_resources import get_distribution

try:
    __version__ = get_distribution("fastafunk").version
except:
    __version__ = "local"

__all__ = ["merge"]

from krakenpy.subcommands import *