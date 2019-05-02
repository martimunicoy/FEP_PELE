# -*- coding: utf-8 -*-


# FEP_PELE imports
from FEP_PELE.FreeEnergy.Command import Command


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Class definitions
class EmptyCommand(Command):
    def __init__(self, settings, name="Standard", label="Standard Command",
                 path=""):
        self._name = name
        self._label = label
        Command.__init__(self, settings)
        self._path = path
