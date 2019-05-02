# -*- coding: utf-8 -*-


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Class definitions
class SamplingMethod(object):
    def __init__(self, settings):
        self.settings = settings

    @property
    def name(self):
        return self._name

    def getShiftedLambdas(self, lambda_):
        return self._getShiftedLambdas(lambda_)
