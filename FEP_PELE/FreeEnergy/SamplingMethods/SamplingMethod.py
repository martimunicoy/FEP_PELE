# -*- coding: utf-8 -*-


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


class SamplingMethod(object):
    def __init__(self, settings):
        self.settings = settings

    def getShiftedLambdas(self, lambda_):
        return self._getShiftedLambdas(lambda_)
