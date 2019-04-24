# -*- coding: utf-8 -*-


# FEP_PELE imports
from .SamplingMethod import SamplingMethod

from FEP_PELE.FreeEnergy import Constants as co

from FEP_PELE.TemplateHandler import Lambda


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Class definitions
class DoubleWideSampling(SamplingMethod):
    def __init__(self, settings):
        self._name = co.SAMPLING_METHODS_DICT["DOUBLE_WIDE"]
        SamplingMethod.__init__(self, settings)

    def _getShiftedLambdas(self, lambda_):
        shifted_lambdas = []

        previous_lambda = lambda_.previous_lambda
        if (previous_lambda is None):
            shifted_value = float(0.0)
        else:
            shifted_value = float(lambda_.value -
                                  (lambda_.value - previous_lambda.value) /
                                  2.0)

        shifted_lambdas.append(Lambda.Lambda(shifted_value,
                                             lambda_type=lambda_.type))

        next_lambda = lambda_.next_lambda
        if (next_lambda is None):
            shifted_value = float(1.0)
        else:
            shifted_value = float(lambda_.value +
                                  (next_lambda.value - lambda_.value) /
                                  2.0)

        shifted_lambdas.append(Lambda.Lambda(shifted_value,
                                             lambda_type=lambda_.type))

        return shifted_lambdas
