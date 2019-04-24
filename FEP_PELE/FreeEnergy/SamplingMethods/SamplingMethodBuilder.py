# -*- coding: utf-8 -*-


# Python imports
from .DoubleWideSampling import DoubleWideSampling
from .DoubleEndedSampling import DoubleEndedSampling
from FEP_PELE.FreeEnergy.Constants import SAMPLING_METHODS_DICT as methods_dict


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


class SamplingMethodBuilder(object):
    def __init__(self, settings):
        self.settings = settings
        self.sampling_method_name = settings.sampling_method

    def createSamplingMethod(self):
        if (self.sampling_method_name == methods_dict["DOUBLE_WIDE"]):
            return DoubleWideSampling(self.settings)
        if (self.sampling_method_name == methods_dict["DOUBLE_ENDED"]):
            return DoubleEndedSampling(self.settings)
        else:
            print("Sampling method name {} not recogniced".format(
                self.sampling_method_name))
            exit(1)
