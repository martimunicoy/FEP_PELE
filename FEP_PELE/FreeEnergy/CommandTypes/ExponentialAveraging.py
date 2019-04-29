# -*- coding: utf-8 -*-


# Python imports
import sys


# FEP_PELE imports
from FEP_PELE.FreeEnergy import Constants as co
from FEP_PELE.FreeEnergy.Command import Command
from FEP_PELE.FreeEnergy.Analysis import Calculators
from FEP_PELE.FreeEnergy.Analysis import FEPAnalysis

from FEP_PELE.Utils.InOut import isThereAPath


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Class definitions
class ExponentialAveraging(Command):
    def __init__(self, settings):
        self._name = co.COMMAND_NAMES_DICT["EXPONENTIAL_AVERAGING"]
        self._label = co.COMMAND_LABELS_DICT["EXPONENTIAL_AVERAGING"]
        Command.__init__(self, settings)
        self._path = self.settings.calculation_path

    @property
    def name(self):
        return self._name

    @property
    def path(self):
        return self._path

    def run(self):
        self._start()

        if (not isThereAPath(self.path)):
            print("Error: no lambda calculation was found in the expected " +
                  "path {}. You need to run ".format(self.path) +
                  "LambdaSimulation and a Sampling command before calling a " +
                  "Calculator command. Check your parameters.")
            sys.exit(1)

        print(" - Sampling method: {}".format(self.settings.sampling_method))

        print(" - Retrieving lambda folders from {}".format(self.path))

        lambda_folders = self._getLambdaFolders()

        print("  - Checking lambda folders")

        self._checkLambdaFolders(lambda_folders)

        print("  - {} lambda folders were found".format(len(lambda_folders)))

        print(" - Calculating Free Energy change")

        analysis = FEPAnalysis.FEPAnalysis(lambda_folders,
                                           self.settings.sampling_method,
                                           divisions=5)

        print(" - Plotting energetic histogram")

        analysis.plotHistogram()

        dEs, stdevs = analysis.getResults()

        self._printResults(dEs, stdevs)

        self._finish()

    def _printResults(self, dEs, stdevs):
        if (len(dEs) == 1):
            results_type = ("", "")
        else:
            results_type = ("(Direct) ", "(Reverse) ")

        i = 0
        for dE, stdev in zip(dEs, stdevs):
            print(" - {}Relative Free Energy results:".format(results_type[i]))

            print("  - Prediction " +
                  "{:.2f} kcal/mol".format(dE))
            print("  - Error " +
                  "{:.3f} kcal/mol".format(stdev))

            i += 1
