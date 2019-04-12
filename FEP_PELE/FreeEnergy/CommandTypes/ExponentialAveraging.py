# -*- coding: utf-8 -*-


# Python imports
import sys


# FEP_PELE imports
from FEP_PELE.FreeEnergy.Command import Command
from FEP_PELE.Utils.InOut import isThereAPath
from FEP_PELE.FreeEnergy import Calculators

# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Class definitions
class ExponentialAveraging(Command):
    def __init__(self, settings):
        Command.__init__(self, settings)

    def run(self):
        path = self.settings.calculation_path

        if (not isThereAPath(path)):
            print("Error: no lambda calculation was found in the expected " +
                  "path {}. You need to run LambdaSimulation ".format(path) +
                  "and a Sampling command before calling a Calculator " +
                  "command. Check your parameters.")
            sys.exit(1)

        lambda_folders = self._getLambdaFolders(path)

        results = []

        for lambda_folder in lambda_folders:
            energies = lambda_folder.getDeltaEnergyValues()
            lamda_result = Calculators.calculateThermodynamicAverage(energies)
            results.append(Calculators.zwanzigEquation(lamda_result))

        result = Calculators.calculateMean(results)

        print("##############")
        print("   Results")
        print("##############")
        print(" - Relative Free Energy prediction " +
              "{:.2f} kcal/mol".format(result))
