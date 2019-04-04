# -*- coding: utf-8 -*-


# Python imports
import sys


# FEP_PELE imports
from FEP_PELE.FreeEnergy.Command import Command
from FEP_PELE.Utils.InOut import isThereAPath
from FEP_PELE.FreeEnergy import Constants as co
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
        path = self.settings.general_path + co.CALCULATION_FOLDER

        if (not isThereAPath(path)):
            print("Error: no lambda calculation was found in the expected " +
                  "path {}. You need to run LambdaSimulation ".format(path) +
                  "and a Sampling command before calling a Calculator " +
                  "command. Check your parameters.")
            sys.exit(1)

        lambda_folders = self._getLambdaFolders(path)

        energies = []

        for lambda_folder in lambda_folders:
            energies += lambda_folder.getDeltaEnergyValues()
            print(energies)

        result = Calculators.calculateThermodynamicAverage(energies)
        print(result)
        result = Calculators.zwanzingEquation(result)
        print(result)
