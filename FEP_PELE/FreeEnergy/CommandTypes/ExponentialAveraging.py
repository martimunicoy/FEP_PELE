# -*- coding: utf-8 -*-


# Python imports
import sys


# FEP_PELE imports
from FEP_PELE.FreeEnergy import Constants as co
from FEP_PELE.FreeEnergy.Command import Command
from FEP_PELE.FreeEnergy import Calculators

from FEP_PELE.Utils.InOut import isThereAPath

from FEP_PELE.TemplateHandler import Lambda

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
        Command.__init__(self, settings)

    @property
    def name(self):
        return self._name

    def run(self):
        print("#######################")
        print(" Exponential Averaging")
        print("#######################")
        path = self.settings.calculation_path

        if (not isThereAPath(path)):
            print("Error: no lambda calculation was found in the expected " +
                  "path {}. You need to run LambdaSimulation ".format(path) +
                  "and a Sampling command before calling a Calculator " +
                  "command. Check your parameters.")
            sys.exit(1)

        if (self.settings.splitted_lambdas):
            lambda_folders = self._getLambdaFolders(path + '?_' +
                                                    Lambda.STERIC_LAMBDA +
                                                    '/')
            lambda_folders += self._getLambdaFolders(path + '?_' +
                                                     Lambda.COULOMBIC_LAMBDA +
                                                     '/')
        else:
            lambda_folders = self._getLambdaFolders(path)

        result = 0.

        for lambda_folder in lambda_folders:
            energies = lambda_folder.getDeltaEnergyValues()
            lamda_average = Calculators.calculateThermodynamicAverage(energies)
            lambda_energy = Calculators.zwanzigEquation(lamda_average)
            print(lambda_folder.lambda_value, lambda_energy)
            result += lambda_energy

        print(" - Relative Free Energy prediction " +
              "{:.2f} kcal/mol".format(result))
