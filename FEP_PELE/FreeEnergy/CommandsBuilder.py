# -*- coding: utf-8 -*-


# FEP_PELE imports
from . import Constants as co


# Python imports
from .CommandTypes.LambdasSimulation import LambdasSimulation
from .CommandTypes.DoubleWideSampling import DoubleWideSampling
from .CommandTypes.ExponentialAveraging import ExponentialAveraging
from .CommandTypes.SolvationFreeEnergyCalculation \
    import SolvationFreeEnergyCalculation


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


class CommandsBuilder(object):
    def __init__(self, settings):
        self.settings = settings
        self.commands_names = settings.command_names

    def createCommands(self):
        commands = []
        for command_name in self.commands_names:
            command = self.getCommandFromName(command_name)
            commands.append(command)

        return commands

    def getCommandFromName(self, command_name):
        if (command_name == co.COMMAND_NAMES_DICT["LAMBDA_SIMULATION"]):
            return LambdasSimulation(self.settings)
        elif (command_name == co.COMMAND_NAMES_DICT["DOUBLE_WIDE_SAMPLING"]):
            return DoubleWideSampling(self.settings)
        elif (command_name == co.COMMAND_NAMES_DICT["EXPONENTIAL_AVERAGING"]):
            return ExponentialAveraging(self.settings)
        elif (command_name == co.COMMAND_NAMES_DICT[
                "SOLVATION_FREE_ENERGY_CALCULATION"]):
            return SolvationFreeEnergyCalculation(self.settings)
        else:
            print("Command {} not recogniced".format(command_name))
            exit(1)
