# -*- coding: utf-8 -*-


# Python imports
import sys


# FEP_PELE imports
from FEP_PELE.FreeEnergy.Command import Command
from FEP_PELE.TemplateHandler.AlchemicalTemplateCreator import \
    AlchemicalTemplateCreator
from FEP_PELE.Utils.InOut import clear_directory
from FEP_PELE.Utils.InOut import full_clear_directory
from FEP_PELE.Utils.InOut import write_lambda_value_to_control_file
from FEP_PELE.Utils.InOut import getFileFromPath
from FEP_PELE.PELETools import PELEConstants as pele_co
from FEP_PELE.PELETools.PELERunner import PELERunner


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Class definitions
class LambdasSimulation(Command):
    def __init__(self, settings):
        Command.__init__(self, settings)

    def run(self):
        alchemicalTemplateCreator = AlchemicalTemplateCreator(
            self.settings.initial_template,
            self.settings.final_template,
            self.settings.atom_links)

        full_clear_directory(self.settings.simulation_path)

        for lambda_value in self.settings.lambdas:
            print("##############")
            print(" Lambda: " + str(lambda_value))
            print("##############")
            print(" - Creating alchemical template")

            alchemicalTemplateCreator.create(
                lambda_value,
                self.settings.general_path +
                pele_co.HETEROATOMS_TEMPLATE_PATH +
                self.settings.final_template_name)

            print("   Done")

            print(" - Running PELE")

            print("  - Initial minimization")

            self._minimize()

            print("  - Simulation")

            self._simulate(lambda_value)

            print("   Done")

    def _minimize(self):
        path = self.settings.minimization_path

        clear_directory(path)

        runner = PELERunner(self.settings.serial_pele,
                            number_of_processors=1)

        try:
            runner.run(self.settings.min_control_file)
        except SystemExit as exception:
            print("LambdasSimulation error: \n" + str(exception))
            sys.exit(1)

    def _simulate(self, lambda_value):
        path = self.settings.simulation_path + str(lambda_value) + "/"

        clear_directory(path)

        runner = PELERunner(
            self.settings.mpi_pele,
            number_of_processors=self.settings.number_of_processors)

        control_file_name = getFileFromPath(self.settings.sim_control_file)

        write_lambda_value_to_control_file(self.settings.sim_control_file,
                                           lambda_value,
                                           path + control_file_name)

        try:
            runner.run(path + control_file_name)
        except SystemExit as exception:
            print("LambdasSimulation error: \n" + str(exception))
            sys.exit(1)
