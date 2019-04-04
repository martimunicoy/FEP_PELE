# -*- coding: utf-8 -*-


# Python imports
from subprocess import check_output


# FEP_PELE imports
from FEP_PELE.FreeEnergy import Constants as co
from FEP_PELE.FreeEnergy.Command import Command
from FEP_PELE.TemplateHandler.AlchemicalTemplateCreator import \
    AlchemicalTemplateCreator
from FEP_PELE.Utils.InOut import clear_directory
from FEP_PELE.Utils.InOut import write_lambda_value_to_control_file
from FEP_PELE.Utils.InOut import getFileFromPath
from FEP_PELE.PELETools import PELEConstants as pele_co


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

        clear_directory(self.settings.general_path + co.SIMULATION_FOLDER)

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

            path = self.settings.general_path + co.MINIMIZATION_FOLDER
            clear_directory(path)

            print("  - Initial minimization")
            check_output([self.settings.serial_pele,
                          self.settings.min_control_file])

            print("  - Simulation")

            path = self.settings.general_path + co.SIMULATION_FOLDER + \
                str(lambda_value) + "/"

            clear_directory(path)

            control_file_name = getFileFromPath(self.settings.sim_control_file)

            write_lambda_value_to_control_file(self.settings.sim_control_file,
                                               lambda_value,
                                               path +
                                               control_file_name)
            check_output(["mpirun", "-n",
                          str(self.settings.number_of_processors),
                          "--oversubscribe", self.settings.mpi_pele,
                          path + control_file_name])

            print("   Done")
