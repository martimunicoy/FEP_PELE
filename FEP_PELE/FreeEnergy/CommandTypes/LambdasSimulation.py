# -*- coding: utf-8 -*-


# Python imports
import sys
import random


# FEP_PELE imports
from FEP_PELE.FreeEnergy.Command import Command
from FEP_PELE.FreeEnergy import Constants as co

from FEP_PELE.TemplateHandler.AlchemicalTemplateCreator import \
    AlchemicalTemplateCreator

from FEP_PELE.Utils.InOut import clear_directory
from FEP_PELE.Utils.InOut import full_clear_directory
from FEP_PELE.Utils.InOut import write_lambda_value_to_control_file
from FEP_PELE.Utils.InOut import getFileFromPath
from FEP_PELE.Utils.InOut import writeLambdaTitle

from FEP_PELE.PELETools import PELEConstants as pele_co
from FEP_PELE.PELETools.PELERunner import PELERunner
from FEP_PELE.PELETools.ControlFileCreator import \
    ControlFileFromTemplateCreator


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
        print("####################")
        print(" Lambda Simulations")
        print("####################")

        alchemicalTemplateCreator = AlchemicalTemplateCreator(
            self.settings.initial_template,
            self.settings.final_template,
            self.settings.atom_links)

        # @TODO: add capacity to restart without losing previous information
        # Clear all directories
        full_clear_directory(self.settings.simulation_path)
        full_clear_directory(self.settings.calculation_path)
        full_clear_directory(self.settings.minimization_path)

        for lambda_value in self.settings.lambdas:
            writeLambdaTitle(lambda_value)

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

        self._writeMinimizationControlFile()

        runner = PELERunner(self.settings.serial_pele,
                            number_of_processors=1)

        try:
            runner.run(self.settings.minimization_path +
                       co.MINIMIZATION_CF_NAME)
        except SystemExit as exception:
            print("LambdasSimulation error: \n" + str(exception))
            sys.exit(1)

    def _simulate(self, lambda_value):
        path = self.settings.simulation_path + str(lambda_value) + "/"

        control_file_name = getFileFromPath(self.settings.sim_control_file)

        clear_directory(path)

        self._writeSimulationControlFile(path, control_file_name)

        runner = PELERunner(
            self.settings.mpi_pele,
            number_of_processors=self.settings.number_of_processors)

        try:
            runner.run(path + control_file_name)
        except SystemExit as exception:
            print("LambdasSimulation error: \n" + str(exception))
            sys.exit(1)

    def _writeMinimizationControlFile(self):
        cf_creator = ControlFileFromTemplateCreator(
            self.settings.min_control_file)

        cf_creator.replaceFlag("INPUT_PDB_NAME", self.settings.input_pdb)
        cf_creator.replaceFlag("SOLVENT_TYPE", self.settings.solvent_type)
        cf_creator.replaceFlag("LOG_PATH", self.settings.minimization_path +
                               co.SINGLE_LOGFILE_NAME)
        cf_creator.replaceFlag("TRAJECTORY_PATH",
                               self.settings.minimization_path +
                               getFileFromPath(self.settings.input_pdb))

        cf_creator.write(self.settings.minimization_path +
                         co.MINIMIZATION_CF_NAME)

    def _writeSimulationControlFile(self, path, name):
        cf_creator = ControlFileFromTemplateCreator(
            self.settings.sim_control_file)

        cf_creator.replaceFlag("INPUT_PDB_NAME",
                               self.settings.minimization_path +
                               getFileFromPath(self.settings.input_pdb))
        cf_creator.replaceFlag("SOLVENT_TYPE", self.settings.solvent_type)
        cf_creator.replaceFlag("LOG_PATH", path + co.SINGLE_LOGFILE_NAME)
        cf_creator.replaceFlag("REPORT_PATH", path + co.SINGLE_REPORT_NAME)
        cf_creator.replaceFlag("TRAJECTORY_PATH", path +
                               co.SINGLE_TRAJECTORY_NAME)
        cf_creator.replaceFlag("SEED", random.randint(0, 999999))

        cf_creator.write(path + name)
