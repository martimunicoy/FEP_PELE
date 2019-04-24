# -*- coding: utf-8 -*-


# Python imports
import sys
import random


# FEP_PELE imports
from FEP_PELE.FreeEnergy.Command import Command
from FEP_PELE.FreeEnergy import Constants as co

from FEP_PELE.TemplateHandler import Lambda
from FEP_PELE.TemplateHandler.AlchemicalTemplateCreator import \
    AlchemicalTemplateCreator

from FEP_PELE.Utils.InOut import clear_directory
from FEP_PELE.Utils.InOut import full_clear_directory
from FEP_PELE.Utils.InOut import getFileFromPath
from FEP_PELE.Utils.InOut import writeLambdaTitle

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
class LambdasSampling(Command):
    def __init__(self, settings):
        self._name = co.COMMAND_NAMES_DICT["LAMBDAS_SAMPLING"]
        self._label = co.COMMAND_LABELS_DICT["LAMBDAS_SAMPLING"]
        Command.__init__(self, settings)

    @property
    def name(self):
        return self._name

    def run(self):
        self._start()

        alchemicalTemplateCreator = AlchemicalTemplateCreator(
            self.settings.initial_template,
            self.settings.final_template,
            self.settings.atom_links)

        if (not self.settings.restart):
            full_clear_directory(self.settings.simulation_path)
            full_clear_directory(self.settings.calculation_path)
            full_clear_directory(self.settings.minimization_path)
        else:
            clear_directory(self.settings.simulation_path)
            clear_directory(self.settings.calculation_path)
            clear_directory(self.settings.minimization_path)

        if (self.settings.splitted_lambdas):
            self._run_with_splitted_lambdas(alchemicalTemplateCreator)
        else:
            lambdas = self.settings.lambdas
            self._lambdasCheckUp(lambdas, num=1)
            self._lambdasCheckUp(lambdas, num=2)
            self._run(alchemicalTemplateCreator, lambdas, Lambda.DUAL_LAMBDA)

        self._finish()

    def _run(self, alchemicalTemplateCreator, lambdas, lambdas_type, num=0,
             constant_lambda=None):
        for lambda_ in Lambda.IterateOverLambdas(lambdas, lambdas_type):
            if (self.checkPoint.check((self.name, str(num) +
                                       str(lambda_.type) +
                                       str(lambda_.value)))):
                continue

            writeLambdaTitle(lambda_)

            print(" - Creating alchemical template")

            self._createAlchemicalTemplate(alchemicalTemplateCreator,
                                           lambda_, constant_lambda)

            print("   Done")

            print(" - Running PELE")

            print("  - Initial minimization")

            self._minimize()

            print("  - Simulation")

            self._simulate(lambda_, num)

            print("   Done")

            self.checkPoint.save((self.name, str(num) + str(lambda_.type) +
                                  str(lambda_.value)))

        return []

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

    def _simulate(self, lambda_, num):
        path = self.settings.simulation_path
        if (lambda_.type != Lambda.DUAL_LAMBDA):
            path += str(num) + '_' + lambda_.type + "/"
        path += str(lambda_.value) + "/"

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

    def starting(self):
        print("#################")
        print(" Lambda Sampling")
        print("#################")

        pr