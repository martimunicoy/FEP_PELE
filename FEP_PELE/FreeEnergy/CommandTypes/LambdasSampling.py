# -*- coding: utf-8 -*-


# Python imports
import sys
import random


# FEP_PELE imports
from FEP_PELE.FreeEnergy.Command import Command
from FEP_PELE.FreeEnergy import Constants as co
from FEP_PELE.FreeEnergy.Constants import SAMPLING_METHODS_DICT as METHODS_DICT

from FEP_PELE.TemplateHandler import Lambda

from FEP_PELE.Utils.InOut import clear_directory
from FEP_PELE.Utils.InOut import create_directory
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
        self._path = self.settings.simulation_path

    def run(self):
        self._start()

        create_directory(self.settings.simulation_path)

        self._run()

        self._finish()

    def _run(self):
        for lmb in self.lambdas:
            if (self.checkPoint.check((self.name, str(lmb.index) +
                                       str(lmb.type) + str(lmb.value)))):
                continue

            writeLambdaTitle(lmb)

            print(" - Creating alchemical template")

            ctt_lmb = self.getConstantLambda(lmb)

            self._createAlchemicalTemplate(lmb, ctt_lmb)

            print(" - Running PELE")

            print("  - Initial minimization")

            self._minimize()

            print("  - Simulation")

            self._simulate(lmb)

            self.checkPoint.save((self.name, str(lmb.index) + str(lmb.type) +
                                  str(lmb.value)))

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

    def _simulate(self, lmb):
        path = self.path
        if (lmb.type != Lambda.DUAL_LAMBDA):
            path += str(lmb.index) + '_' + lmb.type + "/"
        path += str(lmb.value) + "/"

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
        cf_creator.replaceFlag("TOTAL_PELE_STEPS",
                               self.settings.total_PELE_steps)

        cf_creator.write(path + name)
