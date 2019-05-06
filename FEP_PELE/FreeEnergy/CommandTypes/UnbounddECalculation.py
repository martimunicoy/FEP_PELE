# -*- coding: utf-8 -*-


# Python imports
import sys


# FEP_PELE imports
from FEP_PELE.FreeEnergy.Command import Command
from FEP_PELE.FreeEnergy import Constants as co

from FEP_PELE.Utils.InOut import clear_directory

from FEP_PELE.TemplateHandler import Lambda

from FEP_PELE.PELETools.PELERunner import PELERunner
from FEP_PELE.PELETools.ControlFileCreator import \
    ControlFileFromTemplateCreator
from FEP_PELE.PELETools import PELEConstants as pele_co


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Class definitions
class UnbounddECalculation(Command):
    def __init__(self, settings):
        self._name = co.COMMAND_NAMES_DICT["UNBOUND_DE_CALCULATION"]
        self._label = co.COMMAND_LABELS_DICT["UNBOUND_DE_CALCULATION"]
        Command.__init__(self, settings)
        self._path = self.settings.calculation_path + 'unbound/'

    def run(self):
        self._start()

        print(" - Calculating energy difference of both lambda edges")

        clear_directory(self.path)

        c_lambdas = self.settings.c_lambdas
        s_lambdas = self.settings.lj_lambdas
        if (len(c_lambdas) < 1):
            c_lambdas = self.settings.lambdas
        if (len(s_lambdas) < 1):
            s_lambdas = self.settings.lambdas

        first_lambda_value = 0.
        last_lambda_value = 1.

        if (self.settings.splitted_lambdas):
            if (self.alchemicalTemplateCreator.explicit_is_final):
                first_lambda_value = min(s_lambdas)
                last_lambda_value = max(c_lambdas)
            else:
                first_lambda_value = min(c_lambdas)
                last_lambda_value = max(s_lambdas)

        first_lambda = Lambda.Lambda(first_lambda_value,
                                     lambda_type=Lambda.DUAL_LAMBDA)
        last_lambda = Lambda.Lambda(last_lambda_value,
                                    lambda_type=Lambda.DUAL_LAMBDA)

        initial_energy = self._run(first_lambda)
        final_energy = self._run(last_lambda)

        print(" - Relative Unbound Free Energy prediction " +
              "{:.2f} kcal/mol".format(final_energy - initial_energy))

    def _run(self, lambda_):
        runner = PELERunner(self.settings.serial_pele,
                            number_of_processors=1)

        self._createAlchemicalTemplate(lambda_, None)

        return self._minimize(runner, lambda_)

    def _minimize(self, runner, lambda_):

        clear_directory(self.settings.minimization_path)

        self._writeMinimizationControlFile(lambda_, self.path)

        try:
            output = runner.run(self.settings.minimization_path +
                                co.MINIMIZATION_CF_NAME)
        except SystemExit as exception:
            print("UnboundLambdasSimulation error: \n" + str(exception))
            sys.exit(1)

        for line in output.split('\n'):
            if line.startswith(pele_co.ENERGY_RESULT_LINE):
                energy = float(line.strip().split()[-1])
                break
        else:
            print("Error: energy calculation failed")

        return energy

    def _writeMinimizationControlFile(self, lambda_, out_path):
        cf_creator = ControlFileFromTemplateCreator(
            self.settings.min_control_file)

        if (self.alchemicalTemplateCreator.explicit_is_final):
            input_pdb = self.settings.final_ligand_pdb
        else:
            input_pdb = self.settings.initial_ligand_pdb

        cf_creator.replaceFlag("INPUT_PDB_NAME", input_pdb)
        cf_creator.replaceFlag("SOLVENT_TYPE", self.settings.solvent_type)
        cf_creator.replaceFlag("LOG_PATH", self.settings.minimization_path +
                               co.SINGLE_LOGFILE_NAME)
        cf_creator.replaceFlag("TRAJECTORY_PATH",
                               out_path +
                               co.PDB_OUT_NAME.format(str(lambda_.value)))

        cf_creator.write(self.settings.minimization_path +
                         co.MINIMIZATION_CF_NAME)
