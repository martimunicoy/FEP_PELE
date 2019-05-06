# -*- coding: utf-8 -*-


# Python imports
import sys


# FEP_PELE imports
from FEP_PELE.FreeEnergy.Command import Command
from FEP_PELE.FreeEnergy import Constants as co
from FEP_PELE.FreeEnergy.Constants import SAMPLING_METHODS_DICT as METHODS_DICT

from FEP_PELE.Utils.InOut import clear_directory
from FEP_PELE.Utils.InOut import join_splitted_models
from FEP_PELE.Utils.InOut import remove_splitted_models
from FEP_PELE.Utils.InOut import deleteAllFilesWithExtension
from FEP_PELE.Utils.InOut import getFileFromPath
from FEP_PELE.Utils.InOut import create_directory

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

        print(" - Calculating energy differences for each delta lambda")

        clear_directory(self.path)

        if (self.settings.splitted_lambdas):
            delta_energies_info = self._run_with_splitted_lambdas()
        else:
            delta_energies_info = self._run(self.settings.lambdas,
                                       Lambda.DUAL_LAMBDA)

        #join_splitted_models(self.path, co.PDB_OUT_NAME.format('*'))
        #remove_splitted_models(self.path, co.PDB_OUT_NAME.format('*'))
        #deleteAllFilesWithExtension(self.path, 'conf')
        #deleteAllFilesWithExtension(self.path, 'txt')

        if (self.settings.sampling_method == METHODS_DICT["DOUBLE_WIDE"]):
            self._printDWSResults(delta_energies_info)
        elif (self.settings.sampling_method == METHODS_DICT["DOUBLE_ENDED"]):
            self._printDESResults(delta_energies_info)

        self._finish()

    def _run(self, lambdas, lambdas_type, num=0,
             constant_lambda=None):
        delta_energies = []

        atoms_to_minimize = self._getAtomIdsToMinimize()

        runner = PELERunner(self.settings.serial_pele,
                            number_of_processors=1)

        lambdas = self.lambdasBuilder.build(lambdas, lambdas_type)

        for lambda_ in lambdas:
            self._createAlchemicalTemplate(lambda_, constant_lambda)

            initial_energy = self._minimize(runner, lambda_, num)

            final_energies, factors = self._calculateEnergeticDifference(
                lambda_, constant_lambda, num, atoms_to_minimize)

            for final_energy, factor in zip(final_energies, factors):
                delta_energies.append((final_energy - initial_energy, factor))

                print(lambda_, factor, final_energy)

        return delta_energies

    def _minimize(self, runner, lambda_, num):

        clear_directory(self.settings.minimization_path)

        out_path = self._getOutPath(lambda_, num)

        create_directory(out_path)

        self._writeMinimizationControlFile(lambda_, out_path)

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

    def _calculateEnergeticDifference(self, lambda_, constant_lambda, num,
                                      atoms_to_minimize):
        runner = PELERunner(self.settings.serial_pele,
                            number_of_processors=1)

        energies = []
        factors = []

        out_path = self._getOutPath(lambda_, num)

        for shif_lambda in self.sampling_method.getShiftedLambdas(lambda_):
            self._createAlchemicalTemplate(shif_lambda, constant_lambda)

            pdb_path = out_path + co.PDB_OUT_NAME.format(str(lambda_.value) +
                                                         '0')

            logfile_name = out_path + co.LOGFILE_NAME.format(
                str(lambda_.value) + '_' + str(shif_lambda.value))

            lambda_folder = out_path + \
                self._getLambdaFolderName(lambda_, shif_lambda)
            lambda_folder += '/'
            clear_directory(lambda_folder)

            self._preparePDB(pdb_path, lambda_folder, shif_lambda)

            self._writeRecalculationControlFile(
                self.settings.sp_control_file,
                lambda_folder + getFileFromPath(pdb_path),
                logfile_name,
                lambda_folder + co.SINGLE_POINT_CF_NAME.format(''))

            output = runner.run(lambda_folder +
                                co.SINGLE_POINT_CF_NAME.format(''))

            for line in output.split('\n'):
                if line.startswith(pele_co.ENERGY_RESULT_LINE):
                    energies.append(float(line.strip().split()[-1]))
                    break
            else:
                print("Error: energy calculation failed")
                print(output)

            if (lambda_.value < shif_lambda.value):
                factors.append(float(+1.0))
            else:
                factors.append(float(-1.0))

            deleteAllFilesWithExtension(lambda_folder, 'conf')

        deleteAllFilesWithExtension(out_path, 'txt')

        return energies, factors

    def _getOutPath(self, lambda_, num):
        out_path = self.path
        if (lambda_.type != Lambda.DUAL_LAMBDA):
            out_path += str(num) + '_' + lambda_.type + "/"

        return out_path

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
                               co.PDB_OUT_NAME.format(str(lambda_.value) +
                                                      '0'))

        cf_creator.write(self.settings.minimization_path +
                         co.MINIMIZATION_CF_NAME)

    def _writeRecalculationControlFile(self, template_path, pdb_name,
                                       logfile_name,
                                       output_path):
        builder = ControlFileFromTemplateCreator(template_path)

        builder.replaceFlag("INPUT_PDB_NAME", pdb_name)
        builder.replaceFlag("SOLVENT_TYPE", self.settings.solvent_type)
        builder.replaceFlag("LOG_PATH", logfile_name)

        builder.write(output_path)

    def _printDESResults(self, delta_energies_info):
        direct_e = []
        reverse_e = []

        for delta_energy, factor in delta_energies_info:
            if (factor > 0):
                direct_e.append(delta_energy)
            else:
                reverse_e.append(delta_energy)

        print("  - (Direct) Prediction " +
              "{:.2f} kcal/mol".format(sum(direct_e)))

        print("  - (Reverse) Prediction " +
              "{:.2f} kcal/mol".format(sum(reverse_e)))

    def _printDWSResults(self, delta_energies_info):
        result = 0.
        for delta_energy, factor in delta_energies_info:
            result += factor * delta_energy

        print(" - Relative Unbound Free Energy prediction " +
              "{:.2f} kcal/mol".format(result))
