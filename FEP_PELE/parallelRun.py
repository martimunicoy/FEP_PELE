# -*- coding: utf-8 -*-


# Python imports
import os
import sys
import argparse
import copy
# import multiprocessing
# import multiprocessing.pool
from multiprocessing import Pool, current_process
from functools import partial


# FEP_PELE imports
from FreeEnergy.InputFileParser import InputFileParser
from FreeEnergy.CommandsBuilder import CommandsBuilder
from TemplateHandler import Lambda
from Utils.InOut import clear_directory, copyFolder, copySymLink
from Utils.InOut import getLastFolderFromPath
from FreeEnergy import Constants as co


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Class definitions
"""
class NoDaemonProcess(multiprocessing.Process):
    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass

    daemon = property(_get_daemon, _set_daemon)


class NonDaemonPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess
"""


# Function definitions
def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', metavar='PATH', type=str, nargs=1,
                        help='Path to input file')

    args = parser.parse_args()

    path_to_input_file = args.input_file[0]

    return path_to_input_file


def setFile(path_to_file, relative_path):
    try:
        checkFile(path_to_file)
    except NameError:
        path_to_file = relative_path + path_to_file

    return path_to_file


def prepareSettingsForLambda(original_settings, lmb):
    settings = copy.deepcopy(original_settings)

    path = settings.general_path + 'runs/'
    relative_path = '../../'
    if (lmb.type != Lambda.DUAL_LAMBDA):
        path += lmb.type + '/'
        relative_path += '../'
    path += lmb.folder_name + "/"

    clear_directory(path)
    copyFolder(settings.general_path + 'DataLocal',
               path + 'DataLocal')
    copySymLink(settings.general_path + 'Data',
                path + 'Data')
    copySymLink(settings.general_path + 'Documents',
                path + 'Documents')

    settings.setGeneralPath(path)
    settings.setMinimizationPath(
        path + getLastFolderFromPath(original_settings.minimization_path))
    settings.setSimulationPath(
        original_settings.general_path +
        os.path.relpath(original_settings.simulation_path))
    settings.setCalculationPath(
        original_settings.general_path +
        os.path.relpath(original_settings.calculation_path))

    if (settings.initial_template is not None):
        settings.setInitialTemplate(setFile(settings.initial_template,
                                            relative_path))
    if (settings.final_template is not None):
        settings.setFinalTemplate(setFile(settings.final_template,
                                          relative_path))

    if (settings.initial_ligand_pdb is not None):
        settings.setInitialLigandPdb(setFile(settings.initial_ligand_pdb,
                                             relative_path))

    if (settings.final_ligand_pdb is not None):
        settings.setFinalLigandPdb(setFile(settings.final_ligand_pdb,
                                           relative_path))

    if (settings.pp_control_file is not None):
        settings.setPPControlFile(setFile(settings.pp_control_file,
                                          relative_path))

    if (settings.sp_control_file is not None):
        settings.setSPControlFile(setFile(settings.sp_control_file,
                                          relative_path))

    if (settings.min_control_file is not None):
        settings.setMinControlFile(setFile(settings.min_control_file,
                                           relative_path))

    if (settings.sim_control_file is not None):
        settings.setSimControlFile(setFile(settings.sim_control_file,
                                           relative_path))

    if (settings.input_pdb is not None):
        settings.setInputPDB(setFile(settings.input_pdb, relative_path))

    return settings


def parallelCommandRunner(original_settings, lmb):
    pid = current_process().pid

    settings = prepareSettingsForLambda(original_settings,
                                        lmb)

    print(" - Running commands for {} lambda: {:6.4f}".format(lmb.type,
                                                              lmb.value))

    sys.stdout = open(settings.general_path + str(pid) + ".out", "w")
    sys.stderr = open(settings.general_path + str(pid) + ".err", "w")

    commandsBuilder = CommandsBuilder(settings)
    commands = commandsBuilder.createCommands()

    if (len(commands) > 1):
        raise IOError("Only a single command is accepted by parallel runner")

    command = commands[0]

    command.setPID(pid)
    command.setLambdas([lmb, ])

    command.run()


def main():
    path_to_input_file = parseArguments()

    inputFileParser = InputFileParser(path_to_input_file)
    settings = inputFileParser.createSettings()

    lambdasBuilder = Lambda.LambdasBuilder()
    lambdas = lambdasBuilder.buildFromSettings(settings)

    if (len(settings.command_names) > 1):
        raise IOError("Only a single command is accepted by parallel runner")

    max_parallel_commands = settings.number_of_processors
    if (settings.command_names[0] ==
            co.COMMAND_NAMES_DICT["LAMBDAS_SAMPLING"]):
        max_parallel_commands = int(settings.number_of_processors /
                                    settings.parallel_PELE_runs)
        settings.setNumberOfProcessors(settings.parallel_PELE_runs)

    pCommandRunner = partial(parallelCommandRunner, settings)

    with Pool(max_parallel_commands) as pool:
        pool.map(pCommandRunner, lambdas)


if __name__ == '__main__':
    main()
