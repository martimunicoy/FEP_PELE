# -*- coding: utf-8 -*-


# Python imports
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


def prepareSettingsForLambda(original_settings, lambda_):
    settings = copy.deepcopy(original_settings)

    path = settings.general_path + 'runs/'
    relative_path = '../../'
    if (lambda_.type != Lambda.DUAL_LAMBDA):
        path += lambda_.type + '/'
        relative_path += '../'
    path += lambda_.folder_name + "/"

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


def setLambdas(command, lambda_):
    command.settings.setLambdas([])
    command.settings.setStericLambdas([])
    command.settings.setCoulombicLambdas([])

    lambdas = [lambda_.value, ]

    # If command is exponential averaging, make pairs of lambdas
    if (command.name == co.COMMAND_NAMES_DICT["EXPONENTIAL_AVERAGING"]):
        if (lambda_.next_lambda is not None):
            lambdas.append(lambda_.next_lambda.value)

    if (lambda_.type == Lambda.DUAL_LAMBDA):
        command.settings.setLambdas(lambdas)
    elif (lambda_.type == Lambda.STERIC_LAMBDA):
        command.settings.setStericLambdas(lambdas)
    elif (lambda_.type == Lambda.COULOMBIC_LAMBDA):
        command.settings.setCoulombicLambdas(lambdas)


def parallelCommandRunner(original_settings, lambda_):
    pid = current_process().pid

    settings = prepareSettingsForLambda(original_settings,
                                        lambda_)

    print(" - Running commands for {} lambda: {:6.4f}".format(lambda_.type,
                                                              lambda_.value))

    sys.stdout = open(settings.general_path + str(pid) + ".out", "w")
    sys.stderr = open(settings.general_path + str(pid) + ".err", "w")

    commandsBuilder = CommandsBuilder(settings)
    commands = commandsBuilder.createCommands()

    for command in commands:
        command.setPath(original_settings.simulation_path)
        command.setPID(pid)
        setLambdas(command, lambda_)
        print(command.settings.lambdas)
        print(command.settings.lj_lambdas)
        print(command.settings.c_lambdas)

        command.run()


def main():
    path_to_input_file = parseArguments()

    inputFileParser = InputFileParser(path_to_input_file)
    settings = inputFileParser.createSettings()

    lambdas = []
    lambdasBuilder = Lambda.LambdasBuilder()

    if (settings.splitted_lambdas):
        lambdas += lambdasBuilder.build(settings.lj_lambdas,
                                        Lambda.STERIC_LAMBDA)

        lambdas += lambdasBuilder.build(settings.c_lambdas,
                                        Lambda.COULOMBIC_LAMBDA)

    else:
        lambdas += lambdasBuilder.build(settings.lj_lambdas,
                                        Lambda.DUAL_LAMBDA)

    number_of_parallel_simulations = 5
    max_parallel_commands = int(settings.number_of_processors /
                                number_of_parallel_simulations)

    settings.setNumberOfProcessors(number_of_parallel_simulations)

    pCommandRunner = partial(parallelCommandRunner, settings)

    with Pool(max_parallel_commands) as pool:
        pool.map(pCommandRunner, lambdas)


if __name__ == '__main__':
    main()
