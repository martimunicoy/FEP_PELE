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
from Utils.InOut import full_clear_directory, copyFolder, copySymLink


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


def prepareSettingsForLambda(original_settings, lambda_, lambda_type):
    settings = copy.deepcopy(original_settings)

    lambda_ = Lambda.Lambda(lambda_)

    path = settings.general_path + 'runs/'
    relative_path = '../../'
    if (lambda_type != Lambda.DUAL_LAMBDA):
        path += lambda_type + '/'
        relative_path += '../'
    path += lambda_.folder_name + "/"

    full_clear_directory(path)
    copyFolder(settings.general_path + 'DataLocal',
               path + 'DataLocal')
    copySymLink(settings.general_path + 'Data',
                path + 'Data')
    copySymLink(settings.general_path + 'Documents',
                path + 'Documents')

    settings.setGeneralPath(path)

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

    settings.setLambdas([])
    settings.setStericLambdas([])
    settings.setCoulombicLambdas([])

    if (lambda_type == Lambda.DUAL_LAMBDA):
        settings.setLambdas([lambda_.value, ])
    elif (lambda_type == Lambda.STERIC_LAMBDA):
        settings.setStericLambdas([lambda_.value, ])
    elif (lambda_type == Lambda.COULOMBIC_LAMBDA):
        settings.setCoulombicLambdas([lambda_.value, ])

    return settings


def parallelCommandRunner(original_settings, lambda_info):
    lambda_, lambda_type = lambda_info

    pid = current_process().pid

    settings = prepareSettingsForLambda(original_settings,
                                        lambda_, lambda_type)

    print(" - Running commands for lambda: " + str(round(lambda_, 5)))

    sys.stdout = open(settings.general_path + str(pid) + ".out", "w")
    sys.stderr = open(settings.general_path + str(pid) + ".err", "w")

    commandsBuilder = CommandsBuilder(settings)
    commands = commandsBuilder.createCommands()

    for command in commands:
        command.setPath(original_settings.simulation_path)
        command.run()


def main():
    path_to_input_file = parseArguments()

    inputFileParser = InputFileParser(path_to_input_file)
    settings = inputFileParser.createSettings()

    lambdas = []
    lambda_types = []

    if (settings.splitted_lambdas):
        for lambda_ in settings.lj_lambdas:
            lambdas.append(lambda_)
            lambda_types.append(Lambda.STERIC_LAMBDA)
        for lambda_ in settings.c_lambdas:
            lambdas.append(lambda_)
            lambda_types.append(Lambda.COULOMBIC_LAMBDA)

    else:
        for lambda_ in settings.lambdas:
            lambdas.append(lambda_)
            lambda_types.append(Lambda.DUAL_LAMBDA)

    number_of_parallel_simulations = 2
    max_parallel_commands = int(settings.number_of_processors /
                                number_of_parallel_simulations)

    settings.setNumberOfProcessors(number_of_parallel_simulations)

    pCommandRunner = partial(parallelCommandRunner, settings)

    with Pool(max_parallel_commands) as pool:
        pool.map(pCommandRunner, zip(lambdas, lambda_types))


if __name__ == '__main__':
    main()
