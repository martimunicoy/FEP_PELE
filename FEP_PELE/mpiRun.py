# -*- coding: utf-8 -*-


# Python imports
import sys
import os
import argparse
import copy
from subprocess import check_output, CalledProcessError, STDOUT


# External imports
from mpi4py import MPI


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


def prepareSettingsForLambdas(original_settings, lambdas, lambda_types, path,
                              relative_path):
    settings = copy.deepcopy(original_settings)

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

    settings.setSimulationPath(original_settings.simulation_path)
    settings.setCalculationPath(original_settings.calculation_path)
    settings.setMinimizationPath(path + original_settings.minimization_path)

    d_lambdas = []
    lj_lambdas = []
    c_lambdas = []

    for lambda_, lambda_type in zip(lambdas, lambda_types):
        if (lambda_type == Lambda.DUAL_LAMBDA):
            d_lambdas.append(lambda_)
        elif (lambda_type == Lambda.STERIC_LAMBDA):
            lj_lambdas.append(lambda_)
        elif (lambda_type == Lambda.COULOMBIC_LAMBDA):
            c_lambdas.append(lambda_)

    settings.setLambdas(d_lambdas)
    settings.setStericLambdas(lj_lambdas)
    settings.setCoulombicLambdas(c_lambdas)

    return settings


def parallelCommandRunner(original_settings, lambdas, lambda_types):
    lambda_, lambda_type = lambda_info

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


def prepareFolder(original_settings, path):
    full_clear_directory(path)
    copyFolder(original_settings.general_path + 'DataLocal',
               path + 'DataLocal')
    copySymLink(original_settings.general_path + 'Data',
                path + 'Data')
    copySymLink(original_settings.general_path + 'Documents',
                path + 'Documents')


def getLambdas(settings, rank, size):
    lambdas = []
    lambda_types = []

    if (settings.splitted_lambdas):
        for i1, lambda_ in enumerate(settings.lj_lambdas):
            if (i1 % size == rank):
                lambdas.append(lambda_)
                lambda_types.append(Lambda.STERIC_LAMBDA)
        i1 += 1
        for i2, lambda_ in enumerate(settings.c_lambdas):
            if ((i2 + i1) % size == rank):
                lambdas.append(lambda_)
                lambda_types.append(Lambda.COULOMBIC_LAMBDA)

    else:
        for index, lambda_ in enumerate(settings.lambdas):
            if (index % size == rank):
                lambdas.append(lambda_)
                lambda_types.append(Lambda.DUAL_LAMBDA)

    return lambdas, lambda_types


def mpiRun(original_settings):
    comm = MPI.COMM_WORLD

    size = comm.Get_size()
    rank = comm.Get_rank()

    lambdas, lambda_types = getLambdas(original_settings, rank, size)

    path = original_settings.general_path + 'runs/' + str(rank) + '/'
    relative_path = '../../'

    prepareFolder(original_settings, path)

    settings = prepareSettingsForLambdas(original_settings, lambdas,
                                         lambda_types, path, relative_path)

    settings.setNumberOfProcessors(5)

    settings.write(path + 'mpiRun.inp')

    args = [sys.executable, os.path.dirname(sys.argv[0]) + '/' + 'run.py', path + 'mpiRun.inp']
    print(args)

    try:
        check_output(args, stderr=STDOUT)

    except CalledProcessError as exception:
        print("Error:")
        print(exception.output.decode('utf-8').strip())


def main():
    path_to_input_file = parseArguments()

    inputFileParser = InputFileParser(path_to_input_file)
    settings = inputFileParser.createSettings()

    mpiRun(settings)

    """


    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    number_of_parallel_simulations = 5
    max_parallel_commands = int(settings.number_of_processors /
                                number_of_parallel_simulations)

    settings.setNumberOfProcessors(number_of_parallel_simulations)

    arallelCommandRunner(settings, lambdas, lambda_types)
    """


if __name__ == '__main__':
    main()
