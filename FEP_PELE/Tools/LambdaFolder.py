# -*- coding: utf-8 -*-


# FEP_PELE imports
from FEP_PELE.Utils.InOut import checkPath
from FEP_PELE.Utils.InOut import getLastFolderFromPath
from FEP_PELE.Utils.InOut import getFileFromPath
from FEP_PELE.Utils.InOut import getPathFromFile
from FEP_PELE.Utils.InOut import get_all_files_from_with_extension
from FEP_PELE.FreeEnergy import Constants as co
from FEP_PELE.PELETools.SimulationParser import Report

# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


class LambdaFolder(object):
    def __init__(self, path):
        # Setting default values
        try:
            checkPath(path)
        except NameError:
            raise AttributeError("Invalid path to lambda folder: " +
                                 "{}".format(path))

        name = getLastFolderFromPath(path)

        lambda_values = name.split('_')

        if (len(lambda_values) != 2):
            raise AttributeError("Invalid name of lambda folder: " +
                                 "{}".format(name))

        initial_lambda, final_lambda = lambda_values

        try:
            initial_lambda = float(initial_lambda)
            final_lambda = float(final_lambda)
        except ValueError:
            raise AttributeError("Invalid name of lambda folder: " +
                                 "{}".format(name))

        if (initial_lambda > 1) or (initial_lambda < 0):
            raise AttributeError("Invalid name of lambda folder: " +
                                 "{}".format(name))

        if (final_lambda > 1) or (final_lambda < 0):
            raise AttributeError("Invalid name of lambda folder: " +
                                 "{}".format(name))

        if (initial_lambda == final_lambda):
            raise AttributeError("Invalid name of lambda folder: " +
                                 "{}".format(name))

        self.__path = path
        self.__name = name
        self.__initial_lambda = initial_lambda
        self.__final_lambda = final_lambda

        if (initial_lambda < final_lambda):
            direction_name = "FORWARD"
        else:
            direction_name = "BACKWARDS"

        self.__direction = co.DIRECTION_LABELS[direction_name]
        self.__direction_factor = co.DIRECTION_FACTORS[direction_name]

    @property
    def path(self):
        return self.__path

    @property
    def name(self):
        return self.__name

    @property
    def initial_lambda(self):
        return self.__initial_lambda

    @property
    def final_lambda(self):
        return self.__final_lambda

    @property
    def direction(self):
        return self.__direction

    @property
    def direction_factor(self):
        return self.__direction_factor

    def getDeltaEnergyValues(self):
        files = get_all_files_from_with_extension(self.path, 'out')

        energies = []

        for file in files:
            report = Report(getPathFromFile(file), getFileFromPath(file),
                            '_'.join(getFileFromPath(file).split('_')[:-1]) +
                            '_',
                            None)

            energies += report.getMetric(co.PP_DELTA_ENERGIES_CO)

        return energies
