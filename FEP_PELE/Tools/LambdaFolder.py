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

        if (len(name) < 2):
            raise AttributeError("Invalid name of lambda folder: " +
                                 "{}".format(name))

        try:
            lambda_value = float(name[:-1])
        except ValueError:
            raise AttributeError("Invalid name of lambda folder: " +
                                 "{}".format(name))

        if (name[-1] not in co.DIRECTION_CHARS):
            raise AttributeError("Invalid name of lambda folder: " +
                                 "{}".format(name))

        self.__path = path
        self.__name = name
        self.__lambda_value = lambda_value
        self.__direction = co.CHAR_TO_DIRECTION[name[-1]]
        self.__direction_factor = float(co.DIRECTION_FACTORS[self.__direction])

    @property
    def path(self):
        return self.__path

    @property
    def name(self):
        return self.__name

    @property
    def lambda_value(self):
        return self.__lambda_value

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
