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
    def __init__(self, path, lambda_type=None, total_PELE_steps=None):
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
        self.__delta_energies = None
        self.__type = lambda_type

        if (total_PELE_steps is None):
            print("  - LambdaFolder Warning: unknown total number of PELE " +
                  "steps")
        self.__total_PELE_steps = total_PELE_steps

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

    @property
    def type(self):
        return self.__type

    @property
    def total_PELE_steps(self):
        return self.__total_PELE_steps

    def __lt__(self, other):
        if (self.initial_lambda != other.initial_lambda):
            return self.initial_lambda < other.initial_lambda
        else:
            return self.final_lambda < other.final_lambda

    def __eq__(self, other):
        return self.type, self.initial_lambda, self.final_lambda == \
            other.type, other.initial_lambda, other.final_lambda

    def __ne__(self, other):
        return not(self == other)

    def __hash__(self):
        return hash((self.type, self.initial_lambda, self.final_lambda))

    def __str__(self):
        return "Lambda Folder from {} to {}".format(self.initial_lambda,
                                                    self.final_lambda)

    def _getDeltaEnergiesFromReport(self, report):
        energies = []

        report_steps = report.getMetric(co.PP_STEPS_COL)
        report_energies = report.getMetric(co.PP_DELTA_ENERGIES_COL)

        if (len(report_energies) == 0):
            print("  - LambdaFolder Warning: found an empty report file " +
                  "{}".format(report.path))
            return energies

        """
        if (len(report_energies) == 1):
            print("  - LambdaFolder Warning: no steps were accepted in " +
                  "report file {}".format(report.path))
            return energies

        # We discard step 0
        report_steps = report_steps[1:]
        report_energies = report_energies[1:]
        """

        last_accepted_step = int(report_steps[0])
        last_accepted_energy = report_energies[0]

        for step, energy in zip(map(int, report_steps[1:]),
                                report_energies[1:]):
            for i in range(0, step - last_accepted_step):
                energies.append(last_accepted_energy)

            last_accepted_step = int(step)
            last_accepted_energy = energy

        if (self.total_PELE_steps is not None):
            for i in range(0, self.total_PELE_steps - last_accepted_step + 1):
                energies.append(last_accepted_energy)

        return energies

    def getDeltaEnergyValues(self):
        if (self.__delta_energies is not None):
            return self.__delta_energies

        files = get_all_files_from_with_extension(self.path, 'out')

        energies = []

        for file in files:
            report = Report(getPathFromFile(file), getFileFromPath(file),
                            '_'.join(getFileFromPath(file).split('_')[:-1]) +
                            '_',
                            None)

            energies += self._getDeltaEnergiesFromReport(report)

        self.__delta_energies = energies

        return energies

    def getDeltaEnergyValuesByFile(self):
        files = get_all_files_from_with_extension(self.path, 'out')

        energies_by_file = {}

        for file in files:
            report = Report(getPathFromFile(file), getFileFromPath(file),
                            '_'.join(getFileFromPath(file).split('_')[:-1]) +
                            '_',
                            None)

            report_energies = report.getMetric(co.PP_DELTA_ENERGIES_COL)
            if (len(report_energies) == 0):
                print("  - LambdaFolder Warning: found an empty report file " +
                      "{}".format(report.path))
                continue

            # Discard first energy
            report_energies = report_energies[1:]

            energies_by_file[file] = report_energies

        return energies_by_file
