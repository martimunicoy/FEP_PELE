# -*- coding: utf-8 -*-


# FEP_PELE imports
from .Calculators import calculateThermodynamicAverage
from .Calculators import zwanzigEquation
from .Calculators import calculateStandardDeviationOfMean
from .Calculators import calculateMean
from .Calculators import squaredSum
from .Plotters import dEDistributionPlot

from FEP_PELE.FreeEnergy.Constants import SAMPLING_METHODS_DICT as METHODS_DICT
from FEP_PELE.FreeEnergy.Constants import DIRECTION_LABELS


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


class FEPAnalysis(object):
    def __init__(self, lambda_folders, sampling_method, divisions=10):
        self.lambda_folders = lambda_folders
        self.sampling_method = sampling_method
        self.divisions = divisions
        self.averages = self._calculateAverages()

    def _calculateAverages(self):
        averages = {}

        for lambda_folder in self.lambda_folders:
            energies = lambda_folder.getDeltaEnergyValuesByFile()

            sub_energies = [[] for i in range(0, min(self.divisions,
                                                     len(energies)))]

            for file_index, file_energies in enumerate(energies.values()):
                sub_index = file_index % self.divisions
                sub_energies[sub_index] += file_energies

            averages[lambda_folder] = []

            for sub_energy in sub_energies:
                average = zwanzigEquation(
                    calculateThermodynamicAverage(sub_energy))
                averages[lambda_folder].append(average)

        return averages

    def getDeltaEnergies(self):
        energies = {}

        for lambda_folder in self.lambda_folders:
            energies[lambda_folder] = calculateMean(
                self.averages[lambda_folder])

        return energies

    def getStandardDeviations(self):
        stdevs = {}

        for lambda_folder in self.lambda_folders:
            stdevs[lambda_folder] = calculateStandardDeviationOfMean(
                self.averages[lambda_folder])

        return stdevs

    def _DWSResults(self):
        energies = self.getDeltaEnergies()
        stdevs = self.getStandardDeviations()

        return [sum(energies.values()), ], [squaredSum(stdevs.values()), ]

    def _DESResults(self):
        energies = self.getDeltaEnergies()
        stdevs = self.getStandardDeviations()

        forward_e = []
        reverse_e = []
        forward_sd = []
        reverse_sd = []

        for lambda_folder in self.lambda_folders:
            if (lambda_folder.direction ==
                    DIRECTION_LABELS["FORWARD"]):
                forward_e.append(energies[lambda_folder])
                forward_sd.append(stdevs[lambda_folder])

            elif (lambda_folder.direction ==
                    DIRECTION_LABELS["BACKWARDS"]):
                reverse_e.append(energies[lambda_folder])
                reverse_sd.append(stdevs[lambda_folder])

        return [sum(forward_e), sum(reverse_e)], \
            [squaredSum(forward_e), squaredSum(reverse_e)]

    def getResults(self):
        if (self.sampling_method == METHODS_DICT["DOUBLE_WIDE"]):
            return self._DWSResults()
        elif (self.sampling_method == METHODS_DICT["DOUBLE_ENDED"]):
            return self._DESResults()

    def plotHistogram(self):
        energies = {}
        for lambda_folder in self.lambda_folders:
            key = (lambda_folder.initial_lambda, lambda_folder.final_lambda)
            energies[key] = lambda_folder.getDeltaEnergyValues()

        plotter = dEDistributionPlot(energies, self.averages)
        plotter.plotHistogram(20, range=(-2, 2), facecolor='blue', alpha=0.5)
