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
    def __init__(self, lambda_folders, sampling_method, divisions=10,
                 temperature=298.15):
        self.lambda_folders = lambda_folders
        self.sampling_method = sampling_method
        self.divisions = divisions
        self.temperature = temperature
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
                    calculateThermodynamicAverage(sub_energy,
                                                  self.temperature))
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

        for lambda_folder, energy in energies.items():
            energies[lambda_folder] = energy * lambda_folder.direction_factor

        return sum(energies.values()), squaredSum(stdevs.values())

    def _DESResults(self):
        energies = self.getDeltaEnergies()
        stdevs = self.getStandardDeviations()

        direct_e = []
        direct_sd = []
        reverse_e = []
        reverse_sd = []

        for (lambda_folder, energy), stdev in zip(energies.items(),
                                                  stdevs.values()):
            if (lambda_folder.direction == DIRECTION_LABELS["FORWARD"]):
                direct_e.append(energy)
                direct_sd.append(stdev)
            elif (lambda_folder.direction == DIRECTION_LABELS["BACKWARDS"]):
                reverse_e.append(energy)
                reverse_sd.append(stdev)

        return sum(direct_e), squaredSum(direct_sd), \
            sum(reverse_e), squaredSum(reverse_sd),

    def getResults(self):
        if (self.sampling_method == METHODS_DICT["DOUBLE_WIDE"]):
            return self._DWSResults()
        elif (self.sampling_method == METHODS_DICT["DOUBLE_ENDED"]):
            return self._DESResults()

    def printResults(self):
        if (self.sampling_method == METHODS_DICT["DOUBLE_WIDE"]):
            dE, stdev = self.getResults()

            print("  - Prediction " +
                  "{:.2f} kcal/mol".format(dE))
            print("  - Error " +
                  "{:.3f} kcal/mol".format(stdev))

        elif (self.sampling_method == METHODS_DICT["DOUBLE_ENDED"]):
            d_dE, d_stdev, r_dE, r_stdev = self.getResults()

            print("  - (Direct) Prediction " +
                  "{:.2f} kcal/mol".format(d_dE))
            print("  - (Direct) Error " +
                  "{:.3f} kcal/mol".format(d_stdev))

            print("  - (Reverse) Prediction " +
                  "{:.2f} kcal/mol".format(r_dE))
            print("  - (Reverse) Error " +
                  "{:.3f} kcal/mol".format(r_stdev))

    def plotHistogram(self):
        energies = {}

        for lambda_folder in self.lambda_folders:
            key = (lambda_folder.initial_lambda, lambda_folder.final_lambda,
                   lambda_folder.type)
            energies[key] = lambda_folder.getDeltaEnergyValues()

        plotter = dEDistributionPlot(energies, self.averages)
        plotter.plotHistogram(40, range=(-1, 1), facecolor='blue',
                              alpha=0.5)
