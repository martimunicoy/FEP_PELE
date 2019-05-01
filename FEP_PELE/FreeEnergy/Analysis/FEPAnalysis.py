# -*- coding: utf-8 -*-


# FEP_PELE imports
from .Calculators import calculateThermodynamicAverage
from .Calculators import zwanzigEquation
from .Calculators import calculateStandardDeviationOfMean
from .Calculators import calculateMean
from .Calculators import squaredSum
from .Plotters import dEDistributionPlot


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


class FEPAnalysis(object):
    def __init__(self, lambda_folders, divisions=10):
        self.lambda_folders = lambda_folders
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

    def getResults(self):
        energies = self.getDeltaEnergies()
        stdevs = self.getStandardDeviations()

        return sum(energies.values()), squaredSum(stdevs.values())

    def plotHistogram(self):
        energies = {}
        for lambda_folder in self.lambda_folders:
            key = (lambda_folder.initial_lambda, lambda_folder.final_lambda)
            energies[key] = lambda_folder.getDeltaEnergyValues()

        plotter = dEDistributionPlot(energies, self.averages)
        plotter.plotHistogram(20, range=(-2, 2), facecolor='blue', alpha=0.5)
