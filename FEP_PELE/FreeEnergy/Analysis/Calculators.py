# -*- coding: utf-8 -*-


# Python imports
import math


# FEP_PELE imports
from FEP_PELE.FreeEnergy import Constants as co


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


def calculateThermodynamicAverage(energies, temperature=300):
    result = 0.
    beta = float(1 / co.BOLTZMANN_CONSTANT_IN_KCAL_MOL / temperature)

    #print(energies)

    for energy in energies:
        #print(math.exp(energy), ' ', end='')
        result += math.exp(- energy * beta)
    #print()

    result /= len(energies)

    return result


def zwanzigEquation(thermodynamicAverage, temperature=300):
    kBT = co.BOLTZMANN_CONSTANT_IN_KCAL_MOL * temperature
    return - kBT * math.log(thermodynamicAverage)


def calculateMean(values):
    result = 0.

    for value in values:
        result += float(value)

    return result / len(values)


def calculateStandardDeviation(values):
    if (len(values) == 1):
        return 0

    mean = calculateMean(values)

    summatory = 0.

    for value in values:
        summatory += math.pow(value - mean, 2)

    return math.sqrt(summatory / (len(values) - 1))


def calculateStandardDeviationOfMean(values):
    stdev = calculateStandardDeviation(values)

    return stdev / math.sqrt(len(values))


def squaredSum(values):
    summatory = 0.
    for value in values:
        summatory += math.pow(value, 2)

    return math.sqrt(summatory)
