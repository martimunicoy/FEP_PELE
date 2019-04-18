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

    for energy in energies:
        result += math.exp(- energy * beta)

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
