# -*- coding: utf-8 -*-


# Python imports
import os


# FEP_PELE imports
from . import Constants as co
from FEP_PELE.Tools.LambdaFolder import LambdaFolder
from FEP_PELE.Utils.InOut import getFoldersInAPath
from FEP_PELE.Utils.InOut import getLastFolderFromPath


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Class definitions
class Command(object):
    def __init__(self, settings):
        self.settings = settings

        # PELE needs to be working in the general path in order to find the
        # the required Data and Document folders
        os.chdir(self.settings.general_path)

    def _get_delta_lambda(self, position):
        if (position == 0):
            return round(abs(self.settings.lambdas[position]), 3)
        else:
            return round(abs((self.settings.lambdas[position - 1] -
                             self.settings.lambdas[position]) / 2.), 3)

    def _getLambdaFolders(self, path):
        folders = getFoldersInAPath(path)

        selected_folders = []

        for folder in folders:
            name = getLastFolderFromPath(folder)

            if (len(name) < 2):
                continue

            lambda_value = name[:-1]
            direction = name[-1]

            try:
                lambda_value = float(lambda_value)
            except ValueError:
                continue

            if (lambda_value > 1) or (lambda_value < 0):
                continue

            if (direction not in co.DIRECTION_CHARS):
                continue

            selected_folders.append(LambdaFolder(folder))

        return selected_folders
