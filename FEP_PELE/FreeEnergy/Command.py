# -*- coding: utf-8 -*-


# Python imports
import os


# FEP_PELE imports
from . import Constants as co

from FEP_PELE.Tools.LambdaFolder import LambdaFolder

from FEP_PELE.PELETools import PELEConstants as pele_co

from FEP_PELE.Utils.InOut import getFoldersInAPath
from FEP_PELE.Utils.InOut import getLastFolderFromPath

from FEP_PELE.TemplateHandler import Lambda

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

    def _run_with_splitted_lambdas(self, alchemicalTemplateCreator):
        c_lambdas = self.settings.c_lambdas
        s_lambdas = self.settings.lj_lambdas
        if (len(c_lambdas) < 1):
            c_lambdas = self.settings.lambdas
        if (len(s_lambdas) < 1):
            s_lambdas = self.settings.lambdas

        if (alchemicalTemplateCreator.explicit_is_final):
            # Case where the final atomset contains the initial one.
            # Thus, there are atoms that will be created. The best
            # method is to set from a first beggining coulombic lambda
            # to zero to ensure that final atoms will have a null
            # partial charge. Until their Lennard Jones parameters
            # are not the final ones, charges will be zero. Then,
            # progressively will be applied to them.
            c_lambda = Lambda.Lambda(0, lambda_type=Lambda.COULOMBIC_LAMBDA)
            self._run(alchemicalTemplateCreator, s_lambdas,
                      Lambda.STERIC_LAMBDA, num=1, constant_lambda=c_lambda)

            s_lambda = Lambda.Lambda(1, lambda_type=Lambda.STERIC_LAMBDA)
            self._run(alchemicalTemplateCreator, c_lambdas,
                      Lambda.COULOMBIC_LAMBDA, num=2, constant_lambda=s_lambda)
        else:
            # Case where the initial atomset contains the final
            # atomset. So, there are atoms that will disappear.
            # In this way, we need to annihilate first coulombic
            # charges, then, we modify Lennard Jones parameters.
            s_lambda = Lambda.Lambda(0, lambda_type=Lambda.STERIC_LAMBDA)
            self._run(alchemicalTemplateCreator, c_lambdas,
                      Lambda.COULOMBIC_LAMBDA, num=1, constant_lambda=s_lambda)

            c_lambda = Lambda.Lambda(1, lambda_type=Lambda.COULOMBIC_LAMBDA)
            self._run(alchemicalTemplateCreator, s_lambdas,
                      Lambda.STERIC_LAMBDA, num=2, constant_lambda=c_lambda)

    def _createAlchemicalTemplate(self, alchemicalTemplateCreator,
                                  lambda_, constant_lambda):
        path = self.settings.general_path + pele_co.HETEROATOMS_TEMPLATE_PATH

        if (alchemicalTemplateCreator.explicit_is_final):
            path += self.settings.final_template_name
        else:
            path += self.settings.initial_template_name

        if (constant_lambda is not None):
            alchemicalTemplateCreator.applyLambda(constant_lambda)

        alchemicalTemplateCreator.applyLambda(lambda_)
        alchemicalTemplateCreator.writeAlchemicalTemplate(path)
        alchemicalTemplateCreator.reset()

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
