# -*- coding: utf-8 -*-


# Python imports
import os


# FEP_PELE imports
from . import Constants as co
from .CheckPoint import CheckPoint
from .SamplingMethods.SamplingMethodBuilder import SamplingMethodBuilder

from FEP_PELE.Tools.LambdaFolder import LambdaFolder
from FEP_PELE.Tools.PDBTools import PDBParser

from FEP_PELE.PELETools import PELEConstants as pele_co

from FEP_PELE.Utils.InOut import printCommandTitle
from FEP_PELE.Utils.InOut import getFoldersInAPath
from FEP_PELE.Utils.InOut import getLastFolderFromPath

from FEP_PELE.TemplateHandler.AlchemicalTemplateCreator import \
    AlchemicalTemplateCreator
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
        self._settings = settings

        # PELE needs to be working in the general path in order to find the
        # the required Data and Document folders
        os.chdir(self.settings.general_path)

        checkPoint = CheckPoint(settings.general_path +
                                co.CHECKPOINT_NAME,
                                settings)

        try:
            checkPoint.initialize()
        except RuntimeError as e:
            print("  - {} Warning: ".format(self.name) + str(e))

        self._checkPoint = checkPoint

        self._lambdasBuilder = Lambda.LambdasBuilder()

        self._alchemicalTemplateCreator = AlchemicalTemplateCreator(
            settings.initial_template,
            settings.final_template,
            settings.atom_links)

        builder = SamplingMethodBuilder(settings)
        self._s_method = builder.createSamplingMethod()

    @property
    def settings(self):
        return self._settings

    @property
    def checkPoint(self):
        return self._checkPoint

    @property
    def lambdasBuilder(self):
        return self._lambdasBuilder

    @property
    def alchemicalTemplateCreator(self):
        return self._alchemicalTemplateCreator

    @property
    def sampling_method(self):
        return self._s_method

    @property
    def name(self):
        return self._name

    @property
    def label(self):
        return self._label

    @property
    def path(self):
        return self._path

    def _run_with_splitted_lambdas(self):
        output = []

        c_lambdas = self.settings.c_lambdas
        s_lambdas = self.settings.lj_lambdas
        if (len(c_lambdas) < 1):
            c_lambdas = self.settings.lambdas
        if (len(s_lambdas) < 1):
            s_lambdas = self.settings.lambdas

        if (self.alchemicalTemplateCreator.explicit_is_final):
            # Case where the final atomset contains the initial one.
            # Thus, there are atoms that will be created. The best
            # method is to set from a first beggining coulombic lambda
            # to zero to ensure that final atoms will have a null
            # partial charge. Until their Lennard Jones parameters
            # are not the final ones, charges will be zero. Then,
            # progressively will be applied to them.
            #self._lambdasCheckUp(s_lambdas, num=1)
            c_lambda = Lambda.Lambda(0, lambda_type=Lambda.COULOMBIC_LAMBDA)
            output += self._run(s_lambdas, Lambda.STERIC_LAMBDA, num=1,
                                constant_lambda=c_lambda)

            #self._lambdasCheckUp(s_lambdas, num=2)
            s_lambda = Lambda.Lambda(1, lambda_type=Lambda.STERIC_LAMBDA)
            output += self._run(c_lambdas, Lambda.COULOMBIC_LAMBDA, num=2,
                                constant_lambda=s_lambda)
        else:
            # Case where the initial atomset contains the final
            # atomset. So, there are atoms that will disappear.
            # In this way, we need to annihilate first coulombic
            # charges, then, we modify Lennard Jones parameters.
            #self._lambdasCheckUp(s_lambdas, num=1)
            s_lambda = Lambda.Lambda(0, lambda_type=Lambda.STERIC_LAMBDA)
            output += self._run(c_lambdas, Lambda.COULOMBIC_LAMBDA, num=1,
                                constant_lambda=s_lambda)

            #self._lambdasCheckUp(s_lambdas, num=2)
            c_lambda = Lambda.Lambda(1, lambda_type=Lambda.COULOMBIC_LAMBDA)
            output += self._run(s_lambdas, Lambda.STERIC_LAMBDA, num=2,
                                constant_lambda=c_lambda)

        return output

    """
    def _lambdasCheckUp(self, lambdas, num):
        if ((self.settings.sampling_method ==
             co.SAMPLING_METHODS_DICT["OVERLAP"]) or
            (self.settings.sampling_method ==
             co.SAMPLING_METHODS_DICT["DOUBLE_ENDED"])):
            self._lambdasCheckUpWithEdges(lambdas, num)
        else:
            self._lambdasCheckUpWithoutEdges(lambdas, num)

    def _lambdasCheckUpWithEdges(self, lambdas, num):
        edges = (0.0, 1.0)
        indexes = (0, -1)
        positions = (0, len(lambdas))

        lambda_ = lambdas[indexes[num - 1]]
        lambda_value = edges[num - 1]

        if (lambda_ != lambda_value):
            print("  - Warning: adding extra lambda {}".format(lambda_value) +
                  ", required for {} ".format(self.settings.sampling_method) +
                  "sampling")
            lambdas.insert(positions[num - 1], lambda_value)
    """

    def _lambdasCheckUpWithoutEdges(self, lambdas, num):
        edges = (0.0, 1.0)
        indexes = (0, -1)
        positions = (0, len(lambdas) - 1)

        lambda_ = lambdas[indexes[num - 1]]
        lambda_value = edges[num - 1]

        if (lambda_ == lambda_value):
            print("  - Warning: removing extra lambda " +
                  "{}, not compatible with ".format(lambda_value) +
                  "{} sampling".format(self.settings.sampling_method))
            del lambdas[positions[num - 1]]

    def _createAlchemicalTemplate(self, lambda_, constant_lambda,
                                  original_lambda=None):
        path = self.settings.general_path + pele_co.HETEROATOMS_TEMPLATE_PATH

        if (self.alchemicalTemplateCreator.explicit_is_final):
            path += self.settings.final_template_name
        else:
            path += self.settings.initial_template_name

        if (constant_lambda is not None):
            self.alchemicalTemplateCreator.applyLambda(constant_lambda)

        change_bonding_params = True
        if (original_lambda is not None):
            self.alchemicalTemplateCreator.applyLambda(original_lambda,
                                                       change_bonding_params)
            change_bonding_params = False

        self.alchemicalTemplateCreator.applyLambda(lambda_,
                                                   change_bonding_params)
        self.alchemicalTemplateCreator.writeAlchemicalTemplate(path)
        self.alchemicalTemplateCreator.reset()

    def _getLambdaFoldersFrom(self, path, lambda_type=Lambda.DUAL_LAMBDA):
        folders = getFoldersInAPath(path)

        selected_folders = []

        for folder in folders:
            name = getLastFolderFromPath(folder)

            lambda_values = name.split('_')

            if (len(lambda_values) != 2):
                continue

            initial_lambda, final_lambda = lambda_values

            try:
                initial_lambda = float(initial_lambda)
            except ValueError:
                continue

            try:
                final_lambda = float(final_lambda)
            except ValueError:
                continue

            if (initial_lambda > 1) or (initial_lambda < 0):
                continue

            if (final_lambda > 1) or (final_lambda < 0):
                continue

            selected_folders.append(LambdaFolder(folder,
                                                 lambda_type=lambda_type))

        return sorted(selected_folders)

    def _getAtomsToMinimize(self):
        fragment_atoms = self.alchemicalTemplateCreator.getFragmentAtoms()

        atoms_to_minimize = []
        for fragment_atom in fragment_atoms:
            atoms_to_minimize.append(fragment_atom.pdb_atom_name)

        return atoms_to_minimize

    def _getPerturbingLinkId(self):
        if (self.alchemicalTemplateCreator.explicit_is_final):
            pdb_parser = PDBParser(self.settings.final_ligand_pdb)
        else:
            pdb_parser = PDBParser(self.settings.initial_ligand_pdb)

        if (len(pdb_parser.links) == 0):
            print("DoubleWideSampling error: ligand not found in " +
                  "ligand PDB: {}".format(pdb_parser))

        if (len(pdb_parser.links) > 1):
            print("DoubleWideSampling error: found more than one link in " +
                  "ligand PDB: {}".format(pdb_parser))

        ligand_link = pdb_parser.links[0]

        return ligand_link.chain + ':' + str(ligand_link.number)

    def _getAtomIdsToMinimize(self):
        atoms_to_minimize = self._getAtomsToMinimize()

        link_id = self._getPerturbingLinkId()

        atom_ids_to_minimize = [link_id + ':' + i for i in atoms_to_minimize]

        return atom_ids_to_minimize

    def _start(self):
        printCommandTitle(self.label)
        print(" - Sampling method: {}".format(self.sampling_method.name))

    def _finish(self):
        pass

    def _getLambdaFolders(self):
        if (self.settings.splitted_lambdas):
            lambda_folders = self._getLambdaFoldersFrom(
                self.path + '?_' + Lambda.STERIC_LAMBDA + '/',
                lambda_type=Lambda.STERIC_LAMBDA)
            lambda_folders += self._getLambdaFoldersFrom(
                self.path + '?_' + Lambda.COULOMBIC_LAMBDA + '/',
                lambda_type=Lambda.COULOMBIC_LAMBDA)
        else:
            lambda_folders = self._getLambdaFoldersFrom(self.path)

        return lambda_folders

    def _checkLambdaFolders(self, lambda_folders):
        pass
