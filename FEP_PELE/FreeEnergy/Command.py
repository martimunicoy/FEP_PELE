# -*- coding: utf-8 -*-


# Python imports
import os


# FEP_PELE imports
from . import Constants as co
from .CheckPoint import CheckPoint
from .SamplingMethods.SamplingMethodBuilder import SamplingMethodBuilder

from FEP_PELE.Tools.LambdaFolder import LambdaFolder
from FEP_PELE.Tools.PDBTools import PDBParser
from FEP_PELE.Tools.PDBTools import PDBModifier
from FEP_PELE.Tools.Math import norm

from FEP_PELE.PELETools import PELEConstants as pele_co

from FEP_PELE.Utils.InOut import printCommandTitle
from FEP_PELE.Utils.InOut import getFoldersInAPath
from FEP_PELE.Utils.InOut import getLastFolderFromPath
from FEP_PELE.Utils.InOut import copyFile
from FEP_PELE.Utils.InOut import getFileFromPath

from FEP_PELE.TemplateHandler.Templates import TemplateOPLS2005
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

        self._ligand_template = \
            self._alchemicalTemplateCreator.explicit_template

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

    @property
    def ligand_template(self):
        return self._ligand_template

    def setPath(self, path):
        self._path = path

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
            # progressively, charges will be introduced.
            c_lambda = Lambda.Lambda(0, lambda_type=Lambda.COULOMBIC_LAMBDA)
            output += self._run(s_lambdas, Lambda.STERIC_LAMBDA, num=1,
                                constant_lambda=c_lambda)

            s_lambda = Lambda.Lambda(1, lambda_type=Lambda.STERIC_LAMBDA)
            output += self._run(c_lambdas, Lambda.COULOMBIC_LAMBDA, num=2,
                                constant_lambda=s_lambda)
        else:
            # Case where the initial atomset contains the final
            # atomset. So, there are atoms that will disappear.
            # In this way, we need to annihilate first coulombic
            # charges, then, we modify Lennard Jones parameters.
            s_lambda = Lambda.Lambda(0, lambda_type=Lambda.STERIC_LAMBDA)
            output += self._run(c_lambdas, Lambda.COULOMBIC_LAMBDA, num=1,
                                constant_lambda=s_lambda)

            c_lambda = Lambda.Lambda(1, lambda_type=Lambda.COULOMBIC_LAMBDA)
            output += self._run(s_lambdas, Lambda.STERIC_LAMBDA, num=2,
                                constant_lambda=c_lambda)

        return output

    def _createAlchemicalTemplate(self, lambda_, constant_lambda, gap=''):
        print("{} - Creating alchemical template".format(gap))

        path = self.settings.general_path + pele_co.HETEROATOMS_TEMPLATE_PATH

        if (self.alchemicalTemplateCreator.explicit_is_final):
            path += self.settings.final_template_name
        else:
            path += self.settings.initial_template_name

        if (constant_lambda is not None):
            print("{}  - Applying {}".format(gap, str(constant_lambda)))
            self.alchemicalTemplateCreator.applyLambda(constant_lambda)

        print("{}  - Applying {}".format(gap, str(lambda_)))
        self.alchemicalTemplateCreator.applyLambda(lambda_)

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

            selected_folders.append(LambdaFolder(
                folder, lambda_type=lambda_type,
                total_PELE_steps=self.settings.total_PELE_steps))

        return sorted(selected_folders)

    def _getAtomsToMinimize(self):
        fragment_atoms = self.alchemicalTemplateCreator.getFragmentAtomNames()

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

    def _getLambdaFolderName(self, lambda_, shifted_lambda):
        return lambda_.folder_name + '_' + shifted_lambda.folder_name

    def _applyMinimizedDistancesTo(self, target_pdb, minimized_pdb):
        # Read distances on minimized_pdb
        min_pdb = PDBParser(minimized_pdb)
        min_link = min_pdb.getLinkWithId(self._getPerturbingLinkId())

        tar_pdb = PDBParser(target_pdb)
        tar_link = tar_pdb.getLinkWithId(self._getPerturbingLinkId())

        template_atoms = self.ligand_template.list_of_atoms
        fragment_bonds = self.ligand_template.get_list_of_fragment_bonds()

        modifier = PDBModifier(tar_pdb)
        modifier.setLinkToModify(tar_link, self.ligand_template)

        core_atoms = self.alchemicalTemplateCreator.getCoreAtoms()

        for ((atom_id1, atom_id2), bond) in fragment_bonds:
            template_atom1 = template_atoms[atom_id1]
            template_atom2 = template_atoms[atom_id2]

            atom1 = min_link.getAtomWithName(template_atom1.pdb_atom_name)
            atom2 = min_link.getAtomWithName(template_atom2.pdb_atom_name)

            length = norm(atom1.coords - atom2.coords)
            f_index = self._getFixedIndex(template_atom1, template_atom2,
                                          core_atoms)

            bond = (atom1.atom_name, atom2.atom_name)

            modifier.modifyBond(bond, length, f_index)

        modifier.write(minimized_pdb)

    def _preparePDB(self, pdb_path, general_path, lambda_, shif_lambda,
                    constant_lambda):
        if ((shif_lambda.type == Lambda.DUAL_LAMBDA) or
                (shif_lambda.type == Lambda.STERIC_LAMBDA)):
            pdb = PDBParser(pdb_path)
            link = pdb.getLinkWithId(self._getPerturbingLinkId())
            modifier = PDBModifier(pdb)
            modifier.setLinkToModify(link, self.ligand_template)

            bonds, lengths, f_indexes = self._getAllAlchemicalBondsInfo(
                lambda_, constant_lambda, link)

            for bond, length, f_index in zip(bonds, lengths, f_indexes):
                modifier.modifyBond(bond, length, f_index)

            modifier.write(general_path + getFileFromPath(pdb_path))

        else:
            copyFile(pdb_path, general_path)

    def _getAllAlchemicalBondsInfo(self, lambda_, constant_lambda, link):
        bonds = []
        lengths = []
        f_indexes = []

        core_atoms = self.alchemicalTemplateCreator.getCoreAtoms()
        template_atoms = self.ligand_template.list_of_atoms

        # Bonds need to be retrived from the current template
        # (may be modified)
        current_template = TemplateOPLS2005(
            self.settings.general_path + pele_co.HETEROATOMS_TEMPLATE_PATH +
            getFileFromPath(self.ligand_template.path_to_template))

        list_of_bonds = current_template.list_of_bonds

        for ((atom_id1, atom_id2), bond) in \
                self.ligand_template.get_list_of_fragment_bonds():
            atom1 = template_atoms[atom_id1]
            atom2 = template_atoms[atom_id2]

            # Select the bond in the template that contains information about
            # the current state of the bond
            bond = list_of_bonds[(atom_id1, atom_id2)]

            bonds.append((atom1.pdb_atom_name, atom2.pdb_atom_name))
            lengths.append(bond.eq_dist)
            f_indexes.append(self._getFixedIndex(atom1, atom2, core_atoms))

            self.alchemicalTemplateCreator.reset()

        return bonds, lengths, f_indexes

    def _getFixedIndex(self, atom1, atom2, core_atoms):
        dist1 = atom1.calculateMinimumDistanceWithAny(core_atoms)
        dist2 = atom2.calculateMinimumDistanceWithAny(core_atoms)

        if (dist1 < dist2):
            return 0
        else:
            return 1

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
