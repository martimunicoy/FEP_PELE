# -*- coding: utf-8 -*-


# Python imports
import sys


# FEP_PELE imports
from FEP_PELE.FreeEnergy import Constants as co
from FEP_PELE.FreeEnergy.Command import Command
from FEP_PELE.Utils.InOut import clear_directory
from FEP_PELE.Utils.InOut import copyFile
from FEP_PELE.Utils.InOut import getFileFromPath
from FEP_PELE.PELETools import PELEConstants as pele_co
from FEP_PELE.PELETools.PELERunner import PELERunner
from FEP_PELE.PELETools.ControlFileCreator import \
    ControlFileFromTemplateCreator


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Class definitions
class SolvationFreeEnergyCalculation(Command):
    def __init__(self, settings):
        Command.__init__(self, settings)

    def run(self):
        print("###################################")
        print(" Solvation Free Energy Calculation")
        print("###################################")

        clear_directory(self.settings.minimization_path)

        runner = PELERunner(self.settings.serial_pele,
                            number_of_processors=1)

        copyFile(self.settings.initial_template,
                 pele_co.HETEROATOMS_TEMPLATE_PATH)

        print(" - Calculating Free Energy of initial ligand with PELE")

        initial_ligand_energy = self._calculateSolvationFreeEnergy(
            runner, self.settings.initial_ligand_pdb)

        print("   Done")

        copyFile(self.settings.final_template,
                 pele_co.HETEROATOMS_TEMPLATE_PATH)

        print(" - Calculating Free Energy of final ligand with PELE")

        final_ligand_energy = self._calculateSolvationFreeEnergy(
            runner, self.settings.final_ligand_pdb)

        print("   Done")

        result = final_ligand_energy - initial_ligand_energy

        print(" - Relative Solvation Free Energy prediction " +
              "{:.2f} kcal/mol".format(result))

    def _calculateSolvationFreeEnergy(self, runner, pdb_path):
        print("  - Minimizing and calculating energy in vacuum")

        self._writeControlFile(pdb_path, pele_co.VACUUM_TYPE_NAME)

        vacuum_energy = self._calculateEnergy(runner)

        print("  - Minimizing and calculating energy in " +
              "{}".format(self.settings.solvent_type))

        self._writeControlFile(pdb_path, self.settings.solvent_type)

        solvated_energy = self._calculateEnergy(runner)

        return solvated_energy - vacuum_energy

    def _writeControlFile(self, input_pdb_name, solvent_type):
        cf_creator = ControlFileFromTemplateCreator(
            self.settings.pp_control_file)

        cf_creator.replaceFlag("INPUT_PDB_NAME", input_pdb_name)
        cf_creator.replaceFlag("SOLVENT_TYPE", solvent_type)
        cf_creator.replaceFlag("LOG_PATH", self.settings.minimization_path +
                               co.LOGFILE_NAME.format(''))
        cf_creator.replaceFlag("TRAJECTORY_PATH",
                               self.settings.minimization_path +
                               getFileFromPath(input_pdb_name) + '.pdb')

        cf_creator.write(self.settings.minimization_path +
                         co.MINIMIZATION_CF_NAME)

    def _calculateEnergy(self, runner):
        try:
            output = runner.run(self.settings.minimization_path +
                                co.MINIMIZATION_CF_NAME)

        except SystemExit as exception:
            print("SolvationFreeEnergyCalculation error: \n" + str(exception))
            sys.exit(1)

        for line in output.split('\n'):
            if line.startswith(pele_co.ENERGY_RESULT_LINE):
                energy = float(line.strip().split()[-1])
                break
            elif line.startswith(pele_co.ENERGY_RESULT_LINE_IN_VACUUM):
                energy = float(line.strip().split()[-1])
                break
        else:
            print("SolvationFreeEnergyCalculation Error: energy calculation " +
                  "failed")
            print(output)

        return energy
