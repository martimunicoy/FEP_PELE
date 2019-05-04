# -*- coding: utf-8 -*-


# Python imports
import sys
from multiprocessing import Pool, current_process
from functools import partial


# FEP_PELE imports
from FEP_PELE.FreeEnergy import Constants as co
from FEP_PELE.FreeEnergy.Command import Command
from FEP_PELE.FreeEnergy.Analysis.Checkers import checkModelCoords

from FEP_PELE.TemplateHandler import Lambda

from FEP_PELE.PELETools import PELEConstants as pele_co
from FEP_PELE.PELETools.SimulationParser import Simulation
from FEP_PELE.PELETools.PELERunner import PELERunner
from FEP_PELE.PELETools.ControlFileCreator import \
    ControlFileFromTemplateCreator

from FEP_PELE.Utils.InOut import create_directory
from FEP_PELE.Utils.InOut import clear_directory
from FEP_PELE.Utils.InOut import write_energies_report
from FEP_PELE.Utils.InOut import join_splitted_models
from FEP_PELE.Utils.InOut import remove_splitted_models
from FEP_PELE.Utils.InOut import writeLambdaTitle
from FEP_PELE.Utils.InOut import copyFile
from FEP_PELE.Utils.InOut import getFileFromPath

from FEP_PELE.Tools.PDBTools import PDBParser
from FEP_PELE.Tools.PDBTools import PDBModifier

# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Class definitions
class dECalculation(Command):
    def __init__(self, settings):
        self._name = co.COMMAND_NAMES_DICT["DE_CALCULATION"]
        self._label = co.COMMAND_LABELS_DICT["DE_CALCULATION"]
        Command.__init__(self, settings)
        self._path = self.settings.calculation_path

    def run(self):
        self._start()

        clear_directory(self.path)

        if (self.settings.splitted_lambdas):
            self._run_with_splitted_lambdas()
        else:
            self._run(self.settings.lambdas)

        self._finish()

    def _run(self, lambdas, lambdas_type=Lambda.DUAL_LAMBDA, num=0,
             constant_lambda=None):

        lambdas = self.lambdasBuilder.build(lambdas, lambdas_type)

        for lambda_ in lambdas:

            if (self.checkPoint.check((self.name, str(num) +
                                       str(lambda_.type) +
                                       str(lambda_.value)))):
                continue

            writeLambdaTitle(lambda_)

            clear_directory(self.path + co.MODELS_FOLDER)

            input_path = self.settings.simulation_path
            if (lambda_.type != Lambda.DUAL_LAMBDA):
                input_path += str(num) + '_' + lambda_.type + "/"
            input_path += str(lambda_.value) + "/"

            print(" - Splitting PELE models")
            simulation = Simulation(input_path, sim_type="PELE",
                                    report_name="report_",
                                    trajectory_name="trajectory_",
                                    logfile_name="logFile_")
            simulation.getOutputFiles()

            with Pool(self.settings.number_of_processors) as pool:
                pool.map(self._parallelTrajectoryWriterLoop,
                         simulation.iterateOverReports)

            """
            # Inactivate bad models
            for model_info_chunks in models_to_discard:
                if (len(model_info_chunks) == 0):
                    continue
                for model_info in model_info_chunks:
                    print("  - Removing bad model {}".format(model_info[2]))
                    for report in simulation.iterateOverReports:
                        if (report.name == model_info[0]):
                            report.models.inactivate(model_info[1])
            """

            # ---------------------------------------------------------------------

            for shif_lambda in self.sampling_method.getShiftedLambdas(lambda_):
                print(" - Applying delta lambda " +
                      str(round(shif_lambda.value - lambda_.value, 5)))

                print("  - Creating alchemical template")

                self._createAlchemicalTemplate(shif_lambda, constant_lambda)

                general_path = self._getGeneralPath(lambda_, shif_lambda, num)

                print("  - Preparing PDB")

                for report_file in simulation.iterateOverReports:
                    for model_id in range(
                            0, report_file.trajectory.models.number):
                        pdb_path = self.path + co.MODELS_FOLDER + \
                            str(model_id) + '-' + report_file.trajectory.name

                        self._preparePDB(pdb_path, general_path, shif_lambda)

                # -----------------------------------------------------------------

                print("  - Calculating energetic differences")

                atoms_to_minimize = self._getAtomIdsToMinimize()

                parallelLoop = partial(self._parallelPELEMinimizerLoop,
                                       shif_lambda, atoms_to_minimize,
                                       general_path)

                with Pool(self.settings.number_of_processors) as pool:
                    pool.map(parallelLoop, simulation.iterateOverReports)

            self.checkPoint.save((self.name, str(num) + str(lambda_.type) +
                                  str(lambda_.value)))

            clear_directory(self.path)

        return []

    def _parallelTrajectoryWriterLoop(self, report_file):
        """
        models_to_discard = []
        """

        for model_id in range(0, report_file.trajectory.models.number):

            model_name = self.path + co.MODELS_FOLDER + str(model_id) + \
                '-' + report_file.trajectory.name

            report_file.trajectory.writeModel(model_id, model_name)

        """
            if (self.settings.safety_check):
                model_okay = checkModelCoords(
                    model_name,
                    self.alchemicalTemplateCreator.getFragmentAtoms())

                if (not model_okay):
                    models_to_discard.append((report_file.name, model_id,
                                              model_name))

        return models_to_discard
        """

    def _preparePDB(self, pdb_path, general_path, shif_lambda):
        if ((shif_lambda.type == Lambda.DUAL_LAMBDA) or
                (shif_lambda.type == Lambda.STERIC_LAMBDA)):
            pdb = PDBParser(pdb_path)
            link = pdb.getLinkWithId(self._getPerturbingLinkId())
            modifier = PDBModifier(pdb)
            modifier.setLinkToModify(link, self.ligand_template)

            bonds, lengths, f_indexes = self._getAllAlchemicalBondsInfo()

            for bond, length, f_index in zip(bonds, lengths, f_indexes):
                print("Modifying bond:", bond, length, f_index)
                modifier.modifyBond(bond, length, f_index)

            modifier.write(general_path + getFileFromPath(pdb_path))

        else:
            copyFile(pdb_path, general_path)

    def _getAllAlchemicalBondsInfo(self):
        bonds = []
        lengths = []
        f_indexes = []

        core_atoms = self.alchemicalTemplate.getCoreAtoms()
        template_atoms = self.ligand_template.list_of_atoms

        for ((atom_id1, atom_id2), bond) in \
                self.ligand_template.get_list_of_fragment_bonds():
            atom1 = template_atoms[atom_id1]
            atom2 = template_atoms[atom_id2]

            bonds.append((atom1.pdb_atom_name, atom2.pdb_atom_name))
            lengths.append(bond.spring)
            f_indexes.append(self._getFixedAtom(atom1, atom2, core_atoms))

        return bonds, lengths, f_indexes

    def _getFixedIndex(self, atom1, atom2, core_atoms):
        dist1 = atom1.calculateMinimumDistanceWithAny(core_atoms)
        dist2 = atom2.calculateMinimumDistanceWithAny(core_atoms)

        if (dist1 < dist2):
            return 0
        else:
            return 1

    def _parallelPELEMinimizerLoop(self, shifted_lambda, atoms_to_minimize,
                                   general_path, report_file):
        create_directory(general_path)

        pid = current_process().pid

        energies = []
        rmsds = []

        for model_id, active in enumerate(report_file.models):
            # Set initial variables
            file_name = str(model_id) + '-' + report_file.trajectory.name
            original_pdb = self.path + co.MODELS_FOLDER + file_name
            shifted_pdb = general_path + file_name
            logfile_name = self.path + co.LOGFILE_NAME.format(pid)

            # Define new PELERunner
            runner = PELERunner(self.settings.serial_pele,
                                number_of_processors=1)

            # In case bad model was previously removed
            """
            if (not active):
                print("  - Warning: skipping bad model from path: " +
                      "\'{}\'".format(shifted_pdb))
                continue
            """

            # Write recalculation control file
            self._writeRecalculationControlFile(
                self.settings.sp_control_file,
                shifted_pdb,
                logfile_name,
                atoms_to_minimize,
                self.path + co.SINGLE_POINT_CF_NAME.format(pid))

            # Run PELE and extract energy prediction
            energies.append(self._getPELEEnergyPrediction(runner, pid))

            # Calculate RMSD between original pdb and shifted one
            rmsd = self._calculateRMSD(original_pdb, shifted_pdb)
            rmsds.append(rmsd)

        # Write trajectories and reports
        write_energies_report(general_path, report_file, energies, rmsds)
        join_splitted_models(general_path, "*-" + report_file.trajectory.name)

        # Clean temporal files
        remove_splitted_models(general_path,
                               "*-" + report_file.trajectory.name)

    def _writeRecalculationControlFile(self, template_path, pdb_name,
                                       logfile_name, atoms_to_minimize,
                                       output_path):
        builder = ControlFileFromTemplateCreator(template_path)

        builder.replaceFlag("INPUT_PDB_NAME", pdb_name)
        builder.replaceFlag("SOLVENT_TYPE", self.settings.solvent_type)
        builder.replaceFlag("LOG_PATH", logfile_name)
        builder.replaceFlag("ATOMS_TO_MINIMIZE",
                            "\"" + '\", \"'.join(atoms_to_minimize) + "\"")

        builder.write(output_path)

    def _calculateRMSD(self, pdb_name, trajectory_name):
        linkId = self._getPerturbingLinkId()

        initial = PDBParser(pdb_name).getLinkWithId(linkId)
        final = PDBParser(trajectory_name).getLinkWithId(linkId)

        return final.calculateRMSDWith(initial)

    def _getGeneralPath(self, lambda_, shifted_lambda, num):
        general_path = self.path
        if (lambda_.type != Lambda.DUAL_LAMBDA):
            general_path += str(num) + '_' + lambda_.type + "/"
        general_path += self._getLambdaFolderName(lambda_, shifted_lambda)
        general_path += "/"

        return general_path

    def _getLambdaFolderName(self, lambda_, shifted_lambda):
        return lambda_.folder_name + '_' + shifted_lambda.folder_name

    def _getPELEEnergyPrediction(self, runner, pid):
        try:
            output = runner.run(self.path +
                                co.SINGLE_POINT_CF_NAME.format(pid))
        except SystemExit as exception:
            print("DoubleWideSampling error: \n" + str(exception))
            sys.exit(1)

        for line in output.split('\n'):
            if line.startswith(pele_co.ENERGY_RESULT_LINE):
                return float(line.strip().split()[-1])

        print("Error: energy calculation failed")
        print(output)
        sys.exit(1)
