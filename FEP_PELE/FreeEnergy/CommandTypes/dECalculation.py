# -*- coding: utf-8 -*-


# Python imports
import sys
from multiprocessing import Pool, current_process
from functools import partial


# FEP_PELE imports
from FEP_PELE.FreeEnergy import Constants as co
from FEP_PELE.FreeEnergy.Command import Command

from FEP_PELE.TemplateHandler import Lambda

from FEP_PELE.PELETools import PELEConstants as pele_co
from FEP_PELE.PELETools.SimulationParser import Simulation
from FEP_PELE.PELETools.PELERunner import PELERunner
from FEP_PELE.PELETools.ControlFileCreator import \
    ControlFileFromTemplateCreator

from FEP_PELE.Utils.InOut import create_directory
from FEP_PELE.Utils.InOut import clear_directory
from FEP_PELE.Utils.InOut import write_energies_report
from FEP_PELE.Utils.InOut import write_recalculated_energies_report
from FEP_PELE.Utils.InOut import join_splitted_models
from FEP_PELE.Utils.InOut import remove_splitted_models
from FEP_PELE.Utils.InOut import writeLambdaTitle
from FEP_PELE.Utils.InOut import copyFile

from FEP_PELE.Tools.PDBTools import PDBParser

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
        atoms_to_minimize = self._getAtomIdsToMinimize()

        for lambda_ in lambdas:
            if (self.checkPoint.check((self.name, str(num) +
                                       str(lambda_.type) +
                                       str(lambda_.value)))):
                continue

            writeLambdaTitle(lambda_)

            clear_directory(self.path + co.MODELS_FOLDER)

            print(" - Splitting PELE models")
            simulation = self._getSimulation(lambda_, num)

            self._splitModels(simulation)

            self._createAlchemicalTemplate(lambda_, constant_lambda)

            self._calculateOriginalEnergies(simulation, lambda_, num)

            clear_directory(self.path)

            for shif_lambda in self.sampling_method.getShiftedLambdas(lambda_):
                print(" - Applying delta lambda " +
                      str(round(shif_lambda.value - lambda_.value, 5)))

                self._createAlchemicalTemplate(shif_lambda, constant_lambda,
                                               gap=' ')

                general_path = self._getGeneralPath(lambda_, num, shif_lambda)
                clear_directory(general_path)

                self._minimize(simulation, lambdas_type, general_path,
                               atoms_to_minimize, gap=' ')

                self._dECalculation(simulation, lambda_, shif_lambda,
                                    general_path, num, gap=' ')

            self.checkPoint.save((self.name, str(num) + str(lambda_.type) +
                                  str(lambda_.value)))

            clear_directory(self.path)

        return []

    def _getSimulation(self, lambda_, num):
        path = self.settings.simulation_path
        if (lambda_.type != Lambda.DUAL_LAMBDA):
            path += str(num) + '_' + lambda_.type + "/"
        path += str(lambda_.value) + "/"

        simulation = Simulation(path, sim_type="PELE",
                                report_name="report_",
                                trajectory_name="trajectory_",
                                logfile_name="logFile_")
        simulation.getOutputFiles()

        return simulation

    def _splitModels(self, simulation):
        with Pool(self.settings.number_of_processors) as pool:
            pool.map(self._parallelTrajectoryWriterLoop,
                     simulation.iterateOverReports)

    def _calculateOriginalEnergies(self, simulation, lambda_, num, gap=''):
        print("{} - Calculating original energies".format(gap))

        path = self.path
        if (lambda_.type != Lambda.DUAL_LAMBDA):
            path += str(num) + '_' + lambda_.type + "/"
        path += lambda_.folder_name + '/'

        clear_directory(path)

        originalEnergiesCalculator = partial(
            self._parallelOriginalEnergiesCalculator, path)

        with Pool(self.settings.number_of_processors) as pool:
            pool.map(originalEnergiesCalculator,
                     simulation.iterateOverReports)

    def _minimize(self, simulation, lambdas_type, general_path,
                  atoms_to_minimize, gap=''):
        if ((self.settings.reminimize) and
            ((lambdas_type == Lambda.DUAL_LAMBDA) or
             (lambdas_type == Lambda.STERIC_LAMBDA))):
            print("{} - Minimizing distances".format(gap))

            parallelLoop = partial(self._parallelPELEMinimizerLoop,
                                   general_path, atoms_to_minimize)

            with Pool(self.settings.number_of_processors) as pool:
                pool.map(parallelLoop, simulation.iterateOverReports)

        else:
            for report in simulation.iterateOverReports:
                for model_id in range(0, report.trajectory.models.number):
                    file_name = str(model_id) + '-' + report.trajectory.name
                    original_pdb = self.path + co.MODELS_FOLDER + file_name
                    copyFile(original_pdb, general_path)

    def _dECalculation(self, simulation, lambda_, shif_lambda, general_path,
                       num, gap=''):
        print("{} - Calculating energetic differences".format(gap))

        parallelLoop = partial(self._parallelPELERecalculatorLoop,
                               lambda_, shif_lambda, general_path, num)

        with Pool(self.settings.number_of_processors) as pool:
            pool.map(parallelLoop, simulation.iterateOverReports)

    def _parallelTrajectoryWriterLoop(self, report_file):
        for model_id in range(0, report_file.trajectory.models.number):

            model_name = self.path + co.MODELS_FOLDER + str(model_id) + \
                '-' + report_file.trajectory.name

            report_file.trajectory.writeModel(model_id, model_name)

    def _parallelOriginalEnergiesCalculator(self, path, report_file):
        pid = current_process().pid

        # Define new PELERunner
        runner = PELERunner(self.settings.serial_pele,
                            number_of_processors=1)

        # Save original energies
        energies = []

        for model_id in range(0, report_file.trajectory.models.number):
            model_name = self.path + co.MODELS_FOLDER + str(model_id) + \
                '-' + report_file.trajectory.name

            logfile_name = path + co.LOGFILE_NAME.format(pid)

            # Write recalculation control file
            self._writeRecalculationControlFile(
                self.settings.sp_control_file,
                model_name,
                self.path + co.SINGLE_POINT_CF_NAME.format(pid),
                logfile_name=logfile_name)

            # Run PELE and extract energy prediction
            energies.append(self._getPELEEnergyPrediction(runner, pid))

        write_recalculated_energies_report(path + report_file.name, energies)

    def _parallelPELEMinimizerLoop(self, general_path, atoms_to_minimize,
                                   report_file):
        create_directory(general_path)

        pid = current_process().pid

        for model_id, active in enumerate(report_file.models):
            # Set initial variables
            file_name = str(model_id) + '-' + report_file.trajectory.name
            original_pdb = self.path + co.MODELS_FOLDER + file_name
            logfile_name = self.path + co.LOGFILE_NAME.format(pid)
            minimized_pdb = general_path + file_name

            # Define new PELERunner
            runner = PELERunner(self.settings.serial_pele,
                                number_of_processors=1)

            # Write recalculation control file
            self._writeRecalculationControlFile(
                self.settings.pp_control_file,
                original_pdb,
                self.path + co.POST_PROCESSING_CF_NAME.format(pid),
                logfile_name=logfile_name,
                trajectory_name=minimized_pdb,
                atoms_to_minimize=atoms_to_minimize)

            runner.run(self.path + co.POST_PROCESSING_CF_NAME.format(pid))

            self._applyMinimizedDistancesTo(original_pdb, minimized_pdb)

    def _parallelPELERecalculatorLoop(self, lambda_, shifted_lambda,
                                      general_path, num, report_file):
        create_directory(general_path)

        pid = current_process().pid

        original_energies = self._getOriginalEnergies(
            self._getGeneralPath(lambda_, num) + report_file.name)
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
                self.path + co.SINGLE_POINT_CF_NAME.format(pid),
                logfile_name=logfile_name)

            # Run PELE and extract energy prediction
            energies.append(self._getPELEEnergyPrediction(runner, pid))

            # Calculate RMSD between original pdb and shifted one
            rmsd = self._calculateRMSD(original_pdb, shifted_pdb)
            rmsds.append(rmsd)

        # Write trajectories and reports
        write_energies_report(general_path, report_file, energies,
                              original_energies, rmsds)
        join_splitted_models(general_path, "*-" + report_file.trajectory.name)

        # Clean temporal files
        remove_splitted_models(general_path,
                               "*-" + report_file.trajectory.name)

    def _getOriginalEnergies(self, path):
        energies = []

        with open(path, 'r') as file:
            file.readline()
            for line in file:
                energies.append(float(line.strip()))

        return energies

    def _writeRecalculationControlFile(self, template_path, pdb_name,
                                       output_path,
                                       logfile_name=None,
                                       trajectory_name=None,
                                       atoms_to_minimize=None):
            builder = ControlFileFromTemplateCreator(template_path)

            builder.replaceFlag("INPUT_PDB_NAME", pdb_name)
            builder.replaceFlag("SOLVENT_TYPE", self.settings.solvent_type)
            if (logfile_name is not None):
                builder.replaceFlag("LOG_PATH", logfile_name)
            if (trajectory_name is not None):
                builder.replaceFlag("TRAJECTORY_PATH", trajectory_name)
            if (atoms_to_minimize is not None):
                builder.replaceFlag("ATOMS_TO_MINIMIZE",
                                    "\"" + '\", \"'.join(atoms_to_minimize) +
                                    "\"")

            builder.write(output_path)

    def _calculateRMSD(self, pdb_name, trajectory_name):
        linkId = self._getPerturbingLinkId()

        initial = PDBParser(pdb_name).getLinkWithId(linkId)
        final = PDBParser(trajectory_name).getLinkWithId(linkId)

        return final.calculateRMSDWith(initial)

    def _getGeneralPath(self, lambda_, num, shifted_lambda=None):
        general_path = self.path
        if (lambda_.type != Lambda.DUAL_LAMBDA):
            general_path += str(num) + '_' + lambda_.type + "/"
        if (shifted_lambda is not None):
            general_path += self._getLambdaFolderName(lambda_, shifted_lambda)
        else:
            general_path += lambda_.folder_name
        general_path += "/"

        return general_path

    def _getPELEEnergyPrediction(self, runner, pid):
        try:
            output = runner.run(self.path +
                                co.SINGLE_POINT_CF_NAME.format(pid))
        except SystemExit as exception:
            print("dECalculation error: \n" + str(exception))
            sys.exit(1)

        for line in output.split('\n'):
            if line.startswith(pele_co.ENERGY_RESULT_LINE):
                return float(line.strip().split()[-1])

        print("Error: energy calculation failed")
        print(output)
        sys.exit(1)
