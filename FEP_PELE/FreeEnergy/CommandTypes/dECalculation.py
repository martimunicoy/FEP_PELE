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
from FEP_PELE.Utils.InOut import join_splitted_models
from FEP_PELE.Utils.InOut import remove_splitted_models
from FEP_PELE.Utils.InOut import writeLambdaTitle


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

        for lmb in self.lambdas:
            if (self.checkPoint.check((self.name, str(lmb.index) +
                                       str(lmb.type) +
                                       str(lmb.value)))):
                continue

            writeLambdaTitle(lmb)

            clear_directory(self.path + co.MODELS_FOLDER)

            print(" - Splitting PELE models")
            simulation = self._getSimulation(lmb, lmb.index)

            self._splitModels(simulation)

            ctt_lmb = None
            if (lmb.index == 2):
                if (lmb.type == Lambda.COULOMBIC_LAMBDA):
                    ctt_lmb = Lambda.Lambda(
                        1.0, lambda_type=Lambda.STERIC_LAMBDA)
                elif (lmb.type == Lambda.STERIC_LAMBDA):
                    ctt_lmb = Lambda.Lambda(
                        1.0, lambda_type=Lambda.COULOMBIC_LAMBDA)

            self._createAlchemicalTemplate(lmb, ctt_lmb)

            self._calculateOriginalEnergies(simulation, lmb, lmb.index)

            for shf_lmb in self.sampling_method.getShiftedLambdas(lmb):
                print(" - Applying delta lambda " +
                      str(round(shf_lmb.value - lmb.value, 5)))

                self._createAlchemicalTemplate(shf_lmb, ctt_lmb,
                                               gap=' ')

                general_path = self._getGeneralPath(lmb, lmb.index, shf_lmb)
                clear_directory(general_path)

                self._dECalculation(simulation, general_path, gap=' ')

            self.checkPoint.save((self.name, str(lmb.index) + str(lmb.type) +
                                  str(lmb.value)))

            clear_directory(self.path)

        self._finish()

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

    def _dECalculation(self, simulation, general_path, gap=''):
        print("{} - Calculating energetic differences".format(gap))

        parallelLoop = partial(self._parallelPELERecalculatorLoop,
                               general_path)

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
                model_name,
                self.path + co.SINGLE_POINT_CF_NAME.format(pid),
                logfile_name=logfile_name)

            # Run PELE and extract energy prediction
            energies.append(self._getPELEEnergyPrediction(runner, pid))

        write_energies_report(path, report_file, energies)

    def _parallelPELERecalculatorLoop(self, general_path, report_file):
        create_directory(general_path)

        pid = current_process().pid
        energies = []

        for model_id, active in enumerate(report_file.models):
            # Set initial variables
            model_name = self.path + co.MODELS_FOLDER + str(model_id) + \
                '-' + report_file.trajectory.name
            logfile_name = self.path + co.LOGFILE_NAME.format(pid)

            # Define new PELERunner
            runner = PELERunner(self.settings.serial_pele,
                                number_of_processors=1)

            # Write recalculation control file
            self._writeRecalculationControlFile(
                model_name,
                self.path + co.SINGLE_POINT_CF_NAME.format(pid),
                logfile_name=logfile_name)

            # Run PELE and extract energy prediction
            energies.append(self._getPELEEnergyPrediction(runner, pid))

        # Write trajectories and reports
        write_energies_report(general_path, report_file, energies)
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

    def _writeRecalculationControlFile(self, pdb_name,
                                       output_path,
                                       logfile_name=None):
            builder = ControlFileFromTemplateCreator(
                self.settings.sp_control_file,)

            builder.replaceFlag("INPUT_PDB_NAME", pdb_name)
            builder.replaceFlag("SOLVENT_TYPE", self.settings.solvent_type)
            if (logfile_name is not None):
                builder.replaceFlag("LOG_PATH", logfile_name)

            builder.write(output_path)

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
