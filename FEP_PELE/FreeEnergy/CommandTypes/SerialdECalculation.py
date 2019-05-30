# -*- coding: utf-8 -*-


# Python imports
import sys


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
from FEP_PELE.Utils.InOut import remove_directory
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
class SerialdECalculation(Command):
    def __init__(self, settings):
        self._name = co.COMMAND_NAMES_DICT["SERIAL_DE_CALCULATION"]
        self._label = co.COMMAND_LABELS_DICT["SERIAL_DE_CALCULATION"]
        Command.__init__(self, settings)
        self._path = self.settings.calculation_path

    def run(self):
        self._start()

        create_directory(self.path)

        for lmb in self.lambdas:
            writeLambdaTitle(lmb)

            create_directory(self.path + str(self.PID) + '_' +
                             co.MODELS_FOLDER)

            print(" - Splitting PELE models")
            simulation = self._getSimulation(lmb)

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

            self._calculateOriginalEnergies(simulation, lmb)

            for shf_lmb in self.sampling_method.getShiftedLambdas(lmb):
                print(" - Applying delta lambda " +
                      str(round(shf_lmb.value - lmb.value, 5)))

                self._createAlchemicalTemplate(shf_lmb, ctt_lmb,
                                               gap=' ')

                general_path = self._getGeneralPath(lmb, shf_lmb)
                clear_directory(general_path)

                self._dECalculation(simulation, general_path, gap=' ')

            remove_directory(self.path + str(self.PID) + '_' +
                             co.MODELS_FOLDER)

        self._finish()

    def _getSimulation(self, lambda_):
        path = self.settings.simulation_path
        if (lambda_.type != Lambda.DUAL_LAMBDA):
            path += str(lambda_.index) + '_' + lambda_.type + "/"
        path += str(lambda_.value) + "/"

        simulation = Simulation(path, sim_type="PELE",
                                report_name="report_",
                                trajectory_name="trajectory_",
                                logfile_name="logFile_")
        simulation.getOutputFiles()

        return simulation

    def _splitModels(self, simulation):
        for report in simulation.iterateOverReports:
            self._trajectoryWriterLoop(report)

    def _calculateOriginalEnergies(self, simulation, lmb, gap=''):
        print("{} - Calculating original energies".format(gap))

        path = self.path + lmb.path
        clear_directory(path)

        for report in simulation.iterateOverReports:
            self._originalEnergiesCalculator(path, report)

    def _dECalculation(self, simulation, general_path, gap=''):
        print("{} - Calculating energetic differences".format(gap))

        for report in simulation.iterateOverReports:
            self._PELERecalculatorLoop(general_path, report)

    def _trajectoryWriterLoop(self, report_file):
        for model_id in range(0, report_file.trajectory.models.number):

            model_name = self.path + str(self.PID) + '_' + co.MODELS_FOLDER + \
                str(model_id) + '-' + report_file.trajectory.name

            report_file.trajectory.writeModel(model_id, model_name)

    def _originalEnergiesCalculator(self, path, report_file):
        # Define new PELERunner
        runner = PELERunner(self.settings.serial_pele,
                            number_of_processors=1)

        # Save original energies
        energies = []

        for model_id in range(0, report_file.trajectory.models.number):
            model_name = self.path + str(self.PID) + '_' + co.MODELS_FOLDER + \
                str(model_id) + '-' + report_file.trajectory.name

            logfile_name = path + co.LOGFILE_NAME.format(self.PID)

            # Write recalculation control file
            self._writeRecalculationControlFile(
                model_name,
                self.path + str(self.PID) + '_' + co.MODELS_FOLDER +
                co.SINGLE_POINT_CF_NAME.format(self.PID),
                logfile_name=logfile_name)

            # Run PELE and extract energy prediction
            energies.append(self._getPELEEnergyPrediction(runner))

        write_energies_report(path, report_file, energies)

    def _PELERecalculatorLoop(self, general_path, report_file):
        create_directory(general_path)

        energies = []
        path = self.path + str(self.PID) + '_' + co.MODELS_FOLDER

        for model_id, active in enumerate(report_file.models):
            # Set initial variables
            file_name = str(model_id) + '-' + report_file.trajectory.name
            model_name = path + file_name
            logfile_name = path + co.LOGFILE_NAME.format(self.PID)

            # Define new PELERunner
            runner = PELERunner(self.settings.serial_pele,
                                number_of_processors=1)

            # Write recalculation control file
            self._writeRecalculationControlFile(
                model_name,
                path + co.SINGLE_POINT_CF_NAME.format(self.PID),
                logfile_name=logfile_name)

            # Run PELE and extract energy prediction
            energies.append(self._getPELEEnergyPrediction(runner))

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

    def _getGeneralPath(self, lmb, shf_lmb=None):
        general_path = self.path
        if (lmb.type != Lambda.DUAL_LAMBDA):
            general_path += str(lmb.index) + '_' + lmb.type + "/"
        if (shf_lmb is not None):
            general_path += self._getLambdaFolderName(lmb, shf_lmb)
        else:
            general_path += lmb.folder_name
        general_path += "/"

        return general_path

    def _getPELEEnergyPrediction(self, runner):
        path = self.path + str(self.PID) + '_' + co.MODELS_FOLDER

        try:
            output = runner.run(path +
                                co.SINGLE_POINT_CF_NAME.format(self.PID))
        except SystemExit as exception:
            print("dECalculation error: \n" + str(exception))
            sys.exit(1)

        for line in output.split('\n'):
            if line.startswith(pele_co.ENERGY_RESULT_LINE):
                return float(line.strip().split()[-1])

        print("Error: energy calculation failed")
        print(output)
        sys.exit(1)
