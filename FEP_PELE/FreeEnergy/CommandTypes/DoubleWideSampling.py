# -*- coding: utf-8 -*-


# Python imports
import os
import sys
from multiprocessing import Pool, current_process
from functools import partial


# FEP_PELE imports
from FEP_PELE.FreeEnergy import Constants as co
from FEP_PELE.FreeEnergy.Command import Command

from FEP_PELE.TemplateHandler.AlchemicalTemplateCreator import \
    AlchemicalTemplateCreator

from FEP_PELE.PELETools import PELEConstants as pele_co
from FEP_PELE.PELETools.SimulationParser import Simulation
from FEP_PELE.PELETools.PELERunner import PELERunner

from FEP_PELE.Utils.InOut import create_directory
from FEP_PELE.Utils.InOut import clear_directory
from FEP_PELE.Utils.InOut import write_recalculation_control_file
from FEP_PELE.Utils.InOut import write_energies_report
from FEP_PELE.Utils.InOut import join_splitted_models
from FEP_PELE.Utils.InOut import isThereAFile


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Class definitions
class DoubleWideSampling(Command):
    def __init__(self, settings):
        Command.__init__(self, settings)
        self.directions = co.DOUBLE_WIDE_SAMPLING_DIRECTIONS

    def run(self):
        alchemicalTemplateCreator = AlchemicalTemplateCreator(
            self.settings.initial_template,
            self.settings.final_template,
            self.settings.atom_links)

        clear_directory(self.settings.general_path + co.CALCULATION_FOLDER)

        for i, lambda_value in enumerate(self.settings.lambdas):
            print("##############")
            print(" Lambda: " + str(lambda_value))
            print("##############")

            clear_directory(self.settings.general_path +
                            co.CALCULATION_FOLDER + co.MODELS_FOLDER)

            print(" - Splitting PELE models")
            simulation = Simulation(self.settings.general_path +
                                    co.SIMULATION_FOLDER +
                                    str(lambda_value), sim_type="PELE",
                                    report_name="report_",
                                    trajectory_name="trajectory_",
                                    logfile_name="logFile_")
            simulation.getOutputFiles()

            with Pool(self.settings.number_of_processors) as pool:
                pool.map(self._parallelTrajectoryWriterLoop,
                         simulation.iterateOverReports)
            print("   Done")

            # ---------------------------------------------------------------------

            delta_lambda = self._get_delta_lambda(i)

            for direction in self.directions:
                dir_factor = co.DIRECTION_FACTORS[direction]
                print(" - Creating alchemical template")
                print("  - Applying delta lambda " +
                      str(dir_factor * delta_lambda))

                alchemicalTemplateCreator.create(
                    lambda_value + dir_factor * delta_lambda,
                    self.settings.general_path +
                    pele_co.HETEROATOMS_TEMPLATE_PATH +
                    self.settings.final_template_name)
                print("   Done")

                # -----------------------------------------------------------------

                print(" - Minimizing and calculating energetic differences")

                parallelLoop = partial(self._parallelPELEMinimizerLoop,
                                       dir_factor * lambda_value,
                                       co.DIRECTION_TO_CHAR[direction])

                with Pool(self.settings.number_of_processors) as pool:
                    pool.map(parallelLoop, simulation.iterateOverReports)

                print("   Done")

    def _parallelTrajectoryWriterLoop(self, report_file):
        for model_id in range(0, report_file.trajectory.models.number):
            report_file.trajectory.writeModel(model_id,
                                              self.settings.general_path +
                                              co.CALCULATION_FOLDER +
                                              co.MODELS_FOLDER +
                                              str(model_id) + '-' +
                                              report_file.trajectory.name)

            model_okay = checkModelCoords(self.settings.general_path +
                                          co.CALCULATION_FOLDER +
                                          co.MODELS_FOLDER +
                                          str(model_id) + '-' +
                                          report_file.trajectory.name,
                                          self.settings.atom_links)

            if (not model_okay):
                os.remove(self.settings.general_path + co.CALCULATION_FOLDER +
                          co.MODELS_FOLDER + str(model_id) + '-' +
                          report_file.trajectory.name)

    def _parallelPELEMinimizerLoop(self, lambda_value, direction_char,
                                   report_file):

        lambda_value = str(round(lambda_value, 3)) + direction_char

        create_directory(self.settings.general_path + co.CALCULATION_FOLDER +
                         lambda_value + "/")

        if (isThereAFile(self.settings.general_path + co.CALCULATION_FOLDER +
                         lambda_value + "/" + report_file.trajectory.name)):
            return

        pid = current_process().pid

        energies = []

        for model_id in range(0, report_file.models.number):
            pdb_name = self.settings.general_path + co.CALCULATION_FOLDER + \
                co.MODELS_FOLDER + str(model_id) + '-' + \
                report_file.trajectory.name

            logfile_name = self.settings.general_path + \
                co.CALCULATION_FOLDER + co.LOGFILE_NAME.format(pid)

            trajectory_name = self.settings.general_path + \
                co.CALCULATION_FOLDER + lambda_value + "/" + \
                str(model_id) + '-' + \
                report_file.trajectory.name

            runner = PELERunner(self.settings.serial_pele,
                                number_of_processors=1)

            # TODO!!!
            if (not isThereAFile(self.settings.general_path +
                                 co.CALCULATION_FOLDER +
                                 co.SINGLE_POINT_CF_NAME.format(pid))):
                continue

            if (isThereAFile(trajectory_name)):
                write_recalculation_control_file(
                    self.settings.sp_control_file,
                    pdb_name,
                    logfile_name,
                    trajectory_name,
                    self.settings.general_path + co.CALCULATION_FOLDER +
                    co.SINGLE_POINT_CF_NAME.format(pid))

                try:
                    output = runner.run(self.settings.general_path +
                                        co.CALCULATION_FOLDER +
                                        co.SINGLE_POINT_CF_NAME.format(pid))
                except SystemExit as exception:
                    print("DoubleWideSampling error: \n" + str(exception))
                    sys.exit(1)

            else:
                write_recalculation_control_file(
                    self.settings.pp_control_file,
                    pdb_name,
                    logfile_name,
                    trajectory_name,
                    self.settings.general_path + co.CALCULATION_FOLDER +
                    co.POST_PROCESSING_CF_NAME.format(pid))

                try:
                    output = runner.run(self.settings.general_path +
                                        co.CALCULATION_FOLDER +
                                        co.POST_PROCESSING_CF_NAME.format(pid))
                except SystemExit as exception:
                    print("DoubleWideSampling error: \n" + str(exception))
                    sys.exit(1)

            for line in output.split('\n'):
                if line.startswith(pele_co.ENERGY_RESULT_LINE):
                    energies.append(float(line.strip().split()[-1]))
                    break
            else:
                print("Error: energy calculation failed")
                print(output)

        path = self.settings.general_path + co.CALCULATION_FOLDER + \
            lambda_value + "/"

        write_energies_report(path, report_file, energies)
        join_splitted_models(path, report_file.trajectory.name)


def checkModelCoords(path, atom_names):
    coords = set()
    # @TODO
    atom_names = ("_C7_", "_H7_", "_H8_", "_H9_")
    with open(path, 'r') as file:
        for line in file:
            atom_name = line[12:16].replace(' ', '_')
            if (atom_name in atom_names):
                coords.add((line[31:38],
                            line[39:46],
                            line[47:54]))

    if (len(coords) != len(atom_names)):
        return False

    return True
