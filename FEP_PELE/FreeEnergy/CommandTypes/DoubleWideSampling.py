# -*- coding: utf-8 -*-


# Python imports
import sys
from multiprocessing import Pool, current_process
from functools import partial


# FEP_PELE imports
from FEP_PELE.FreeEnergy import Constants as co
from FEP_PELE.FreeEnergy.Command import Command
from FEP_PELE.FreeEnergy.Checkers import checkModelCoords

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
from FEP_PELE.Utils.InOut import remove_splitted_models

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
        print("######################")
        print(" Double Wide Sampling")
        print("######################")

        alchemicalTemplateCreator = AlchemicalTemplateCreator(
            self.settings.initial_template,
            self.settings.final_template,
            self.settings.atom_links)

        clear_directory(self.settings.calculation_path)

        for i, lambda_value in enumerate(self.settings.lambdas):
            print("\n# Lambda: {} #\n".format(lambda_value))

            clear_directory(self.settings.calculation_path + co.MODELS_FOLDER)

            print(" - Splitting PELE models")
            simulation = Simulation(self.settings.simulation_path +
                                    str(lambda_value), sim_type="PELE",
                                    report_name="report_",
                                    trajectory_name="trajectory_",
                                    logfile_name="logFile_")
            simulation.getOutputFiles()

            parallelLoop = partial(self._parallelTrajectoryWriterLoop,
                                   alchemicalTemplateCreator)

            with Pool(self.settings.number_of_processors) as pool:
                models_to_discard = pool.map(parallelLoop,
                                             simulation.iterateOverReports)

            # Inactivate bad models
            for model_info in models_to_discard:
                if (len(model_info) == 0):
                    continue
                print("  - Removing bad model {}".format(model_info[0][2]))
                for report in simulation.iterateOverReports:
                    if (report.name == model_info[0][0]):
                        report.models.inactivate(model_info[0][1])

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
                                       lambda_value,
                                       co.DIRECTION_TO_CHAR[direction])

                with Pool(self.settings.number_of_processors) as pool:
                    pool.map(parallelLoop, simulation.iterateOverReports)

                print("   Done")

        clear_directory(self.settings.calculation_path)

    def _parallelTrajectoryWriterLoop(self, alchemicalTemplateCreator,
                                      report_file):

        models_to_discard = []

        for model_id in range(0, report_file.trajectory.models.number):

            model_name = self.settings.calculation_path + co.MODELS_FOLDER + \
                str(model_id) + '-' + report_file.trajectory.name

            report_file.trajectory.writeModel(model_id, model_name)

            if (self.settings.safety_check):
                model_okay = checkModelCoords(
                    model_name, alchemicalTemplateCreator.getFragmentAtoms())

                if (not model_okay):
                    models_to_discard.append((report_file.name, model_id,
                                              model_name))

        return models_to_discard

    def _parallelPELEMinimizerLoop(self, lambda_value, direction_char,
                                   report_file):

        lambda_value = str(round(lambda_value, 3)) + direction_char

        dir_name = self.settings.calculation_path + lambda_value + "/"

        create_directory(dir_name)

        if (isThereAFile(dir_name + report_file.trajectory.name)):
            return

        pid = current_process().pid

        energies = []

        for model_id, active in enumerate(report_file.models):
            pdb_name = self.settings.calculation_path + co.MODELS_FOLDER + \
                str(model_id) + '-' + report_file.trajectory.name

            logfile_name = self.settings.calculation_path \
                + co.LOGFILE_NAME.format(pid)

            trajectory_name = self.settings.calculation_path + lambda_value + \
                "/" + str(model_id) + '-' + report_file.trajectory.name

            runner = PELERunner(self.settings.serial_pele,
                                number_of_processors=1)

            # In case bad model was previously removed
            if (not active):
                print("  - Warning: skipping bad model from path: " +
                      "\'{}\'".format(pdb_name))
                continue

            if (isThereAFile(trajectory_name)):
                write_recalculation_control_file(
                    self.settings.sp_control_file,
                    pdb_name,
                    logfile_name,
                    trajectory_name,
                    self.settings.calculation_path +
                    co.SINGLE_POINT_CF_NAME.format(pid))

                try:
                    output = runner.run(self.settings.calculation_path +
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
                    self.settings.calculation_path +
                    co.POST_PROCESSING_CF_NAME.format(pid))

                try:
                    output = runner.run(self.settings.calculation_path +
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

        path = self.settings.calculation_path + lambda_value + "/"

        # Write trajectories and reports
        write_energies_report(path, report_file, energies)
        join_splitted_models(path, report_file.trajectory.name)

        # Clean temporal files
        remove_splitted_models(path, report_file.trajectory.name)
