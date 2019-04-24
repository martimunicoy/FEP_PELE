# -*- coding: utf-8 -*-


# Python imports
import sys
from multiprocessing import Pool, current_process
from functools import partial


# FEP_PELE imports
from FEP_PELE.FreeEnergy import Constants as co
from FEP_PELE.FreeEnergy.Command import Command
from FEP_PELE.FreeEnergy.Checkers import checkModelCoords
from FEP_PELE.FreeEnergy.SamplingMethods.SamplingMethodBuilder import \
    SamplingMethodBuilder

from FEP_PELE.TemplateHandler import Lambda
from FEP_PELE.TemplateHandler.AlchemicalTemplateCreator import \
    AlchemicalTemplateCreator

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
from FEP_PELE.Utils.InOut import isThereAFile


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

        builder = SamplingMethodBuilder(self.settings)
        self._s_method = builder.createSamplingMethod()

    @property
    def name(self):
        return self._name

    @property
    def path(self):
        return self._path

    @property
    def s_method(self):
        return self._s_method

    def run(self):
        self._start()

        alchemicalTemplateCreator = AlchemicalTemplateCreator(
            self.settings.initial_template,
            self.settings.final_template,
            self.settings.atom_links)

        clear_directory(self.path)

        if (self.settings.splitted_lambdas):
            self._run_with_splitted_lambdas(alchemicalTemplateCreator)
        else:
            self._run(alchemicalTemplateCreator, self.settings.lambdas,
                      Lambda.DUAL_LAMBDA)

        self._finish()

    def _run(self, alchemicalTemplateCreator, lambdas, lambdas_type, num=0,
             constant_lambda=None):

        for lambda_ in Lambda.IterateOverLambdas(lambdas, lambdas_type):
            if (self.checkPoint.check((self.name, str(num) +
                                       str(lambda_.type) +
                                       str(lambda_.value)))):
                continue

            writeLambdaTitle(lambda_)

            clear_directory(self.path + co.MODELS_FOLDER)

            path = self.settings.simulation_path
            if (lambda_.type != Lambda.DUAL_LAMBDA):
                path += str(num) + '_' + lambda_.type + "/"
            path += str(lambda_.value) + "/"

            print(" - Splitting PELE models")
            simulation = Simulation(path, sim_type="PELE",
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
            for model_info_chunks in models_to_discard:
                if (len(model_info_chunks) == 0):
                    continue
                for model_info in model_info_chunks:
                    print("  - Removing bad model {}".format(model_info[2]))
                    for report in simulation.iterateOverReports:
                        if (report.name == model_info[0]):
                            report.models.inactivate(model_info[1])

            print("   Done")

            # ---------------------------------------------------------------------

            for shifted_lambda in self.s_method.getShiftedLambdas(lambda_):
                print(" - Creating alchemical template")
                print("  - Applying delta lambda " +
                      str(round(shifted_lambda.value - lambda_.value, 3)))

                self._createAlchemicalTemplate(alchemicalTemplateCreator,
                                               shifted_lambda, constant_lambda)

                print("   Done")

                # -----------------------------------------------------------------

                print(" - Minimizing and calculating energetic differences")

                atoms_to_minimize = self._getAtomIdsToMinimize(
                    alchemicalTemplateCreator)

                parallelLoop = partial(self._parallelPELEMinimizerLoop,
                                       lambda_, shifted_lambda,
                                       atoms_to_minimize, num)

                with Pool(self.settings.number_of_processors) as pool:
                    pool.map(parallelLoop, simulation.iterateOverReports)

                print("   Done")

            self.checkPoint.save((self.name, str(num) + str(lambda_.type) +
                                  str(lambda_.value)))

            clear_directory(self.path)

        return []

    def _parallelTrajectoryWriterLoop(self, alchemicalTemplateCreator,
                                      report_file):

        models_to_discard = []

        for model_id in range(0, report_file.trajectory.models.number):

            model_name = self.path + co.MODELS_FOLDER + str(model_id) + \
                '-' + report_file.trajectory.name

            report_file.trajectory.writeModel(model_id, model_name)

            if (self.settings.safety_check):
                model_okay = checkModelCoords(
                    model_name, alchemicalTemplateCreator.getFragmentAtoms())

                if (not model_okay):
                    models_to_discard.append((report_file.name, model_id,
                                              model_name))

        return models_to_discard

    def _parallelPELEMinimizerLoop(self, lambda_, shifted_lambda,
                                   atoms_to_minimize, num, report_file):

        lambda_folder = str(round(lambda_.value, 3)) + '_' + \
            str(round(shifted_lambda.value, 3))

        dir_name = self.path
        if (lambda_.type != Lambda.DUAL_LAMBDA):
            dir_name += str(num) + '_' + lambda_.type + "/"
        dir_name += lambda_folder + "/"

        create_directory(dir_name)

        if (isThereAFile(dir_name + report_file.trajectory.name)):
            return

        pid = current_process().pid

        energies = []

        for model_id, active in enumerate(report_file.models):
            pdb_name = self.path + co.MODELS_FOLDER + str(model_id) + \
                '-' + report_file.trajectory.name

            logfile_name = self.path + co.LOGFILE_NAME.format(pid)

            trajectory_name = self.path
            if (lambda_.type != Lambda.DUAL_LAMBDA):
                trajectory_name += str(num) + '_' + lambda_.type + "/"
            trajectory_name += lambda_folder + "/" + \
                str(model_id) + '-' + report_file.trajectory.name

            runner = PELERunner(self.settings.serial_pele,
                                number_of_processors=1)

            # In case bad model was previously removed
            if (not active):
                print("  - Warning: skipping bad model from path: " +
                      "\'{}\'".format(pdb_name))
                continue

            if (isThereAFile(trajectory_name)):
                self._writeRecalculationControlFile(
                    self.settings.sp_control_file,
                    pdb_name,
                    logfile_name,
                    trajectory_name,
                    atoms_to_minimize,
                    self.path + co.SINGLE_POINT_CF_NAME.format(pid))

                try:
                    output = runner.run(self.path +
                                        co.SINGLE_POINT_CF_NAME.format(pid))
                except SystemExit as exception:
                    print("DoubleWideSampling error: \n" + str(exception))
                    sys.exit(1)

            else:
                self._writeRecalculationControlFile(
                    self.settings.pp_control_file,
                    pdb_name,
                    logfile_name,
                    trajectory_name,
                    atoms_to_minimize,
                    self.path + co.POST_PROCESSING_CF_NAME.format(pid))

                try:
                    output = runner.run(self.path +
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

        path = self.path
        if (lambda_.type != Lambda.DUAL_LAMBDA):
            path += str(num) + '_' + lambda_.type + "/"
        path += lambda_folder + "/"

        # Write trajectories and reports
        write_energies_report(path, report_file, energies)
        join_splitted_models(path, "*-" + report_file.trajectory.name)

        # Clean temporal files
        remove_splitted_models(path, "*-" + report_file.trajectory.name)

    def _writeRecalculationControlFile(self, template_path, pdb_name,
                                       logfile_name,
                                       trajectory_name, atoms_to_minimize,
                                       output_path):
        builder = ControlFileFromTemplateCreator(template_path)

        builder.replaceFlag("INPUT_PDB_NAME", pdb_name)
        builder.replaceFlag("SOLVENT_TYPE", self.settings.solvent_type)
        builder.replaceFlag("LOG_PATH", logfile_name)
        builder.replaceFlag("TRAJECTORY_PATH", trajectory_name)
        builder.replaceFlag("ATOMS_TO_MINIMIZE",
                            "\"" + '\", \"'.join(atoms_to_minimize) + "\"")

        builder.write(output_path)
