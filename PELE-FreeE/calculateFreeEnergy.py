from PELETools.SimulationParser import Simulation
from PELETools.Utils import isThereAFile
from TemplateHandler.AlchemicalTemplateCreator import AlchemicalTemplateCreator

from Utils.InOut import create_directory
from Utils.InOut import clear_directory
from Utils.InOut import write_recalculation_control_file
from Utils.InOut import write_energies_report
from Utils.InOut import join_splitted_models

from subprocess import check_output
import multiprocessing
from functools import partial

TEST_GLOBAL_PATH = "/home/municoy/repos/PELE-FreeE/tests/LYS_BNZ-TOL/"
PELE_SERIAL_GLOBAL_PATH = "/home/municoy/builds/PELE/v1.5/PELE-1.5_serial"
ENERGY_RESULT_LINE = "ENERGY VACUUM + SGB + CONSTRAINTS + SELF + NON POLAR:"
TOTAL_ENERGY_COLUMN = 4
PROCESSORS_NUMBER = 36

LAMBDA_VALUES = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]


def parallelTrajectoryWriterLoop(report_file):
    for model_id in range(0, report_file.trajectory.models.number):
        report_file.trajectory.writeModel(model_id, TEST_GLOBAL_PATH +
                                          "calculation/models/" +
                                          str(model_id) + '-' +
                                          report_file.trajectory.name)


def parallelPELEMinimizerLoop(lambda_value, report_file):
    
    if (lambda_value > 0):
        lambda_value = str(round(lambda_value, 3)) + 'f'
    else:
        lambda_value = str(round(lambda_value, 3))[1:] + 'r'

    create_directory(TEST_GLOBAL_PATH + "calculation/" +
                     lambda_value + "/")

    if (isThereAFile(TEST_GLOBAL_PATH + "calculation/" +
                     lambda_value + "/" +
                     report_file.trajectory.name)):
        return

    pid = multiprocessing.current_process().pid

    energies = []

    for model_id in range(0, report_file.models.number):
        pdb_name = TEST_GLOBAL_PATH + "calculation/models/" + \
            str(model_id) + '-' + report_file.trajectory.name

        logfile_name = TEST_GLOBAL_PATH + "calculation/" + "logfile_" + \
            str(pid) + ".txt"

        trajectory_name = TEST_GLOBAL_PATH + "calculation/" + \
            lambda_value + "/" + str(model_id) + '-' + \
            report_file.trajectory.name

        if (isThereAFile(trajectory_name)):
            write_recalculation_control_file(TEST_GLOBAL_PATH +
                                             "pele_sp.conf",
                                             pdb_name,
                                             logfile_name,
                                             trajectory_name,
                                             TEST_GLOBAL_PATH +
                                             "calculation/" +
                                             "pele_sp_" +
                                             str(pid) + ".conf")

            output = check_output([PELE_SERIAL_GLOBAL_PATH,
                                   TEST_GLOBAL_PATH +
                                   "calculation/pele_sp_" +
                                   str(pid) + ".conf"])
        else:
            write_recalculation_control_file(TEST_GLOBAL_PATH +
                                             "pele_recal.conf",
                                             pdb_name,
                                             logfile_name,
                                             trajectory_name,
                                             TEST_GLOBAL_PATH +
                                             "calculation/" +
                                             "pele_recal_" +
                                             str(pid) + ".conf")

            output = check_output([PELE_SERIAL_GLOBAL_PATH,
                                   TEST_GLOBAL_PATH +
                                   "calculation/pele_recal_" +
                                   str(pid) + ".conf"])

        output = output.decode('utf-8').strip()

        for line in output.split('\n'):
            if line.startswith(ENERGY_RESULT_LINE):
                energies.append(float(line.strip().split()[-1]))
                break
        else:
            print("Error: energy calculation failed")
            print(output)

    path = TEST_GLOBAL_PATH + "calculation/" + lambda_value + "/"

    write_energies_report(path, report_file, energies)
    join_splitted_models(path, report_file.trajectory.name)


def main():
    alchemicalTemplateCreator = AlchemicalTemplateCreator(
        TEST_GLOBAL_PATH + "initialTemplates/bnzz",
        TEST_GLOBAL_PATH + "initialTemplates/tolz",
        [('_H1_', '_C7_')])

    clear_directory(TEST_GLOBAL_PATH + "calculation/")

    for i, lambda_value in enumerate(LAMBDA_VALUES):
        print("##############")
        print(" Lambda: " + str(lambda_value))
        print("##############")

        clear_directory(TEST_GLOBAL_PATH + "calculation/models/")

        print(" - Splitting PELE models")
        simulation = Simulation(TEST_GLOBAL_PATH + "simulation/" +
                                str(lambda_value), sim_type="PELE",
                                report_name="report_",
                                trajectory_name="trajectory_",
                                logfile_name="logFile_")
        simulation.getOutputFiles()

        with multiprocessing.Pool(PROCESSORS_NUMBER) as pool:
            pool.map(parallelTrajectoryWriterLoop,
                     simulation.iterateOverReports)
        print("   Done")

        # ---------------------------------------------------------------------

        delta_lambda = get_delta_lambda(i, LAMBDA_VALUES)

        operations = [-1, 1]

        for op in operations:
            print(" - Creating alchemical template")
            print("  - Applying delta lambda " + str(op * delta_lambda))

            alchemicalTemplateCreator.create(
                lambda_value + op * delta_lambda,
                TEST_GLOBAL_PATH +
                "DataLocal/Templates/OPLS2005/HeteroAtoms/tolz")
            print("   Done")

            # -----------------------------------------------------------------

            print(" - Minimizing and calculating energetic differences")

            parallelLool = partial(parallelPELEMinimizerLoop,
                                   op * lambda_value)

            with multiprocessing.Pool(PROCESSORS_NUMBER) as pool:
                pool.map(parallelLool, simulation.iterateOverReports)

            print("   Done")

            """

            with open(TEST_GLOBAL_PATH +
                      "calculation/energies.out", 'a') as file:
                file.write("# From {:.3f}".format(lambda_value) +
                           " to {:.3f}".format(round(lambda_value + op *
                                                     delta_lambda, 3)) +
                           "\n")
                file.write(str(energies).strip('[]'))
                file.write('\n')
            """


def get_delta_lambda(position, lambda_values):
    if (position == 0):
        return round(abs(lambda_values[position]), 3)
    else:
        return round(abs(lambda_values[position - 1] -
                         lambda_values[position]), 3)


if __name__ == '__main__':
    main()
