from __future__ import absolute_import
import subprocess

from TemplateHandler.AlchemicalTemplateCreator import AlchemicalTemplateCreator
from Utils.InOut import clear_directory
from Utils.InOut import write_lambda_value_to_control_file

TEST_GLOBAL_PATH = "/home/municoy/repos/PELE-FreeE/tests/LYS_BNZ-TOL/"
PELE_SERIAL_GLOBAL_PATH = "/home/municoy/builds/PELE/v1.5/PELE-1.5_serial"
PELE_MPI_GLOBAL_PATH = "/home/municoy/builds/PELE/v1.5/PELE-1.5_mpi"
LAMBDA_VALUES = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]


def main():

    alchemicalTemplateCreator = AlchemicalTemplateCreator(
        TEST_GLOBAL_PATH + "initialTemplates/bnzz",
        TEST_GLOBAL_PATH + "initialTemplates/tolz",
        [('_H1_', '_C7_')])

    clear_directory(TEST_GLOBAL_PATH + "simulation/")
    for lambda_value in LAMBDA_VALUES:
        print("##############")
        print(" Lambda: " + str(lambda_value))
        print("##############")
        print(" - Creating alchemical template")
        alchemicalTemplateCreator.create(
            lambda_value,
            TEST_GLOBAL_PATH + "DataLocal/Templates/OPLS2005/HeteroAtoms/tolz")
        print("   Done")

        print(" - Running PELE")

        path = TEST_GLOBAL_PATH + "minimization" + "/"
        clear_directory(path)

        print("  - Initial minimization")
        subprocess.run([PELE_SERIAL_GLOBAL_PATH,
                        TEST_GLOBAL_PATH + "pele_min.conf"])

        print("  - Simulation")

        path = TEST_GLOBAL_PATH + "simulation/" + str(lambda_value) + "/"
        clear_directory(path)
        write_lambda_value_to_control_file(TEST_GLOBAL_PATH + "pele_sim.conf",
                                           lambda_value,
                                           path + "pele_sim.conf")
        subprocess.run(["mpirun", "-n", "36", "--oversubscribe",
                        PELE_MPI_GLOBAL_PATH, path + "/pele_sim.conf"])


if __name__ == '__main__':
    main()
