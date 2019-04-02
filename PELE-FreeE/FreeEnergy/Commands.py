
class Command(object):
    def __init__(self, settings):
        self.settings = settings


class LambdasSampling(Command):
    def __init__(self, settings):
        Command.__init__(self, settings)

    def run():
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

