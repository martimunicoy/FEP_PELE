import os
import glob

from .StringTools import natural_sort

LAMBDA_FLAG = "$LAMBDA$"
PDB_NAME_FLAG = "$INPUT_PDB_NAME$"
LOGFILE_NAME_FLAG = "$LOGFILE_NAME$"
TRAJECTORY_NAME_FLAG = "$OUTPUT_TRAJECTORY_NAME$"
REPORT_FIRST_LINE = "#Task    Step    " + \
    "numberOfAcceptedPeleSteps    currentEnergy    " + \
    "energyDifference\n"


# Functions
def isThereAFile(file_path):
    if os.path.exists(file_path):
        return True
    else:
        return False


def isAPath(path):
    path = str(path)
    if (not os.path.isdir(path)):
        raise NameError("Invalid path")


def clear_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        files = glob.glob(path + "*")
        for f in files:
            try:
                os.remove(f)
            except OSError:
                pass


def clear_file(path):
    if os.path.exists(path):
        os.remove(path)


def create_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)


def write_lambda_value_to_control_file(input_path, lambda_value,
                                       output_path=None):
    if (output_path is None):
        output_path = input_path

    # Read in the file
    with open(input_path, 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace(LAMBDA_FLAG, str(lambda_value))

    # Write the file out again
    with open(output_path, 'w') as file:
        file.write(filedata)


def get_all_files_from_with_extension(path, extension):
    files = glob.glob(path + "*." + extension)
    return files


def write_recalculation_control_file(input_path, pdb_name, logfile_name,
                                     trajectory_name, output_path=None):
    if (output_path is None):
        output_path = input_path

    # Read in the file
    with open(input_path, 'r') as file:
        filedata = file.read()

    # Replace the target string
    filedata = filedata.replace(PDB_NAME_FLAG, str(pdb_name))
    filedata = filedata.replace(LOGFILE_NAME_FLAG, str(logfile_name))
    filedata = filedata.replace(TRAJECTORY_NAME_FLAG, str(trajectory_name))

    # Write the file out again
    with open(output_path, 'w') as file:
        file.write(filedata)


def write_energies_report(output_path, report_file, energies):
    tasks = report_file.getMetric(1)
    steps = report_file.getMetric(2)
    accepted_steps = report_file.getMetric(3)
    original_energies = report_file.getMetric(4)

    with open(output_path + report_file.name, 'w') as file:
        file.write(REPORT_FIRST_LINE)
        for i, energy in enumerate(energies):
            file.write(str(tasks[i]) + "    " +
                       str(steps[i]) + "    " +
                       str(accepted_steps[i]) + "    " +
                       str(energy) + "    " +
                       str(energy - original_energies[i]) +
                       "\n")


def join_splitted_models(path, trajectory_name):
    with open(path + trajectory_name, 'w') as trajectory_file:
        models = glob.glob(path + '*-' + trajectory_name)
        models = natural_sort(models)
        for i, model in enumerate(models):
            trajectory_file.write("MODEL " + str(i + 1) + '\n')
            with open(model) as model_file:
                trajectory_file.writelines(model_file.readlines()[:-1])
            trajectory_file.write("ENDMDL" + '\n')
