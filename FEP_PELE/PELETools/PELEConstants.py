

# Constants regarding PELE standard output
ENERGY_RESULT_LINE = "ENERGY VACUUM + SGB + CONSTRAINTS + SELF + NON POLAR:"

# Constants regarding PELE report file
REPORT_TOTAL_ENERGY_COLUMN = 4

# Hardcoded PELE paths
HETEROATOMS_TEMPLATE_PATH = "DataLocal/Templates/OPLS2005/HeteroAtoms/"

# PELE executable types
PELE_SERIAL_EXEC_TYPE = "SERIAL"
PELE_MPI_EXEC_TYPE = "MPI"
PELE_EXECUTABLE_TYPES = (PELE_SERIAL_EXEC_TYPE, PELE_MPI_EXEC_TYPE)

# PELE Control File flag names
CONTROL_FILE_FLAG_NAMES = ["LICENSE_PATH",
                           "LOG_PATH",
                           "REPORT_PATH",
                           "TRAJECTORY_PATH",
                           "SOLVENT_TYPE",
                           "INPUT_PDB_NAME"]
