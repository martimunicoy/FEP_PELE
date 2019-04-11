

# Default settings
DEF_SERIAL_PELE = ""
DEF_MPI_PELE = ""
DEF_INITIAL_TEMPLATE = ""
DEF_FINAL_TEMPLATE = ""
DEF_ATOM_LINKS = []
DEF_LAMBDAS = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95}
DEF_NUMBER_OF_PROCESSORDS = 4
DEF_COMMANDS = []
DEF_MIN_CONTROL_FILE = "min_pele.conf"
DEF_SIM_CONTROL_FILE = "sim_pele.conf"
DEF_PP_CONTROL_FILE = "pele_recal.conf"
DEF_SP_CONTROL_FILE = "pele_sp.conf"
DEF_SAFETY_CHECK = False

# Input file keys
INPUT_FILE_KEYS = [
    # General path
    "GeneralPath",
    # PELE paths
    "SerialPelePath",
    "MPIPelePath",
    # Template paths
    "InitialTemplate",
    "FinalTemplate",
    # Atom links
    "AtomLink",
    # Free Energy Settings
    "Lambdas",
    "NumberOfProcessors",
    "SafetyCheck",
    # List of commands
    "Commands",
    # PELE control files
    "MinimizationControlFile",
    "SimulationControlFile",
    "PostProcessingControlFile",
    "SinglePointControlFile"]


# Input file dict
CONTROL_FILE_DICT = {
    # General path
    "GENERAL_PATH": INPUT_FILE_KEYS[0],
    # PELE paths
    "SERIAL_PELE_PATH": INPUT_FILE_KEYS[1],
    "MPI_PELE_PATH": INPUT_FILE_KEYS[2],
    # Template paths
    "INITIAL_TEMPLATE": INPUT_FILE_KEYS[3],
    "FINAL_TEMPLATE": INPUT_FILE_KEYS[4],
    # Atom links
    "ATOM_LINK": INPUT_FILE_KEYS[5],
    # Free Energy Settings
    "LAMBDAS": INPUT_FILE_KEYS[6],
    "NUMBER_OF_PROCESSORS": INPUT_FILE_KEYS[7],
    "SAFETY_CHECK": INPUT_FILE_KEYS[8],
    # List of commands
    "COMMANDS": INPUT_FILE_KEYS[9],
    # PELE control files
    "MIN_CONTROL_FILE": INPUT_FILE_KEYS[10],
    "SIM_CONTROL_FILE": INPUT_FILE_KEYS[11],
    "PP_CONTROL_FILE": INPUT_FILE_KEYS[12],
    "SP_CONTROL_FILE": INPUT_FILE_KEYS[13]}

# List of Command names
COMMAND_NAMES_LIST = [
    # Lambda Simulation-related commands
    "LambdaSimulation",
    # Samping-related commands
    "DoubleWideSampling",
    # Free Energy calculation-related commands
    "ExponentialAveraging"]

# Dictionary of Command names
COMMAND_NAMES_DICT = {
    # Lambda Simulation-related commands
    "LAMBDA_SIMULATION": COMMAND_NAMES_LIST[0],
    # Samping-related commands
    "DOUBLE_WIDE_SAMPLING": COMMAND_NAMES_LIST[1],
    # Free Energy calculation-related commands
    "EXPONENTIAL_AVERAGING": COMMAND_NAMES_LIST[2]}

# Folder names
SIMULATION_FOLDER = "simulation/"
MINIMIZATION_FOLDER = "minimization/"
CALCULATION_FOLDER = "calculation/"
MODELS_FOLDER = "models/"

# File names
LOGFILE_NAME = "logfile_{}.txt"
SINGLE_POINT_CF_NAME = "pele_sp_{}.conf"
POST_PROCESSING_CF_NAME = "pele_recal_{}.conf"

# Direction definitions
DIRECTION_NAMES = ['BACKWARDS', 'FORWARD']
DIRECTION_CHARS = ['b', 'f']
DIRECTION_TO_CHAR = {DIRECTION_NAMES[0]: DIRECTION_CHARS[0],
                     DIRECTION_NAMES[1]: DIRECTION_CHARS[1]}
CHAR_TO_DIRECTION = {DIRECTION_CHARS[0]: DIRECTION_NAMES[0],
                     DIRECTION_CHARS[1]: DIRECTION_NAMES[1]}
DIRECTION_FACTORS = {DIRECTION_NAMES[0]: -1,
                     DIRECTION_NAMES[1]: +1}
DOUBLE_WIDE_SAMPLING_DIRECTIONS = DIRECTION_NAMES

# Physical constants
BOLTZMANN_CONSTANT_IN_KCAL_MOL = 0.0019872041

# Report file constants
PP_ABSOLUTE_ENERGIES_COL = 4
PP_DELTA_ENERGIES_CO = 5
