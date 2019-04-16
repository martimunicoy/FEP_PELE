

from FEP_PELE.PELETools import PELEConstants as pele_co

# Default settings
DEF_SERIAL_PELE = ""
DEF_MPI_PELE = ""
DEF_INITIAL_TEMPLATE = ""
DEF_FINAL_TEMPLATE = ""
DEF_ATOM_LINKS = []
DEF_LAMBDAS = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]
DEF_LJ_LAMBDAS = []
DEF_C_LAMBDAS = []
DEF_NUMBER_OF_PROCESSORS = 4
DEF_COMMANDS = []
DEF_MIN_CONTROL_FILE = "min_pele.conf"
DEF_SIM_CONTROL_FILE = "sim_pele.conf"
DEF_PP_CONTROL_FILE = "pele_recal.conf"
DEF_SP_CONTROL_FILE = "pele_sp.conf"
DEF_SAFETY_CHECK = False
DEF_MIN_FOLDER = "minimization/"
DEF_SIM_FOLDER = "simulation/"
DEF_CAL_FOLDER = "calculation/"
DEF_INPUT_PDB = ""
DEF_INITIAL_LIGAND_PDB = ""
DEF_FINAL_LIGAND_PDB = ""
DEF_SOLVENT_TYPE = pele_co.SGBNP_TYPE_NAME
DEF_LAMBDA_SPLITTING = False
DEF_RESTART = False

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
    "StericLambdas",
    "CoulombicLambdas",
    "NumberOfProcessors",
    "SafetyCheck",
    "SolventType",
    "Restart",
    # List of commands
    "Commands",
    # PELE control files
    "MinimizationControlFile",
    "SimulationControlFile",
    "PostProcessingControlFile",
    "SinglePointControlFile",
    # Folder names
    "MinimizationFolder",
    "SimulationFolder",
    "CalculationFolder",
    # Input PDBs
    "InputPDB",
    "InitialLigandPDB",
    "FinalLigandPDB"]

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
    "LJ_LAMBDAS": INPUT_FILE_KEYS[7],
    "C_LAMBDAS": INPUT_FILE_KEYS[8],
    "NUMBER_OF_PROCESSORS": INPUT_FILE_KEYS[9],
    "SAFETY_CHECK": INPUT_FILE_KEYS[10],
    "SOLVENT_TYPE": INPUT_FILE_KEYS[11],
    "RESTART": INPUT_FILE_KEYS[12],
    # List of commands
    "COMMANDS": INPUT_FILE_KEYS[13],
    # PELE control files
    "MIN_CONTROL_FILE": INPUT_FILE_KEYS[14],
    "SIM_CONTROL_FILE": INPUT_FILE_KEYS[15],
    "PP_CONTROL_FILE": INPUT_FILE_KEYS[16],
    "SP_CONTROL_FILE": INPUT_FILE_KEYS[17],
    # Folder names
    "MIN_FOLDER": INPUT_FILE_KEYS[18],
    "SIM_FOLDER": INPUT_FILE_KEYS[19],
    "CAL_FOLDER": INPUT_FILE_KEYS[20],
    # Input PDB path
    "INPUT_PDB": INPUT_FILE_KEYS[21],
    "INITIAL_LIGAND_PDB": INPUT_FILE_KEYS[22],
    "FINAL_LIGAND_PDB": INPUT_FILE_KEYS[23]}

# List of Command names
COMMAND_NAMES_LIST = [
    # Lambda Simulation-related commands
    "LambdaSimulation",
    # Samping-related commands
    "DoubleWideSampling",
    # Free Energy calculation-related commands
    "ExponentialAveraging",
    "SolvationFreeEnergyCalculation"]

# Dictionary of Command names
COMMAND_NAMES_DICT = {
    # Lambda Simulation-related commands
    "LAMBDA_SIMULATION": COMMAND_NAMES_LIST[0],
    # Samping-related commands
    "DOUBLE_WIDE_SAMPLING": COMMAND_NAMES_LIST[1],
    # Free Energy calculation-related commands
    "EXPONENTIAL_AVERAGING": COMMAND_NAMES_LIST[2],
    "SOLVATION_FREE_ENERGY_CALCULATION": COMMAND_NAMES_LIST[3]}

# Folder names
MODELS_FOLDER = "models/"

# File names
LOGFILE_NAME = "logfile_{}.txt"
PDB_OUT_NAME = "pele_out.pdb"
SINGLE_POINT_CF_NAME = "pele_sp_{}.conf"
POST_PROCESSING_CF_NAME = "pele_recal_{}.conf"
MINIMIZATION_CF_NAME = "pele_min.conf"
SINGLE_LOGFILE_NAME = "logfile.txt"
SINGLE_REPORT_NAME = "report.out"
SINGLE_TRAJECTORY_NAME = "trajectory.pdb"
CHECKPOINT_NAME = ".FEP_PELE.ckp"

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
