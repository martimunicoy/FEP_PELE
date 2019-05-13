# -*- coding: utf-8 -*-


# FEP_PELE imports
from FEP_PELE.PELETools import PELEConstants as pele_co


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


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
    "SamplingMethod",
    "NumberOfProcessors",
    "TotalPELESteps",
    "SafetyCheck",
    "SolventType",
    "Reminimize",
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
    "SAMPLING_METHOD": INPUT_FILE_KEYS[9],
    "NUMBER_OF_PROCESSORS": INPUT_FILE_KEYS[10],
    "TOTAL_PELE_STEPS": INPUT_FILE_KEYS[11],
    "SAFETY_CHECK": INPUT_FILE_KEYS[12],
    "SOLVENT_TYPE": INPUT_FILE_KEYS[13],
    "REMINIMIZE": INPUT_FILE_KEYS[14],
    "RESTART": INPUT_FILE_KEYS[15],
    # List of commands
    "COMMANDS": INPUT_FILE_KEYS[16],
    # PELE control files
    "MIN_CONTROL_FILE": INPUT_FILE_KEYS[17],
    "SIM_CONTROL_FILE": INPUT_FILE_KEYS[18],
    "PP_CONTROL_FILE": INPUT_FILE_KEYS[19],
    "SP_CONTROL_FILE": INPUT_FILE_KEYS[20],
    # Folder names
    "MIN_FOLDER": INPUT_FILE_KEYS[21],
    "SIM_FOLDER": INPUT_FILE_KEYS[22],
    "CAL_FOLDER": INPUT_FILE_KEYS[23],
    # Input PDB path
    "INPUT_PDB": INPUT_FILE_KEYS[24],
    "INITIAL_LIGAND_PDB": INPUT_FILE_KEYS[25],
    "FINAL_LIGAND_PDB": INPUT_FILE_KEYS[26]}

# List of Command names
COMMAND_NAMES_LIST = [
    # Bound state-related commands
    "LambdasSampling",
    "dECalculation",
    # Analysis commands
    "ExponentialAveraging",
    # Unbound state-related commands
    "Unbound-dECalculation",
    # Others
    "SolvationFreeEnergyCalculation"]

# Dictionary of Command names
COMMAND_NAMES_DICT = {
    # Bound state-related commands
    "LAMBDAS_SAMPLING": COMMAND_NAMES_LIST[0],
    "DE_CALCULATION": COMMAND_NAMES_LIST[1],
    # Analysis commands
    "EXPONENTIAL_AVERAGING": COMMAND_NAMES_LIST[2],
    # Unbound state-related commands
    "UNBOUND_DE_CALCULATION": COMMAND_NAMES_LIST[3],
    # Others
    "SOLVATION_FREE_ENERGY_CALCULATION": COMMAND_NAMES_LIST[4]}

# Dictionary of Command labels
COMMAND_LABELS_DICT = {
    # Bound state-related commands
    "LAMBDAS_SAMPLING": "Lambda Sampling",
    "DE_CALCULATION": "dE Calculation",
    # Analysis commands
    "EXPONENTIAL_AVERAGING": "Exponential Averaging",
    # Unbound state-related commands
    "UNBOUND_DE_CALCULATION": "Unbound dE Calculation",
    # Others
    "SOLVATION_FREE_ENERGY_CALCULATION": "Solvation Free Energy Calculation"}

# List of sampling methods
SAMPLING_METHODS_LIST = [
    "DoubleWide",
    "Overlap",
    "DoubleEnded"]

# Dictionary of sampling methods names
SAMPLING_METHODS_DICT = {
    "DOUBLE_WIDE": SAMPLING_METHODS_LIST[0],
    "OVERLAP": SAMPLING_METHODS_LIST[1],
    "DOUBLE_ENDED": SAMPLING_METHODS_LIST[2]}

# Dictionary of sampling methods labels
SAMPLING_METHODS_NAMES = {
    "DOUBLE_WIDE": "double-wide sampling [0 <-- M --> 1]",
    "OVERLAP": "overlap sampling [0 --> M <-- 1]",
    "DOUBLE_ENDED": "double-ended sampling [0 <--> 1]"}

# Default settings
DEF_SERIAL_PELE = None
DEF_MPI_PELE = None
DEF_INITIAL_TEMPLATE = None
DEF_FINAL_TEMPLATE = None
DEF_ATOM_LINKS = []
DEF_LAMBDAS = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]
DEF_LJ_LAMBDAS = []
DEF_C_LAMBDAS = []
DEF_SAMPLING_METHOD = SAMPLING_METHODS_DICT["DOUBLE_WIDE"]
DEF_NUMBER_OF_PROCESSORS = 4
DEF_TOTAL_PELE_STEPS = 200
DEF_COMMANDS = []
DEF_MIN_CONTROL_FILE = None
DEF_SIM_CONTROL_FILE = None
DEF_PP_CONTROL_FILE = None
DEF_SP_CONTROL_FILE = None
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
DEF_REMINIMIZE = True

# Folder names
MODELS_FOLDER = "models/"

# File names
LOGFILE_NAME = "logfile_{}.txt"
PDB_OUT_NAME = "pele_out_{}.pdb"
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
DIRECTION_LABELS = {DIRECTION_NAMES[0]: "backwards",
                    DIRECTION_NAMES[1]: "forward"}
DIRECTION_TO_CHAR = {DIRECTION_NAMES[0]: DIRECTION_CHARS[0],
                     DIRECTION_NAMES[1]: DIRECTION_CHARS[1]}
CHAR_TO_DIRECTION = {DIRECTION_CHARS[0]: DIRECTION_NAMES[0],
                     DIRECTION_CHARS[1]: DIRECTION_NAMES[1]}
DIRECTION_FACTORS = {DIRECTION_NAMES[0]: -1,
                     DIRECTION_NAMES[1]: +1}

# Physical constants
BOLTZMANN_CONSTANT_IN_KCAL_MOL = 0.0019872041

# Report file constants
PP_STEPS_COL = 2
PP_ABSOLUTE_ENERGIES_COL = 4
PP_DELTA_ENERGIES_COL = 5
