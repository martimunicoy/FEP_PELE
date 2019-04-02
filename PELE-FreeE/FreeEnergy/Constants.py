# Default settings
DEF_GENERAL_PATH = ""
DEF_INITIAL_TEMPLATE = ""
DEF_FINAL_TEMPLATE = ""
DEF_ATOM_LINKS = []
DEF_LAMBDAS = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]
DEF_NUMBER_OF_PROCESSORDS = 4
DEF_COMMANDS = []

# Input file keys
INPUT_FILE_KEYS = [
    # General path
    "GeneralPath",
    # Template paths
    "InitialTemplate",
    "FinalTemplate",
    # Atom links
    "AtomLink",
    # Free Energy Settings
    "Lambdas",
    "NumberOfProcessors"]

# Input file dict
CONTROL_FILE_DICT = {
    # General path
    "GENERAL_PATH": INPUT_FILE_KEYS[0],
    # Template paths
    "INITIAL_TEMPLATE": INPUT_FILE_KEYS[1],
    "FINAL_TEMPLATE": INPUT_FILE_KEYS[2],
    # Atom links
    "ATOM_LINK": INPUT_FILE_KEYS[3],
    # Free Energy Settings
    "LAMBDAS": INPUT_FILE_KEYS[4],
    "NUMBER_OF_PROCESSORS": INPUT_FILE_KEYS[5]}
