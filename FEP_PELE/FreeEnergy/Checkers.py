# -*- coding: utf-8 -*-


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Function definitions
def checkModelCoords(path, atoms):
    coords = set()

    atom_names = []
    for atom in atoms:
        atom_names.append(atom.pdb_atom_name)

    with open(path, 'r') as file:
        for line in file:
            atom_name = line[12:16].replace(' ', '_')
            if (atom_name in atom_names):
                coords.add((line[31:38],
                            line[39:46],
                            line[47:54]))

    if (len(coords) != len(atom_names)):
        return False

    return True
