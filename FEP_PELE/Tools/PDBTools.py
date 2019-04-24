# -*- coding: utf-8 -*-


# Python imports
import sys


# FEP_PELE imports
from FEP_PELE.Utils.InOut import checkFile
from .Molecules import atomBuilder, linkBuilder, chainBuilder

# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Classes
class PDBParser:
    def __init__(self, pdb_path):
        try:
            checkFile(pdb_path)
        except NameError:
            raise NameError("PDBParser Error: no PDB file found in path " +
                            "{}".format(pdb_path))
        self._path = pdb_path

        self._atoms = []
        self._links = []
        self._chains = []

        self.__temporary_atoms_chunk = []
        self.__temporary_links_chunk = []

        self._processPDB()

    @property
    def all_links(self):
        return self._all_links

    @property
    def atoms(self):
        return self._atoms

    @property
    def links(self):
        return self._links

    @property
    def chains(self):
        return self._chains

    def _processPDB(self):
        with open(self._path, 'r') as file:
            for line in file:
                if (self._foundTER(line)):
                    self._processTER()
                    continue

                if (len(line) <= 6):
                    continue

                line_type = line[0:6]

                if (line_type == "ATOM  "):
                    self._processATOM(line)
                elif (line_type == "HETATM"):
                    self._processHETATM(line)
                elif (line_type == "CONECT"):
                    self._processCONECT()
                elif (line_type == "SEQRES"):
                    self._processSEQRES()
                elif (line_type == "MODEL "):
                    self._processMODEL()
                elif (line_type == "ENDMDL"):
                    self._processENDMDL()
                elif (line_type == "REMARK"):
                    self._processREMARK()
                else:
                    print("PDBParser Warning: unknown line type " +
                          "{}".format(line))

    def _foundTER(self, line):
        if (len(line) < 3):
            return False

        line_type = line[0:3]

        return line_type == "TER"

    def _processTER(self):
        # Add last link
        if (len(self.__temporary_atoms_chunk) == 0):
            print("PDBParser Error: invalid TER location was found")
            sys.exit(1)
        self._links.append(linkBuilder(self.__temporary_atoms_chunk))
        self.__temporary_atoms_chunk = []

        # Add last chain
        self.__temporary_links_chunk.append(self.links[-1])
        self._chains.append(chainBuilder(self.__temporary_links_chunk))
        self.__temporary_links_chunk = []

    def _processATOM(self, line):
        try:
            checkPDBLine(line)
        except NameError as e:
            raise NameError("PDBParser Error: invalid PBD_line, " +
                            str(e) + ': ' + str(line))

        self._atoms.append(atomBuilder(line))

        self._processLink()

    def _processLink(self):
        if (len(self.__temporary_atoms_chunk) == 0):
            self.__temporary_atoms_chunk.append(self.atoms[-1])
            return

        atom1 = self.atoms[-1]
        atom2 = self.atoms[-2]

        if ((atom1.chain != atom2.chain) or
                (atom1.residue_name != atom2.residue_name) or
                (atom1.residue_number != atom2.residue_number)):
            self._links.append(linkBuilder(self.__temporary_atoms_chunk))
            self.__temporary_atoms_chunk = [atom1, ]
            self._processChain()
        else:
            self.__temporary_atoms_chunk.append(atom1)

    def _processChain(self):
        if (len(self.__temporary_links_chunk) == 0):
            self.__temporary_links_chunk.append(self.links[-1])
            return

        link1 = self.links[-1]
        link2 = self.links[-2]

        if (link1.chain != link2.chain):
            self._chain.append(chainBuilder(self._temporary_links_chunk))
            self.__temporary_links_chunk = [link1, ]
        else:
            self.__temporary_links_chunk.append(link1)

    def _processHETATM(self, line):
        self._processATOM(line)
        self._atoms[-1].is_heteroatom = True

    def _processCONECT(self, line):
        pass

    def _processSEQRES(self, line):
        pass

    def _processMODEL(self, line):
        pass

    def _processENDMDL(self, line):
        pass

    def _processREMARK(self, line):
        pass


def checkPDBLine(PDB_line):
    okay = True
    messages = []

    # Check type
    if (type(PDB_line) != str):
        okay = False
        messages.append("invalid type")

    # Check length
    if (len(PDB_line) < 77):
        okay = False
        messages.append("invalid length")

    if (not okay):
        raise NameError(', '.join(messages))


def getLineType(PDB_line):
    try:
        checkPDBLine(PDB_line)
    except NameError as e:
        raise NameError("PDBTools.getLineType Error: invalid PBD_line, " +
                        str(e))

    return PDB_line[0:6]


def getAtomName(PDB_line):
    try:
        checkPDBLine(PDB_line)
    except NameError as e:
        raise NameError("PDBTools.getAtomName Error: invalid PBD_line, " +
                        str(e))

    return PDB_line[12:16]


class PDBHandler:
    def __init__(self, simulation):
        self.simulation = simulation
        self.trajectory = None
        self.are_atoms_indexed = False
        self.indexedAtoms = None
        self.system_size = None

    def getSystemSize(self):
        if (self.system_size is not None):
            return self.system_size

        if (self.trajectory is None):
            self.trajectory = self.simulation[0].trajectory

        path = self.trajectory.path + '/' + self.trajectory.name

        with open(path) as trajectory_file:
            for size, line in enumerate(trajectory_file):
                if line.startswith("ENDMDL"):
                    break

        return size + 1

    def indexAtoms(self):
        self.indexedAtoms = {}

        if (self.trajectory is None):
            self.trajectory = self.simulation[0].trajectory

        if self.system_size is None:
            self.system_size = self.getSystemSize()

        path = self.trajectory.path + '/' + self.trajectory.name

        with open(path) as trajectory_file:
            for i, line in enumerate(trajectory_file):
                if len(line) < 80:
                    continue
                if (int(i / (self.system_size + 1)) > 0):
                    break

                linetype = line.split()[0]
                chain = line[21]
                number = int(line[23:26])
                name = line[12:16]

                if linetype in ["HETATM", "ATOM"]:
                    if chain in self.indexedAtoms:
                        if number in self.indexedAtoms[chain]:
                            if type(self.indexedAtoms[chain][number]) != dict:
                                self.indexedAtoms[chain][number] = {}
                            self.indexedAtoms[chain][number][name] = i
                        else:
                            if type(self.indexedAtoms[chain]) != dict:
                                self.indexedAtoms[chain] = {}
                            self.indexedAtoms[chain][number] = {}
                            self.indexedAtoms[chain][number][name] = {}
                    else:
                        self.indexedAtoms[chain] = {}
                        self.indexedAtoms[chain][number] = {}
                        self.indexedAtoms[chain][number][name] = {}

                    self.indexedAtoms[chain][number][name] = i

        self.are_atoms_indexed = True

    def getAtomLineInPDB(self, atom_data):
        if (not self.are_atoms_indexed):
            print("PDBHandler:getAtomLineInPDB: Error, PDB atoms are not " +
                  "indexed")
            sys.exit(1)

        chain, number, name = atom_data
        number = int(number)

        try:
            line_number = self.indexedAtoms[chain][number][name]
        except KeyError:
            print("PDBHandler:getAtomLineInPDB: Error, atom " +
                  "{} not found in PDBHandler indexed ".format(atom_data) +
                  "information")
            sys.exit(1)
        return line_number

    def getLinkLinesInPDB(self, link_data):
        if (not self.are_atoms_indexed):
            print("PDBHandler:getAtomLineInPDB: Error, PDB atoms are not " +
                  "indexed")
            sys.exit(1)

        chain, number = link_data
        number = int(number)

        try:
            link_lines = self.indexedAtoms[chain][number]
        except KeyError:
            print("PDBHandler:getLinkLinesInPDB: Error, link " +
                  "{} not found in PDBHandler indexed ".format(link_data) +
                  "information")
            sys.exit(1)
        return link_lines

    def currentModel(self, line):
        return int(line / (self.system_size + 1))
