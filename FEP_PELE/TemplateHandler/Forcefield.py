from .Patterns import PATTERN_OPLS2005_RESX_LINE
from .Patterns import PATTERN_OPLS2005_NBON
from .Patterns import PATTERN_OPLS2005_BOND
from .Patterns import PATTERN_OPLS2005_THETA
from .Patterns import PATTERN_OPLS2005_PHI
from .Patterns import PATTERN_OPLS2005_IPHI


class Atom:
    def __init__(self, atom_id, parent_id, location, atom_type, pdb_atom_name,
                 unknown, x_zmatrix=0, y_zmatrix=0, z_zmatrix=0, sigma=0,
                 epsilon=0, charge=0, radnpSGB=0, radnpType=0, sgbnpGamma=0,
                 sgbnpType=0, is_linker=False, is_unique=False):
        self.atom_id = int(atom_id)
        self.parent_id = int(parent_id)
        self.location = str(location)
        self.atom_type = str(atom_type)
        self.pdb_atom_name = str(pdb_atom_name)
        self.unknown = int(unknown)
        self.x_zmatrix = float(x_zmatrix)
        self.y_zmatrix = float(y_zmatrix)
        self.z_zmatrix = float(z_zmatrix)
        self.sigma = float(sigma)
        self.epsilon = float(epsilon)
        self.charge = float(charge)
        self.radnpSGB = float(radnpSGB)
        self.radnpType = float(radnpType)
        self.sgbnpGamma = float(sgbnpGamma)
        self.sgbnpType = float(sgbnpType)
        self.bonds = []
        self.thetas = []
        self.phis = []
        self.iphis = []
        self.is_unique = bool(is_unique)
        self.is_linker = bool(is_linker)
        self.is_root = False
        self.parents = []
        self.childs = []

    def __eq__(self, other):
        return bool((self.atom_id, self.pdb_atom_name) ==
                    (other.atom_id, other.pdb_atom_name))

    def __ne__(self, other):
        return not(self == other)

    def __hash__(self):
        return hash((self.atom_id, self.pdb_atom_name))

    def write_resx(self):
        return PATTERN_OPLS2005_RESX_LINE.format(self.atom_id, self.parent_id,
                                                 self.location, self.atom_type,
                                                 self.pdb_atom_name,
                                                 self.unknown, self.x_zmatrix,
                                                 self.y_zmatrix,
                                                 self.z_zmatrix)

    def write_nbon(self):
        return PATTERN_OPLS2005_NBON.format(self.atom_id, self.sigma,
                                            self.epsilon, self.charge,
                                            self.radnpSGB, self.radnpType,
                                            self.sgbnpGamma, self.sgbnpType)

    def calculateMinimumDistanceWithAny(self, list_of_atoms):
        return self._recursiveDistanceFollowerToAny(list_of_atoms)

    def _recursiveDistanceFollowerToAny(self, list_of_atoms, current_atom=None,
                                        distance=0, already_visited=None):
        if (current_atom is None):
            current_atom = self
        if (already_visited is None):
            already_visited = []

        if (current_atom in list_of_atoms):
            return distance

        already_visited.append(current_atom)

        distance += 1

        # Calculate distance for each child
        other_distances = []

        for child in current_atom.childs:
            if (child not in already_visited):
                other_distances.append(
                    self._recursiveDistanceFollowerToAny(list_of_atoms,
                                                         child, distance,
                                                         already_visited))

        # Calculate distance for each parent
        for parent in current_atom.parents:
            if (parent not in already_visited):
                other_distances.append(
                    self._recursiveDistanceFollowerToAny(list_of_atoms,
                                                         parent, distance,
                                                         already_visited))

        # Filter None values
        other_distances = list(filter(None.__ne__, other_distances))

        # Return minimum distance
        if (len(other_distances) > 0):
            return min(other_distances)


class Bond:
    def __init__(self, atom1, atom2, spring, eq_dist, is_unique=False,
                 is_linker=False):
        self.atom1 = int(atom1)
        self.atom2 = int(atom2)
        self.spring = float(spring)
        self.eq_dist = float(eq_dist)
        self.is_unique = bool(is_unique)
        self.is_linker = bool(is_linker)

    def write_bond(self):
        return PATTERN_OPLS2005_BOND.format(self.atom1, self.atom2,
                                            self.spring, self.eq_dist)


class Theta:
    def __init__(self, atom1, atom2, atom3, spring, eq_angle,
                 is_unique=False):
        self.atom1 = int(atom1)
        self.atom2 = int(atom2)
        self.atom3 = int(atom3)
        self.spring = float(spring)
        self.eq_angle = float(eq_angle)
        self.is_unique = bool(is_unique)

    def write_theta(self):
        return PATTERN_OPLS2005_THETA.format(self.atom1, self.atom2,
                                             self.atom3, self.spring,
                                             self.eq_angle)


class Phi:
    def __init__(self, atom1, atom2, atom3, atom4, constant, prefactor, nterm,
                 improper, is_unique=False):
        self.atom1 = int(atom1)
        self.atom2 = int(atom2)
        self.atom3 = int(atom3)
        self.atom4 = int(atom4)
        self.constant = float(constant)
        self.prefactor = float(prefactor)
        self.nterm = float(nterm)
        self.improper = bool(improper)
        self.is_unique = bool(is_unique)

    def write_phi(self):
        if not self.improper:
            return PATTERN_OPLS2005_PHI.format(self.atom1, self.atom2,
                                               self.atom3, self.atom4,
                                               self.constant, self.prefactor,
                                               self.nterm)

    def write_iphi(self):
        if self.improper:
            return PATTERN_OPLS2005_IPHI.format(self.atom1, self.atom2,
                                                self.atom3, self.atom4,
                                                self.constant, self.prefactor,
                                                self.nterm)
