# Python imports
import copy


# Script information
__author__ = "Marti Municoy"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Marti Municoy"
__email__ = "marti.municoy@bsc.es"


# Constant definitions
INITIAL_DUMMY_BOND_EQDIST = float(0.1)
INITIAL_DUMMY_BOND_SPRING = float(0.0)
INITIAL_DUMMY_THETA_SPRING = float(0.0)


# Class definitions
class Combiner:
    def __init__(self, explicit_template, lambda_parameter,
                 combiner_function, atoms_pairs, bonds_pairs,
                 thetas_pairs, explicit_is_final):
        self.explicit_template = explicit_template
        self.lambda_parameter = lambda_parameter
        self.combiner_function = combiner_function
        self.atoms_pairs = atoms_pairs
        self.bonds_pairs = bonds_pairs
        self.thetas_pairs = thetas_pairs
        self.new_template = copy.deepcopy(explicit_template)
        self.explicit_is_final = explicit_is_final

    def combine_epsilons(self):
        # Modify paired atoms
        for atom_pair in self.atoms_pairs:
            result = self.combiner_function(atom_pair[0].epsilon,
                                            atom_pair[1].epsilon)
            index = self._getAtomPairIndex(atom_pair)
            self.new_template.list_of_atoms[index].epsilon = result

        # Modify rest of atoms
        atoms = self.explicit_template.get_list_of_fragment_atoms()
        all_paired_atoms = [i for sub in self.atoms_pairs for i in sub]

        for key, atom in atoms:
            if (atom in all_paired_atoms):
                continue
            epsilon = atom.epsilon
            result = self.combiner_function(0, epsilon)

            self.new_template.list_of_atoms[key].epsilon = result

    def combine_sigmas(self):
        # Modify paired atoms
        for atom_pair in self.atoms_pairs:
            result = self.combiner_function(atom_pair[0].sigma,
                                            atom_pair[1].sigma)
            index = self._getAtomPairIndex(atom_pair)
            self.new_template.list_of_atoms[index].sigma = result

        # Modify rest of atoms
        atoms = self.explicit_template.get_list_of_fragment_atoms()
        all_paired_atoms = [i for sub in self.atoms_pairs for i in sub]

        for key, atom in atoms:
            if (atom in all_paired_atoms):
                continue
            sigma = atom.sigma
            result = self.combiner_function(0, sigma)

            self.new_template.list_of_atoms[key].sigma = result

    def combine_charges(self):
        # Modify paired atoms
        for atom_pair in self.atoms_pairs:
            result = self.combiner_function(atom_pair[0].charge,
                                            atom_pair[1].charge)
            index = self._getAtomPairIndex(atom_pair)
            self.new_template.list_of_atoms[index].charge = result

        # Modify rest of atoms
        atoms = self.explicit_template.get_list_of_fragment_atoms()
        all_paired_atoms = [i for sub in self.atoms_pairs for i in sub]

        for key, atom in atoms:
            if (atom in all_paired_atoms):
                continue
            charge = atom.charge
            result = self.combiner_function(0, charge)

            self.new_template.list_of_atoms[key].charge = result

    def combine_radnpSGB(self):
        # Modify paired atoms
        for atom_pair in self.atoms_pairs:
            result = self.combiner_function(atom_pair[0].radnpSGB,
                                            atom_pair[1].radnpSGB)
            index = self._getAtomPairIndex(atom_pair)
            self.new_template.list_of_atoms[index].radnpSGB = result

        # Modify rest of atoms
        atoms = self.explicit_template.get_list_of_fragment_atoms()
        all_paired_atoms = [i for sub in self.atoms_pairs for i in sub]

        for key, atom in atoms:
            if (atom in all_paired_atoms):
                continue
            radnpSGB = atom.radnpSGB
            result = self.combiner_function(0, radnpSGB)

            self.new_template.list_of_atoms[key].radnpSGB = result

    def combine_radnpType(self):
        # Modify paired atoms
        for atom_pair in self.atoms_pairs:
            result = self.combiner_function(atom_pair[0].radnpType,
                                            atom_pair[1].radnpType)
            index = self._getAtomPairIndex(atom_pair)
            self.new_template.list_of_atoms[index].radnpType = result

        # Modify rest of atoms
        atoms = self.explicit_template.get_list_of_fragment_atoms()
        all_paired_atoms = [i for sub in self.atoms_pairs for i in sub]

        for key, atom in atoms:
            if (atom in all_paired_atoms):
                continue
            radnpType = atom.radnpType
            result = self.combiner_function(0, radnpType)

            self.new_template.list_of_atoms[key].radnpType = result

    def combine_SGBNPGamma(self):
        # Modify paired atoms
        for atom_pair in self.atoms_pairs:
            result = self.combiner_function(atom_pair[0].sgbnpGamma,
                                            atom_pair[1].sgbnpGamma)
            index = self._getAtomPairIndex(atom_pair)
            self.new_template.list_of_atoms[index].sgbnpGamma = result

        # Modify rest of atoms
        atoms = self.explicit_template.get_list_of_fragment_atoms()
        all_paired_atoms = [i for sub in self.atoms_pairs for i in sub]

        for key, atom in atoms:
            if (atom in all_paired_atoms):
                continue
            sgbnpGamma = atom.sgbnpGamma
            result = self.combiner_function(0, sgbnpGamma)

            self.new_template.list_of_atoms[key].sgbnpGamma = result

    def combine_SGBNPType(self):
        # Modify paired atoms
        for atom_pair in self.atoms_pairs:
            result = self.combiner_function(atom_pair[0].sgbnpType,
                                            atom_pair[1].sgbnpType)
            index = self._getAtomPairIndex(atom_pair)
            self.new_template.list_of_atoms[index].sgbnpType = result

        # Modify rest of atoms
        atoms = self.explicit_template.get_list_of_fragment_atoms()
        all_paired_atoms = [i for sub in self.atoms_pairs for i in sub]

        for key, atom in atoms:
            if (atom in all_paired_atoms):
                continue
            sgbnpType = atom.sgbnpType
            result = self.combiner_function(0, sgbnpType)

            self.new_template.list_of_atoms[key].sgbnpType = result

    def combine_BondEqDist(self):
        # Modify paired bonds
        for bond_pair in self.bonds_pairs:
            result = self.combiner_function(bond_pair[0].eq_dist,
                                            bond_pair[1].eq_dist)
            index = self._getBondPairIndex(bond_pair)
            self.new_template.list_of_bonds[index].eq_dist = result

        # Modify rest of bonds
        bonds = self.explicit_template.get_list_of_fragment_bonds()
        all_paired_bonds = [i for sub in self.bonds_pairs for i in sub]

        for key, bond in bonds:
            if (bond in all_paired_bonds):
                continue
            eq_dist = bond.eq_dist
            result = self.combiner_function(INITIAL_DUMMY_BOND_EQDIST, eq_dist)

            self.new_template.list_of_bonds[key].eq_dist = result

    def combine_BondSprings(self):
        # Modify paired bonds
        for bond_pair in self.bonds_pairs:
            result = self.combiner_function(bond_pair[0].spring,
                                            bond_pair[1].spring)
            index = self._getBondPairIndex(bond_pair)
            self.new_template.list_of_bonds[index].spring = result

        # Modify rest of bonds
        bonds = self.explicit_template.get_list_of_fragment_bonds()
        all_paired_bonds = [i for sub in self.bonds_pairs for i in sub]

        for key, bond in bonds:
            if (bond in all_paired_bonds):
                continue
            spring = bond.spring
            result = self.combiner_function(INITIAL_DUMMY_BOND_SPRING, spring)

            self.new_template.list_of_bonds[key].spring = result

    def combine_ThetaSprings(self):
        # Modify paired thetas
        for thetas_pair in self.thetas_pairs:
            result = self.combiner_function(thetas_pair[0].spring,
                                            thetas_pair[1].spring)
            index = self._getThetaPairIndex(thetas_pair)
            self.new_template.list_of_thetas[index].spring = result

        # Modify rest of thetas
        thetas = self.explicit_template.get_list_of_fragment_thetas()
        all_paired_thetas = [i for sub in self.thetas_pairs for i in sub]

        for key, theta in thetas:
            if (theta in all_paired_thetas):
                continue
            spring = theta.spring
            result = self.combiner_function(INITIAL_DUMMY_THETA_SPRING, spring)

            self.new_template.list_of_thetas[key].spring = result

    def add_connecting_atoms_pair(self, atoms_pair):
        self.atoms_pairs.append(atoms_pair)

    def add_connecting_bonds_pair(self, bonds_pair):
        self.bonds_pairs.append(bonds_pair)

    def get_resulting_template(self):
        return self.new_template

    def _getAtomPairIndex(self, atom_pair):
        if (self.explicit_is_final):
            index = atom_pair[1].atom_id
        else:
            index = atom_pair[0].atom_id
        return index

    def _getBondPairIndex(self, bond_pair):
        if (self.explicit_is_final):
            index = (bond_pair[1].atom1, bond_pair[1].atom2)
        else:
            index = (bond_pair[0].atom1, bond_pair[0].atom2)
        return index

    def _getThetaPairIndex(self, theta_pair):
        if (self.explicit_is_final):
            index = (theta_pair[1].atom1, theta_pair[1].atom2,
                     theta_pair[1].atom3)
        else:
            index = (theta_pair[0].atom1, theta_pair[0].atom2,
                     theta_pair[0].atom3)
        return index


class CombineLinearly(Combiner):
    def __init__(self, explicit_template, lambda_parameter,
                 atoms_pairs=[], bonds_pairs=[], thetas_pairs=[],
                 explicit_is_final=True):
        Combiner.__init__(self, explicit_template, lambda_parameter,
                          self.combiner_function, atoms_pairs, bonds_pairs,
                          thetas_pairs, explicit_is_final)
        self.explicit_is_final = explicit_is_final

    def combiner_function(self, initial_value, final_value):
        if (self.explicit_is_final):
            result = float(initial_value * (1 - self.lambda_parameter) +
                           final_value * self.lambda_parameter)
        else:
            result = float(final_value * (1 - self.lambda_parameter) +
                           initial_value * self.lambda_parameter)

        return result
