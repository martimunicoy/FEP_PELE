import copy


class Combiner:
    def __init__(self, explicit_template, lambda_parameter,
                 combiner_function, atoms_pairs, bonds_pairs,
                 explicit_is_final):
        self.explicit_template = explicit_template
        self.lambda_parameter = lambda_parameter
        self.combiner_function = combiner_function
        self.atoms_pairs = atoms_pairs
        self.bonds_pairs = bonds_pairs
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
            radnpType = atom.radnpType
            result = self.combiner_function(0, radnpType)

            self.new_template.list_of_atoms[key].radnpType = result

    def combine_bond_eq_dist(self):
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
            result = self.combiner_function(0, eq_dist)

            self.new_template.list_of_bonds[key].eq_dist = result

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


class CombineLinearly(Combiner):
    def __init__(self, explicit_template, lambda_parameter,
                 atoms_pairs=[], bonds_pairs=[], explicit_is_final=True):
        Combiner.__init__(self, explicit_template, lambda_parameter,
                          self.combiner_function, atoms_pairs, bonds_pairs,
                          explicit_is_final)
        self.explicit_is_final = explicit_is_final

    def combiner_function(self, initial_value, final_value):
        if (self.explicit_is_final):
            result = float(initial_value * (1 - self.lambda_parameter) +
                           final_value * self.lambda_parameter)
        else:
            result = float(final_value * (1 - self.lambda_parameter) +
                           initial_value * self.lambda_parameter)

        return result
