class Reductor:
    def __init__(self, template, lambda_to_reduce):
        self.template = template
        self.lambda_to_reduce = lambda_to_reduce

    def reduce_epsilons(self, function):
        atoms = self.template.get_list_of_fragment_atoms()
        for key, atom in atoms:
            epsilon = atom.epsilon
            result = function(epsilon)
            self.template.list_of_atoms[key].epsilon = result

    def reduce_sigmas(self, function):
        atoms = self.template.get_list_of_fragment_atoms()
        for key, atom in atoms:
            sigma = atom.sigma
            result = function(sigma)
            self.template.list_of_atoms[key].sigma = result

    def reduce_charges(self, function):
        atoms = self.template.get_list_of_fragment_atoms()
        for key, atom in atoms:
            charge = atom.charge
            result = function(charge)
            self.template.list_of_atoms[key].charge = result

    def reduce_radnpSGB(self, function):
        atoms = self.template.get_list_of_fragment_atoms()
        for key, atom in atoms:
            radnpSGB = atom.radnpSGB
            result = function(radnpSGB)
            self.template.list_of_atoms[key].radnpSGB = result

    def reduce_radnpType(self, function):
        atoms = self.template.get_list_of_fragment_atoms()
        for key, atom in atoms:
            radnpType = atom.radnpType
            result = function(radnpType)
            self.template.list_of_atoms[key].radnpType = result

    def reduce_bond_eq_dist(self, function):
        bonds = self.template.get_list_of_fragment_bonds()
        for key, bond in bonds:
            eq_dist = bond.eq_dist
            result = function(eq_dist)
            self.template.list_of_bonds[key].eq_dist = result


class ReduceLinearly(Reductor):
    def __init__(self, template, lambda_to_reduce):
        Reductor.__init__(self, template, lambda_to_reduce)

    def reduce_value(self, value):
        result = value * self.lambda_to_reduce
        return result


class ReduceExponentially(Reductor):
    def __init__(self, template, lambda_to_reduce):
        Reductor.__init__(self, template, lambda_to_reduce)

    def reduce_value(self, value):
        result = value * self.lambda_to_reduce ** 2
        return result
