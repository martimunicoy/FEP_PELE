
import sys

from .Templates import TemplateOPLS2005
from .Combiner import CombineLinearly

from .Lambda import DUAL_LAMBDA, STERIC_LAMBDA, COULOMBIC_LAMBDA


# Class definitions
class AlchemicalTemplateCreator:
    def __init__(self, initial_template_path, final_template_path,
                 pdb_atom_name_pairs=[]):
        self.initial_template = TemplateOPLS2005(initial_template_path)
        self.final_template = TemplateOPLS2005(final_template_path)
        self.pdb_atom_name_pairs = pdb_atom_name_pairs
        self.explicit_template, self.implicit_template = \
            self.dectectExplicitAndImplicitTemplates()

        self.alchemicalTemplate = self.explicit_template

    def dectectExplicitAndImplicitTemplates(self, first_guess=True):
        explicit_guess = self.initial_template
        implicit_guess = self.final_template
        if (not first_guess):
            explicit_guess = self.final_template
            implicit_guess = self.initial_template

        for atom in implicit_guess.list_of_atoms.values():
            if (atom.pdb_atom_name in
                    [i.pdb_atom_name for i in
                     explicit_guess.list_of_atoms.values()]):
                continue
            elif ((atom.pdb_atom_name in
                   [i[0] for i in self.pdb_atom_name_pairs]) or
                  (atom.pdb_atom_name in
                   [i[1] for i in self.pdb_atom_name_pairs])):
                continue
            else:
                if (first_guess):
                    return self.dectectExplicitAndImplicitTemplates(False)
                else:
                    print("Error: there are unique atoms in both states")
                    sys.exit(1)

        return explicit_guess, implicit_guess

    def getFragmentAtoms(self):
        return detect_fragment_atoms(self.explicit_template,
                                     self.implicit_template)

    def reset(self):
        self.alchemicalTemplate = self.explicit_template

    def getFragmentAtomsAndBonds(self):
        fragment_atoms = detect_fragment_atoms(self.alchemicalTemplate,
                                               self.implicit_template)

        fragment_bonds = detect_fragment_bonds(fragment_atoms,
                                               self.alchemicalTemplate)

        set_fragment_atoms(list_of_fragment_atoms=fragment_atoms)
        set_fragment_bonds(list_of_fragment_bonds=fragment_bonds)

        atoms_pairs = []

        for name1, name2 in self.pdb_atom_name_pairs:
            if (not self.explicit_is_final):
                t = name1
                name1 = name2
                name2 = t

            atoms_pairs.append(set_connecting_atoms(self.implicit_template,
                                                    name1,
                                                    self.alchemicalTemplate,
                                                    name2))

        bonds_pairs = []

        for atoms_pair in atoms_pairs:
            bonds_pairs.append(set_connecting_bonds(atoms_pair,
                                                    self.implicit_template,
                                                    self.alchemicalTemplate))

        return atoms_pairs, bonds_pairs

    def applyLambda(self, _lambda):
        atoms_pairs, bonds_pairs = self.getFragmentAtomsAndBonds()

        combiner = CombineLinearly(self.alchemicalTemplate,
                                   _lambda.value, atoms_pairs,
                                   bonds_pairs, self.explicit_is_final)

        if ((_lambda.type == DUAL_LAMBDA) or
                (_lambda.type == STERIC_LAMBDA)):
            combiner.combine_sigmas()
            combiner.combine_bond_eq_dist()
            combiner.combine_radnpSGB()

        if ((_lambda.type == DUAL_LAMBDA) or
                (_lambda.type == COULOMBIC_LAMBDA)):
            combiner.combine_charges()

        self.alchemicalTemplate = combiner.get_resulting_template()

    def writeAlchemicalTemplate(self, output_path):
        self.alchemicalTemplate.write_template_to_file(
            template_new_name=output_path)

    @property
    def explicit_is_final(self):
        return self.explicit_template == self.final_template


def detect_fragment_atoms(explicit_template, implicit_template):
    fragment_atoms = []
    core_atoms = find_equal_pdb_atom_names(explicit_template,
                                           implicit_template)
    for key, atom in explicit_template.list_of_atoms.items():
        pdb_atom_name = atom.pdb_atom_name
        if pdb_atom_name not in core_atoms:
            fragment_atoms.append(atom)
    return fragment_atoms


def find_equal_pdb_atom_names(template1, template2):
    pdb_atom_names_tmpl_1 = [template1.list_of_atoms[n].pdb_atom_name
                             for n in range(1,
                                            len(template1.list_of_atoms) + 1)]

    pdb_atom_names_tmpl_2 = [template2.list_of_atoms[n].pdb_atom_name
                             for n in range(1,
                                            len(template2.list_of_atoms) + 1)]

    return list(set(pdb_atom_names_tmpl_1).intersection(pdb_atom_names_tmpl_2))


def find_unique_pdb_atom_names(template1, template2):
    names1 = [template1.list_of_atoms[n].pdb_atom_name
              for n in range(1, len(template1.list_of_atoms) + 1)]

    names2 = [template2.list_of_atoms[n].pdb_atom_name
              for n in range(1, len(template2.list_of_atoms) + 1)]

    common_atoms = set(names1).intersection(names2)

    unique_atoms = set(names1).union(names2) - common_atoms

    return list(unique_atoms)


def set_fragment_atoms(list_of_fragment_atoms):
    for atom in list_of_fragment_atoms:
        atom.is_unique = True


def detect_fragment_bonds(list_of_fragment_atoms, exclusive_template):
    fragment_bonds = []
    fragment_indexes = []
    for atom in list_of_fragment_atoms:
        fragment_indexes.append(atom.atom_id)
    fragment_indexes = list(set(fragment_indexes))
    for key, bond in exclusive_template.list_of_bonds.items():
        if key[0] in fragment_indexes or key[1] in fragment_indexes:
            fragment_bonds.append(bond)
    return fragment_bonds


def set_fragment_bonds(list_of_fragment_bonds):
    for bond in list_of_fragment_bonds:
        bond.is_unique = True


def set_connecting_atoms(template1, pdb_atom_name1, template2, pdb_atom_name2):
    atom1 = template1.get_atom_by_pdb_atom_name(pdb_atom_name1)
    atom2 = template2.get_atom_by_pdb_atom_name(pdb_atom_name2)

    atom1.is_linker = True
    atom2.is_linker = True

    return (atom1, atom2)


def set_connecting_bonds(atoms_pair, template1, template2):
    b1 = None
    b2 = None

    for key, bond in template1.list_of_bonds.items():
        if (atoms_pair[0].atom_id in key):
            index = int(key.index(atoms_pair[0].atom_id) == 0)
            if (not template1.list_of_atoms[key[index]].is_unique):
                b1 = bond

    for key, bond in template2.list_of_bonds.items():
        if (atoms_pair[1].atom_id in key):
            index = int(key.index(atoms_pair[1].atom_id) == 0)
            if (not template2.list_of_atoms[key[index]].is_unique):
                b2 = bond

    if (b1 is None or b2 is None):
        raise NameError("No bond was found between atoms " +
                        "{} and ".format(atoms_pair[0].pdb_atom_name) +
                        "{}".format(atoms_pair[1].pdb_atom_name))

    return (b1, b2)
