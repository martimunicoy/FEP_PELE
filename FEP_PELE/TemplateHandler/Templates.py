from .Headers import HEADER_OPLS2005
from .Patterns import PATTERN_OPLS2005_RESX_HEADER
from .Forcefield import Atom, Bond, Theta, Phi


class TemplateOPLS2005:

    def __init__(self, path_to_template):
        self.path_to_template = path_to_template
        self.template_name = ""
        self.num_nbon_params = 0
        self.num_bond_params = 0
        self.num_angle_params = 0
        self.num_dihedr_params = 0
        self.num_nonnull = 0
        self.list_of_atoms = {}
        self.list_of_bonds = {}
        self.list_of_thetas = {}
        self.list_of_phis = []
        self.list_of_iphis = []
        self.unique_atoms = []
        self.read_template()

    def read_template(self):
        template = file_to_list_of_lines(self.path_to_template)
        for line in template[2:3]:
            self.template_name = get_string_from_line(line=line,
                                                      index_initial=0,
                                                      index_final=5)
            self.num_nbon_params = int(get_string_from_line(line=line,
                                                            index_initial=6,
                                                            index_final=11))
            self.num_bond_params = int(get_string_from_line(line=line,
                                                            index_initial=13,
                                                            index_final=17))
            self.num_angle_params = int(get_string_from_line(line=line,
                                                             index_initial=18,
                                                             index_final=24))
            self.num_dihedr_params = int(get_string_from_line(line=line,
                                                              index_initial=25,
                                                              index_final=31))
            self.num_nonnull = int(get_string_from_line(line=line,
                                                        index_initial=32,
                                                        index_final=39))
        for line in template[3:]:
            if line.startswith("NBON"):
                index = template.index(line)
                break
            try:
                atom_id = get_string_from_line(line=line,
                                               index_initial=0,
                                               index_final=6)
                parent_id = get_string_from_line(line=line,
                                                 index_initial=6,
                                                 index_final=11)
                location = get_string_from_line(line=line,
                                                index_initial=12,
                                                index_final=13)
                atom_type = get_string_from_line(line=line,
                                                 index_initial=16,
                                                 index_final=20)
                pdb_atom_name = get_string_from_line(line=line,
                                                     index_initial=21,
                                                     index_final=25)
                unknown = get_string_from_line(line=line,
                                               index_initial=26,
                                               index_final=31)
                x_zmatrix = get_string_from_line(line=line,
                                                 index_initial=32,
                                                 index_final=43)
                y_zmatrix = get_string_from_line(line=line,
                                                 index_initial=44,
                                                 index_final=55)
                z_zmatrix = get_string_from_line(line=line,
                                                 index_initial=56,
                                                 index_final=67)
                atom = Atom(atom_id=atom_id, parent_id=parent_id,
                            location=location, atom_type=atom_type,
                            pdb_atom_name=pdb_atom_name, unknown=unknown,
                            x_zmatrix=x_zmatrix, y_zmatrix=y_zmatrix,
                            z_zmatrix=z_zmatrix)
                self.list_of_atoms.setdefault(atom.atom_id, atom)
                if pdb_atom_name not in self.unique_atoms:
                    self.unique_atoms.append(pdb_atom_name)
                else:
                    raise ValueError("ERROR: PDB ATOM NAME " +
                                     "{} ".format(pdb_atom_name) +
                                     "ALREADY EXISTS in the template" +
                                     " {}!".format(self.path_to_template))
            except ValueError:
                raise ValueError("Unexpected type in line " +
                                 "{}".format(template.index(line)) +
                                 " of {}".format(self.path_to_template) +
                                 "\n{}".format(line))

        for line in template[index + 1:]:
            if line.startswith("BOND"):
                index = template.index(line)
                break
            try:
                id = int(get_string_from_line(line=line,
                                              index_initial=0,
                                              index_final=6))
                self.list_of_atoms[id].sigma = float(
                    get_string_from_line(line=line,
                                         index_initial=7,
                                         index_final=14))
                self.list_of_atoms[id].epsilon = float(
                    get_string_from_line(line=line,
                                         index_initial=15,
                                         index_final=23))
                self.list_of_atoms[id].charge = float(
                    get_string_from_line(line=line,
                                         index_initial=24,
                                         index_final=34))
                self.list_of_atoms[id].radnpSGB = float(
                    get_string_from_line(line=line,
                                         index_initial=35,
                                         index_final=43))
                self.list_of_atoms[id].radnpType = float(
                    get_string_from_line(line=line,
                                         index_initial=44,
                                         index_final=52))
                self.list_of_atoms[id].sgbnpGamma = float(
                    get_string_from_line(line=line,
                                         index_initial=53,
                                         index_final=66))
                self.list_of_atoms[id].sgbnpType = float(
                    get_string_from_line(line=line,
                                         index_initial=67,
                                         index_final=80))
            except ValueError:
                raise ValueError(
                    "Unexpected type in line {}".format(template.index(line)) +
                    " of {}\n{}".format(self.path_to_template, line))

        for line in template[index + 1:]:
            if line.startswith("THET"):
                index = template.index(line)
                break
            try:
                id_atom1 = int(get_string_from_line(line=line,
                                                    index_initial=0,
                                                    index_final=6))
                id_atom2 = int(get_string_from_line(line=line,
                                                    index_initial=6,
                                                    index_final=12))
                spring = get_string_from_line(line=line,
                                              index_initial=13,
                                              index_final=21)
                eq_dist = get_string_from_line(line=line,
                                               index_initial=23,
                                               index_final=28)
                # Create bond instance
                bond = Bond(atom1=id_atom1, atom2=id_atom2,
                            spring=spring, eq_dist=eq_dist)
                self.list_of_bonds.setdefault((id_atom1, id_atom2), bond)
                # Set which atom is bonded with
                self.list_of_atoms[id_atom1].bonds.append(bond)

            except ValueError:
                raise ValueError(
                    "Unexpected type in line {}".format(template.index(line)) +
                    " of {}\n{}".format(self.path_to_template, line))
        for line in template[index + 1:]:
            if line.startswith("PHI"):
                index = template.index(line)
                break
            try:
                id_atom1 = int(get_string_from_line(line=line,
                                                    index_initial=0,
                                                    index_final=6))
                id_atom2 = int(get_string_from_line(line=line,
                                                    index_initial=6,
                                                    index_final=12))
                id_atom3 = int(get_string_from_line(line=line,
                                                    index_initial=13,
                                                    index_final=18))
                spring = get_string_from_line(line=line,
                                              index_initial=19,
                                              index_final=29)
                eq_angle = get_string_from_line(line=line,
                                                index_initial=31,
                                                index_final=40)
                # Create bond instance
                theta = Theta(atom1=id_atom1, atom2=id_atom2, atom3=id_atom3,
                              spring=spring, eq_angle=eq_angle)
                self.list_of_thetas.setdefault((id_atom1, id_atom2, id_atom3),
                                               theta)
                self.list_of_atoms[id_atom1].thetas.append(theta)

            except ValueError:
                raise ValueError(
                    "Unexpected type in line {}".format(template.index(line)) +
                    " of {}\n{}".format(self.path_to_template, line))
        for line in template[index + 1:]:
            if line.startswith("IPHI"):
                index = template.index(line)
                break
            try:
                id_atom1 = int(get_string_from_line(line=line,
                                                    index_initial=0,
                                                    index_final=5))
                id_atom2 = int(get_string_from_line(line=line,
                                                    index_initial=6,
                                                    index_final=11))
                id_atom3 = int(get_string_from_line(line=line,
                                                    index_initial=12,
                                                    index_final=17))
                id_atom4 = int(get_string_from_line(line=line,
                                                    index_initial=18,
                                                    index_final=23))
                constant = get_string_from_line(line=line,
                                                index_initial=26,
                                                index_final=32)
                preafactor = get_string_from_line(line=line,
                                                  index_initial=33,
                                                  index_final=38)
                nterm = get_string_from_line(line=line,
                                             index_initial=39,
                                             index_final=42)
                # Create bond instance
                phi = Phi(atom1=id_atom1, atom2=id_atom2, atom3=id_atom3,
                          atom4=id_atom4, constant=constant,
                          prefactor=preafactor, nterm=nterm, improper=False)
                self.list_of_phis.append(phi)
                self.list_of_atoms[id_atom1].phis.append(phi)

            except ValueError:
                raise ValueError(
                    "Unexpected type in line {}".format(template.index(line)) +
                    " of {}\n{}".format(self.path_to_template, line))
        for line in template[index + 1:]:
            if line.startswith("END"):
                break
            try:
                id_atom1 = int(get_string_from_line(line=line,
                                                    index_initial=0,
                                                    index_final=6))
                id_atom2 = int(get_string_from_line(line=line,
                                                    index_initial=7,
                                                    index_final=12))
                id_atom3 = int(get_string_from_line(line=line,
                                                    index_initial=13,
                                                    index_final=18))
                id_atom4 = int(get_string_from_line(line=line,
                                                    index_initial=19,
                                                    index_final=24))
                constant = get_string_from_line(line=line,
                                                index_initial=26,
                                                index_final=34)
                preafactor = get_string_from_line(line=line,
                                                  index_initial=34,
                                                  index_final=39)
                nterm = get_string_from_line(line=line,
                                             index_initial=40,
                                             index_final=43)
                # Create bond instance
                phi = Phi(atom1=id_atom1, atom2=id_atom2, atom3=id_atom3,
                          atom4=id_atom4, constant=constant,
                          prefactor=preafactor, nterm=nterm, improper=True)
                self.list_of_iphis.append(phi)

            except ValueError:
                raise ValueError(
                    "Unexpected type in line {}".format(template.index(line)) +
                    " of {}\n{}".format(self.path_to_template, line))

    def write_header(self):
        return HEADER_OPLS2005 + PATTERN_OPLS2005_RESX_HEADER.format(
            self.template_name, self.num_nbon_params,
            self.num_bond_params, self.num_angle_params,
            self.num_dihedr_params, self.num_nonnull)

    def write_xres(self):
        content = []
        for n in range(1, len(self.list_of_atoms) + 1):
            line = self.list_of_atoms[n].write_resx()
            content.append(line)
        return "".join(content)

    def write_nbon(self):
        content = []
        for n in range(1, len(self.list_of_atoms) + 1):
            line = self.list_of_atoms[n].write_nbon()
            content.append(line)
        return "".join(content)

    def write_bond(self):
        content = []
        for key in self.list_of_bonds.keys():
            line = self.list_of_bonds[key].write_bond()
            content.append(line)
        return "".join(content)

    def write_theta(self):
        content = []
        for key in self.list_of_thetas.keys():
            line = self.list_of_thetas[key].write_theta()
            content.append(line)
        return "".join(content)

    def write_phis(self):
        content = []
        for phi in self.list_of_phis:
            line = phi.write_phi()
            content.append(line)
        return "".join(content)

    def write_iphis(self):
        content = []
        for phi in self.list_of_iphis:
            line = phi.write_iphi()
            content.append(line)
        return "".join(content)

    def write_template(self):
        header = self.write_header()
        xres = self.write_xres()
        nbon_header = "NBON\n"
        nbon = self.write_nbon()
        bond_header = "BOND\n"
        bond = self.write_bond()
        theta_header = "THET\n"
        theta = self.write_theta()
        phis_header = "PHI\n"
        phis = self.write_phis()
        iphis_header = "IPHI\n"
        iphis = self.write_iphis()
        ending = "END"

        return header + xres + nbon_header + nbon + bond_header + bond + \
            theta_header + theta + phis_header + phis + iphis_header + \
            iphis + ending

    def write_template_to_file(self, template_new_name=None):
        if not template_new_name:
            name = self.template_name.lower() + "z"
        else:
            name = template_new_name
        with open(name, "w") as template:
            template.write(self.write_template())

    def get_list_of_fragment_atoms(self):
        atoms = []
        for key, atom in self.list_of_atoms.items():
            if atom.is_unique:
                atoms.append((key, atom))
        return atoms

    def get_list_of_fragment_bonds(self):
        bonds = []
        for key, bond in self.list_of_bonds.items():
            if bond.is_unique:
                bonds.append((key, bond))
        return bonds

    def get_list_of_fragment_thetas(self):
        thetas = []
        for key, theta in self.list_of_thetas.items():
            if theta.is_unique:
                thetas.append((key, theta))
        return thetas

    def get_list_of_fragment_phis(self):
        phis = []
        for phi in self.list_of_phis:
            if phi.is_unique:
                phis.append(phi)
        return phis

    def get_list_of_fragment_iphis(self):
        iphis = []
        for iphi in self.list_of_iphis:
            if iphi.is_unique:
                iphis.append(iphi)
        return iphis

    def get_atom_by_pdb_atom_name(self, atom_name):
        for key, atom in self.list_of_atoms.items():
            if atom.pdb_atom_name == atom_name:
                return atom
        else:
            raise NameError("PDB atom name {}".format(atom_name) +
                            "not found in template " +
                            "{}".format(self.template_name))

    def getBondBetweenAtomNamesPair(self, atom_names_pair):
        if (len(atom_names_pair) != 2):
            raise NameError("Bond {} of incompatible type".format(
                atom_names_pair))

        for (index1, index2), bond in self.list_of_bonds.items():
            name1 = self.list_of_atoms[index1].pdb_atom_name
            name2 = self.list_of_atoms[index2].pdb_atom_name

            if ((name1 in atom_names_pair) and (name2 in atom_names_pair)):
                return bond
        else:
            raise NameError("Bond {} not found in template".format(bond))


def file_to_list_of_lines(file_path):
    with open(file_path, "r") as template:
        content = template.readlines()
    return content


def get_string_from_line(line, index_initial, index_final):
    string = line[index_initial:index_final]
    return string.strip()
