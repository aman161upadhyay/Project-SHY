from atom import Atom


class Protein:
    """
    Model for a protein.
    Reads a PDB file.
    Info https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    """

    def __init__(self, name, data):
        self.name = name
        self.lines = []
        for line in data.readlines():
            if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("CONECT"):
                self.lines.append(line)

    @staticmethod
    def parse_line_atom_details(line):
        """ Get atom object from line string"""
        atom_name = line[12:16].strip()
        residue_number = int(line[7:11].strip())
        aminoacid = int(line[22:26].strip())
        b_factor = float(line[60:66].strip())
        occupancy = float(line[54:60].strip())
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        return Atom(atom_name, aminoacid, b_factor, occupancy, residue_number, x, y, z)

    def contains_amino_acid(self, amino_acid_name):
        """ Check if a protein contains an amino acid """
        for line in self.lines:
            if line[17:20].strip() == amino_acid_name:
                return True
        return False

    def get_atoms(self, amino_acid_name, atom_name):
        """ Returns list of atom objects only for SE and MSE """
        results = []
        for line in self.lines:
            if line[12:16].strip() == atom_name and line[17:20].strip() == amino_acid_name:
                atom = Protein.parse_line_atom_details(line)
                results.append(atom)
        return results

    def get_proximal_atoms(self, center_atom, neighbour_atom_type, max_neighbour_distance):
        """
        Get atoms inside a particular radii's sphere
        (exception) Atoms belonging to water residues are excluded
        """
        results = []
        for line in self.lines:
            if line[12:16].strip()[0] == neighbour_atom_type and line[17:20].strip() != "HOH":
                atom = Protein.parse_line_atom_details(line)
                if Atom.distance(center_atom, atom) <= max_neighbour_distance:
                    results.append(atom)
        return results

    def get_connected_hydrogens(self, center_atom):
        """ Return list of connected hydrogen atoms connected to center atom """
        hydrogen_identifier = center_atom.atom_name[1:]

        if hydrogen_identifier == '':
            return []

        results = []
        for line in self.lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):

                if center_atom.aminoacid == int(line[22:26].strip()):
                    line_atom_name = line[12:16].strip()
                    if line_atom_name[0] == 'H' and line_atom_name[1:].startswith(hydrogen_identifier):
                        results.append(Protein.parse_line_atom_details(line))
        return results

    def get_connected_atoms(self, center_atom):
        """ Return list of atoms connected to center atom """
        results = []
        for line in self.lines:
            if line.startswith("CONECT"):
                connections = [int(line[start_index: start_index + 5].strip()) for start_index in
                               range(6, len(line) - 1, 5)]
                if connections[0] == center_atom.residue_number:
                    results = [self.get_atom(connection) for connection in connections[1:]]
        return results

    def get_atom(self, residue_number):
        """ Get atom with a specific residue number """
        if residue_number == 0:
            raise ValueError
        for line in self.lines:
            if int(line[7:11].strip()) == residue_number:
                return Protein.parse_line_atom_details(line)
