from atom import Atom


class Protein:
    """
    Model for a protein.
    Reads a PDB file.
    Info https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    """

    def __init__(self, name, lines):
        self.name = name
        self.lines = lines

    @staticmethod
    def parse_line_atom_details(line):
        residue_number = line[7:11].strip()
        x = line[30:38].strip()
        y = line[38:46].strip()
        z = line[46:54].strip()
        return Atom(residue_number, x, y, z)

    def contains_amino_acid(self, amino_acid_name):
        """ Check if a protein contains an amino acid """
        for line in self.lines:
            if line[16:20] == amino_acid_name:
                return True
        return False

    def get_atoms(self, amino_acid_name, atom_name):
        """ Returns details of a type of atom """
        results = []
        for line in self.lines:
            if line[12:16] == atom_name and line[16:20] == amino_acid_name:
                atom = Protein.parse_line_atom_details(line)
                results.extend(atom)
        return results

    def get_atom(self, residue_number):
        if residue_number == 0:
            raise ValueError
        for line in self.lines:
            if int(line[7:11]) == residue_number:
                return Protein.parse_line_atom_details(line)

    def get_neighbouring_atoms_of_type(self, center_atom, neighbour_atom_type, max_neighbour_distance):
        """
        Get atoms' info inside a particular radii's sphere
        (except) Atoms belonging to water residues are excluded
        """
        results = []
        for line in self.lines:
            if line[12:16] == neighbour_atom_type and line[17:20] != "HOH":
                atom = Protein.parse_line_atom_details(line)
                if Atom.distance(center_atom, atom) <= max_neighbour_distance:
                    results.extend(atom)
        return results

    def get_connected_atoms(self, center_atom):
        """ Return list of atoms connected to center atom """
        results = []
        for line in self.lines:
            if line.startswith("CONECT"):
                connections = line.split()[1:]
                if connections[0] == center_atom.residue_number:
                    results = [self.get_atom(connection) for connection in connections[1:]]
        return results
