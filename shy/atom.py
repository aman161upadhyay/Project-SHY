from numpy import array, dot, arccos
from numpy.linalg import norm


class Atom:
    """ Represents location and residue number of an atom """

    def __init__(self, residue_number=None, x=None, y=None, z=None):
        self.residue_number = residue_number
        self.coordinates = array([x, y, z])

    @staticmethod
    def distance(atom1, atom2):
        """ Distance between two atoms """
        dist = norm(atom1.coordinates - atom2.coordinates)
        return dist

    @staticmethod
    def get_bond_angles(center_atom, neighbour_atom1, neighbour_atom2):
        """ Get the angle between the two bonds of center to neighbours """
        bond1 = neighbour_atom1.coordinates - center_atom.coordinates
        bond2 = neighbour_atom2.coordinates - center_atom.coordinates
        angle = arccos(dot(bond1, bond2) / (norm(bond1) * norm(bond2)))
        return angle
