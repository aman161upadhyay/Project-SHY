import numpy


class Atom:
    """ Represents location and residue number of an atom """

    def __init__(self, residue_number=None, x=None, y=None, z=None):
        self.residue_number = residue_number
        self.coordinates = numpy.array([x, y, z])

    @staticmethod
    def distance(atom1, atom2):
        dist = numpy.linalg.norm(atom1.coordinates-atom2.coordinates)
        return dist
