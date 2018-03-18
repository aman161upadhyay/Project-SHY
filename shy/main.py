import numpy
import os
import sys

from shy.atom import Atom
from shy.protein import Protein


class ContinueLabelled(Exception):
    pass


if len(sys.argv) < 2:
    raise RuntimeError("Expected Usage: main.py <path-to-pdb-files>")
else:
    dir_path = sys.argv[1]

if not os.path.exists(dir_path):
    raise RuntimeError("Folder path does not exist")

proteins = []
for filename in os.listdir(dir_path):
    if filename.endswith(".pdb"):
        protein_name = filename[0:4]
        protein_data = open(os.path.join(dir_path, filename), 'r')
        proteins.append(Protein(protein_name, protein_data))
        protein_data.close()


count = 0
probable_hbond = []

for protein in proteins:
    print(count)
    count+=1
    print("Analyzing " + protein.name)
    se_atoms = protein.get_atoms("MSE", "SE")

    for se in se_atoms:
        connected_atoms = protein.get_connected_atoms(se)
        proximal_hydrogen = protein.get_proximal_atoms(se, "H", 2.5)

        for hydrogen in proximal_hydrogen:
            try:
                for anti_antecedent in connected_atoms:
                    bond_angle = Atom.get_bond_angles(se, hydrogen, anti_antecedent)
                    if bond_angle < numpy.pi/2:
                        raise ContinueLabelled

            except ContinueLabelled:
                continue

            proximal_oxygen = protein.get_proximal_atoms(hydrogen, "O", 1)
            proximal_nitrogen = protein.get_proximal_atoms(hydrogen, "N", 1.05)

            for oxygen in proximal_oxygen:
                bond_angle = Atom.get_bond_angles(hydrogen, se, oxygen)
                if bond_angle > numpy.pi / 2:
                    probable_hbond.append((protein.name, se, hydrogen, oxygen))

            for nitrogen in proximal_nitrogen:
                bond_angle = Atom.get_bond_angles(hydrogen, se, nitrogen)
                if bond_angle > numpy.pi / 2:
                    probable_hbond.append((protein.name, se, hydrogen, nitrogen))

    print(len(probable_hbond))
    print("--------------------")

for bonds in probable_hbond:
    print(bonds)


