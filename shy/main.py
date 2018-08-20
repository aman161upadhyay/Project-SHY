import numpy
import os
import sys
import matplotlib.pyplot as plt
from matplotlib import interactive
from math import degrees

from atom import Atom
from protein import Protein


class ContinueLabelled(Exception):
    pass


if len(sys.argv) < 2:
    raise RuntimeError("Expected Usage: main.py <folder-path-of-pdb-files>")
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

count = 1
probable_hbond = []
distances = []
angles = []
# nitrogen_distances = []
# nitrogen_angles = []

for protein in proteins:

    se_atoms = protein.get_atoms("MSE", "SE")

    for se in se_atoms:
        connected_atoms = protein.get_connected_atoms(se)
        proximal_oxygen = protein.get_proximal_atoms(se, "O", 7.0)
        proximal_nitrogen = protein.get_proximal_atoms(se, "N", 7.0)
        proximal_hydrogen = protein.get_proximal_hydrogens(se, "H", 5.0)

        for hydrogen in proximal_hydrogen:
            try:
                for anti_antecedent in connected_atoms:
                    bond_angle = Atom.get_bond_angles(se, hydrogen, anti_antecedent)
                    if bond_angle < numpy.pi / 2:
                        raise ContinueLabelled

            except ContinueLabelled:
                continue

            for oxygen in proximal_oxygen:
                if Atom.distance(hydrogen, oxygen) <= 1.1:
                    bond_angle = Atom.get_bond_angles(hydrogen, se, oxygen)
                    if bond_angle >= numpy.pi / 9:
                        distances.append(Atom.distance(se, oxygen))
                        angles.append(degrees(bond_angle))
                        probable_hbond.append((protein.name, se.residue_number, hydrogen.residue_number,
                                               oxygen.residue_number,
                                               se.coordinates, hydrogen.coordinates, oxygen.coordinates))

            for nitrogen in proximal_nitrogen:
                if Atom.distance(hydrogen, nitrogen) <= 1.1:
                    bond_angle = Atom.get_bond_angles(hydrogen, se, nitrogen)
                    if bond_angle >= numpy.pi / 9:
                        distances.append(Atom.distance(se, nitrogen))
                        angles.append(degrees(bond_angle))
                        probable_hbond.append((protein.name, se.residue_number, hydrogen.residue_number,
                                               nitrogen.residue_number,
                                               se.coordinates, hydrogen.coordinates, nitrogen.coordinates))

colors = '#205520'
area = numpy.pi * 3

# Plot

plt.scatter(distances, angles, s=area, c=colors, alpha=0.5)
plt.title('Scatter plot')
plt.xlabel('distance')
plt.ylabel('bond_angle')
plt.show()

# plt.scatter(nitrogen_distances, nitrogen_angles, s=area, c=colors, alpha=0.5)
# plt.title('Scatter plot N')
# plt.xlabel('distance')
# plt.ylabel('bond_angle')
# plt.show()

#  sys.stdout = open('output_03.txt', 'a')
# print("file_no", count)
# count += 1
# print("Analyzing " + protein.name)
# print(len(probable_hbond))

# for bonds in probable_hbond:
#   print(bonds)

# sys.stdout.flush()
