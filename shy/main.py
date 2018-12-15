import numpy
import os
import sys
import matplotlib.pyplot as plt
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
count_mse = 0
probable_hbond = []
b_factors = []
distances1 = []
angles1 = []
distances2 = []
angles2 = []
distances3 = []
angles3 = []
distances4 = []
angles4 = []

for protein in proteins:

    se_atoms = protein.get_atoms("MSE", "SE")
    count_mse += len(se_atoms)

    for se in se_atoms:
        connected_atoms = protein.get_connected_atoms(se)
        proximal_oxygen = protein.get_proximal_atoms(se, "O", 7.0)
        proximal_nitrogen = protein.get_proximal_atoms(se, "N", 7.0)

        for oxygen in proximal_oxygen:
            connected_hydrogen = []
            connected_hydrogen.extend(protein.get_connected_hydrogens(oxygen))
            for hydrogen in connected_hydrogen:
                try:
                    for anti_antecedent in connected_atoms:
                        bond_angle = Atom.get_bond_angles(se, hydrogen, anti_antecedent)
                        if bond_angle < numpy.pi / 2:
                            raise ContinueLabelled

                except ContinueLabelled:
                    continue

                if Atom.distance(hydrogen, oxygen) <= 1.1 and Atom.distance(se, oxygen) >= 1.0:
                    bond_angle = Atom.get_bond_angles(hydrogen, se, oxygen)
                    if bond_angle >= numpy.pi / 9:
                        distances1.append(Atom.distance(se, oxygen))
                        angles1.append(degrees(bond_angle))

                if (Atom.distance(hydrogen, oxygen) <= 1.1 and Atom.distance(se, oxygen) <= 4.0
                        and degrees(Atom.get_bond_angles(hydrogen, se, oxygen)) >= 130
                        and Atom.distance(se, oxygen) >= 1.0):
                    distances2.append(Atom.distance(se, oxygen))
                    angles2.append(degrees(Atom.get_bond_angles(hydrogen, se, oxygen)))
                    probable_hbond.append((protein.name, se.residue_number, hydrogen.residue_number,
                                           oxygen.residue_number, se.occupancy, se.b_factor,
                                           se.coordinates, hydrogen.coordinates, oxygen.coordinates))
                    b_factors.append(se.b_factor)

        for nitrogen in proximal_nitrogen:
            connected_hydrogen = []
            connected_hydrogen.extend(protein.get_connected_hydrogens(nitrogen))

            for hydrogen in connected_hydrogen:
                try:
                    for anti_antecedent in connected_atoms:
                        bond_angle = Atom.get_bond_angles(se, hydrogen, anti_antecedent)
                        if bond_angle < numpy.pi / 2:
                            raise ContinueLabelled

                except ContinueLabelled:
                    continue

                if Atom.distance(hydrogen, nitrogen) <= 1.1 and Atom.distance(se, nitrogen) >= 3.0:
                    bond_angle = Atom.get_bond_angles(hydrogen, se, nitrogen)
                    if bond_angle >= numpy.pi / 9:
                        distances3.append(Atom.distance(se, nitrogen))
                        angles3.append(degrees(bond_angle))
                if (Atom.distance(hydrogen, nitrogen) <= 1.1 and Atom.distance(se, nitrogen) <= 4.0
                        and degrees(Atom.get_bond_angles(hydrogen, se, nitrogen)) >= 130
                        and Atom.distance(se, nitrogen) >= 1.0):
                    distances4.append(Atom.distance(se, nitrogen))
                    angles4.append(degrees(Atom.get_bond_angles(hydrogen, se, nitrogen)))
                    probable_hbond.append((protein.name, se.residue_number, hydrogen.residue_number,
                                           nitrogen.residue_number, se.occupancy, se.b_factor,
                                           se.coordinates, hydrogen.coordinates, nitrogen.coordinates))
                    b_factors.append(se.b_factor)

    sys.stdout = open('1213_01.txt', 'a')
    print("file_no", count)
    count += 1
    print("Analyzing " + protein.name)
    print(len(probable_hbond))

number_hydrogen_bonds = len(probable_hbond)

colors = '#205520'
colors2 = '#202060'
area = numpy.pi * 2
area2 = numpy.pi * 2

plt.scatter(distances1, angles1, s=area, c=colors, alpha=0.2)
plt.scatter(distances2, angles2, s=area2, c=colors2, alpha=0.6)
plt.xlim(0.0, 7.0)
plt.title('Scatter_Plot_Oxygen')
plt.xlabel('distance')
plt.ylabel('bond_angle')
plt.show()

plt.scatter(distances3, angles3, s=area, c=colors, alpha=0.2)
plt.scatter(distances4, angles4, s=area2, c=colors2, alpha=0.6)
plt.xlim(0.0, 7.0)
plt.title('Scatter_Plot_Nitrogen')
plt.xlabel('distance')
plt.ylabel('bond_angle')
plt.show()

for bonds in probable_hbond:
    print(bonds)
for factors in b_factors:
    print(factors)

sum_b_factors = sum(b_factors)
print("Sum of B factors is: ", sum_b_factors)

avg_b_factors = numpy.mean(b_factors)
std_dev_b_factors = numpy.std(b_factors)
var_b_factors = numpy.var(b_factors)
print("mean of B factors is: ", avg_b_factors)
print("SD of B factors is: ", std_dev_b_factors)
print("Variance of B factors is: ", var_b_factors)

plt.hist(b_factors, bins=100)
plt.title("histogram_b_factors")
plt.show()

print("Sum of MSE is: ", count_mse)
print("ratio is:", count_mse/number_hydrogen_bonds)

sys.stdout.flush()

print("The concept of waiting bewilders me. There are always deadlines. There are always ticking clocks.")