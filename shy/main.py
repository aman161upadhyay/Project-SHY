import numpy
import os
import sys
import matplotlib.pyplot as plt
from math import degrees
from atom import Atom
from protein import Protein
from collections import Counter
import pandas


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
        protein_data = open(os.path.join(dir_path, filename), 'r', encoding="utf8", errors='ignore')
        proteins.append(Protein(protein_name, protein_data))
        protein_data.close()

# -----------------------------------------------------------------------------------------------------------------
# DECLARATION CODES START HERE
# -----------------------------------------------------------------------------------------------------------------


count = 1
count_mse = 0
probable_h_bond = []
b_factors = []
distances1 = []
angles1 = []
distances2 = []
angles2 = []
distances3 = []
angles3 = []
distances4 = []
angles4 = []
pdb_list = []
bifurcated = []
normalized_b_factors = []
control_b_factors = []
normalized_control_b_factors = []
i_plus_1 = []
i_plus_2 = []
i_plus_3 = []
i_plus_4 = []
i_minus_1 = []
i_minus_2 = []
i_minus_3 = []
i_minus_4 = []
acceptor_n_ss = []
acceptor_o_ss = []
donor_ss = []
primary = []
close = []
closer = []

temp_se = 999999
temp_ox = 999999
temp_ni = 999999
temp_b_factor = 99.99

# -----------------------------------------------------------------------------------------------------------------
# PROTEINS PROCESSING START HERE
# -----------------------------------------------------------------------------------------------------------------

for protein in proteins:

    se_atoms = protein.get_atoms("MSE", "SE")
    count_mse += len(se_atoms)
    name = protein.name

    primary = protein.neighbors()

    # -----------------------------------------------------------------------------------------------------------------
    # SE ATOMS' PROCESSING AND OXYGEN'S START HERE
    # -----------------------------------------------------------------------------------------------------------------

    for se in se_atoms:
        connected_atoms = protein.get_connected_atoms(se)
        proximal_oxygen = protein.get_proximal_atoms(se, "O", 7.0)
        proximal_nitrogen = protein.get_proximal_atoms(se, "N", 7.0)

        control_b_factors.append(se.b_factor)

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
                    probable_h_bond.append((protein.name, se.atom_number, hydrogen.atom_number,
                                            oxygen.atom_number, se.occupancy, se.b_factor,
                                            se.coordinates, hydrogen.coordinates, oxygen.coordinates,
                                            se.residue_number, oxygen.residue_number))
                    res = se.residue_number
                    res1 = oxygen.residue_number
                    acceptor_o_ss.append((protein.get_ss(name, res)))
                    donor_ss.append((protein.get_ss(name, res1)))

                    if temp_b_factor != se.b_factor:
                        b_factors.append(se.b_factor)
                    temp_b_factor = se.b_factor

                    pdb_list.append(protein.name)

                    if se.atom_number == temp_se:
                        bifurcated.append((protein.name, se.atom_number, oxygen.atom_number))
                    if oxygen.atom_number == temp_ox:
                        bifurcated.append((protein.name, se.atom_number, oxygen.atom_number))
                    temp_ox = oxygen.atom_number

                    try:
                        close.append(primary[se.residue_number - 1])
                    except IndexError:
                        continue
                    try:
                        close.append(primary[se.residue_number - 2])
                    except IndexError:
                        continue
                    try:
                        close.append(primary[se.residue_number - 3])
                    except IndexError:
                        continue
                    try:
                        close.append(primary[se.residue_number - 4])
                    except IndexError:
                        continue
                    try:
                        close.append(primary[se.residue_number + 1])
                    except IndexError:
                        continue
                    try:
                        close.append(primary[se.residue_number + 2])
                    except IndexError:
                        continue
                    try:
                        close.append(primary[se.residue_number + 3])
                    except IndexError:
                        continue
                    try:
                        close.append(primary[se.residue_number + 4])
                    except IndexError:
                        continue
                    try:
                        closer.append(primary[se.residue_number - 1])
                    except IndexError:
                        continue
                    try:
                        closer.append(primary[se.residue_number - 2])
                    except IndexError:
                        continue
                    try:
                        closer.append(primary[se.residue_number + 1])
                    except IndexError:
                        continue
                    try:
                        closer.append(primary[se.residue_number + 2])
                    except IndexError:
                        continue

        # -----------------------------------------------------------------------------------------------------------------
        # SE ATOMS' PROCESSING AND NITROGEN'S START HERE
        # -----------------------------------------------------------------------------------------------------------------

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
                    probable_h_bond.append((protein.name, se.atom_number, hydrogen.atom_number,
                                            nitrogen.atom_number, se.occupancy, se.b_factor,
                                            se.coordinates, hydrogen.coordinates, nitrogen.coordinates,
                                            se.residue_number, nitrogen.residue_number))

                    res = se.residue_number
                    res1 = nitrogen.residue_number
                    acceptor_n_ss.append((protein.get_ss(name, res)))
                    donor_ss.append((protein.get_ss(name, res1)))

                    if temp_b_factor != se.b_factor:
                        b_factors.append(se.b_factor)
                    temp_b_factor = se.b_factor

                    pdb_list.append(protein.name)
                    if se.atom_number == temp_se:
                        bifurcated.append((protein.name, se.atom_number, nitrogen.atom_number))
                    if nitrogen.atom_number == temp_ni:
                        bifurcated.append((protein.name, se.atom_number, nitrogen.atom_number))
                    temp_ni = nitrogen.atom_number

                    try:
                        close.append(primary[se.residue_number - 1])
                    except IndexError:
                        continue
                    try:
                        close.append(primary[se.residue_number - 2])
                    except IndexError:
                        continue
                    try:
                        close.append(primary[se.residue_number - 3])
                    except IndexError:
                        continue
                    try:
                        close.append(primary[se.residue_number - 4])
                    except IndexError:
                        continue
                    try:
                        close.append(primary[se.residue_number + 1])
                    except IndexError:
                        continue
                    try:
                        close.append(primary[se.residue_number + 2])
                    except IndexError:
                        continue
                    try:
                        close.append(primary[se.residue_number + 3])
                    except IndexError:
                        continue
                    try:
                        close.append(primary[se.residue_number + 4])
                    except IndexError:
                        continue
                    try:
                        closer.append(primary[se.residue_number - 1])
                    except IndexError:
                        continue
                    try:
                        closer.append(primary[se.residue_number - 2])
                    except IndexError:
                        continue
                    try:
                        closer.append(primary[se.residue_number + 1])
                    except IndexError:
                        continue
                    try:
                        closer.append(primary[se.residue_number + 2])
                    except IndexError:
                        continue

        temp_se = se.atom_number

    sys.stdout = open('0114_03.txt', 'a')
    print("file_no", count)
    count += 1
    print("Analyzing " + protein.name)
    print(len(probable_h_bond))

# -----------------------------------------------------------------------------------------------------------------
# SCATTER PLOTS START HERE
# -----------------------------------------------------------------------------------------------------------------

number_hydrogen_bonds = len(probable_h_bond)

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

for bonds in probable_h_bond:
    print(bonds)

# -----------------------------------------------------------------------------------------------------------------
# B_FACTORS CODES START HERE
# -----------------------------------------------------------------------------------------------------------------

print("Number of unique residues is: ", len(b_factors))

avg_b_factors = numpy.mean(b_factors)
std_dev_b_factors = numpy.std(b_factors)
var_b_factors = numpy.var(b_factors)
print("mean of B factors is: ", avg_b_factors)
print("SD of B factors is: ", std_dev_b_factors)
print("Variance of B factors is: ", var_b_factors)

for factors in b_factors:
    normalized_b_factor = (factors - avg_b_factors) / std_dev_b_factors
    normalized_b_factors.append(normalized_b_factor)

for factors in normalized_b_factors:
    print(factors)

avg_normalized_b_factors = numpy.mean(normalized_b_factors)
std_dev_normalized_b_factors = numpy.std(normalized_b_factors)
var_normalized_b_factors = numpy.var(normalized_b_factors)
print("Mean of normalized_B factors is: ", avg_normalized_b_factors)
print("SD of normalized_B factors is: ", std_dev_normalized_b_factors)
print("Variance of normalized_B factors is: ", var_normalized_b_factors)

avg_control_b_factors = numpy.mean(control_b_factors)
std_dev_control_b_factors = numpy.std(control_b_factors)
var_control_b_factors = numpy.var(control_b_factors)
print("mean of control_B factors is: ", avg_control_b_factors)
print("SD of control_B factors is: ", std_dev_control_b_factors)
print("Variance of control_B factors is: ", var_control_b_factors)

for factors in control_b_factors:
    normalized_control_b_factor = (factors - avg_control_b_factors) / std_dev_control_b_factors
    normalized_control_b_factors.append(normalized_control_b_factor)

avg_normalized_control_b_factors = numpy.mean(normalized_control_b_factors)
std_dev_normalized_control_b_factors = numpy.std(normalized_control_b_factors)
var_normalized_control_b_factors = numpy.var(normalized_control_b_factors)
print("Mean of normalized_control_B factors is: ", avg_normalized_control_b_factors)
print("SD of normalized_control_B factors is: ", std_dev_normalized_control_b_factors)
print("Variance of normalized_control_B factors is: ", var_normalized_control_b_factors)

plt.hist(b_factors, bins=100)
plt.title("histogram_b_factors")
plt.show()

plt.hist(normalized_b_factors, bins=100)
plt.title("histogram_normalized_b_factors")
plt.show()

plt.hist(control_b_factors, bins=100)
plt.title("histogram_control_b_factors")
plt.show()

plt.hist(normalized_control_b_factors, bins=100)
plt.title("histogram_normalized_control_b_factors")
plt.show()

# -----------------------------------------------------------------------------------------------------------------
# RANDOM CODES START HERE
# -----------------------------------------------------------------------------------------------------------------

print("Sum of MSE is: ", count_mse)
print("ratio is:", count_mse / number_hydrogen_bonds)

for bi in bifurcated:
    print(bi)
print("Number of bifurcated HBs:", len(bifurcated))

for donors in donor_ss:
    print("Se profile", donors)
for acceptors in acceptor_o_ss:
    print("Ox profile", acceptors)
for acceptors in acceptor_n_ss:
    print("Ni profile", acceptors)

temp_pdb = "temp_pdb"
for pdb_name in pdb_list:
    if pdb_name == temp_pdb:
        continue
    print(pdb_name + ".pdb.hadded.pdb")
    temp_pdb = pdb_name

letter_counts = Counter(close)
df = pandas.DataFrame.from_dict(letter_counts, orient='index')
df.plot(kind='bar')
plt.title("I+-4")
plt.show()

letter_counts_ = Counter(closer)
df1 = pandas.DataFrame.from_dict(letter_counts_, orient='index')
df1.plot(kind='bar')
plt.title("I+-2")
plt.show()

sys.stdout.flush()

print("The concept of waiting bewilders me. There are always deadlines. There are always ticking clocks.")
