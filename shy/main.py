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
donor_n_ss = []
donor_o_ss = []
acceptor_ss = []
primary = []
close = []
closer = []
donor = []
donor_residues = []
test_acceptor_sasa = []
test_donor_sasa = []
control_acceptor_sasa = []
control_donor_sasa = []
normalized_test_donor_sasa = []
normalized_test_acceptor_sasa = []
normalized_control_acceptor_sasa = []
normalized_control_donor_sasa = []
d_a_difference = []
secondary_structure = []

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

            temp_path = 'Downloads/sasa_All/' + protein.name
            sasa_file = open(temp_path + ".rsa", 'r')
            for line in sasa_file.readlines():
                if (line.startswith("HEM") and line[8:9].strip() == se.chain and
                        int(line[9:13].strip()) == se.residue_number and
                        line[4:7].strip() == se.residue_name):
                    sasa = float(line[14:22].strip())
                    control_acceptor_sasa.append(sasa)
                if (line.startswith("RES") and line[8:9].strip() == oxygen.chain and
                        int(line[9:13].strip()) == oxygen.residue_number and
                        line[4:7].strip() == oxygen.residue_name):
                    sasa = float(line[14:22].strip())
                    control_donor_sasa.append(sasa)
            sasa_file.close()

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
                                            se.chain, se.residue_number, oxygen.residue_number))

                    res = se.residue_number
                    res1 = oxygen.residue_number
                    donor_o_ss.append((protein.get_ss(name, se.chain, res)))
                    acceptor_ss.append((protein.get_ss(name, oxygen.chain, res1)))
                    secondary_structure.append((protein.get_ss(name, se.chain, res)))
                    secondary_structure.append((protein.get_ss(name, oxygen.chain, res1)))

                    if temp_b_factor != se.b_factor:
                        b_factors.append(se.b_factor)
                    temp_b_factor = se.b_factor

                    control_b_factors.pop()

                    pdb_list.append(protein.name)

                    donor.append((oxygen.atom_name, oxygen.residue_name))
                    donor_residues.append(oxygen.residue_name)

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

                    d_a_difference.append((abs(se.residue_number - oxygen.residue_number)))

                    temp_path = 'Downloads/sasa_All/' + protein.name
                    sasa_file = open(temp_path + ".rsa", 'r')
                    for line in sasa_file.readlines():
                        if (line.startswith("HEM") and line[8:9].strip() == se.chain and
                                int(line[9:13].strip()) == se.residue_number and
                                line[4:7].strip() == se.residue_name):
                            sasa = float(line[14:22].strip())
                            test_acceptor_sasa.append(sasa)
                        if (line.startswith("RES") and line[8:9].strip() == oxygen.chain and
                                int(line[9:13].strip()) == oxygen.residue_number and
                                line[4:7].strip() == oxygen.residue_name):
                            sasa = float(line[14:22].strip())
                            test_donor_sasa.append(sasa)
                    sasa_file.close()

        # -----------------------------------------------------------------------------------------------------------------
        # SE ATOMS' PROCESSING AND NITROGEN'S START HERE
        # -----------------------------------------------------------------------------------------------------------------

        for nitrogen in proximal_nitrogen:

            temp_path = 'Downloads/sasa_All/' + protein.name
            sasa_file = open(temp_path + ".rsa", 'r')
            for line in sasa_file.readlines():
                if (line.startswith("HEM") and line[8:9].strip() == se.chain and
                        int(line[9:13].strip()) == se.residue_number and
                        line[4:7].strip() == se.residue_name):
                    sasa = float(line[14:22].strip())
                    control_acceptor_sasa.append(sasa)
                if (line.startswith("RES") and line[8:9].strip() == nitrogen.chain and
                        int(line[9:13].strip()) == nitrogen.residue_number and
                        line[4:7].strip() == nitrogen.residue_name):
                    sasa = float(line[14:22].strip())
                    control_donor_sasa.append(sasa)
            sasa_file.close()

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
                                            se.chain, se.residue_number, nitrogen.residue_number))

                    res = se.residue_number
                    res1 = nitrogen.residue_number
                    donor_n_ss.append((protein.get_ss(name, se.chain, res)))
                    acceptor_ss.append((protein.get_ss(name, nitrogen.chain, res1)))
                    secondary_structure.append((protein.get_ss(name, se.chain, res)))
                    secondary_structure.append((protein.get_ss(name, nitrogen.chain, res1)))

                    if temp_b_factor != se.b_factor:
                        b_factors.append(se.b_factor)
                    temp_b_factor = se.b_factor

                    control_b_factors.pop()

                    donor.append((nitrogen.atom_name, nitrogen.residue_name))
                    donor_residues.append(nitrogen.residue_name)

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

                    d_a_difference.append((abs(se.residue_number - nitrogen.residue_number)))

                    temp_path = 'Downloads/sasa_All/' + protein.name
                    sasa_file = open(temp_path + ".rsa", 'r')
                    for line in sasa_file.readlines():
                        if (line.startswith("HEM") and line[8:9].strip() == se.chain and
                                int(line[9:13].strip()) == se.residue_number and line[4:7].strip() == se.residue_name):
                            sasa = float(line[14:22].strip())
                            test_acceptor_sasa.append(sasa)
                        if (line.startswith("RES") and line[8:9].strip() == nitrogen.chain and
                                int(line[9:13].strip()) == nitrogen.residue_number and
                                line[4:7].strip() == nitrogen.residue_name):
                            sasa = float(line[14:22].strip())
                            test_donor_sasa.append(sasa)
                    sasa_file.close()

        temp_se = se.atom_number

    sys.stdout = open('0301.txt', 'a')
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

x = numpy.mean(b_factors)
avg_b_factors = float("{0:.5f}".format(x))
x = numpy.std(b_factors)
std_dev_b_factors = float("{0:.5f}".format(x))
x = numpy.var(b_factors)
var_b_factors = float("{0:.5f}".format(x))
print("mean of B factors is: ", avg_b_factors)
print("SD of B factors is: ", std_dev_b_factors)
print("Variance of B factors is: ", var_b_factors)

for factors in b_factors:
    x = (factors - avg_b_factors) / std_dev_b_factors
    normalized_b_factor = float("{0:.4f}".format(x))
    normalized_b_factors.append(normalized_b_factor)

x = numpy.mean(normalized_b_factors)
avg_normalized_b_factors = float("{0:.4f}".format(x))
x = numpy.std(normalized_b_factors)
std_dev_normalized_b_factors = float("{0:.4f}".format(x))
x = numpy.var(normalized_b_factors)
var_normalized_b_factors = float("{0:.4f}".format(x))
print("Mean of normalized_B factors is: ", avg_normalized_b_factors)
print("SD of normalized_B factors is: ", std_dev_normalized_b_factors)
print("Variance of normalized_B factors is: ", var_normalized_b_factors)

x = numpy.mean(control_b_factors)
avg_control_b_factors = float("{0:.4f}".format(x))
x = numpy.std(control_b_factors)
std_dev_control_b_factors = float("{0:.4f}".format(x))
x = numpy.var(control_b_factors)
var_control_b_factors = float("{0:.4f}".format(x))
print("mean of control_B factors is: ", avg_control_b_factors)
print("SD of control_B factors is: ", std_dev_control_b_factors)
print("Variance of control_B factors is: ", var_control_b_factors)

for factors in control_b_factors:
    x = (factors - avg_control_b_factors) / std_dev_control_b_factors
    normalized_control_b_factor = float("{0:.4f}".format(x))
    normalized_control_b_factors.append(normalized_control_b_factor)

x = numpy.mean(normalized_control_b_factors)
avg_normalized_control_b_factors = float("{0:.4f}".format(x))
x = numpy.std(normalized_control_b_factors)
std_dev_normalized_control_b_factors = float("{0:.4f}".format(x))
x = numpy.var(normalized_control_b_factors)
var_normalized_control_b_factors = float("{0:.4f}".format(x))
print("Mean of normalized_control_B factors is: ", avg_normalized_control_b_factors)
print("SD of normalized_control_B factors is: ", std_dev_normalized_control_b_factors)
print("Variance of normalized_control_B factors is: ", var_normalized_control_b_factors)

plt.hist(b_factors, bins=100)
plt.title("histogram_b_factors")
plt.show()

plt.hist(control_b_factors, bins=100)
plt.title("histogram_control_b_factors")
plt.show()

plt.hist(normalized_b_factors, bins=100, alpha=0.3, color='b', label='histogram_normalized_b_factors', normed=True)
plt.hist(normalized_control_b_factors, bins=150, alpha=0.2, color='g',
         label='histogram_normalized_control_b_factors', normed=True)
plt.legend(loc='upper right')
plt.show()

# -----------------------------------------------------------------------------------------------------------------
# RANDOM CODES START HERE
# -----------------------------------------------------------------------------------------------------------------

print("Sum of MSE is: ", count_mse)
print("ratio is:", count_mse / number_hydrogen_bonds)

for bi in bifurcated:
    print(bi)
print("Number of bifurcated HBs:", len(bifurcated))

for acceptors in acceptor_ss:
    print("Se profile", acceptors)
for donors in donor_o_ss:
    print("Ox profile", donors)
for donors in donor_n_ss:
    print("Ni profile", donors)
#
# temp_pdb = "temp_pdb"
# for pdb_name in pdb_list:
#     if pdb_name == temp_pdb:
#         continue
#     print(pdb_name + ".pdb.hadded.pdb")
#     temp_pdb = pdb_name
#
# letter_counts = Counter(close)
# df0 = pandas.DataFrame.from_dict(letter_counts, orient='index')
# df0.plot(kind='bar')
# plt.title("I+-4")
# plt.show()
#
# letter_counts_ = Counter(closer)
# df1 = pandas.DataFrame.from_dict(letter_counts_, orient='index')
# df1.plot(kind='bar')
# plt.title("I+-2")
# plt.show()

letter_counts_ = Counter(donor_residues)
df2 = pandas.DataFrame.from_dict(letter_counts_, orient='index')
df2.plot(kind='bar')
plt.title("Donor_Residues")
plt.show()

# letter_counts_ = Counter(acceptor_ss)
# df3 = pandas.DataFrame.from_dict(letter_counts_, orient='index')
# df3.plot(kind='bar')
# plt.title("Acceptor_SS")
# plt.show()
#
# donor_ss = donor_o_ss + donor_n_ss
#
# letter_counts_ = Counter(donor_ss)
# df4 = pandas.DataFrame.from_dict(letter_counts_, orient='index')
# df4.plot(kind='bar')
# plt.title("Donor_SS")
# plt.show()

#
# NEIGHBOURING RESIDUES
#

s = pandas.Series(d_a_difference)
values = s.value_counts()
print(values)

plt.hist(s, bins=380)
plt.title("I+-x")
plt.show()

#
# Secondary Structures to delete redundant data
#

second = []
ss = 0
sl = 0
sh = 0
hs = 0
hh = 0
hl = 0
ls = 0
lh = 0
ll = 0
flag = 0
last_sec = None
for sec in secondary_structure:
    flag += 1
    if flag % 2 == 0:
        if last_sec == "Sheet" and sec == "Sheet":
            ss += 1
            second.append("ss")
        if last_sec == "Sheet" and sec == "Helix":
            sh += 1
            second.append("sh")
        if last_sec == "Sheet" and sec == "Loop":
            sl += 1
            second.append("sl")
        if last_sec == "Helix" and sec == "Sheet":
            hs += 1
            second.append("hs")
        if last_sec == "Helix" and sec == "Helix":
            hh += 1
            second.append("hh")
        if last_sec == "Helix" and sec == "Loop":
            hl += 1
            second.append("hl")
        if last_sec == "Loop" and sec == "Helix":
            lh += 1
            second.append("lh")
        if last_sec == "Loop" and sec == "Loop":
            ll += 1
            second.append("ll")
        if last_sec == "Loop" and sec == "Sheet":
            ls += 1
            second.append("ls")

    last_sec = sec

print(ss, sh, sl, ll, ls, lh, hh, hs, hl)

letter_counts_ = Counter(second)
df9 = pandas.DataFrame.from_dict(letter_counts_, orient='index')
df9.plot(kind='bar')
plt.title("Secondary Structures Preference")
plt.show()

# -----------------------------------------------------------------------------------------------------------------
# SASA CODES START HERE
# -----------------------------------------------------------------------------------------------------------------


avg_test_donor_sasa = numpy.mean(test_donor_sasa)
std_dev_test_donor_sasa = numpy.std(test_donor_sasa)
var_test_donor_sasa = numpy.var(test_donor_sasa)
print("mean of donor_sasa is: ", avg_test_donor_sasa)
print("SD of donor_sasa is: ", std_dev_test_donor_sasa)
print("Variance of donor_sasa is: ", var_test_donor_sasa)

avg_control_donor_sasa = numpy.mean(control_donor_sasa)
std_dev_control_donor_sasa = numpy.std(control_donor_sasa)
var_control_donor_sasa = numpy.var(control_donor_sasa)
print("mean of control_donor_sasa is: ", avg_control_donor_sasa)
print("SD of control_donor_sasa is: ", std_dev_control_donor_sasa)
print("Variance of control_donor_sasa is: ", var_control_donor_sasa)

avg_test_acceptor_sasa = numpy.mean(test_acceptor_sasa)
std_dev_test_acceptor_sasa = numpy.std(test_acceptor_sasa)
var_test_acceptor_sasa = numpy.var(test_acceptor_sasa)
print("mean of acceptor_sasa is: ", avg_test_acceptor_sasa)
print("SD of acceptor_sasa is: ", std_dev_test_acceptor_sasa)
print("Variance of acceptor_sasa is: ", var_test_acceptor_sasa)

avg_control_acceptor_sasa = numpy.mean(control_acceptor_sasa)
std_dev_control_acceptor_sasa = numpy.std(control_acceptor_sasa)
var_control_acceptor_sasa = numpy.var(control_acceptor_sasa)
print("mean of control_acceptor_sasa is: ", avg_control_acceptor_sasa)
print("SD of control_acceptor_sasa is: ", std_dev_control_acceptor_sasa)
print("Variance of control_acceptor_sasa is: ", var_control_acceptor_sasa)

# -----------------------------------------------------------------------------------------------------------------
# END STARTs HERE
# -----------------------------------------------------------------------------------------------------------------

sys.stdout.flush()

print("The concept of waiting bewilders me. There are always deadlines. There are always ticking clocks.")
