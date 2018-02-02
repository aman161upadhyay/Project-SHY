import os
import sys

from shy.protein import Protein

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

for protein in proteins:
    if not protein.contains_amino_acid("MSE"):
        print(protein.name+".pdb.hadded.pdb")
