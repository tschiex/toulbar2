#!/usr/bin/env python3
import argparse
AA3to1 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
          'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
          'TYR': 'Y', 'VAL': 'V', 'HSE': 'H', 'HSD': 'H'}  # Rosetta does delta/epsilon adaptative protonation

parser = argparse.ArgumentParser(
    description='Translates a simple Proteus constraints.list file in a resfile.')
parser.add_argument('-i', '--input', help='Input resfile')
parser.add_argument('-o', '--output', help='Output resfile')

args = parser.parse_args()
constraints = {}
with open(args.input, 'r') as f:
    for line in f:
        pos, AA = line.split()
        constraints[pos] = AA3to1[AA]

with open(args.output, 'w') as f:
    f.write("ALLAA\n")
    f.write("start\n")
    for pos, AA in constraints.items():
        f.write(pos+" A PIKAA "+AA+'\n')
