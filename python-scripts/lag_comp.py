#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os, sys
import matplotlib.pyplot as plt
from utils import dissim, read_cfn_gzip, read_sim_mat

python_path = "python3 /home/mruffini/softs/toulbar2-diverse/python-scripts/"
tb2 = "/home/mruffini/softs/toulbar2-diverse/build/bin/Linux/toulbar2"
AAs = "ARNDCQEGHILKMFPSTWYV"
n_aa = len(AAs)
AA3to1 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
          'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
          'TYR': 'Y', 'VAL': 'V'}

parser = argparse.ArgumentParser()

parser.add_argument("--name", required=True,
                    help="Problem name")
parser.add_argument("--niter", default=20, type=int,
                    help="Number of lagrange iterations")
parser.add_argument("--divmin", default=1, type=int,
                    help="Minimum diversity between two solutions")
parser.add_argument("--cpd", action="store_true", default=False,
                    help="Computational Protein Design")
parser.add_argument("--msim", default=None,
                    help="Similarity matrix (cpd)")
args = parser.parse_args()
name = args.name
divmin = args.divmin
if args.cpd:
    cpd_str = " --cpd "
else:
    cpd_str = ""

cfn_filename = name + ".cfn.gz"
cfn = read_cfn_gzip(cfn_filename)
sols_mdd_filename = name + "_divmin" + str(divmin) + "_nsols.sols"
sol_filename = name + ".gmec"

mult_div_cmd = python_path + "mult_div_regular.py  -i " + cfn_filename + " -o " + sols_mdd_filename + \
               " --divmin " + str(divmin) + " --nsols 2 --type mdd" + cpd_str

os.system(mult_div_cmd)

if args.msim:
    msim = read_sim_mat(args.msim)
else:
    msim = None

# Recover cstar and gmec from sols_mdd_file

with open(sols_mdd_filename, 'r') as sols_mdd:
    lines = sols_mdd.readlines()
    sol_file = open(sol_filename, 'w')
    sol = lines[1]
    sol_file.write(sol)
    sol = [int(i) for i in sol.split(" ")]
    sol_file.close()
    cstar_line = 5
    xstar_line = 4
    if (cpd_str != ""):
        cstar_line = 7
        xstar_line = 5
    cstar = float(lines[cstar_line][:-1])
    xstar = [int(xi) for xi in lines[xstar_line][:-1].split(' ')]

print("cstar " + str(cstar))
"""
# Compute qbest
ql_filename = name + "_ql.txt"
qplot_cmd = python_path + "qplot.py -i " + cfn_filename + " -o " + ql_filename + " -s " + \
            sol_filename + " --divmin 1" + cpd_str

os.system(qplot_cmd)

with open(ql_filename, 'r') as ql_file:
    lines = ql_file.readlines()
    qbest = float(lines[-1].split(' ')[1])

print("qbest " + str(qbest))
"""


#######################################
############ Supergradient ############
#######################################

def read_ql_list(output_filename):
    with open(output_filename, 'r') as f:
        lines = f.readlines()
        ql_line = lines[-1]
        xbest = lines[-3][1:-2]
        lbest = float(lines[-4].split(" ")[3])
    ql_list = [float(ql) for ql in ql_line[1:-1].split(', ')]
    xbest = [int(xi) for xi in xbest.split(', ')]
    return (xbest, ql_list, lbest)


vars = list(cfn['variables'].keys())


def step_plot(step, params, l, divmin, niter):
    for h in params:
        output_filename = name + "_" + step + "_h" + str(h) + "2.lag"
        cmd = python_path + "divmin_lagrangian.py -i " + cfn_filename + " -o " + output_filename + \
              " -s " + sol_filename + " -l " + str(l) + " --divmins " + str(divmin) + \
              " --niter " + str(niter) + " --step " + step + " --stepparam " \
              + str(h) + cpd_str
        os.system(cmd)
        (xbest, ql_list, lbest) = read_ql_list(output_filename)
        # Compute diversity measure between the first solution and xbest
        div = 0
        for var_index, v in enumerate(vars):
            if args.cpd:
                div += dissim(cfn['variables'][v][sol[var_index]][0], cfn['variables'][v][xbest[var_index]][0], AAs,
                              msim)
            else:
                div += dissim(sol[var_index], xbest[var_index], None, None)
        E = ql_list[-1] + lbest * (div - divmin)
        plt.plot(ql_list, label=f'{step} {h}\n(D,E)= ({div}, {E:.4})')
    plt.plot([0, niter], [cstar, cstar], label="cstar " + str(cstar))
    plt.legend()
    plt.title(f'{name}, {step}, divmin {divmin}')
    plt.xlabel('Number of iterations t')
    plt.ylabel("Best dual value qbest_t")
    plt.savefig(name + "_" + step + "_divmin" + str(divmin))
    plt.close()


# Constant step size
"""
step_plot("cst_stepsize", [0.05, 0.01, 0.005, 0.001], 0, divmin, args.niter)
step_plot("cst_steplength", [0.05, 0.01, 0.005, 0.001], 0, divmin, args.niter)
step_plot("squaresum_stepsize", [0.1, 1, 10], 0, divmin, args.niter)
step_plot("nonsum_stepsize", [0.1, 1, 10], 0, divmin, args.niter)
step_plot("nonsum_steplength", [0.1, 1, 10], 0, divmin, args.niter)
step_plot("polyak", [0.1, 1, 10], 0, divmin, args.niter)
"""

step_plot("squaresum_stepsize", [0.1], 0, divmin, args.niter)