#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os, sys
import matplotlib.pyplot as plt

python_path = "python3 /home/mruffini/softs/toulbar2-diverse/python-scripts/"
tb2 = "/home/mruffini/softs/toulbar2-diverse/build/bin/Linux/toulbar2"


parser = argparse.ArgumentParser()

parser.add_argument("--name", required=True,
                    help="Problem name")
parser.add_argument("--niter", default=20, type=int,
                    help="Number of lagrange iterations")
parser.add_argument("--cpd", action="store_true", default=False,
                    help="Computational Protein Design")
args = parser.parse_args()
name = args.name
if args.cpd:
    cpd_str = " --cpd "
else:
    cpd_str = ""

cfn_filename = name + ".cfn.gz"
sols_mdd_filename = name + "_divmin10_nsols2.sols"
sol_filename = name + ".gmec"

mult_div_cmd = python_path + "mult_div_regular.py  -i " + cfn_filename + " -o " + sols_mdd_filename + \
               " --divmin 10 --nsols 2 --type mdd" + cpd_str

os.system(mult_div_cmd)

# Recover cstar and gmec from sols_mdd_file

with open(sols_mdd_filename, 'r') as sols_mdd:
    lines = sols_mdd.readlines()
    sol_file = open(sol_filename, 'w')
    sol_file.write(lines[1])
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
    ql_list = [float(ql) for ql in ql_line[1:-1].split(', ')]
    xbest = [int(xi) for xi in xbest.split(', ')]
    return (xbest, ql_list)


def step_plot(step, params, l, divmin, niter):
    for h in params:
        output_filename = name + "_" + step + "_h" + str(h) + ".lag"
        cmd = python_path + "divmin_lagrangian.py -i " + cfn_filename + " -o " + output_filename + \
              " -s " + sol_filename + " -l " + str(l) + " --divmins " + str(divmin) + \
              " --niter " + str(niter) + " --step " + step + " --stepparam " \
              + str(h) + cpd_str
        os.system(cmd)
        (xbest, ql_list) = read_ql_list(output_filename)
        dist = 0
        for i, xbesti in enumerate(xbest):
            if xbesti != xstar[i]:
                dist += 1
        plt.plot(ql_list, label=step + " " + str(h) + '\n' + "d(xbest, xstar):" + str(dist))
    plt.plot([0, niter], [cstar, cstar], label="cstar")
    plt.legend()
    plt.title(f'{name}, {step}')
    plt.xlabel('Number of iterations t')
    plt.ylabel("Best dual value qbest_t")
    plt.savefig(name + "_" + step + "_divmin" + str(divmin))
    plt.close()


# Constant step size

step_plot("cst_stepsize", [0.05, 0.01, 0.005, 0.001], 0, 10, args.niter)
step_plot("cst_steplength", [0.05, 0.01, 0.005, 0.001], 0, 10, args.niter)
step_plot("squaresum_stepsize", [0.1, 1, 10], 0, 10, args.niter)
step_plot("nonsum_stepsize", [0.1, 1, 10], 0, 10, args.niter)
step_plot("nonsum_steplength", [0.1, 1, 10], 0, 10, args.niter)
step_plot("polyak", [0.1, 1, 10], 0, 10, args.niter)
