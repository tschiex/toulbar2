#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from collections import OrderedDict
import sys
import os
from decimal import Decimal
import argparse
import copy
import matplotlib.pyplot as plt

from utils import read_sols, write_cfn, read_sim_mat, dissim, get_domain, sols_to_cpd_sols, read_cfn_gzip, get_optimum

AAs = "ARNDCQEGHILKMFPSTWYV"
n_aa = len(AAs)
AA3to1 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
          'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
          'TYR': 'Y', 'VAL': 'V'}
toulbar2 = "/home/tschiex/toulbar2-diverse/build/bin/Linux/toulbar2 "

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", required=True,
                    help="The input file to process (.cfn.gz)", metavar="FILE")
parser.add_argument("-o", "--output", required=True,
                    help="Output file")
parser.add_argument("-s", "--sol", default=None,
                    help="The file containing the solution already computed")
parser.add_argument("--divmin", required=True, nargs="*", type=float,
                    help="Minimum diversity")

# For now we assume that all variables have the same domain
#  (union of all domains - and that one similarity matrix is given)
parser.add_argument("--msim", default=None,
                    help="Similarity matrix - If None: Hamming distance")
# Except in the CPD case: a similarity matrix on amino acids is given.
parser.add_argument("--cpd", action="store_true", default=False,
                    help="Computational Protein Design (for sim matrix)")

args = parser.parse_args()

############
cfn_filename = args.input
cfn = read_cfn_gzip(cfn_filename)
if cfn_filename[-7:] == '.cfn.gz':
    name = cfn_filename[:-7]
else:
    name = cfn_filename

sol_filename = args.sol

sols = read_sols(sol_filename)
nsols = len(sols)

if args.cpd:
    for i, soli in enumerate(sols):
        sols[i] = soli[:int(len(soli) / 2)]
    sols_to_cpd_sols(sols, cfn)

if sols == []:
    print("No solution given")
    sys.exit(0)

cfn_tmp_file = name + "_tmp.cfn"

output_file = open(args.output, 'w')

vars_todo = list(cfn['variables'].keys())
print('nvars:' + str(len(vars_todo)))

if args.msim != None:
    msim = read_sim_mat(args.msim)
else:
    msim = None

divmin = args.divmin


#######
def add_unary_costs(cfninit, vars_todo, sols, l):
    cfn = copy.deepcopy(cfninit)
    for f_name in cfn['functions'].keys():
        if (len(cfn['functions'][f_name]['scope']) == 1):
            var = cfn['functions'][f_name]['scope'][0]
            if isinstance(var, int):
                var = list(cfn['variables'].keys())[var]
            if (var in vars_todo):
                vars_todo.remove(var)
                var_index = list(cfn['variables'].keys()).index(var)
                domain = get_domain(var, cfn)
                n_vals = len(domain)
                if 'defaultcost' in cfn['functions'][f_name].keys():
                    costs = cfn['functions'][f_name]['costs']
                    cfn['functions'][f_name]['costs'] = cfn['functions'][f_name]['defaultcost'] * n_vals
                    cfn['functions'][f_name].remove('defaultcost')
                    for i in range(0, len(costs), 2):
                        cfn['functions'][f_name]['costs'][costs[i]] = costs[i + 1]
                if args.cpd:
                    for d in range(n_vals):
                        for sol_idx, sol in enumerate(sols):
                            cfn['functions'][f_name]['costs'][d] = cfn['functions'][f_name]['costs'][d] +\
                                          Decimal(-l * dissim(domain[d][0], sol[var_index][0], AAs, msim))
                else:
                    for d in range(n_vals):
                        for sol_idx, sol in enumerate(sols):
                            cfn['functions'][f_name]['costs'][d] = cfn['functions'][f_name]['costs'][d] + \
                                                     Decimal(-l * dissim(d, sol[var_index], domain, msim))
    while (vars_todo != []):
        var = vars_todo.pop()
        domain = get_domain(var, cfn)
        n_vals = len(domain)
        f_name = "diverse_unary_" + var
        cfn['functions'][f_name] = OrderedDict()
        cfn['functions'][f_name]['scope'] = [var]
        cfn['functions'][f_name]['costs'] = [0] * n_vals
        var_index = list(cfn['variables'].keys()).index(var)
        if args.cpd:
            for d in range(n_vals):
                for sol in sols:
                    cfn['functions'][f_name]['costs'][d] = cfn['functions'][f_name]['costs'][d] + \
                                Decimal(-l * dissim(domain[d][0], sol[var_index][0], AAs, msim))
        else:
            for d in range(n_vals):
                for sol_idx, sol in enumerate(sols):
                    cfn['functions'][f_name]['costs'][d] = cfn['functions'][f_name]['costs'][d] + \
                        Decimal(-l * dissim(d, sol[var_index], domain, msim))
    return cfn


#######

output_file.write("# " + args.input + "\t" + args.sol + "\tdivmin" + str(args.divmin) + '\n')

tb2log = "tmp.tb2"
tb2_cmd = toulbar2 + cfn_tmp_file + ' -s -w="tmp.sol" | tee > ' + tb2log

lambdas = np.linspace(0, 1, 11)

qbest = None
lbest = None
qlist = []
vars = list(cfn['variables'].keys())

for l in lambdas:
    print("l=" + str(l))
    vars_todo = copy.deepcopy(vars)
    cfn_tmp = add_unary_costs(cfn, vars_todo, sols, l)
    write_cfn(cfn_tmp, cfn_tmp_file)
    print(tb2_cmd)
    os.system(tb2_cmd)
    (ql, xl) = get_optimum(tb2log)
    print("ql -l*k: " + str(ql))
    ql = ql + l * divmin[0]
    print("ql: " + str(ql))
    if qbest == None:
        qbest = ql
        lbest = l
        xstar = xl
    else:
        if (ql >= qbest):
            qbest = ql
            lbest = l
            xstar = xl

    qlist.append(ql)

output_file.write("lambdas: " + str(lambdas) + '\n')
output_file.write("qlist: " + str(qlist) + '\n')
output_file.write("lbest: " + str(lbest) + '\n')
output_file.write(str(xstar) + '\n')
output_file.write("qbest: " + str(qbest))
output_file.close()

plt.plot(lambdas, qlist)
plt.title("q(lambda) depending on lambda, for divmin = 1, cfn " + args.input + '\n' + "qbest= " + str(qbest))
plt.xlabel("lambda")
plt.ylabel("q(lambda)")

plt.savefig(name + "_plot")
plt.close()
