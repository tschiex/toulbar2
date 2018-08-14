#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 12 16:53:50 2018

@author: mruffini
"""

"""
Input : CFN file (.cnf.gz) describing an minimization problem
        A file describing a set of solutions for the problem
        A list of parameters lamdba_i (one per solution)
        Similarity matrices for the variable domains - or dissimilarity matrices
        From a similarity S, dissimilarity is computed as follows: D(x,y) = S(x,x) - S(x,y)
        What we use is the dissimilarity measure. With that defintion, it remains positive \o/
        
Output : A CFN file describing the problem where extra unary costs are added
        to the values used in the given solutions
            (the resulting CFN should solve min_n (C(x) -lambda * (sum_i dist(x, sol_i))))
            (see https://www.overleaf.com/read/xdwwgmsfvkbc for details)
"""

import numpy as np
from collections import OrderedDict
import sys
from decimal import Decimal
import argparse
from utils import read_sols, write_cfn, read_sim_mat, dissim, get_domain, sols_to_cpd_sols, read_cfn

AAs = "ARNDCQEGHILKMFPSTWYV"
n_aa = len(AAs)
AA3to1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
           'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
           'TYR':'Y','VAL':'V'}

    
########################


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", required = True,
                  help="The input file to process (.cfn.gz)", metavar="FILE")
parser.add_argument("-o", "--output", default= None,
                  help="The name of the output file for the modified network (.cfn.gz)")
parser.add_argument("-s", "--sols", default = None,
                  help="The file containing the solutions already computed (1 solution per line)")
parser.add_argument("-l","--lambdas", required = True, nargs="*", type=float, 
                  help="Per solution Lagrangian multipliers")
#For now we assume that all variables have the same domain (union of all domains - and that one similarity matrix is given)
parser.add_argument("--msim", default = None,
                  help="Similarity matrix - If None: Hamming distance")
#Except in the CPD case: a similarity matrix on amino acids is given.
parser.add_argument("--cpd",action = "store_true" ,default = False,
                  help="Computational Protein Design (for sim matrix)")

args = parser.parse_args()


cfn_filename = args.input
cfn = read_cfn(cfn_filename)
if (cfn_filename[-7:] == '.cfn.gz'):
    name = cfn_filename[:-7]

sol_filename = args.sols

cfn_output = args.output if (args.output) else name + '_' + os.basename(sol_filename) + '_l' + str(args.lambda) + '.cfn.gz'

##########################

sols = read_sols(sol_filename)
if (len(sols) != len(args.lambdas)):
    print("We need one multiplier per solution.")
    sys.exit(1)

if args.cpd:
    sols_to_cpd_sols(sols,cfn)

if sols == []:
    write_cfn(cfn, cfn_output)
    sys.exit(0)

vars_todo = list(cfn['variables'].keys())

if args.msim != None:
    msim = read_sim_mat(args.msim)
else:
    msim = None

for f_name in cfn['functions'].keys():
    
    if (len(cfn['functions'][f_name]['scope']) == 1):
        var = cfn['functions'][f_name]['scope'][0]
        if isinstance(var, int):
            var = list(cfn['variables'].keys())[var]
        if (var in vars_todo) :
            vars_todo.remove(var)
            var_index = list(cfn['variables'].keys()).index(var)
            domain = get_domain(var, cfn)
            n_vals = len(domain)
            if 'defaultcost' in cfn['functions'][f_name].keys():
                costs = cfn['functions'][f_name]['costs']
                cfn['functions'][f_name]['costs'] = cfn['functions'][f_name]['defaultcost'] * n_vals
                cfn['functions'][f_name].remove('defaultcost')
                for i in range(0, len(costs),2):
                    cfn['functions'][f_name]['costs'][costs[i]] = costs[i+1]
            if args.cpd:
                for d in range(n_vals):
                    for sol_idx, sol in enumerate(sols):
                        cfn['functions'][f_name]['costs'][d] += Decimal(-args.lambdas[sol_idx] * dissim(domain[d][0], sol[var_index][0], AAs, msim))
            else:
                for d in range(n_vals):
                    for sol_idx, sol in enumerate(sols):
                        cfn['functions'][f_name]['costs'][d] += Decimal(-args.lambdas[sol_idx] * dissim(d, sol[var_index], domain, msim))  ## OR m_dissim[sol[var_index]][d] ??????????

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
                cfn['functions'][f_name]['costs'][d] += Decimal(-l * dissim(d[0], sol[var_index[0]], AAs, msim))
    else:
        for d in range(n_vals):
            for sol in sols:
                cfn['functions'][f_name]['costs'][d] += Decimal(-l * dissim(d, sol[var_index], domain, msim))  ## OR m_dissim[sol[var_index]][d] ??????????

    
write_cfn(cfn, cfn_output)



