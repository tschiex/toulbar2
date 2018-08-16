#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 17:36:49 2018

@author: mruffini
"""

"""
Input: CFN file (.cnf.gz) describing an minimization problem
        A file describing a set of solutions for the problem
        Minimum distance from each solution
Output: Next best solution that respects the diversity constraints

Method: To do so, we use the lagrangian relaxation and solve the dual problem. It provides a lower bound
on the optimal solution.

!!! Duality gap ??

"""


import numpy as np
from collections import OrderedDict
import sys
from decimal import Decimal
import argparse
from utils import execute, read_sols, write_cfn, read_sim_mat, dissim, get_domain, sols_to_cpd_sols, read_cfn_gzip, get_optimum

AAs = "ARNDCQEGHILKMFPSTWYV"
n_aa = len(AAs)
AA3to1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
           'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
           'TYR':'Y','VAL':'V'}
toulbar2 = " ../build/bin/Linux/toulbar2 "


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", required = True,
                  help="The input file to process (.cfn.gz)", metavar="FILE")
parser.add_argument("-o", "--output", required = True,
                  help="Output file - for each iteration: lambda values + solution + cost")
parser.add_argument("-s", "--sols", default = None,
                  help="The file containing the solutions already computed (1 solution per line)")
parser.add_argument("-l","--lambdas", required = True, nargs="*", type=float, 
                  help="Per solution initial Lagrangian multipliers")
parser.add_argument("--divmins", required = True, nargs="*", type=float, 
                  help="Per solution minimum diversity")
parser.add_argument("--niter", default = 1, type = int,
                    help="Number of iterations")
#For now we assume that all variables have the same domain (union of all domains - and that one similarity matrix is given)
parser.add_argument("--msim", default = None,
                  help="Similarity matrix - If None: Hamming distance")
#Except in the CPD case: a similarity matrix on amino acids is given.
parser.add_argument("--cpd",action = "store_true" ,default = False,
                  help="Computational Protein Design (for sim matrix)")

args = parser.parse_args()



############
cfn_filename = args.input
cfn = read_cfn_gzip(cfn_filename)
if (cfn_filename[-7:] == '.cfn.gz'):
    name = cfn_filename[:-7]

sol_filename = args.sols

sols = read_sols(sol_filename)
nsols = len(sols)
if (len(sols) != len(args.lambdas)):
    print("We need one multiplier per solution.")
    sys.exit(1)
    
if args.cpd:
    sols_to_cpd_sols(sols,cfn)
    
    
if sols == []:
    print("No solution given")
    sys.exit(0)
    
cfn_tmp_file = name + "_tmp.cfn"

output_file = open(args.output, 'w')

vars_todo = list(cfn['variables'].keys())

if args.msim != None:
    msim = read_sim_mat(args.msim)
else:
    msim = None
    
lambdas = np.array(args.lambdas)
divmins = np.array(args.divmins)

#######

def add_unary_costs(cfn,vars_todo, sols, lambdas):
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
                    cfn['functions'][f_name]['costs'][d] += Decimal(-args.lambdas[sol_idx] * dissim(d[0], sol[var_index[0]], AAs, msim))
        else:
            for d in range(n_vals):
                for sol_idx, sol in enumerate(sols):
                    cfn['functions'][f_name]['costs'][d] += Decimal(-args.lambdas[sol_idx] * dissim(d, sol[var_index], domain, msim))  ## OR m_dissim[sol[var_index]][d] ??????????
    return cfn

output_file.write("# " + args.input + "\t" + args.sols + "\tdivmin" + str(args.divmins) + '\n')

tb2log="tmp.tb2"
tb2_cmd = toulbar2 + cfn_tmp_file + ' -w="tmp.sol" | tee > ' + tb2log
####### TODO
dual_lb = None
dual_best = None
b=1
#######

for t in range(args.niter):
    output_file.write("# Iteration " + str(t+1) + " with lambdas="+ str(lambdas) + '\n')
    cfn_tmp = add_unary_costs(cfn, ["Q1","Q2","Q3","Q4"], sols, lambdas)
    write_cfn(cfn_tmp, cfn_tmp_file)
    execute("Iteration " + str(t+1), tb2_cmd)
    (dual_cost, sol)=get_optimum(tb2log)
    output_file.write(str(sol) + '\n')
    dual_cost += np.sum(lambdas*divmins)
    if dual_best ==None:
        dual_best=dual_cost
    else:
        dual_best = max(dual_best, dual_cost)
    output_file.write(dual_cost + '\n')
    alpha = (1+nsols)/(t+nsols)
    subgrad =  divmins
    for j, solj in enumerate(sols):
        for var in vars_todo:
            var_index = list(cfn['variables'].keys()).index(var)
            domain = get_domain(var, cfn)
            subgrad[j] -= dissim(sol[var_index], solj[var_index], domain, msim)
    dual_approx = (1+b)*dual_best
    if(dual_best == dual_cost):
        b *= 1.5
    else:
        b /= 1.5        
    stepsize = alpha*(dual_approx - dual_cost)/np.sum(subgrad**2)
    for j in range(len(lambdas)):
        lambdas[j] = max(0, lambdas[j]+stepsize*subgrad[j])
    



