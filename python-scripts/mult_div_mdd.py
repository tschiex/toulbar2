#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 11:37:53 2018

@author: mruffini

"""

from collections import OrderedDict
import sys
import time
import argparse

import utils

AAs = "ARNDCQEGHILKMFPSTWYV"
toulbar2 = "/home/mruffini/softs/toulbar2-diverse/build/bin/Linux/toulbar2 "

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", required=True,
                    help="The input file to process (.cfn.gz)", metavar="FILE")
parser.add_argument("-o", "--output", default=None,
                    help="The name of the output file where the solutions are written")
parser.add_argument("--divmin", default=1, type=int,
                    help="Hard constraint - minimum dissimilarity distance from each previous solution\n (divmin "
                         "differences is accepted ; =0 if no constraint)")
parser.add_argument("--nsols", default=1, type=int,
                    help="Number of diverse good solutions to compute")
parser.add_argument("--type", default='smdd',
                    help="method: emdd (one exact mdd for all solutions - default), smdd (one exact mdd per solution) "
                         "or mmdd (one relaxed mdd with all solutions + one exact mdd per solution)")
parser.add_argument("--relax", default=1, type = int,
                   help="if type == mmdd, node selection method for relaxation (1=random ; 2=largest alphap)")
parser.add_argument("--maxwidth", default=None, type=int,
                    help="Maximum layer width in relaxed mdd")

# CPD
parser.add_argument("--cpd", action="store_true", default=False,
                    help="Computational Protein Design - addition of sequence variables + possiblity of using a "
                         "similarity matrix")
parser.add_argument("--msim", default=None,
                    help="Similarity matrix - If None: Hamming distance")

args = parser.parse_args()

cfn_filename = args.input
cfn = utils.read_cfn_gzip(cfn_filename)
if cfn_filename[-7:] == '.cfn.gz':
    name = cfn_filename[:-7]
else:
    name = cfn_filename

sols_filename = args.output if args.output else \
    name + '_divmin' + str(args.divmin) + "_k" + str(args.nsols) + ".sol"

cfn_tmp = name + '_divmin' + str(args.divmin) + '_' + args.type + '.tmp.cfn'
n_vars = len(cfn['variables'])

# Sequence variables (cpd)

if args.cpd:
    # Adding sequence variables
    utils.add_seq_vars(cfn)


############################

# Functions to create a constraint from a previous solution

def mdd_dissim(vars_list, sols, name, val_list=None, msim=None, max_width=None):
    # n_vars = len(vars_list)
    fname = 'mdd_' + name
    cfn['functions'][fname] = OrderedDict()
    cfn['functions'][fname]['scope'] = []
    for var in vars_list:
        cfn['functions'][fname]['scope'].append(list(cfn["variables"].keys()).index(var))
    cfn['functions'][fname]['type'] = 'mdd'
    cfn['functions'][fname]['params'] = OrderedDict()
    if (msim != None):
        print("mdd constraint not available with similarity matrix for now")
        sys.exit()
    # Comment on choisit le coût ? Pour l'instant : on prend \top (contrainte dure) mais on peut peut régler ?
    # pénalité si la contrainte dure n'est pas satisfiable ? 
    cfn['functions'][fname]['params']['cost'] = float(cfn['problem']['mustbe'][1:]) + 1
    if args.msim == None:
        cfn['functions'][fname]['params']['matrix'] = "identity"
    else:
        cfn['functions'][fname]['params']['matrix'] = args.msim
    cfn["functions"][fname]["params"]["above"] = 1
    cfn["functions"][fname]["params"]["distance"] = args.divmin - 1
    if max_width == None:
        cfn["functions"][fname]["params"]["relax"] = 0
        cfn["functions"][fname]["params"]["max_width"] = (args.divmin + 1)**(len(sols))  # to keep exact mdd
    else:
        cfn["functions"][fname]["params"]["relax"] = args.relax
        cfn["functions"][fname]["params"]["max_width"] = max_width
    cfn["functions"][fname]["params"]["labels"] = sols


## Compute k solutions
# vars_list = list of variables on which the diversity constraint applies
# if cpd: sequence variables
if args.cpd:
    vars_list = list(cfn['variables'].keys())[n_vars:]
    val_list = AAs
else:
    vars_list = list(cfn['variables'].keys())
    # val_list : only if there is a similarity matrix
    # All variables must have the same domain ; all values in the order corresponding to the one in the similarity matrix
    val_list = list(cfn['variables'].values())[0]

if args.msim != None:
    msim = utils.read_sim_mat(args.msim)
else:
    msim = None
    val_list = None

utils.write_cfn(cfn, cfn_tmp)
sols_file = open(sols_filename, 'w')
tb2log = "tmp.tb2"
tb2_cmd = toulbar2 + cfn_tmp + ' -s -w="tmp.sol" | tee > ' + tb2log
sols = []
start_time = time.clock()
for k in range(args.nsols):
    utils.execute("Looking for solution " + str(k + 1) + " with toulbar2", tb2_cmd)
    t = time.clock() - start_time
    (opt, sol) = utils.get_optimum(tb2log)
    sols.append(sol)
    # print solution & write it in solution file
    sols_file.write("Solution " + str(k + 1) + '\n')
    if (not args.cpd):
        sols_file.write(' '.join([str(val) for val in sols[k]]))
        print(' '.join([str(val) for val in sols[k]]))
        sols_file.write('\n')
    else:
        sols_file.write(' '.join([str(val) for val in sols[k]]))
        print(' '.join([str(val) for val in sols[k]]))
        sols_file.write('\n')
        utils.sols_to_cpd_sols([sols[k]], cfn)
        sols_file.write(' '.join([str(val) for val in sols[k][:n_vars]]))
        print(' '.join([str(val) for val in sols[k][:n_vars]]))
        sols_file.write('\n')
    sols_file.write(str(opt) + '\n')
    sols_file.write(str(t) + " seconds\n")
    if args.type == "emdd":
        if args.cpd:
            mdd_dissim(vars_list,
                           [[utils.get_domain(vars_list[var], cfn).index(sol[var]) for var in range(n_vars)] for sol in
                            sols],
                           "all_sols", val_list, msim, None)
        else:
            mdd_dissim(vars_list, sols, "all_sols", val_list, msim, None)
    elif args.type == "smdd":
        if args.cpd:
            mdd_dissim(vars_list,
                           [[utils.get_domain(vars_list[var], cfn).index(sols[k][var]) for var in range(n_vars)]],
                           "sol_" + str(k + 1), val_list, msim, None)
        else:
            mdd_dissim(vars_list, [sols[k]], "sol_" + str(k + 1), val_list, msim, None)
    elif args.type == "mmdd":
        if args.cpd:
            mdd_dissim(vars_list,
                           [[utils.get_domain(vars_list[var], cfn).index(sol[var]) for var in range(n_vars)] for sol in
                            sols],
                           "all_sols", val_list, msim, args.maxwidth)
            mdd_dissim(vars_list,
                           [[utils.get_domain(vars_list[var], cfn).index(sols[k][var]) for var in range(n_vars)]],
                           "sol_" + str(k + 1), val_list, msim, None)
        else:
            mdd_dissim(vars_list, sols, "all_sols", val_list, msim, args.maxwidth)
            mdd_dissim(vars_list, [sols[k]], "sol_" + str(k + 1), val_list, msim, None)
    else:
        print("--type must be emdd, smdd or mmdd")
        sys.exit()
    utils.write_cfn(cfn, cfn_tmp)

sols_file.close()
