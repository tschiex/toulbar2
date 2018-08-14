#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 11:37:53 2018

@author: mruffini
"""

import numpy as np
import gzip
import json
from collections import OrderedDict
import sys
import argparse
from utils import read_sols, write_cfn, read_sim_mat, get_domain, add_seq_vars, sols_to_cpd_sols, dissim, read_cfn


AAs = "ARNDCQEGHILKMFPSTWYV"

############################


def add_wreg_dissim(vars_list, sol, sol_name, val_list=None, msim=None):
    n_vars = len(vars_list)
    fname = 'wregular_' + sol_name
    cfn['functions'][fname] = OrderedDict()
    cfn['functions'][fname]['scope'] = vars_list
    cfn['functions'][fname]['type'] = 'wregular'
    cfn['functions'][fname]['params'] = OrderedDict()

    end_state = n_vars * (args.divmin + 1)
    cfn['functions'][fname]['params']['nb_state'] = end_state + \
        1  # States from 0 to end
    cfn['functions'][fname]['params']['starts'] = [[0, 0.0]]
    cfn['functions'][fname]['params']['ends'] = [[end_state, 0.0]]
    cfn['functions'][fname]['params']['transitions'] = []

    # from 0
    var = vars_list[0]
    domain = get_domain(var, cfn)
    for val_index in range(len(domain)):
        val = domain[val_index]
        delta = min(dissim(val, sol[0], val_list, msim), args.divmin)
        cfn['functions'][fname]['params']['transitions'].append(
            [0, val_index, (delta * n_vars + 1), 0.0])

    for d in range(args.divmin):  # line in the diversity counting automaton
        for pos in range(1, n_vars):
            var = vars_list[pos]
            domain = get_domain(var, cfn)
            for val_index in range(len(domain)):
                val = domain[val_index]
                delta = dissim(val, sol[pos], val_list, msim)
                current_state = d*n_vars + pos
                next_state = (min(d+delta, args.divmin)) * n_vars + pos + 1
                cfn['functions'][fname]['params']['transitions'].append(
                    [current_state, val_index, next_state, 0.0])
    # last line:
    d = args.divmin
    for pos in range(1, n_vars):
        var = vars_list[pos]
        domain = get_domain(var, cfn)
        for val_index in range(len(domain)):
            current_state = d*n_vars + pos
            next_state = current_state + 1
            cfn['functions'][fname]['params']['transitions'].append(
                [current_state, val_index, next_state, 0.0])


def add_sreg_dissim(vars_list, sol, sol_name, type_string, val_list, msim):
    n_vars = len(vars_list)
    fname = type_string + '_' + sol_name
    cfn['functions'][fname] = OrderedDict()
    cfn['functions'][fname]['scope'] = vars_list
    cfn['functions'][fname]['type'] = type_string
    cfn['functions'][fname]['params'] = OrderedDict()

    cfn['functions'][fname]['params']['metric'] = 'var'  # or edit ?
    cfn['functions'][fname]['params']['cost'] = float(
        cfn['problem']['mustbe'][1:])
    end_state = n_vars * (args.divmin + 1)
    cfn['functions'][fname]['params']['nb_states'] = end_state + \
        1  # States from 0 to end
    cfn['functions'][fname]['params']['starts'] = [0]
    cfn['functions'][fname]['params']['ends'] = [end_state]
    cfn['functions'][fname]['params']['transitions'] = []

    # from 0
    var = vars_list[0]
    domain = get_domain(var, cfn)
    for val_index in range(len(domain)):
        val = domain[val_index]
        delta = min(dissim(val, sol[0], val_list, msim), args.divmin)
        cfn['functions'][fname]['params']['transitions'].append(
            [0, val_index, (delta * n_vars + 1)])

    for d in range(args.divmin):  # line in the diversity counting automaton
        for pos in range(1, n_vars):
            var = vars_list[pos]
            domain = get_domain(var, cfn)
            for val_index in range(len(domain)):
                val = domain[val_index]
                delta = dissim(val, sol[pos], val_list, msim)
                current_state = d*n_vars + pos
                next_state = (min(d+delta, args.divmin)) * n_vars + pos + 1
                cfn['functions'][fname]['params']['transitions'].append(
                    [current_state, val_index, next_state])
    # last line:
    d = args.divmin
    for pos in range(1, n_vars):
        var = vars_list[pos]
        domain = get_domain(var, cfn)
        for val_index in range(len(domain)):
            current_state = d*n_vars + pos
            next_state = current_state + 1
            cfn['functions'][fname]['params']['transitions'].append(
                [current_state, val_index, next_state])


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", required=True,
                    help="The input file to process (.cfn.gz)", metavar="FILE")
parser.add_argument("-o", "--output", default=None,
                    help="The name of the output file for the modified network (.cfn.gz)")
parser.add_argument("-s", "--sols", default=None,
                    help="The file containing the solutions already computed (1 solution per line)")
# For now we assume that all variables have the same domain (union of all domains - and that one similarity matrix is given)
parser.add_argument("--msim", default=None,
                    help="Similarity matrix - If None: Hamming distance")
# Except in the CPD case: a similarity matrix on amino acids is given.
parser.add_argument("--cpd", action="store_true", default=False,
                    help="Computational Protein Design (for sim matrix)")

parser.add_argument("--divmin", default=1, type=int,
                    help="Hard constraint - minimum dissimilarity distance from each solution\n (=0 if no constraint)")

parser.add_argument("--regular", default='wregular',
                    help="Implementation of regular constraint to use (wregular (default), sregular or sregulardp)")

args = parser.parse_args()

cfn_filename = args.input
cfn = read_cfn(cfn_filename)
if (cfn_filename[-7:] == '.cfn.gz'):
    name = cfn_filename[:-7]

sol_filename = args.sols

cfn_output = args.output if (args.output) else name + '_' + str(sol_filename.split(
    '/')[-1]) + '_divmin' + str(args.divmin) + '_' + args.regular + '.cfn.gz'
n_vars = len(cfn['variables'])

##########################

sols = read_sols(sol_filename)
if args.cpd:
    sols_to_cpd_sols(sols, cfn)

if sols == []:
    write_cfn(cfn, cfn_output)
    sys.exit(0)

#############################
###       end_idem       ####
#############################
# Sequence variables (cpd)

if args.cpd:
    # Adding sequence variables
    add_seq_vars(cfn)
    vars_list = list(cfn['variables'].keys())[n_vars:]
    val_list = AAs
else:
    vars_list = list(cfn['variables'].keys())
    val_list = list(cfn['variables'].values())[0]

if args.msim != None:
    msim = read_sim_mat(args.msim)
else:
    msim = None
    val_list = None

if args.regular == "wregular":
    for sol_index in range(len(sols)):
        sol = sols[sol_index]
        sol_name = 'sol'+str(sol_index)
        add_wreg_dissim(vars_list, sol, sol_name, val_list, msim)
elif args.regular == "sregular":
    for sol_index in range(len(sols)):
        sol = sols[sol_index]
        sol_name = 'sol'+str(sol_index)
        add_sreg_dissim(vars_list, sol, sol_name, 'sregular', val_list, msim)
elif args.regular == "sregulardp":
    for sol_index in range(len(sols)):
        sol = sols[sol_index]
        sol_name = 'sol'+str(sol_index)
        add_sreg_dissim(vars_list, sol, sol_name, 'sregulardp', val_list, msim)
else:
    print("--regular must be wregular, sregular or sregulardp")
    sys.exit()

write_cfn(cfn, cfn_output)
