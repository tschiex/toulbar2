#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 16:52:23 2018

@author: mruffini
"""
#############################
### idem divsum_unary.py ####
#############################
# uses inhomogeneous couting automata decomposed in ternary constraints

import numpy as np
import gzip
import json
from collections import OrderedDict
import sys
import argparse
from utils import read_sols, write_cfn, read_sim_mat, get_domain, add_seq_vars, sols_to_cpd_sols, dissim, read_cfn

AAs = "ARNDCQEGHILKMFPSTWYV"

# Dissimilarity counting variables
"""
!!! Here, we need the (dis)similarity values between pairs of values to be integers
"""
# all vars have their domain included in val_list - msim and val_list = same order


def counting_dissim(vars_list, sol, sol_name, val_list=None, msim=None):
    n_vars = len(vars_list)  # = len(sol)
    # Q0
    q = 'q_' + sol_name + '_'
    cfn['variables'][q + '-1'] = 1
    for i in range(n_vars):
        var = vars_list[i]
        soli = sol[i]
        qi = q+str(i)
        qip = q + str(i-1)
        fi = 'f' + qi
        cfn['functions'][fi] = OrderedDict()
        cfn['functions'][fi]['scope'] = [qip, var, qi]
        cfn['functions'][fi]['defaultcost'] = float(
            cfn['problem']['mustbe'][1:])
        cfn['functions'][fi]['costs'] = []
        for qip_val in get_domain(qip, cfn):
            for val in get_domain(var, cfn):
                delta = dissim(val, soli, val_list, msim)
                qi_val = min(qip_val + delta, args.divmin)
                cfn['functions'][fi]['costs'] += [qip_val, val, qi_val, 0]
        cfn['variables'][qi] = args.divmin+1
    # Adding unary cost on Qn
    cfn['functions']['unary_' + qi] = OrderedDict()
    cfn['functions']['unary_' + qi]['scope'] = [qi]
    cfn['functions']['unary_' + qi]['costs'] = []
    for j in range(args.divmin):
        cfn['functions']['unary_' +
                         qi]['costs'].append(float(cfn['problem']['mustbe'][1:]))
    cfn['functions']['unary_' + qi]['costs'].append(0)


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

args = parser.parse_args()

cfn_filename = args.input
cfn = read_cfn(cfn_filename)
if (cfn_filename[-7:] == '.cfn.gz'):
    name = cfn_filename[:-7]

sol_filename = args.sols

cfn_output = args.output if (args.output) else name + '_' + str(
    sol_filename.split('/')[-1]) + '_divmin' + str(args.divmin) + '.cfn.gz'

n_vars = len(cfn['variables'])

##########################

sols = read_sols(sol_filename)
if args.cpd:
    sols_to_cpd_sols(sols, cfn)

if sols == []:
    write_cfn(cfn, cfn_output)
    sys.exit(0)

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

for sol_index in range(len(sols)):
    sol = sols[sol_index]
    sol_name = 'sol'+str(sol_index)
    counting_dissim(vars_list, sol, sol_name, val_list, msim)

write_cfn(cfn, cfn_output)
