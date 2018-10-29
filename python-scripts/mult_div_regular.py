#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 11:37:53 2018

@author: mruffini

"""

from collections import OrderedDict
import sys

import argparse
import time
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
parser.add_argument("--type", default='wregular',
                    help="Implementation of regular constraint (wregular (default), compt, sregular or sregulardp")

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


def add_wreg_dissim(vars_list, sol, sol_name, val_list=None, msim=None):
    n_vars = len(vars_list)
    fname = 'wregular_' + sol_name
    cfn['functions'][fname] = OrderedDict()
    cfn['functions'][fname]['scope'] = vars_list
    cfn['functions'][fname]['type'] = 'wregular'
    cfn['functions'][fname]['params'] = OrderedDict()

    end_state = n_vars * (args.divmin + 1)
    cfn['functions'][fname]['params']['nb_states'] = end_state + 1  # States from 0 to end
    cfn['functions'][fname]['params']['starts'] = [[0, 0.0]]
    cfn['functions'][fname]['params']['ends'] = [[end_state, 0.0]]
    cfn['functions'][fname]['params']['transitions'] = []

    # from 0 
    var = vars_list[0]
    domain = utils.get_domain(var, cfn)
    for val_index in range(len(domain)):
        val = domain[val_index]
        delta = min(utils.dissim(val, sol[0], val_list, msim), args.divmin)
        cfn['functions'][fname]['params']['transitions'].append([0, val_index, (delta * n_vars + 1), 0.0])

    for d in range(args.divmin):  # line in the diversity counting automaton
        for pos in range(1, n_vars):
            var = vars_list[pos]
            domain = utils.get_domain(var, cfn)
            for val_index in range(len(domain)):
                val = domain[val_index]
                delta = utils.dissim(val, sol[pos], val_list, msim)
                current_state = d * n_vars + pos
                next_state = (min(d + delta, args.divmin)) * n_vars + pos + 1
                cfn['functions'][fname]['params']['transitions'].append([current_state, val_index, next_state, 0.0])
    # last line:
    d = args.divmin
    for pos in range(1, n_vars):
        var = vars_list[pos]
        domain = utils.get_domain(var, cfn)
        for val_index in range(len(domain)):
            current_state = d * n_vars + pos
            next_state = current_state + 1
            cfn['functions'][fname]['params']['transitions'].append([current_state, val_index, next_state, 0.0])


def counting_dissim(vars_list, sol, sol_name, val_list=None,
                    msim=None):  # all vars have their domain included in val_list - msim and val_list = same order
    n_vars = len(vars_list)  # = len(sol)
    # Q0
    q = 'q_' + sol_name + '_'
    cfn['variables'][q + '-1'] = 1
    for i in range(n_vars):
        var = vars_list[i]
        soli = sol[i]
        qi = q + str(i)
        qip = q + str(i - 1)
        fi = 'f' + qi
        cfn['functions'][fi] = OrderedDict()
        cfn['functions'][fi]['scope'] = [qip, var, qi]
        cfn['functions'][fi]['defaultcost'] = float(cfn['problem']['mustbe'][1:])
        cfn['functions'][fi]['costs'] = []
        for qip_val in utils.get_domain(qip, cfn):
            for val in utils.get_domain(var, cfn):
                delta = utils.dissim(val, soli, val_list, msim)
                qi_val = min(qip_val + delta, args.divmin)
                cfn['functions'][fi]['costs'] += [qip_val, val, qi_val, 0]
        cfn['variables'][qi] = args.divmin + 1
    # Adding unary cost on Qn
    cfn['functions']['unary_' + qi] = OrderedDict()
    cfn['functions']['unary_' + qi]['scope'] = [qi]
    cfn['functions']['unary_' + qi]['costs'] = []
    for j in range(args.divmin):
        cfn['functions']['unary_' + qi]['costs'].append(float(cfn['problem']['mustbe'][1:]))
    cfn['functions']['unary_' + qi]['costs'].append(0)


def add_sreg_dissim(vars_list, sol, sol_name, type_string, val_list, msim):
    n_vars = len(vars_list)
    fname = type_string + '_' + sol_name
    cfn['functions'][fname] = OrderedDict()
    cfn['functions'][fname]['scope'] = vars_list
    cfn['functions'][fname]['type'] = type_string
    cfn['functions'][fname]['params'] = OrderedDict()

    cfn['functions'][fname]['params']['metric'] = 'var'  # or edit ?
    cfn['functions'][fname]['params']['cost'] = float(cfn['problem']['mustbe'][1:])
    end_state = n_vars * (args.divmin + 1)
    cfn['functions'][fname]['params']['nb_states'] = end_state + 1  # States from 0 to end
    cfn['functions'][fname]['params']['starts'] = [0]
    cfn['functions'][fname]['params']['ends'] = [end_state]
    cfn['functions'][fname]['params']['transitions'] = []

    # from 0 
    var = vars_list[0]
    domain = utils.get_domain(var, cfn)
    for val_index in range(len(domain)):
        val = domain[val_index]
        delta = min(utils.dissim(val, sol[0], val_list, msim), args.divmin)
        cfn['functions'][fname]['params']['transitions'].append([0, val_index, (delta * n_vars + 1)])

    for d in range(args.divmin):  # line in the diversity counting automaton
        for pos in range(1, n_vars):
            var = vars_list[pos]
            domain = utils.get_domain(var, cfn)
            for val_index in range(len(domain)):
                val = domain[val_index]
                delta = utils.dissim(val, sol[pos], val_list, msim)
                current_state = d * n_vars + pos
                next_state = (min(d + delta, args.divmin)) * n_vars + pos + 1
                cfn['functions'][fname]['params']['transitions'].append([current_state, val_index, next_state])
    # last line:
    d = args.divmin
    for pos in range(1, n_vars):
        var = vars_list[pos]
        domain = utils.get_domain(var, cfn)
        for val_index in range(len(domain)):
            current_state = d * n_vars + pos
            next_state = current_state + 1
            cfn['functions'][fname]['params']['transitions'].append([current_state, val_index, next_state])


## Compute k solutions
# vars_list = list of variables on which the diversity contraint applies
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
start_time=time.clock()
for k in range(args.nsols):
    utils.execute("Looking for solution " + str(k + 1) + " with toulbar2", tb2_cmd)
    (opt, sol) = utils.get_optimum(tb2log)
    sols = [sol]
    t = time.clock() - start_time
    if (k > 0 and (args.type == "wregular" or args.type == "compt")):
        sols = [sols[0][:-(k) * (n_vars + 1)]]
    sol_name = "sol" + str(k + 1)
    # print solution and add it to the solution file
    sols_file.write("Solution " + str(k + 1) + '\n')
    if (not args.cpd):
        sols_file.write(' '.join([str(val) for val in sols[0]]))
        print(' '.join([str(val) for val in sols[0]]))
        sols_file.write('\n')
    else:
        sols_file.write(' '.join([str(val) for val in sols[0]]))
        print(' '.join([str(val) for val in sols[0]]))
        sols_file.write('\n')
        utils.sols_to_cpd_sols(sols, cfn)
        sols_file.write(' '.join([str(val) for val in sols[0][:n_vars]]))
        print(' '.join([str(val) for val in sols[0][:n_vars]]))
        sols_file.write('\n')
    sols_file.write(str(opt) + '\n')
    sols_file.write(str(t) + " seconds\n")
    if (k < args.nsols - 1):
        # Add constraint corresponding to the solution
        if args.type == "wregular":
            add_wreg_dissim(vars_list, sols[0], sol_name, val_list, msim)
        elif args.type == "compt":
            counting_dissim(vars_list, sols[0], sol_name, val_list, msim)
        elif args.type == "sregular":
            add_sreg_dissim(vars_list, sols[0], sol_name, 'sregular', val_list, msim)
        elif args.type == "sregulardp":
            add_sreg_dissim(vars_list, sols[0], sol_name, 'sregulardp', val_list, msim)
        else:
            print("--regular must be wregular, compt, sregular or sregulardp")
            sys.exit()
    utils.write_cfn(cfn, cfn_tmp)

sols_file.close()
