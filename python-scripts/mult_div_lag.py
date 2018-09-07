#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
                    help="Hard constraint - minimum dissimilarity distance from each previous solution\n (divmin differences is accepted ; =0 if no constraint)")
parser.add_argument("--nsols", default=1, type=int,
                    help="Number of diverse good solutions to compute")
parser.add_argument("--niter", default=1, type=int,
                    help="Maximum number of iterations to perform at each solution search")
parser.add_argument("--step", default="cst_stepsize",
                    help="Stepsize strategy: cst_stepsize, cst_steplength, squaresum_stepsize, nonsum_stepsize, nonsum_steplength, polyak")
parser.add_argument("--stepparam", default=0.1, type=float,
                    help="Parameter for the step size")

# CPD
parser.add_argument("--cpd", action="store_true", default=False,
                    help="Computational Protein Design - addition of sequence variables + possiblity of using a similarity matrix")
parser.add_argument("--msim", default=None,
                    help="Similarity matrix - If None: Hamming distance")

args = parser.parse_args()

cfn_filename = args.input
name = cfn_filename[:-7]

sols_filename = args.output if args.output else name + '_divmin' + str(args.divmin) + "_k" + str(args.nsols) + ".sol"

# Compute first solution

tb2_cmd = f'{toulbar2} {cfn_filename} -s -w={sols_filename}'
utils.execute("Looking for solution " + str(1), tb2_cmd)
# Init
lag_filename = 'tmp.lag'
l = ""
divmins = ""
min_divs = []


def read_newsol(lag_filename):
    with open(lag_filename, 'r') as f:
        lines = f.readlines()
        xbest = lines[-3][1:-2]
        mindiv = float(lines[-4].split(" ")[7])
    xbest = [int(xi) for xi in xbest.split(', ')]
    return (xbest, mindiv)


python_cmd = "python3 /home/mruffini/softs/toulbar2-diverse/python-scripts/divmin_lagrangian.py"

# Loop
for k in range(1, args.nsols):
    l += " 0"
    divmins += " " + str(args.divmin)
    lag_cmd = f'{python_cmd} -i {cfn_filename} -o {lag_filename} -s {sols_filename} -l{l} --divmin{divmins}' \
              f' --niter {args.niter} --step {args.step} --stepparam {args.stepparam}'
    if args.msim:
        lag_cmd += f' --msim {args.msim}'
    if args.cpd:
        lag_cmd += ' --cpd'
    utils.execute("Looking for solution " + str(k + 1), lag_cmd)
    # Read solution
    (xbest, mindiv) = read_newsol(lag_filename)
    newsol = ' '.join(str(xi) for xi in xbest)
    with open(sols_filename, 'a') as f:
        f.write(str(newsol) + '\n')
    min_divs.append(mindiv)

with open(sols_filename, 'a') as f:
    f.write("min_divs= " + str(min_divs))
    f.close()
