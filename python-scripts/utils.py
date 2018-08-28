#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 16:52:50 2018

@author: mruffini
"""

import numpy as np
import gzip
import json
import re
from collections import OrderedDict
from json import encoder
import decimal
import subprocess


class DecimalEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, decimal.Decimal):
            return '{:.12f}'.format(o)
        return super(DecimalEncoder, self).default(o)


AAs = "ARNDCQEGHILKMFPSTWYV"
n_aa = len(AAs)
AA3to1 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H',
          'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
          'TYR': 'Y', 'VAL': 'V'}


########################

def read_cfn_gzip(cfn_filename):
    file = gzip.open(cfn_filename, 'r')
    first_byte = file.peek(1)
    if first_byte.decode("utf-8")[0] == '#':
        file.readline()
    cfn = json.load(file, object_pairs_hook=OrderedDict, parse_float=decimal.Decimal)
    for f_name in cfn['functions'].keys():
        for i, costi in enumerate(cfn['functions'][f_name]['costs']):
            cfn['functions'][f_name]['costs'][i] = decimal.Decimal(costi)
    return cfn


def read_cfn(cfn_filename):
    file = open(cfn_filename, 'r')
    first_byte = file.peek(1)
    if first_byte.decode("utf-8")[0] == '#':
        file.readline()
    cfn = json.load(file, object_pairs_hook=OrderedDict, parse_float=decimal.Decimal)
    for f_name in cfn['functions'].keys():
        for i, costi in cfn['functions'][f_name]['costs'].enumerate():
            cfn['functions'][f_name]['costs'][i] = decimal.Decimal(costi)
    return cfn


def read_sols(sol_filename):
    sols = []
    if sol_filename == None:
        return sols
    else:
        sol_file = open(sol_filename, 'r')
        for sol in sol_file:
            if (sol[-1] == '\n'):
                sol = sol[:-1]
            sol = list(map(int, sol.split(' ')))
            sols.append(sol)
        return sols


def sols_to_cpd_sols(sols, cfn):
    print("Sol length:" + str(len(sols[0])))
    print("nvars: " + str(len(list(cfn['variables'].values()))))
    for sol in sols:
        for pos in range(len(sol)):
            sol[pos] = list(cfn['variables'].values())[pos][sol[pos]][0]


def write_cfn(cfn, output_filename):
    cfn_str = json.dumps(cfn, indent=2, cls=DecimalEncoder)
    # cfn_bytes = cfn_str.encode('utf-8')
    with open(output_filename, 'w') as fout:
        fout.write(cfn_str)


def write_cfn_gzip(cfn, output_filename):
    cfn_str = json.dumps(cfn, indent=2, cls=DecimalEncoder)
    cfn_bytes = cfn_str.encode('utf-8')
    with gzip.GzipFile(output_filename, 'w') as fout:
        fout.write(cfn_bytes)


def read_sim_mat(msim_filename):
    msim = []
    msim_file = open(msim_filename, 'r')
    for line in msim_file:
        line = line.split(' ')
        line = [int(a) for a in line if (a != '' and a != '\n')]
        msim.append(line)
    return msim


# msim = read_sim_mat("/home/mruffini/softs/msd/positive/pars/store/blosum62")

def dissim(val1, val2, val_list=None, msim=None):
    if val1 == val2:
        return 0
    else:
        if msim == None:
            return 1
        else:
            i1 = val_list.index(val1)
            i2 = val_list.index(val2)
            return (msim[i1][i1] - msim[i1][i2])


def get_domain(var, cfn):
    if isinstance(cfn['variables'][var], int):
        domain = range(0, abs(cfn['variables'][var]))
    else:
        domain = list(cfn['variables'][var])
    return domain


def get_aa_domain(var, cfn):  # for sequence variables in cpd
    aa_domain = []

    for rot in cfn['variables'][var]:
        if rot[0] not in aa_domain:
            aa_domain.append(rot[0])
    return aa_domain


def add_seq_vars(cfn):
    """
    Adds sequence variables with domain = possible amino acids at the corresponding position
    And the necessary binary cost functions to ensure that rotamer and sequence variables correspond to the same aa
    """
    rot_vars = list(cfn['variables'].keys())
    for var in rot_vars:
        cfn['variables'][var + '_aa'] = get_aa_domain(var, cfn)
    # Extra binary cost functions to ensure that rotamer and sequence variables correspond to the same aa
    for rot_var in rot_vars:
        aa_var = rot_var + '_aa'
        fname = 'fseq' + str(rot_var) + '_' + str(aa_var)
        cfn['functions'][fname] = OrderedDict()
        cfn['functions'][fname]['scope'] = [rot_var, aa_var]
        cfn['functions'][fname]['defaultcost'] = float(cfn['problem']['mustbe'][1:])
        cfn['functions'][fname]['costs'] = []
        for rot in cfn['variables'][rot_var]:
            for aa in cfn['variables'][aa_var]:
                if rot[0] == aa[0]:
                    cfn['functions'][fname]['costs'] += [rot, aa, 0]


def execute(message, commandline):
    print(message)
    print(commandline)
    proc = subprocess.call(commandline, shell=True)


def get_optimum(tb2log):
    lines = open(tb2log).readlines()
    for i in range(len(lines)):
        if "New solution:" in lines[i]:
            seq = re.split(' ', lines[i + 1])[1:]
            for pos in range(len(seq)):
                seq[pos] = int(seq[pos])
        if "Optimum:" in lines[i]:
            s = re.split(' ', lines[i])
            return (float(s[1]), seq)
    return (None, None)




    """
def msim_to_mdissim(msim):
    if msim == None:
        return None
    dim = np.shape(msim)[0]
    mdissim = np.zeros((dim,dim))
    for i in range(dim):
        for j in range(dim):
            mdissim[i][j] = msim[i][i] - msim[i][j]
    return mdissim

def dissim_hamming(index_val1, index_val2):
    if index_val1 == index_val2:
        return 0
    else:
        return 1

def dissim_msim(index_val1, index_val2, msim):

    #We assume that all variables have the same domain, and domain values are always in the same order
    #and that same order applies to the similarity matrix msim

    delta = msim[index_val1][index_val1] - msim[index_val1][index_val2]
    return delta
    """


"""
def dissim_cpd_hamming(val1, val2):
    if val1[0] == val2[0]:
        return 0
    else:
        return 1

def dissim_cpd_msim(val1, val2, msim):
    """
# msim similarity matrix - var order cf AAs
"""
    aa_i1 = AAs.index(val1[0])
    aa_i2 = AAs.index(val2[0])
    delta = msim[aa_i1][aa_i1] - msim[aa_i1][aa_i2]
    return delta
   

#mdissim = msim_to_mdissim(read_sim_mat("/home/mruffini/softs/msd/positive/pars/store/blosum62"))
    

def hamming_dissim_matrix(n_vals):
    return(np.ones((n_vals, n_vals))- np.identity(n_vals))

#print(hamming_dissim_matrix(5))


def dissim_cpd(msim, domain):
    n_vals = len(domain)
    mdissim = np.zeros((n_vals,n_vals))
    for val1 in range(n_vals):
        for val2 in range(n_vals): # We need to go through all ordered pairs because the dissimilarity measure is not symetric, even if the similarity is
            if(domain[val1][0] != domain[val2][0]): # if same aa type, dissim remains 0
                index_val1 = AAs.index(domain[val1][0])
                index_val2 = AAs.index(domain[val2][0])
                mdissim[val1][val2] = msim[index_val1][index_val1] - msim[index_val1][index_val2]
    return(mdissim)
 """
