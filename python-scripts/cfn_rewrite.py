#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from utils import read_cfn_gzip, write_cfn_gzip

cfn = read_cfn_gzip("water.cfn.gz")

for var in cfn["variables"].keys():
    cfn["variables"][var] = ["v" + val for val in cfn["variables"][var]]

write_cfn_gzip(cfn, "water.cfn.gz")