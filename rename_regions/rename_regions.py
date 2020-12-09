#!/usr/bin/env python

""""
Renames gbk files

Only for antiSMASH 5+ output, which uses the string "region"

Input: inputfolder, new string
Output: (none)

"""

import os
import sys
from pathlib import Path
import argparse
from shutil import move

__author__ = "Jorge Navarro"
__version__ = "1"
__maintainer__ = "Jorge Navarro"
__email__ = "j.navarro@wi.knaw.nl"

def parameter_parser():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--inputfolder", type=Path, \
        help="Path to fungiSMASH results folder", required=True)
    parser.add_argument("-s", "--string", type=str, required=True, 
        help="All found .gbk files will be renamed '[string].region[XXX].gbk'",)
    
    return parser.parse_args()


if __name__ == "__main__":
    parameters = parameter_parser()
    
    i = parameters.inputfolder
    if not i.is_dir():
        sys.exit("Error: {} is not a valid folder".format(i))

    s = parameters.string
    bad_chars= set("/|\:<>*")
    if len(set(s) & bad_chars) > 0:
        sys.exit("Error: illegal characters used with parameter --string")
    
    gbks = list(sorted(i.glob("*region*.gbk")))
    num_gbks = len(gbks)
    
    padding_space = 0
    if num_gbks < 10:
        padding_space = 1
    elif num_gbks < 100:
        padding_space = 2
    else:
        padding_space = 3 # no genome has more than 999 BGCs
        
    for num, gbk in enumerate(gbks):
        p = gbk.parent
        old = gbk.name
        new = "{0}.region{1:0{2}d}.gbk".format(s, num+1, padding_space)
        move(p/old, p/new)
    print("Renamed {} files".format(num_gbks))
        
