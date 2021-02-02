#!/usr/bin/env python

""""
Copy GenBank files, filtering by (annotated) taxonomic lineage

Only for antiSMASH 5+ output.

Requires: biopython, Python 3

Input: GenBank file(s) (e.g. bundle from MIBiG)
Output: files that pass specified taxonomic filter are copied to outputfolder

"""

import os
import sys
from pathlib import Path
import argparse
from Bio import SeqIO
from shutil import copyfile

__author__ = "Jorge Navarro"
__version__ = "1"
__maintainer__ = "Jorge Navarro"
__email__ = "j.navarro@wi.knaw.nl"


def parameter_parser():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=Path,
        help="Input. Either a GenBank file or a folder with fungiSMASH results", 
        required=True)
    parser.add_argument("-o", "--outputfolder", type=Path,
        help="Output folder where GenBank file(s) will be copied to",
        required=True)
    parser.add_argument("-t", "--taxonomy_filter", type=str, help="Filter \
        extraction using taxonomy within the GenBank file(s) (if annotated)",
        required=True)
    parser.add_argument("--include", default=["region", "BGC"], type=str, 
        help="Optional list of strings used to filter GenBank files to be \
        analyzed. Default: 'region', 'BGC'")
    
    return parser.parse_args()



if __name__ == "__main__":
    parameters = parameter_parser()
    
    i = parameters.input
    o = parameters.outputfolder
    tax_filter = parameters.taxonomy_filter
    
    # Validate parameters
    if not (i.is_dir() or i.is_file()):
        sys.exit("Error: {} is not a valid folder or file".format(i))

    
    # Get a list of all GenBank files
    if i.is_dir():
        gbks = list()
        for gbk in sorted(i.glob("*.gbk")):
            if any([word in gbk.name for word in parameters.include]):
                gbks.append(gbk)
        if len(gbks) == 0:
            sys.exit("Error: input is a folder, but no files with the specified parameters in the '--include' argument were found")
    else:
        gbks = list([i])
        if i.suffix != ".gbk":
            sys.exit("Error: input is a file, but does not contain the .gbk extension")
        if not any([word in i.name for word in parameters.include]):
            sys.exit("Error: input is a file, but does not match the strings specified with the '--include' argument")
         
    if len(gbks) == 0:
        sys.exit("No valid GenBank files were found")
    
    os.makedirs(o, exist_ok=True)
            
    # Extract sequences
    copied = 0
    for gbk in gbks:
        try:
            records = list(SeqIO.parse(gbk, "genbank"))
        except ValueError as e:
            print("Error: cannot parse {}".format(gbk))
            continue
        else:
            if tax_filter.lower() in [t.lower() for t in records[0].annotations["taxonomy"]]:
                copied += 1
                copyfile(gbk, o / gbk.name)
    
    
    # Finish
    print("\nDone. Copied {} files".format(copied))
