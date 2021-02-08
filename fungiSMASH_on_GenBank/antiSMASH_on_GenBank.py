#! /usr/bin/env python

"""
Launches fungiSMASH in a folder full of gbk files and extracts their 'final' result
"""

import warnings
warnings.simplefilter(action='ignore', category='FutureWarning')
import argparse
import os
import sys
from pathlib import Path
import subprocess
from subprocess import STDOUT
from multiprocessing import Pool
import shutil


def parameter_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputfolder", help="Place with files", 
        required=True, type=Path)
    parser.add_argument("-r", "--results", help="Base output folder for \
        antiSMASH intermediate results", required=True, type=Path)
    parser.add_argument("-o", "--outputfolder", help="Base output directory \
        used for antiSMASH _final_ GenBank files", required=True, type=Path)
    parser.add_argument("-p", "--processes", type=int, help="Number of \
        simultaneous antiSMASH processes. Default: 2", default=2)
    parser.add_argument("-c", "--cpus", type=int, help="CPUs used for each \
        antiSMASH run. Default: 2", default=2)

    return parser.parse_args()

def launch_antismash(args, params, cpus, gbk):
    cmd = []
    cmd.append("antismash")
    cmd.append("--cpus")
    cmd.append(str(cpus))
    if len(params) > 0:
        cmd.extend(params)
    cmd.append("--output-dir")
    cmd.append(str(args.results / gbk.stem))
    cmd.append(str(gbk))
    
    final_gbk = args.results / gbk.stem / gbk.name
    # Check if results already exist
    if final_gbk.is_file():
        #print("{} already exists, skipping antiSMASH for it".format(final_gbk))
        # If so, check if "final gbk" is already in the output folder
        if not (args.outputfolder / gbk.name).is_file():
            shutil.copy(final_gbk, args.outputfolder)
        return False
    
    proc = subprocess.run(cmd, stderr=STDOUT, encoding="utf-8")
    try:
        proc.check_returncode()
    except subprocess.CalledProcessError:
        print("Error: {}".format(gbk))
        print(proc.stderr)
        return False
    else:
        shutil.copy(final_gbk, args.outputfolder)
        return True

if __name__ == '__main__':
    args = parameter_parser()
    
    cpus = args.cpus
    procs = args.processes
    
    if cpus < 1 or procs < 1:
        sys.exit("Error. Use positive numbers for --cpus and --processes")
    
    if not args.results.is_dir():
        os.makedirs(args.results)
    if not args.outputfolder.is_dir():
        os.makedirs(args.outputfolder)
    
    # Read parameter file. Each line should contain parameter-space-value
    params = list()
    parameter_file = Path(__file__).parent / "antismash5_parameters.tsv"
    if parameter_file.is_file():
        with open(parameter_file) as f:
            for line in f:
                x = line.strip()
                if x == "" or x[0] == "#":
                    continue
                params.extend(x.split(" "))
                
    with Pool(procs) as pool:
        for gbk in args.inputfolder.glob("*.gbk"):
            pool.apply_async(launch_antismash(args, params, cpus, gbk))
        pool.close()
        pool.join()
    
    
