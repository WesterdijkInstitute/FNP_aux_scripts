#!/usr/bin/env python

""""
Extract protein sequences corresponding to core biosynthetic genes annotated
by fungiSMASH.

Only for antiSMASH 5+ output.

Requires: biopython

Input: inputfolder with fungiSMASH results or a single GenBank file
Output: a fasta file

"""

import os
import sys
from pathlib import Path
import argparse
from Bio import SeqIO

__author__ = "Jorge Navarro"
__version__ = "1"
__maintainer__ = "Jorge Navarro"
__email__ = "j.navarro@wi.knaw.nl"


def parameter_parser():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=Path,
        help="Input. Either a GenBank file or a folder with fungiSMASH results", 
        required=True)
    parser.add_argument("-n", "--name", type=str, required=True, 
        help="Name of fasta file.")
    
    return parser.parse_args()


def extract_cbp(gbk: list, fasta: dict) -> None:
    try:
        records = list(SeqIO.parse(gbk, "genbank"))
    except ValueError as e:
        sys.exit("Error: cannot parse {}".format(gbk))
    else:
        for locus_num, record in enumerate(records):
            cds_num = 0
            for feature in record.features:
                if feature.type == "CDS":
                    CDS = feature
                    cds_num += 1
                    
                    protein_type = ""
                    # catch fungiSMASH 5 annotation
                    role = ""
                    if "gene_kind" in CDS.qualifiers:
                        role = CDS.qualifiers["gene_kind"][0]
                        if role == "biosynthetic" and "gene_functions" in CDS.qualifiers:
                            # get any extra labels
                            name = ""
                            if "name" in CDS.qualifiers:
                                name = CDS.qualifiers["name"][0]
                            gene = ""
                            if "gene" in CDS.qualifiers:
                                gene = CDS.qualifiers["gene"][0]
                            protein_id = ""
                            if "proteinId" in CDS.qualifiers:
                                protein_id = CDS.qualifiers["proteinId"][0]
                            elif "protein_id" in CDS.qualifiers:
                                protein_id = CDS.qualifiers["protein_id"][0]
                            
                            
                            protein_types = []
                            for x in CDS.qualifiers["gene_functions"]:
                                if x.startswith("biosynthetic (rule-based-clusters)"):
                                    protein_types.append(x.split("biosynthetic (rule-based-clusters) ")[1].split(":")[0])
                            if len(set(protein_types)) == 1:
                                protein_type = protein_types[0]
                            elif "NRPS_PKS" in CDS.qualifiers:
                                protein_type = "PKS-NRPS_hybrid"
                    
                            try:
                                sequence = CDS.qualifiers["translation"][0]
                            except KeyError:
                                sys.exit("Error: file {} missing 'translation' qualifier".format(gbk))
                            
                            identifier = "{}~L{}+CDS{}".format(gbk.stem, locus_num, cds_num)
                            # write up extra labels
                            identifier = "{} BiosyntheticType:{}".format(identifier, protein_type)
                            if name != "":
                                identifier = "{} name:{}".format(identifier, name)
                            if gene != "":
                                identifier = "{} gene:{}".format(identifier, gene)
                            if protein_id != "":
                                identifier = "{} protein:{}".format(identifier, protein_id)
                            
                            fasta[identifier] = sequence
    return


def write_fasta_file(fasta: dict, n: str) -> None:
    with open("{}.fasta".format(n), "w") as f:
        for header, sequence in fasta.items():
            # format sequence
            l = len(sequence)
            seq80_a = "\n".join([sequence[row*80:(row+1)*80] for row in range(l // 80)])
            seq80_b = ""
            remainder = l % 80
            if remainder > 0:
                seq80_b = "{}\n".format(sequence[-remainder:])
            sequence80 = "{}\n{}".format(seq80_a, seq80_b)
            
            f.write(">{}\n{}".format(header, sequence80))


if __name__ == "__main__":
    parameters = parameter_parser()
    
    i = parameters.input
    if not (i.is_dir() or i.is_file()):
        sys.exit("Error: {} is not a valid folder or file".format(i))

    n = parameters.name
    bad_chars= set("/|\:<>*")
    if len(set(n) & bad_chars) > 0:
        sys.exit("Error: illegal characters used with parameter --name")
    
    if i.is_dir():
        gbks = list(sorted(i.glob("*region*.gbk")))
        num_gbks = len(gbks)
        if num_gbks == 0:
            sys.exit("Error: input is a folder, but no files with the *region*.gbk structure were found")
    else:
        gbks = list([i])
        if i.suffix != ".gbk":
            sys.exit("Error: input is a file, but does not contain the .gbk extension")
            
            
    fasta = dict()
    for gbk in gbks:
        extract_cbp(gbk, fasta)
    
    print("Got {} sequences".format(len(fasta)))
    
    write_fasta_file(fasta, n)
