#!/usr/bin/env python

""""
Extract protein sequences corresponding to core biosynthetic genes annotated
by fungiSMASH.

Only for antiSMASH 5+ output.

Requires: biopython, Python 3

Input: inputfolder with fungiSMASH results or a single GenBank file
Output: a fasta file

Can filter by:
    - Taxonomy (if GenBank is annotated)
    - Name (using a filename substring)
    - Domain. Recognizes a limited set of domains annotated in the GenBank 
    by antiSMASH during BGC identification:
        * KS: "PKS_KS(Iterative-KS)"
        * AT: "PKS_AT"
        * A: "AMP-binding"

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
    parser.add_argument("-b", "--biosynthetic", default=False, 
        action="store_true", help="If activated, only extract proteins flagged \
        fungiSMASH as 'biosynthetic'")
    parser.add_argument("--include", default=["region"], type=str, 
        help="Optional list of strings used to filter GenBank files to be \
        analyzed")
    parser.add_argument("-t", "--taxonomy_filter", type=str, help="Filter \
        extraction using taxonomy within the GenBank file(s) (if annotated)")
    parser.add_argument("-d", "--domain", type=str, help="Extract domain \
        sub-sequences from core biosynthetic proteins identified by antiSMASH. \
        Valid arguments: KS, AT, A")
    
    return parser.parse_args()


def extract_cbp(gbk: list, fasta: dict, biosynthetic: bool, tax_filter: str, domain: str) -> None:
    try:
        records = list(SeqIO.parse(gbk, "genbank"))
    except ValueError as e:
        sys.exit("Error: cannot parse {}".format(gbk))
    else:
        if tax_filter:
            if tax_filter.lower() not in [t.lower() for t in records[0].annotations["taxonomy"]]:
                return
        
        for locus_num, record in enumerate(records):
            cds_num = 0
            for feature in record.features:
                if feature.type == "CDS":
                    CDS = feature
                    cds_num += 1
                    
                    # user specified 'biosynthetic'. Try to catch fungiSMASH's annotation
                    if biosynthetic:
                        role = ""
                        if "gene_kind" in CDS.qualifiers:
                            role = CDS.qualifiers["gene_kind"][0]
                            # CDS has annotation but it's not 'biosynthetic'
                            if role != "biosynthetic":
                                continue
                        # CDS doesn't have annotation
                        else:
                            continue
                        
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
                    
                    # Get protein type for header
                    protein_types = []
                    protein_type = ""
                    if "gene_functions" in CDS.qualifiers:
                        for x in CDS.qualifiers["gene_functions"]:
                            if x.startswith("biosynthetic (rule-based-clusters)"):
                                protein_types.append(x.split("biosynthetic (rule-based-clusters) ")[1].split(":")[0])
                        if len(set(protein_types)) == 1:
                            protein_type = protein_types[0]
                        elif "NRPS_PKS" in CDS.qualifiers:
                            protein_type = "PKS-NRPS_hybrid"
                    
                    # Get sequence
                    try:
                        sequence = CDS.qualifiers["translation"][0]
                    except KeyError:
                        sys.exit("Error: file {} missing 'translation' qualifier".format(gbk))
                        
                    basic_label = "{}~L{}+CDS{}".format(gbk.stem, locus_num, cds_num)
                    label_extra = ""
                    # write up extra labels
                    if protein_type != "":
                        label_extra = "{} BiosyntheticType:{}".format(label_extra, protein_type)
                    if name != "":
                        label_extra = "{} name:{}".format(label_extra, name)
                    if gene != "":
                        label_extra = "{} gene:{}".format(label_extra, gene)
                    if protein_id != "":
                        label_extra = "{} protein:{}".format(label_extra, protein_id)
                        
                    # Extract subsequence(s) if user needs specific annotated domain
                    if domain:
                        if "NRPS_PKS" in CDS.qualifiers:
                            domains = dict()
                            for q in CDS.qualifiers["NRPS_PKS"]:
                                if q.startswith("Domain"):
                                    name = q.replace("Domain: ", "").split(" ")
                                    dom_acc = name[0]
                                    
                                    if domain != dom_acc:
                                        continue
                                    
                                    coords = name[1][1:-2]
                                    start = int(coords.split("-")[0])
                                    end = int(coords.split("-")[1])
                                    domains[dom_acc] = (start, end)
                                    
                                    fasta["{}@{}_{}{}".format(basic_label, coords, domain, label_extra)] = sequence[start:end]
                            
                    else:
                        fasta["{}{}".format(basic_label, label_extra)] = sequence
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
    
    # Validate parameters
    i = parameters.input
    if not (i.is_dir() or i.is_file()):
        sys.exit("Error: {} is not a valid folder or file".format(i))

    n = parameters.name
    bad_chars= set("/|\:<>*")
    if len(set(n) & bad_chars) > 0:
        sys.exit("Error: illegal characters used with parameter --name")
        
    if parameters.domain:
        d_alias = {"KS": "PKS_KS(Iterative-KS)",
                   "AT": "PKS_AT",
                   "A": "AMP-binding"}
        try:
            domain = d_alias[parameters.domain.upper()]
        except KeyError:
            sys.exit("Error: found unsupported argument for --domain parameter ({})".format(parameters.domain))
    
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
            
            
    # Extract sequences
    fasta = dict()
    for gbk in gbks:
        extract_cbp(gbk, fasta, parameters.biosynthetic, parameters.taxonomy_filter, domain)
    
    
    # Write output
    print("Got {} sequence(s)".format(len(fasta)))
    write_fasta_file(fasta, n)
