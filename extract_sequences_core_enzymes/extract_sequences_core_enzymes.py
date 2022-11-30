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

import sys
from pathlib import Path
import argparse
from Bio import SeqIO

__author__ = "Jorge Navarro"
__version__ = "1.2"
__maintainer__ = "Jorge Navarro"
__email__ = "jorge.navarromunoz@wur.nl"


def parameter_parser():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-i", "--input", type=Path,
        help="Input. Either a GenBank file or a folder with fungiSMASH results", 
        required=True)
    parser.add_argument("-o", "--output", type=str, help="Name of fasta file. \
        Default: same as input")
    parser.add_argument("-b", "--biosynthetic", default=False, 
        action="store_true", help="If activated, only extract sequences \
        flagged by fungiSMASH as 'biosynthetic'")
    parser.add_argument("--include", default=["region"], nargs='+', 
        help="Optional list of strings used to filter GenBank files to be \
        analyzed. Default: region")
    parser.add_argument("-t", "--taxonomy_filter", type=str, help="Filter \
        extraction using taxonomy within the GenBank file(s) (if annotated)")
    parser.add_argument("-d", "--domain", type=str, help="Extract domain \
        sub-sequences from core biosynthetic proteins identified by antiSMASH. \
        Valid arguments, choose one from: KS, AT, A")
    parser.add_argument("-w", "--whole", default=False, action='store_true',
        help="If activated and --domain was used, extract the whole sequence, \
        instead of only the domain subsequence.")
    
    return parser.parse_args()


def extract_cbp(gbk: list, fasta: dict, biosynthetic: bool, tax_filter: str, 
    domain: str, whole: bool) -> None:
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
                        exit(f"Error: file {gbk} missing 'translation' qualifier")
                        
                    basic_label = f"{gbk.stem}~L{locus_num}+CDS{cds_num}"
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
                                    
                                    if whole:
                                        fasta[f"{basic_label}{label_extra}"] = sequence
                                    else:
                                        coords = name[1][1:-2]
                                        start = int(coords.split("-")[0])
                                        end = int(coords.split("-")[1])
                                        domains[dom_acc] = (start, end)
                                        
                                        fasta[f"{basic_label}@{coords}_{domain}{label_extra}"] = sequence[start:end]
                            
                    else:
                        fasta[f"{basic_label}{label_extra}"] = sequence
    return


def write_fasta_file(fasta: dict, filename: str) -> None:
    with open(f"{filename}.fasta", "w") as f:
        for header, sequence in fasta.items():
            # format sequence
            l = len(sequence)
            seq80_a = "\n".join([sequence[row*80:(row+1)*80] for row in 
                range(l // 80)])
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
        exit(f"Error: {i} is not a valid folder or file")

    if parameters.output:
        bad_chars= set("/|\:<>*")
        if len(set(parameters.output) & bad_chars) > 0:
            exit("Error: illegal characters used with parameter --output")
        output_filename = parameters.output
        
    if parameters.domain:
        d_alias = {"KS": "PKS_KS(Iterative-KS)",
                   "AT": "PKS_AT",
                   "A": "AMP-binding"}
        try:
            domain = d_alias[parameters.domain.upper()]
        except KeyError:
            exit(f"Error: found unsupported argument for --domain " \
                f"parameter ({parameters.domain})")
    else:
        domain = ""
    
    # Get a list of all GenBank files
    if i.is_dir():
        gbks = list()
        for gbk in sorted(i.glob("*.gbk")):
            if any([word in gbk.name for word in parameters.include]):
                gbks.append(gbk)
        if len(gbks) == 0:
            exit("Error: input is a folder, but no files with the specified \
                parameters in the '--include' argument were found")
    else:
        gbks = list([i])
        if i.suffix != ".gbk":
            exit("Error: input is a file, but does not contain the .gbk \
                extension")
        if not any([word in i.name for word in parameters.include]):
            exit("Error: input is a file, but does not match the strings \
                specified with the '--include' argument")
            
            
    # Extract sequences
    fasta = dict()
    for gbk in gbks:
        extract_cbp(gbk, fasta, parameters.biosynthetic, 
            parameters.taxonomy_filter, domain, parameters.whole)
    
    
    # Write output
    print(f"Got {len(fasta)} sequence(s)")
    write_fasta_file(fasta, output_filename)
