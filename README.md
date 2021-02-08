# FNP auxiliary scripts

Small Python 3 scripts to help with fungal biosynthetic gene cluster mining.


## Rename regions

Use this script to easily rename all the 'region' GenBank files from the fungiSMASH results. This is specially important for results with generic names (e.g. 'scaffold'), if using tools like [BiG-SCAPE](https://git.wur.nl/medema-group/BiG-SCAPE/).

```
usage: rename_regions.py [-h] -i INPUTFOLDER -s STRING

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFOLDER, --inputfolder INPUTFOLDER
                        Path to fungiSMASH results folder
  -s STRING, --string STRING
                        All found .gbk files will be renamed
                        '[string].region[XXX].gbk'
```

Note that this script will write consecutive region numbers (instead of fungiSMASH's default, which uses the scaffold number).


## Filter copy

Use this script to copy GenBank files according to a taxonomical filter. It depends on Biopython.

```
usage: filtercopy_genbank.py [-h] -i INPUT -o OUTPUTFOLDER -t TAXONOMY_FILTER [--include INCLUDE]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input. Either a GenBank file or a folder with fungiSMASH results
  -o OUTPUTFOLDER, --outputfolder OUTPUTFOLDER
                        Output folder where GenBank file(s) will be copied to
  -t TAXONOMY_FILTER, --taxonomy_filter TAXONOMY_FILTER
                        Filter extraction using taxonomy within the GenBank file(s) (if annotated)
  --include INCLUDE     Optional list of strings used to filter GenBank files to be analyzed. Default: 'region', 'BGC'
```

For this filtering to work, original GenBank files must be properly annotated, e.g.:
```
LOCUS       BGC0000003             37168 bp    DNA     linear   PLN 03-MAY-2005
DEFINITION  Alternaria alternata AF-toxin biosynthesis gene cluster (Aft9-1,
            Aft10-1, Aft11-1, Aft12-1, AftR-2, Aft3-2), complete cds.
ACCESSION   BGC0000003
VERSION     BGC0000003.1
KEYWORDS    .
SOURCE      Alternaria alternata
  ORGANISM  Alternaria alternata
            Eukaryota; Fungi; Dikarya; Ascomycota; Pezizomycotina;
            Dothideomycetes; Pleosporomycetidae; Pleosporales; Pleosporineae;
            Pleosporaceae; Alternaria; Alternaria alternata group.
...
```


## Extract sequences

Use this script to extract the amino acid sequences corresponding to genes marked by fungiSMASH 5 as 'biosynthetic'. It depends on Biopython.

```
usage: extract_sequences_core_enzymes.py [-h] -i INPUT -n NAME [-b] [--include INCLUDE] [-t TAXONOMY_FILTER] [-d DOMAIN]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input. Either a GenBank file or a folder with fungiSMASH results
  -n NAME, --name NAME  Name of fasta file.
  -b, --biosynthetic    If activated, only extract proteins flagged fungiSMASH as 'biosynthetic'
  --include INCLUDE     Optional list of strings used to filter GenBank files to be analyzed
  -t TAXONOMY_FILTER, --taxonomy_filter TAXONOMY_FILTER
                        Filter extraction using taxonomy within the GenBank file(s) (if annotated)
  -d DOMAIN, --domain DOMAIN
                        Extract domain sub-sequences from core biosynthetic proteins identified by antiSMASH. Valid arguments: KS, AT, A
```

The headers will have the format `[GenBank base name]~L[Locus number]+CDS[CDS number] BiosyntheticType:[annotated fungiSMASH type] name:[name if annotated] gene:[gene ID if annotated] protein:[protein ID if anntoatedl]`. If extracting a particular domain subsequence, additional information will be added. For example:
```
>BGC0000003.1~L0+CDS1@0-339_PKS_KS(Iterative-KS) BiosyntheticType:T1PKS gene:AFT9-1 protein:BAD97694.1
MDPQQRLLLETTYEALENAGIPQANTNGSNTSVHVAMFTRDYDRNVYKDTVGIPKYQVTGTGEAIMSNRISHIFNLHGPS
MTIDTGCSGAMTAVSQACMSLRSGDCDIALAGAVNLIMSPDHHISMSNLHMLNAEGKSYAFDSRGAGYGRGEGVATIVMK
RLDDAVRCHDPIRAVILDAVINQDGYTAGITLPSSEAQAQLERKALNRVGLKPQEVAYIEAHGTGTAAGDAAELDALSSV
FCVDRDLPLYVGSVKSNIGHLEAASGMAALIKATLMLENEAIPPSINFSRPKENLRIDERNIKIPTALQPWPKGASARIC
VNSFGYGGTNAHAILERAP
```

## Process GenBank files with fungiSMASH

Use this script to launch instances of fungiSMASH on a folder with GenBank files. 

```
usage: antiSMASH_on_GenBank.py [-h] -i INPUTFOLDER -r RESULTS -o OUTPUTFOLDER
                               [-p PROCESSES] [-c CPUS]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFOLDER, --inputfolder INPUTFOLDER
                        Place with files
  -r RESULTS, --results RESULTS
                        Base output folder for antiSMASH intermediate results
  -o OUTPUTFOLDER, --outputfolder OUTPUTFOLDER
                        Base output directory used for antiSMASH _final_
                        GenBank files
  -p PROCESSES, --processes PROCESSES
                        Number of simultaneous antiSMASH processes. Default: 2
  -c CPUS, --cpus CPUS  CPUs used for each antiSMASH run. Default: 2
```
