# FNP auxiliary scripts

Small Python scripts to help with fungal biosynthetic gene cluster mining

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

## Extract sequences

Use this script to extract the amino acid sequences corresponding to genes marked by fungiSMASH 5 as 'biosynthetic'.

```
usage: extract_sequences.py [-h] -i INPUT -n NAME

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input. Either a GenBank file or a folder with
                        fungiSMASH results
  -n NAME, --name NAME  Name of fasta file.
```

The headers will have the format `[GenBank base name]~L[Locus number]+CDS[CDS number] BiosyntheticType:[annotated fungiSMASH type] name:[name, optional] gene:[gene ID, optional] protein:[protein ID, optional]`
