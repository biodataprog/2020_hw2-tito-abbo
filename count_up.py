#!/usr/bin/env python3

# this is a python script template
# this next line will download the file using curl

gff="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
fasta="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

import os,gzip,itertools,csv,re

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence



if not os.path.exists(gff):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")
genecount = 0  
genelengths = []  
with gzip.open(gff,"rt") as fh:
    # now add code to process this
    gff = csv.reader(fh,delimiter="\t")
    for row in gff:
        if row[0].startswith("#"):
            continue
        if row[2] == "gene":
            genecount += 1
            gl = int(row[4]) - int(row[3])
            genelengths.append(gl)
totalcoding = sum(genelengths)
print("The GFF file reports ",genecount,"genes in the E. coli genome.")
print("Their total lengths add up to ",totalcoding,"Bp.")
with gzip.open(fasta,"rt") as f:
    seqs = dict(aspairs(f))
n=0
for k,v in seqs.items():
   seq_id = k
   seq = v
bp = 0
for line in seq:
    bp += len(line)
print ("The FASTA file reports ",bp,"Bp in the E. coli genome")
percode = float(100 * (totalcoding / bp))
percode = round(percode, 1)
print ("Therefore, according to the two files ",percode,"% of the E. coli genome is coding")


