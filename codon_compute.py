#!/usr/bin/env python3

import os, gzip, itertools

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
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence

url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))

DictCG = {'A':0, 'T':0, 'C':0, 'G':0}
#to automate making the codon dictionaries I found this idea from stackflow
#https://stackoverflow.com/questions/20001045/how-to-generate-the-keys-of-a-dictionary-using-permutations
from itertools import permutations
codondict1 = {"".join(key):0 for key in permutations('ATGCATGCATGC', 3)}
codondict2 = {"".join(key):0 for key in permutations('ATGCATGCATGC', 3)}
with gzip.open(file1,"rt") as fh:
    seqs = aspairs(fh)
    n = 0
    lens1 = 0
    for seq in seqs:
        n += 1 
        bp1 = 0
        for line in seq[1]:
            bp1 += len(line)
            for base in line:        
                DictCG[base] += 1
        codons = []
        for codon in range(0, len(seq[1]), 3):
            codons.append(seq[1][codon : codon + 3])
        for codon in range(0, len(codons)):
            if codons[codon] in codondict1:
                codondict1[codons[codon]] += 1 
        lens1 += bp1   
print("according to the suplied FASTA files:")
print ("The Salmonella genonome has ",n,"genes and their total lengths sum to ",lens1, "Bp") 

with gzip.open(file2,"rt") as fh:
    seqs = aspairs(fh)
    n2 = 0
    lens2 = 0
    for seq in seqs:
        n2 += 1 
        bp2 = 0
        for line in seq[1]:
            bp2 += len(line)
            for base in line:
                DictCG[base] += 1
        codons = []
        for codon in range(0, len(seq[1]), 3):
            codons.append(seq[1][codon : codon + 3])
        for codon in range(0, len(codons)):
            if codons[codon] in codondict2:
                codondict2[codons[codon]] += 1 
        lens2 += bp2   
perGC = round(float(100 * ((DictCG['G'] + DictCG['C']) / (lens1+lens2))), 1)
print ("The Tuberculosis genonome has ",n2,"genes and their total lengths sum to ",lens2, "Bp") 
print("The G+C percentage in the dataset is ",perGC,"%") 
# the following two comented lines are useful to checking if the dictionary entries are consistant with lengths.
#print(3*(sum(codondict1.values())))
#print(3*(sum(codondict2.values())))
for key in codondict1:
    codondict1[key] = float(codondict1[key] / (lens1 / 3))
    codondict1[key] = round(codondict1[key], 4) 
for key in codondict2:
    codondict2[key] = float(codondict2[key] / (lens2 / 3))
    codondict2[key] = round(codondict2[key], 4) 
# insures that frequencies sum to 1 
#print(sum(codondict1.values()))
#print(sum(codondict2.values())) 
print("-----------------------------------------------")
print("codon statistics table (Sp1= Salmonella; Sp2= Tuberculosis):")
print("\t".join(['Codon', 'Frequency in Sp1', 'Frequency in Sp2']))
for key in codondict1:
    print("\t".join([key, str(codondict1[key]), "\t", str(codondict2[key])]))