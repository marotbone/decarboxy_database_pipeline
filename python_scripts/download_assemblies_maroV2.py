#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 13:22:59 2020
This script takes as input sys.argv[1] a mapping_table created by create_mapping_table.py.
A file with 2 columns: 1) Protein accession number 2) corresponding ncbi contig id (accession?) 
on which corresponding protein was found.
As sys.argv[2] a complete-mapping output file is generated which is more readable and contains additional information
about contigs.
sys.argv[3] needs to be the path were contigs should be stored. Files in there will be
named after ncbi contig id. To quickly find interesting files complete-mapping output sys.argv[2] is usefull.

@author: maro
"""
import sys
import os
import zipfile
from Bio import ExPASy
from Bio import SwissProt

from Bio import SeqIO
from Bio import Entrez
Entrez.api_key = '6201482c3419ae91449619dccb91fe413e08'
#Entrez.api_key = '9a0fbb3042bcab68f6eb9fb0ad1b52d3d309'
Entrez.email = "widmer.maro@gmail.com"
#Entrez.email = 'thomas.roder@bioinformatics.unibe.ch'

assembly_ids=[]


with open(sys.argv[1], "r") as f:
    for line in f:
        assembly_id=line.split()[1]#take assembly id
        prot_accession=line.split()[0]
        if assembly_id=="no_contig_found":
            continue
        else:
            assembly_ids.append((prot_accession,assembly_id))
f=open(sys.argv[2],"w")
f.close()
# print(len(assembly_ids))
# assembly_ids=list(set(assembly_ids)) #remove doubled entries
print(len(assembly_ids))
cnt=len(assembly_ids)
for i in assembly_ids: 
    attempt=0
    while attempt<=5:
        try:
            if os.path.isfile('{}/{}.gb'.format(sys.argv[3],i[1])):
                print("file for {} already exists".format(i[1]))
                record = list(SeqIO.parse('{}/{}.gb'.format(sys.argv[3],i[1]), "genbank"))[0]
                with open(sys.argv[2],"a") as a:
                    a.write("{}\t{}\t{}\t{}\n".format(i[0],i[1],record.id,record.description))
                cnt-=1
                attempt=6
                continue
            else:  
                handle = Entrez.efetch(db="nucleotide", id=i[1], rettype="gb", retmode="text")
                record = SeqIO.read(handle, "genbank")  
                with open(sys.argv[2],"a") as a:
                    a.write("{}\t{}\t{}\t{}\n".format(i[0],i[1],record.id,record.description))
                SeqIO.write(record, '{}/{}.gb'.format(sys.argv[3],i[1]), "genbank")
                cnt-=1
                print("{} files to download".format(cnt))
                attempt=6
        except:
            print("There was an error, trying again for {}".format(i[1]))
            attempt+=1
#zipfile.ZipFile('genbank_repository/test.gb.gzip', mode='w').write("test.gb")
