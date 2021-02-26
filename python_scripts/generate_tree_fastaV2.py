#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 13:40:10 2020
This script generates a FASTA file out of a genbank file with nice headers
which can later be used by phyml.(After conversion with fasta_2_phylip.py) 
sys.argv[1]=filename.gb
sys.argv[2]=int(minimal length)
sys.argv[3]=outputname.fasta
sys.artv[4]=genename for fasta header

@author: maro
"""
from Bio import SeqIO
import re
import sys
all_lact=[]
cnt=0

def clean_string(string):
    string=string.replace(" ","_")
    string=string.replace(",","")
    string=string.replace("(","")
    string=string.replace(")","")
    string=string.replace(":","|")
    string=string.replace("[","")
    string=string.replace("]","")    
    return string


# argv=["","filtered_neighbouring_transporters.gb",300]
file_name=sys.argv[1][:-3]
# file_name=argv[1][:-3]
for seq_record in SeqIO.parse(sys.argv[1], "genbank"):
# for seq_record in SeqIO.parse(argv[1], "genbank"):    
    if len(seq_record)<int(sys.argv[2]):
    # if len(seq_record)<int(argv[2]):
        print("{} filtered with length {}".format(seq_record.id,len(seq_record)))
        continue
    record=seq_record
    name=seq_record.annotations["organism"].split(" ")
    try:
    	species_name="{}_{}".format(name[0],name[1])
    except IndexError:
    	species_name="{}".format(name[0])
    	print("no genus name for {} found".format(seq_record.annotations["organism"]))
    	continue
    gene_id=seq_record.id.split(" ")[0]
    record.description=""
    if "score" in seq_record.annotations["source"]:
        
#        substring = seq_record.annotations["source"]
#        substring=substring.replace(" ","_")
#        substring=substring.replace(",","")
#        substring=substring.replace("(","")
#        substring=substring.replace(")","")
#        substring=substring.replace(":","|")
        record.id="{}_{}_{}_{}".format(gene_id,species_name,sys.argv[4],clean_string(seq_record.annotations["source"]))
        cnt+=1
    else:
        
        record.id="{}_{}_{}".format(gene_id,species_name,sys.argv[4])
    all_lact.append(record)
print("{} transporters found".format(cnt))
SeqIO.write(all_lact, sys.argv[3], "fasta")
