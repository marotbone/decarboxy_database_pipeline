#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 10:02:12 2020
sys.argv[1]="ncbi search term"
sys.argv[2]=outputfilename
@author: maro
"""

import sys
from Bio import ExPASy
from Bio import Entrez
from Bio import SeqIO
import time
from socket import error as SocketError
import errno

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

Entrez.api_key = '6201482c3419ae91449619dccb91fe413e08'
Entrez.email = "widmer.maro@gmail.com"

records=[]
try:
    for seq_record in SeqIO.parse(sys.argv[2], "genbank"):
        records.append(seq_record)
        stored_sequences=len(records)
except FileNotFoundError:
    stored_sequences=0




attempt=0
while attempt<=5:    
    try:
        search_handle = Entrez.esearch(db="protein", term=sys.argv[1], idtype="acc", usehistory="y")
        # handle = Entrez.esearch(db="protein", term="(((Burkholderiales[Organism]) AND ornithine decarboxylase[Protein Name])) ", idtype="acc",retmax=100000)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        attempt=6
        time.sleep(1)
    except RuntimeError:
        print("Runtime error occured...trying again...")
        attempt+=1
        time.sleep(5)


# total_records=int(record["Count"])
# print(record["IdList"])
# accessions=record["IdList"]
acc_list = search_results["IdList"]
count = int(search_results["Count"])
print("{} genes found..".format(count))

if count!=stored_sequences:
    

    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]
    
    
    batch_size = 50
    out_handle = open(sys.argv[2], "w")
    for start in range(0, count, batch_size):
        attempt=0
        end = min(count, start + batch_size)
        print("Going to download record %i to %i" % (start + 1, end))
        fetch_handle = Entrez.efetch(
            db="protein",
            rettype="gb",
            retmode="text",
            retstart=start,
            retmax=batch_size,
            webenv=webenv,
            query_key=query_key,
            idtype="acc",
        )
        time.sleep(1)
        while attempt<=5:
            try:
                data = fetch_handle.read()
                time.sleep(1)
                attempt=6
            except SocketError as e:
                if e.errno != errno.ECONNRESET:
                    raise # Not error we are looking for
                attempt+=1
        fetch_handle.close()
        out_handle.write(data)
    out_handle.close()
else:
    print("File with same number of sequences already exists")