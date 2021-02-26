#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 10:43:08 2020
This script reads a protein FASTA file and converts it into nucleotide FASTA file.
The ncbi protein accession number needs to follow after ">". Example: >EQC57918.1 blablalalba
sys.argv[1]=protein FASTA file
sys.argv[2]=output nucleotide filename.fasta

@author: maro
"""

import time
import re
import sys
from Bio import SeqIO
from Bio import Entrez
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from http.client import IncompleteRead
Entrez.api_key = '6201482c3419ae91449619dccb91fe413e08'
Entrez.email = "widmer.maro@gmail.com"


records=[]
mapping_dict={}
nuc_records=[]
contig_dict={}

# pattern1 = re.compile(r'[A-Z]{4}\d{8}')
# pattern2= re.compile(r'gap\(\d+\)')
# pattern = re.compile(r'join\(')

def clean_string(string):
    string=string.replace(" ","_")
    string=string.replace(",","")
    string=string.replace("(","")
    string=string.replace(")","")
    string=string.replace(":","|")
    string=string.replace("[","")
    string=string.replace("]","")    
    return string


def search_assembly_id(prot_accession):
    attempt=0
    while attempt < 5:
        try:
            fetch_handle = Entrez.elink(dbfrom="protein",db = 'nuccore',
                                        id=prot_accession,
                                        linkname="protein_nuccore")
            entrez_record = Entrez.read(fetch_handle)
            assembly_id=entrez_record[0]["LinkSetDb"][0]['Link'][0]['Id']
            attempt=5   
        except RuntimeError:
            print("Runtime error searching for {} ...trying again...".format(prot_accession))
            time.sleep(3)
            attempt+=1
        except IndexError:
            print("Index error {} ...trying again...".format(prot_accession))
            time.sleep(3)
            if attempt == 4:
                assembly_id=False
            attempt+=1
        except UnboundLocalError:
            print("This protein record was suppressed because it is no longer annotated on any genome".format(prot_accession))
            assembly_id=False
            attempt+=6
    return assembly_id


def download_assembly_genbank(assembly_id):
    attempt=0
    while attempt < 5:
        try:
            handle = Entrez.efetch(db="nucleotide", id=assembly_id, rettype="gb", retmode="text")
            assembly = SeqIO.read(handle, "genbank")
            handle.close()
            attempt=5
        except RuntimeError:
            print("Runtime error downloading {}...trying again...".format(assembly_id))
            time.sleep(3)
            attempt+=1
        except IncompleteRead:
            print("Incomplete read for {}...trying again...".format(assembly_id))
            time.sleep(3)
            attempt+=1
    return assembly    

def download_assembly_fasta(assembly_id):
    attempt=0
    while attempt < 5:
        try:
            handle = Entrez.efetch(db="nucleotide", id=assembly_id, rettype="fasta", retmode="text")
            fasta = SeqIO.read(handle, "fasta")
            handle.close()
            attempt=5
        except RuntimeError:
            print("Runtime error downloading {}...trying again...".format(assembly_id))
            time.sleep(3)
            attempt+=1
    return fasta

def extract_CDS_location(assembly_SeqRecord,protein_accession):
    flag=0
    for feature in assembly.features:
        try:                        
                    # print(assembly.id)   
            if feature.qualifiers["protein_id"][0] == protein_accession: #if a feature in assembly has the proteins id
                f=SeqFeature(feature.location) #save location as SeqFeature
                flag=1
                print("{} in {} found".format(record.id,assembly.id))
                continue
                    
        except KeyError:
                continue
    if flag==1:
        return f
    elif flag==0:
        return False
  
def build_scaffold_sequence(contig_info,contig_dict):
    
    pattern1 = re.compile(r'[A-Z]{4}\d{8}')
    pattern2= re.compile(r'gap\(\d+\)')
    
    if contig_info==False:
        print("No contig info found for available".format())
        return "N"
    else:        
        match1 = pattern1.findall(contig_info) #regex for contig accessions
        match2 = pattern2.findall(contig_info) #regex for gaps
        if match1 and match2:
            print("contigs = {}".format(match1))
            print("gaps = {}".format(match2))
            scaffold_seq=Seq("")
            for i in range(len(match1)-1):
                
                if match1[i] in contig_dict:
                    print("{} was allready downloaded")
                    scaffold_seq+=contig_dict[match1[i]]
                else:
                    print("downloading {}...".format(match1[i]))
                    assembly=download_assembly_genbank(match1[i])
                    contig_dict[match1[i]]=assembly.seq
                    scaffold_seq+=assembly.seq
                gap="N"*int(match2[i][4:len(match2[i])-1])
                scaffold_seq+=gap
            if match1[-1] in contig_dict:
                print("{} was allready downloaded")
                scaffold_seq+=contig_dict[match1[i]]
            else:
                
                print("downloading {}...".format(match1[-1]))
                assembly=download_assembly_genbank(match1[-1])
                scaffold_seq+=assembly.seq
            return scaffold_seq, contig_dict
def sequence_counter(file_name,file_extension):    
    protein_ids=[]
    try:
        for seq_record in SeqIO.parse(file_name, file_extension):
            protein_ids.append(seq_record.id)
        count=len(protein_ids)
    except FileNotFoundError:
        count=0
    return count
       
if sequence_counter(sys.argv[1], "fasta") != sequence_counter(sys.argv[2],"fasta"):
    
    nuc_seq=""
    for record in SeqIO.parse(sys.argv[1], "fasta"):
    # for record in SeqIO.parse("filtered_neighbouring_transporters.fasta", "fasta"):
        old_locus_tag=False
        attempt=0
        while attempt<=5:        
            try:
                
                assembly_id=search_assembly_id(record.id)
                if assembly_id==False:
                    print("No assembly found...obsolet entry..skipping record...\n")
                    break
                else:	
                    assembly = download_assembly_genbank(assembly_id)
                try:
                    contig_info=assembly.annotations["contig"]
                except KeyError:
                    contig_info=False
                flag=0
                f=extract_CDS_location(assembly,record.id)
                
                if f==False:
                    print("Could not find protein accession {} in assembly {}".format(record.id,assembly.id))
                    print("Downloading all FASTA sequences of {}".format(assembly.id))
                    handle = Entrez.efetch(db="nucleotide", id=assembly_id, rettype="fasta_cds_na", retmode="text")
                    for entry in SeqIO.parse(handle, "fasta"):
                        if record.id in entry.id:
                            nuc_seq=entry.seq
                            continue
                    
                    # if assembly.id.startswith("NZ_"):                 
                    #     ID=assembly.id[3:]
                    #     print("Trying with {}".format(ID))
                    #     assembly=download_assembly_genbank(ID)
                    #     f=extract_CDS_location(assembly,record.id)
                    #     if f == False:
                    #         print("Could not find protein accession {} in assembly {}".format(record.id,assembly.id))
                            
                            # attempt=6
                            # continue
                if assembly.seq.startswith("N") and nuc_seq=="":
                    print("no sequence found in {}".format(assembly.id))
                    print("downloading FASTA sequence")
                    assembly_fasta=download_assembly_fasta(assembly_id)
                    assembly.seq=assembly_fasta.seq
                # elif assembly.seq.startswith("N") and nuc_seq=="" and f==False:
                #     print("Possibly obsolete entry...skipping...\n")
                    # if assembly.id.startswith("NZ_"):
                    #     print("no sequence found in {}".format(assembly.id))
                    #     ID=assembly.id[3:]
                    #     print("trying with {}".format(ID))
                    #     assembly=download_assembly_genbank(ID)
                    #     if assembly.seq.startswith("N"):
                    #         print("still no sequence found...building scaffold")
                    #         if contig_info==False:
                    #             print("No contig info for {}".format(assembly.id))
                    #         else:
                    #             assembly.seq,contig_dict=build_scaffold_sequence(contig_info,contig_dict)
                    # elif contig_info==True:
                    #     print("no sequence found in {}...building scaffold".format(assembly.id))
                    #     assembly.seq,contig_dict=build_scaffold_sequence(contig_info,contig_dict)
                    # else:
                    #     print("no nucleotide sequence available for {}\n".format(record.id))
                if nuc_seq=="":        
                    try:
                        nuc_seq=f.extract(assembly.seq)   
                    except AttributeError:
                        print("Possibly obsolete entry...skipping...\n")
                        break
                elif len(nuc_seq)>50:
                    nuc_seq=nuc_seq
                if nuc_seq.startswith("N"):
                    print("no CDS for {}\n".format(record.id))
                else:
                    print("CDS OK for {}\n".format(record.id))
               
                record.id="{}_{}_{}".format(record.id,sys.argv[3],clean_string(assembly.annotations["organism"]))
                nuc_record=SeqRecord(nuc_seq,id=record.id,description="")
                nuc_records.append(nuc_record)
                nuc_seq=""
                attempt=6
            except RuntimeError:
                print("Blob error...trying again...")
                time.sleep(10)
                attempt+=1
    
    # SeqIO.write(nuc_records, "nucleotides_neighbouring_transporters.fasta", "fasta")
    SeqIO.write(nuc_records, sys.argv[2], "fasta")       
else:
    print("A nucleotide FASTA file with same number of sequences has already been generated")
