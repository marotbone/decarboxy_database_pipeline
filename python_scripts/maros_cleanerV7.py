#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 19:54:04 2020
sys.argv[1] = input .gb file
sys.argv[2] = min length (int)
sys.argv[3] = filtering similarity threshold
sys.argv[4] = summary output file name
sys.argv[5] = filtered output file name (without extension) .fasta and .gb will be created
@author: maro
"""
import sys
from Bio import SeqIO
from Bio import Align
#from Bio.Align import 
from itertools import compress
from Bio.Align import substitution_matrices

trans_flag={}

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
# aligner.match_score = 1
# aligner.mismatch_score = -1
# aligner.open_gap_score = -1
# aligner.extend_gap_score = -1

name=sys.argv[1][:-3]#remove .gb ending and use as name

# argv=["","all_ncbi_tyrosine_decarboxylases.gb",500,0.7]
# name=argv[1][:-3]
records=[]
#read genbank file into records
for seq_record in SeqIO.parse(name+".gb", "genbank"):
    records.append(seq_record)

f = open(sys.argv[4], "w")
f.write("Filtering log for {} total sequences ={}\n".format(name,len(records)))
f.write("Filename: {}, max_length: {}, max_similarity: {}\n\n".format(sys.argv[1],sys.argv[2],sys.argv[3]))
# f.write("Filename: {}, max_length: {}, max_similarity: {}\n\n".format(argv[1],argv[2],argv[3]))
f.close()
#print(range(len(records)))
boolean_list=[True]*len(records) #index list with booleans
for i in range(len(records)-1):
    try:
        ref_score=aligner.score(records[i].seq, records[i].seq)      
        i_genus=records[i].annotations["organism"].split(" ")[0]
        i_name=records[i].annotations["organism"].split(" ")[1]
        i_species="{}.{}".format(i_genus[0],i_name)
        if "score" in records[i].annotations["source"]:
            trans_flag[i_species]=True
    except IndexError:
        print(records[i].annotations["organism"])
        print("no species name found for {}".format(records[i].id))
        f = open(sys.argv[4], "a")
        f.write("{} no species name found\n".format(records[i].id))
        f.close()
        boolean_list[i]=False
        continue
    except ValueError:
        f = open(sys.argv[4], "a")
        f.write("{} strange character in sequence\n".format(records[i].id))
        f.close()
        boolean_list[i]=False
        continue
    # print(i_species)
    if len(records[i])<int(sys.argv[2]): #check length of record
    # if len(records[i])<int(argv[2]): #check length of record
        i_lng=len(records[i])
        boolean_list[i]=False
        print("{} has been removed with length: {}".format(records[i].id,i_lng))
        f = open(sys.argv[4], "a")
        f.write("{} has been removed with length: {}\n".format(records[i].id,i_lng))
        f.close()
        continue
    
    # elif "sp." in i_species:
    #     boolean_list[i]=False
    #     print("{} has been removed ->sp. in oranism name".format(records[i].id))
    #     f = open(sys.argv[4], "a")
    #     f.write("{} has been removed -> sp. in organism name\n".format(records[i].id))
    #     f.close()
    #     continue
    elif boolean_list[i]==True:
        for j in range(i+1,len(records)):
            if boolean_list[j]==True:
                try:
                    j_genus=records[j].annotations["organism"].split(" ")[0]
                    j_name=records[j].annotations["organism"].split(" ")[1]
                    
                    j_species="{}.{}".format(j_genus[0],j_name)
                except IndexError:
                    print(records[j].annotations["organism"])
                    print("no species name found for {}".format(records[j].id))
                    f = open(sys.argv[4], "a")
                    f.write("{} no species name found\n".format(records[j].id))
                    f.close()
                    boolean_list[j]=False
                    continue
                i_id=records[i].id
                j_id=records[j].id
                # alignment = aligner.align(records[i].seq, records[j].seq)
                #divide align score by query length to receive 1 for exact equal sequences
                try:
                    align_score=aligner.score(records[i].seq, records[j].seq)/ref_score
                except ValueError:
                    f = open(sys.argv[4], "a")
                    f.write("{} strange character in sequence\n".format(records[j].id))
                    f.close()
                    boolean_list[j]=False
                    continue
                # if i_species==j_species:
                    # print(alignment[0])
                    # print(align_score/len(records[i]))
                    # print(i_species,i_id)
                    # print(j_species,j_id)
                    
                    
                min=len(records[i].seq)-20
                max=len(records[i].seq)+20
                #check score and species name before deleting
                if align_score>float(sys.argv[3]) and i_species==j_species and min<len(records[j].seq)<max:
                # if align_score>float(argv[3]) and i_species==j_species and min<len(records[j].seq)<max:
    
                    
                    if "score" in records[j].annotations["source"]: #if there is a transporter
                        try:
                            if trans_flag[j_species]==True:
                                print("{} has been removed by {} with score: {}".format(records[j].id,records[i].id,align_score))
                                f = open(sys.argv[4], "a")
                                f.write("{} has been removed by {} with score: {} (transporter found)\n".format(records[j].id,records[i].id,align_score))
                                f.close()
                                boolean_list[j]=False
                            else:
                                trans_flag[j_species]=True
                        except KeyError:
                            trans_flag[j_species]=True
                            
                    # print(records[j].annotations["source"])
                    elif "score" not in records[j].annotations["source"]: #only remove entry when no transporter found
                        print("{} has been removed by {} with score: {}".format(records[j].id,records[i].id,align_score))
                        f = open(sys.argv[4], "a")
                        f.write("{} has been removed by {} with score: {}\n".format(records[j].id,records[i].id,align_score))
                        f.close()
                        boolean_list[j]=False
                
filtered=list(compress(records, boolean_list))    
seq_removed=len(records)-sum(boolean_list)
print("{} of {} sequences have been removed".format(seq_removed,len(records)))
f= open(sys.argv[4], "a")
f.write("{} sequences have been removed".format(seq_removed))
f.close()
SeqIO.write(filtered, "{}.fasta".format(sys.argv[5]), "fasta")
SeqIO.write(filtered, "{}.gb".format(sys.argv[5]), "genbank")
