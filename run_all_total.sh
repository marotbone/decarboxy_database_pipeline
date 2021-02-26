#!/bin/bash

#SBATCH --mail-user=maro.widmer@students.unibe.ch
#SBATCH --mail-type=ALL
#SBATCH --job-name="dec_search"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#SBATCH --output=dec_search-%j.out
#SBATCH --error=dec_search-%j.error

#~ bash run_all_bacilli_hdc.sh
#~ ##bash run_all_bacilli_ldc.sh
#~ bash run_all_bacilli_odc.sh
#~ bash run_all_bacilli_tdc.sh
#~ bash run_all_burkholderiales_hcd.sh
#~ bash run_all_burkholderiales_ldc.sh
#~ bash run_all_burkholderiales_odc.sh
#~ bash run_all_burkholderiales_tdc.sh
#~ bash run_all_enterobacteriaceae_hdc.sh
#~ bash run_all_enterobacteriaceae_ldc.sh
#~ bash run_all_enterobacteriaceae_odc.sh
#~ bash run_all_enterobacteriaceae_tdc.sh
#~ bash run_all_halomonadacae_morganellaceae_hdc.sh
#~ bash run_all_halomonadacae_morganellaceae_ldc.sh
#~ bash run_all_halomonadaceae_morganellaceae_odc.sh
#~ bash run_all_halomonadaceae_morganellaceae_tdc.sh
#~ bash run_all_micrococcales_hcd.sh
#~ bash run_all_micrococcales_lcd.sh
#~ bash run_all_micrococcales_odc.sh
#~ bash run_all_micrococcales_tdc.sh
#~ ##bash run_all_pseudomonadales_hcd.sh
#~ bash run_all_pseudomonadales_ldc.sh
#~ bash run_all_pseudomonadales_odc.sh
#~ bash run_all_pseudomonadales_tdc.sh
bash run_all_clostridia_ldc.sh
bash run_all_clostridia_tdc.sh
bash run_all_clostridia_hdc.sh
bash run_all_clostridia_odc.sh


#concatenating all nucleotide sequences
cat data/all_*/transporters/nt_sequences/*.fasta >> data/all_nt_sequences/all_transporters_nt_sequences.fasta
cat data/all_*/decarboxylases/nt_sequences/*.fasta >> data/all_nt_sequences/all_decarboxylases_nt_sequences.fasta

#filter for similar sequences
python data/all_nt_sequences/sequence_cleanerV2.py all_transporters_nt_sequences.fasta
python data/all_nt_sequences/sequence_cleanerV2.py all_decarboxylases_nt_sequences.fasta

