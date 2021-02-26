#!/bin/bash

#SBATCH --mail-user=maro.widmer@students.unibe.ch
#SBATCH --mail-type=ALL
#SBATCH --job-name="dec_search"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --mem=16G
#SBATCH --output=dec_search-%j.out
#SBATCH --error=dec_search-%j.error

IFS=
sample_name=all_enterobacteriaceae_ornithine_decarboxylases
work_dir=data/${sample_name}
search_term="((ornithine decarboxylase[protein] AND enterobacteriaceae[organism] NOT Klebsiella pneumoniae[organism] NOT Escherichia coli[organism] NOT Salmonella enterica[organism] NOT Enterobacter hormaechei[organism] NOT Shigella sonnei[organism])) OR RLJ10911.1 OR SSN07774 OR STM40075.1 OR RDT55966.1 OR WP_141011571.1 OR GCE72201.1 OR TNX72561.1 OR CBK86600.1 OR ALA01421.1 OR TJX75091 OR RUK92436.1"
reference_transporter_path=transporter_references/potE_reference_transporter_L.Saerimneri.fasta
decarboxylase_name=ODC
transporter_name=potE
decarboxylase_sim_threshold=0.8
transporter_sim_threshold=0.8
decarboxylase_minlength=500
transporter_minlength=200


bash decarboxylase_searcher.sh ${sample_name} ${work_dir} ${search_term} ${reference_transporter_path} ${decarboxylase_name} ${transporter_name} ${decarboxylase_sim_threshold} ${transporter_sim_threshold} ${decarboxylase_minlength} ${transporter_minlength} 2>&1| tee -a ${sample_name}.log

