#!/bin/bash

IFS=
sample_name=all_entero_histidine_decarboxylases
work_dir=data/${sample_name}
search_term="Enterobacteriaceae[Orgn] AND histidine decarboxylase[Protein Name]"
reference_transporter_path=transporter_references/hdcP_reference_transporter_histidine_histamin_enterobacter.fasta
decarboxylase_name=HDC
transporter_name=hdcP
decarboxylase_sim_threshold=0.8
transporter_sim_threshold=0.8
decarboxylase_minlength="a"
transporter_minlength=200


bash decarboxylase_searcher.sh ${sample_name} ${work_dir} ${search_term} ${reference_transporter_path} ${decarboxylase_name} ${transporter_name} ${decarboxylase_sim_threshold} ${transporter_sim_threshold} ${decarboxylase_minlength} ${transporter_minlength} 2>&1| tee -a ${sample_name}.log

