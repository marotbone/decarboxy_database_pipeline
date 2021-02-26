#!/bin/bash


IFS=
sample_name=all_micrococcales_lysine_decarboxylases
work_dir=data/${sample_name}
search_term="Micrococcales[Organism] AND lysine decarboxylase[Protein Name] "
reference_transporter_path=transporter_references/cadB_reference_transporter_E.coli.fasta
decarboxylase_name=LDC
transporter_name=cadB
decarboxylase_sim_threshold=0.8
transporter_sim_threshold=0.8
decarboxylase_minlength="a"
transporter_minlength=200


bash decarboxylase_searcher.sh ${sample_name} ${work_dir} ${search_term} ${reference_transporter_path} ${decarboxylase_name} ${transporter_name} ${decarboxylase_sim_threshold} ${transporter_sim_threshold} ${decarboxylase_minlength} ${transporter_minlength} 2>&1| tee -a ${sample_name}.log

