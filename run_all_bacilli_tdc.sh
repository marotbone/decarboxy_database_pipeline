#!/bin/bash

IFS=
sample_name=all_bacilli_tyrosine_decarboxylases
work_dir=data/${sample_name}
search_term="Bacilli [Organism] AND tyrosine decarboxylase[Protein Name] "
reference_transporter_path=transporter_references/tyrP_reference_transporter_l.brevis.fasta
decarboxylase_name=TDC
transporter_name=tyrP
decarboxylase_sim_threshold=0.8
transporter_sim_threshold=0.8
decarboxylase_minlength="a"
transporter_minlength=200


bash decarboxylase_searcher.sh ${sample_name} ${work_dir} ${search_term} ${reference_transporter_path} ${decarboxylase_name} ${transporter_name} ${decarboxylase_sim_threshold} ${transporter_sim_threshold} ${decarboxylase_minlength} ${transporter_minlength} 2>&1| tee -a ${sample_name}.log

