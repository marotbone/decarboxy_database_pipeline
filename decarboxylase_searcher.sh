IFS=
sample_name=$1
work_dir=$2
search_term=$3
reference_transporter_path=$4
decarboxylase_name=$5
transporter_name=$6
decarboxylase_sim_threshold=$7
transporter_sim_threshold=$8
decarboxylase_minlength=$9
transporter_minlength=${10}

contig_repository=data/contig_repository

printf "${search_term}\n"

printf "****************************************************************************************************\n"
mkdir  -p ${contig_repository}
printf "DOWNLOAD DECARBOXYLASE AA SEQUENCES for ${sample_name}...\n"
mkdir -p ${work_dir}/decarboxylases/aa_sequences
python python_scripts/download_all_gb_from_ncbiV3.py ${search_term} ${work_dir}/decarboxylases/aa_sequences/${sample_name}.gb

printf "____________________________________________________________________________________________________\n\n"
number_of_decarboxylases="$(grep -c LOCUS ${work_dir}/decarboxylases/aa_sequences/${sample_name}.gb)"
printf "${number_of_decarboxylases} DECARBOXYLASES DOWNLOADED!\n\n"

if [ ${number_of_decarboxylases} -gt 2 ]
then

	printf "CREATE MAPPING TABLE...\n"
	python python_scripts/create_mapping_tableV3.py ${work_dir}/decarboxylases/aa_sequences/${sample_name}.gb ${work_dir}/mapping_table.txt
	printf "____________________________________________________________________________________________________\n\n"


	printf "DOWNLOAD ASSEMBLIES...\n"
	python python_scripts/download_assemblies_maroV2.py ${work_dir}/mapping_table.txt ${work_dir}/assembly_mapping_table.txt ${contig_repository}
	printf "____________________________________________________________________________________________________\n\n"


	printf "SEARCH FOR NEIGHBOURING TRANSPORTERS...\n"
	mkdir -p ${work_dir}/transporters/aa_sequences
	python python_scripts/find_neighbour_ncbiV3.py ${work_dir}/decarboxylases/aa_sequences/${sample_name}.gb ${work_dir}/mapping_table.txt ${contig_repository} ${reference_transporter_path} ${work_dir}/neartrans.txt
	mv ${work_dir}/neighbouring_transporters.* ${work_dir}/transporters/aa_sequences
	printf "____________________________________________________________________________________________________\n\n"


	printf "FILTERING DECARBOXYLASES...\n" 
	python python_scripts/maros_cleanerV7.1.py ${work_dir}/decarboxylases/aa_sequences/${sample_name}.gb ${decarboxylase_minlength} ${decarboxylase_sim_threshold} ${work_dir}/decarboxylases/aa_sequences/filtering_${sample_name}.txt ${work_dir}/decarboxylases/aa_sequences/filtered_${sample_name}
	printf "____________________________________________________________________________________________________\n\n"

	
	printf "GENERATE NICE FASTA HEADER FOR DECARBOXYLASES...\n"
	python python_scripts/generate_tree_fastaV2.py ${work_dir}/decarboxylases/aa_sequences/filtered_${sample_name}.gb 0 ${work_dir}/decarboxylases/aa_sequences/taxa_filtered_${sample_name}.fasta ${decarboxylase_name}
	printf "____________________________________________________________________________________________________\n\n"
	
	
	printf "CONVERT DECARBOXYLASE AA SEQUENCES INTO NUCLEOTIDES...\n"
	mkdir -p ${work_dir}/decarboxylases/nt_sequences
	python python_scripts/get_nucleotides_from_prot_accessionV6.py ${work_dir}/decarboxylases/aa_sequences/filtered_${sample_name}.fasta ${work_dir}/decarboxylases/nt_sequences/nucleotides_filtered_${sample_name}.fasta ${decarboxylase_name}
	printf "____________________________________________________________________________________________________\n\n"

	printf "ALIGN DECARBOXYLASE SEQUENCES..\n"
	clustalo -i ${work_dir}/decarboxylases/aa_sequences/taxa_filtered_${sample_name}.fasta -o ${work_dir}/decarboxylases/aa_sequences/aligned_taxa_filtered_${sample_name}.fasta -t Protein --force	
	printf "____________________________________________________________________________________________________\n\n"

	
	#~ printf "FIND CONSERVED BLOCKS IN ALIGNED DECARBOXYLASES...\n"
	#~ Gblocks  ${work_dir}/decarboxylases/aa_sequences/aligned_taxa_filtered_${sample_name}.fasta -g
	
	#~ printf "____________________________________________________________________________________________________\n\n"
	
	printf "CONVERT INTO .PHYLIP ...\n"
	python python_scripts/fasta_2_phylip.py -i ${work_dir}/decarboxylases/aa_sequences/aligned_taxa_filtered_${sample_name}.fasta -o ${work_dir}/decarboxylases/aa_sequences/aligned_taxa_filtered_${sample_name}.phylip -r
	printf "____________________________________________________________________________________________________\n\n"

	
	printf "CALCULATE DECARBOXYLASE TREE...\n\n"
	phyml-mpi -d aa -m Blosum62 -b -4 -v 0.0 -c 4 -a e -f m --no_memory_check -i ${work_dir}/decarboxylases/aa_sequences/aligned_taxa_filtered_${sample_name}.phylip 
	printf "____________________________________________________________________________________________________\n\n"

	
	number_of_transporters="$(grep -c LOCUS ${work_dir}/transporters/aa_sequences/neighbouring_transporters.gb)"
	printf "${number_of_transporters} $6 TRANSPORTERS FOUND!\n"
	if [ ${number_of_transporters} -gt 2 ]
	then
	
		printf "FILTERING TRANSPORTERS...\n"
		python python_scripts/maros_cleanerV7.1.py ${work_dir}/transporters/aa_sequences/neighbouring_transporters.gb ${transporter_minlength} ${transporter_sim_threshold} ${work_dir}/transporters/aa_sequences/filtering_neighbouring_transporters.txt ${work_dir}/transporters/aa_sequences/filtered_neighbouring_transporters
		printf "____________________________________________________________________________________________________\n\n"

		
		printf "GENERATE NICE FASTA HEADER FOR TRANSPORTERS...\n"
		python python_scripts/generate_tree_fastaV2.py ${work_dir}/transporters/aa_sequences/filtered_neighbouring_transporters.gb 0 ${work_dir}/transporters/aa_sequences/taxa_filtered_neighbouring_transporters.fasta ${transporter_name}
		printf "____________________________________________________________________________________________________\n\n"



		printf "CONVERT TRANSPORTER AA SEQUENCES INTO NUCLEOTIDES..\n"
		mkdir -p ${work_dir}/transporters/nt_sequences
		python python_scripts/get_nucleotides_from_prot_accessionV6.py ${work_dir}/transporters/aa_sequences/filtered_neighbouring_transporters.fasta ${work_dir}/transporters/nt_sequences/nucleotides_filtered_neighbouring_transporters.fasta ${transporter_name}
		printf "____________________________________________________________________________________________________\n\n"

	
		printf "ALIGN TRANSPORTER SEQUENCES..\n"
		clustalo -i ${work_dir}/transporters/aa_sequences/taxa_filtered_neighbouring_transporters.fasta -o ${work_dir}/transporters/aa_sequences/aligned_taxa_filtered_neighbouring_transporters.fasta -t Protein --force
		printf "____________________________________________________________________________________________________\n\n"


		#~ printf "FIND CONSERVED BLOCKS IN ALIGNED TRANSPORTERS...\n"
		#~ Gblocks  ${work_dir}/transporters/aa_sequences/aligned_taxa_filtered_neighbouring_transporters.fasta -g
		#~ printf "____________________________________________________________________________________________________\n\n"


		printf "CONVERT INTO .PHYLIP ...\n"
		python python_scripts/fasta_2_phylip.py -i ${work_dir}/transporters/aa_sequences/aligned_taxa_filtered_neighbouring_transporters.fasta -o ${work_dir}/transporters/aa_sequences/aligned_taxa_filtered_neighbouring_transporters.phylip -r
		printf "____________________________________________________________________________________________________\n\n"

		printf "CALCULATE TRANSPORTER TREE...\n"
		phyml-mpi -d aa -m Blosum62 -b -4 -v 0.0 -c 4 -a e -f m --no_memory_check  -i ${work_dir}/transporters/aa_sequences/aligned_taxa_filtered_neighbouring_transporters.phylip 
		printf "____________________________________________________________________________________________________\n\n"
	else
		printf "NO OR LESS THAN TWO TRANSPORTERS FOUND!\n\n!!!"
		printf "____________________________________________________________________________________________________\n\n"
	fi
else 
	printf "NO DECARBOXYLASES FOUND!\n\n"
	printf "____________________________________________________________________________________________________\n\n"
fi
printf "****************************************************************************************************\n\n"
