#!/bin/bash

# Define the directory containing your input files
input_dir="non_rRNA_reads"

# Define the path to your Kaiju database
db_path="./kaiju_db_nr_euk.fmi"

# Define the path to your nodes file
nodes_path="./nodes.dmp"
names_path=".names.dmp"

# Create an array to store the input file names
input_files=()

# Loop over each non_rRNA.fq.fq file in the input directory
for input_file in "${input_dir}"/*_non_rRNA.fq.fq; do
    # Add the file name to the input_files array
    input_files+=("${input_file}")
done

# Convert the array to a comma-separated string
input_files_str=$(IFS=','; echo "${input_files[*]}")

# Run kaiju-multi on all input files
kaiju-multi -v -z 20 -t "${nodes_path}" -f "${db_path}" -i "${input_files_str}" > all_samples.out
kaiju2krona -t "${nodes_path}" -n "${names_path}" -i all_samples.out -o all_samples.out.krona
kaiju2table -t "${nodes_path}" -n "${names_path}" -r species -l superkingdom,phylum,class,order,family,genus,species -o all_samples.tsv all_samples.out
ktImportText -o all_sample_prg.html all_samples.out.krona