#!/bin/bash

# Unzipping, checking fastq files qualities, quality trimming, dereplicating, and host removing
dedup_dir="deduplicated_reads"
non_host_dir="non_host_reads"
non_rRNA_dir="non_rRNA_reads"
mkdir -p "${dedup_dir}"
mkdir -p "${non_host_dir}"
mkdir -p "${non_rRNA_dir}"

for dir in R*_*; do # Loop through all directories starting with R and ending with an underscore
    for fq in "${dir}"/*.fq; do # Loop through all fastq files in the directory

        # Quality trimming and processing of the paired-end reads

        if [[ "${fq}" =~ _L[0-9]+_1\.fq$ ]]; then # Check if the file is a forward read
            forward_read="${fq}"
            reverse_read="${forward_read/_1.fq/_2.fq}" # Define the reverse read file path
            # Check for inconsistency between forward and reverse reads
            if [[ ! -f "${reverse_read}" ]]; then
                echo "Error: Reverse read file not found for ${forward_read}"
                continue
            fi
            # Define output file paths
            forward_fastp_read_output="${dir}/$(basename "${forward_read}" .fq)_fastp.fq"
            reverse_fastp_read_output="${dir}/$(basename "${reverse_read}" .fq)_fastp.fq"

            # Quality trimming of the paired-end reads
            fastp -i "${fq}" -I "${reverse_read}" -o "${forward_fastp_read_output}" -O "${reverse_fastp_read_output}" \
            -h "${dir}/$(basename "${fq}" .fq)_fastp.html" -j "${dir}/$(basename "${fq}" .fq)_fastp.json"

            # Dereplicating the reads
            ## Run vsearch to merge the paired-end reads
            ### Define output file path
            merged_fastq="${dir}/$(basename "${fq}" .fq)_merged.fq"
            vsearch --fastq_mergepairs "${forward_fastp_read_output}" --reverse "${reverse_fastp_read_output}" --fastqout "${merged_fastq}"
            ## Run czid-dedup for deduplication. This command requires identical matches to be considered duplicates.
            ### Define output file path
            dedup_fastq="${dedup_dir}/$(basename "${fq}" .fq)_dedup.fq"
            cd-hit-dup -i "${merged_fastq}" -o "${dedup_fastq}"

            # Remove vector and host contamination
            base_name=$(basename "${fq}" .fq)
            ## Remove vector contamination using UniVec database
            ### Define output file path and index
            non_vector_reads="${non_host_dir}/${base_name}_non_vector.fq"
            UniVec_index="/data2/MetaPro_databases/univec_core/bowtie/UniVec_Core"
            bowtie2 -U "${dedup_fastq}" -x "${UniVec_index}" --un "${non_vector_reads}" -S /dev/null -p 16 --very-sensitive-local
            ## Remove host reads using Bowtie2
            ### Define output file path and index
            non_host_reads="${non_host_dir}/${base_name}_non_host.fq"
            host_index="/data2/Cricket_Greg_MetaTranscriptomic/X204SC24030899-Z01-F001/Adomesticus_genome/ado_genome_index"
            ### Save host reads using Bowtie2
            host_reads="${non_host_dir}/${base_name}_host.fq"
            bowtie2 -U "${non_vector_reads}" -x "${host_index}" --al "${host_reads}" --un "${non_host_reads}" -S /dev/null -p 16 --very-sensitive-local

            # Remove rRNA contamination
            ## Define output file path and index
            non_rRNA_reads="${non_rRNA_dir}/${base_name}_non_rRNA.fq"
            rRNA_reads="${non_rRNA_dir}/${base_name}_rRNA.fq"
            db_path="/home/edouard/sortmerna/rRNA_databases_v4/smr_v4.3_sensitive_db.fasta"

            # Define the base directory for SortMeRNA's working directories
            base_wdir="/data2/Cricket_Greg_MetaTranscriptomic/X204SC24030899-Z01-F001/01.RawData/sortmerna_wdir"

            # Create a new working directory for this run
            wdir="${base_wdir}/${base_name}"
            mkdir -p "${wdir}"

            ## Run SortMeRNA
            sortmerna --ref "${db_path}" --reads "${non_host_reads}" --other "${non_rRNA_reads}" --aligned "${rRNA_reads}" --workdir "${wdir}" --threads 16 --fastx
        fi
    done
done
echo "Processing complete."

# Kaiju

## Define the directory containing your input files
input_dir="non_rRNA_reads"
## Define the path to your Kaiju database
db_path="/data2/kaiju_db_refseq_nr/kaiju_db_refseq_nr.fmi"
## Define the path to your nodes file
nodes_path="/data2/kaiju_db_refseq_nr/nodes.dmp"
## Create an array to store the input file names
input_files=()
## Loop over each non_rRNA.fq.fq file in the input directory
for input_file in "${input_dir}"/*_non_rRNA.fq.fq; do
    # Add the file name to the input_files array
    input_files+=("${input_file}")
done
## Convert the array to a comma-separated string
input_files_str=$(IFS=','; echo "${input_files[*]}")
## Run kaiju-multi on all input files
kaiju-multi -v -z 19 -t "${nodes_path}" -f "${db_path}" -i "${input_files_str}" > all_samples.out