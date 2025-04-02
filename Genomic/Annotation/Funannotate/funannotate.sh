#!/bin/bash

# Gregarine Genome Analysis Pipeline
# Author: Edouard Bessette
# Date: April 2025
# Description: A reproducible pipeline for genome annotation and comparative analysis of *Leidyana gryllorum*.

set -e  # Exit on error
set -o pipefail  # Exit on pipe failure

# Define paths and parameters
GENOME="input_genome.fasta"
MASKED_GENOME="masked_genome.fasta"
ANNOTATION_DIR="annotation_results"
THREADS=20
SPECIES="Leidyana_gryllorum"
RNA_FASTQ="../putative_gregarineRNA_mapped.fastq"
TRINITY_TRANSCRIPTS="/data2/Cricket_Greg_MetaTranscriptomic/X204SC24030899-Z01-F001/01.RawData/LgryllorumAlign/putative_gregarineRNA_trinity_transcripts/Trinity-GG.nr.fasta"
FUNANNOTATE_DB="/data2/funannotate_db"
EGGNOG_ANNOTATIONS="Lgryllorum.emapper.annotations"
IPRSCAN_XML="interproscan.xml"

mkdir -p $ANNOTATION_DIR

# Step 1: Mask the genome
echo "Masking the genome..."
windowmasker -mk_counts -infmt fasta -sformat oascii -in $GENOME -out genome_counts.txt
windowmasker -ustat genome_counts.txt -infmt fasta -in $GENOME -outfmt 'fasta' -dust true -dust_level 20 -out $MASKED_GENOME

# Step 2: Train Funannotate models
echo "Training Funannotate models..."
conda activate funannotate
nohup funannotate train \
    -i $MASKED_GENOME \
    -o $ANNOTATION_DIR \
    -s $RNA_FASTQ \
    --trinity $TRINITY_TRANSCRIPTS \
    --species "$SPECIES" \
    --max_intronlen 150 \
    --cpus $THREADS > funannotate_train.log 2>&1 &

# Step 3: Predict genes with Funannotate
echo "Predicting genes with Funannotate..."
nohup funannotate predict \
    -i $MASKED_GENOME \
    -o $ANNOTATION_DIR \
    -s "$SPECIES" \
    --busco_db protists \
    --busco_seed_species toxoplasma \
    --optimize_augustus \
    --organism other \
    --min_intronlen 10 \
    --max_intronlen 150 \
    --protein_evidence /home/edouard/GregarineStuff/Lgryllorum_Genome_analysis/Funannotate/GCF_019968955.1/protein.faa \
    /home/edouard/GregarineStuff/Lgryllorum_Genome_analysis/Funannotate/GCF_000223845.1/protein.faa \
    /data2/SwissProt/uniprot_sprot.fasta \
    --trnascan ../tRNAscan.out \
    --cpus $THREADS --header_length 20 > funannotate_predict.log 2>&1 &

# Step 4: Update annotation with RNAseq data
echo "Updating annotation with RNAseq data..."
nohup funannotate update \
    -i $ANNOTATION_DIR \
    --cpus $THREADS \
    --max_intronlen 150 \
    --species "$SPECIES" \
    --no_trimmomatic \
    --no_normalize_reads \
    -s $RNA_FASTQ > funannotate_update.log 2>&1 &

# Step 5: Functional annotation
echo "Running functional annotation..."
nohup funannotate annotate \
    -i $ANNOTATION_DIR \
    --force \
    --header_length 25 \
    --sbt ~/GregarineStuff/Lgryllorum_Genome_analysis/Funannotate/template.sbt \
    --species "$SPECIES" \
    --out funannotate_annotate/ \
    --database $FUNANNOTATE_DB \
    --busco_db protists \
    --eggnog $EGGNOG_ANNOTATIONS \
    --iprscan $IPRSCAN_XML \
    --rename ACNH67_ \
    --cpus 5 > funannotate_annotate.log 2>&1 &

# Step 6: Evaluate annotation with BUSCO
echo "Evaluating annotation with BUSCO..."
conda activate busco
busco -i funannotate_annotate/predicted_proteins.fasta \
    -l apicomplexa_odb12 \
    -o busco \
    -m prot \
    -f

# Step 7: Visualize annotation with JBrowse
echo "Visualizing annotation with JBrowse..."
# Instructions for JBrowse visualization can be added here if needed.

# Step 8: KEGGaNOG visualization
echo "Running KEGGaNOG visualization..."
conda activate kegganog
cat > files.txt << EOF
Run_v3/Lgryllorum.emapper.annotations
GCF_000223845.1/Gni_emapper.annotations
EOF
KEGGaNOG -i files.txt -M -o kegganog_output -dpi 900 -g

# Step 9: Extract GO terms for WEGO
echo "Extracting GO terms for WEGO..."
python <<EOF
import re
from collections import defaultdict

gff_path = "$ANNOTATION_DIR/predicted_genes.gff3"
output_file = "go_terms_for_wego.txt"

gene_go_terms = defaultdict(list)
with open(gff_path, 'r') as gff_file:
    for line in gff_file:
        if 'Ontology_term=' in line:
            match = re.search(r'ID=([^;]+);.*Ontology_term=([^;]+)', line)
            if match:
                gene_id = match.group(1)
                go_terms = match.group(2).split(',')
                gene_go_terms[gene_id].extend(go_terms)

with open(output_file, 'w') as out:
    for gene_id, go_terms in gene_go_terms.items():
        for go_term in go_terms:
            out.write(f"{gene_id}\t{go_term}\n")
EOF

echo "Pipeline completed successfully!"
