# Import necessary modules
import re
from collections import defaultdict

# Define the path to your GFF3 file
gff_path = 'Leidyana_gryllorum_copy.gff3'

# Dictionary to store GO terms for each gene
gene_go_terms = defaultdict(list)

# Set to store all gene IDs
all_gene_ids = set()

# Open the GFF3 file and read the lines
with open(gff_path, 'r') as gff_file:
    for line in gff_file:
        # Check if the line contains a gene ID
        if 'ID=' in line and 'gene' in line:
            # Extract the gene ID using regex
            id_match = re.search(r'ID=([^;]+)', line)
            if id_match:
                gene_id = id_match.group(1)
                all_gene_ids.add(gene_id)
                
        # Check if the line contains GO terms
        if 'Ontology_term=' in line:
            # Extract the gene ID and GO terms using regex
            match = re.search(r'ID=([^;]+);.*Ontology_term=([^;]+)', line)
            if match:
                gene_id = match.group(1)
                go_terms = match.group(2).split(',')
                gene_go_terms[gene_id].extend(go_terms)

# Open the output file and write the GO terms in the desired format
with open('go_terms_for_wego.txt', 'w') as output_file:
    for gene_id, go_terms in gene_go_terms.items():
        for go_term in go_terms:
            output_file.write(f"{gene_id}\t{go_term}\n")

# Calculate metrics
total_genes = len(all_gene_ids)
genes_with_go_terms = len(gene_go_terms)
genes_without_go_terms = total_genes - genes_with_go_terms

# Print metrics
print(f"Total number of genes: {total_genes}")
print(f"Number of genes with GO terms: {genes_with_go_terms}")
print(f"Number of genes without GO terms: {genes_without_go_terms}")
