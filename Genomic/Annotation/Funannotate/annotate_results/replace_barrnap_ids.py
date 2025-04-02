import sys

# Input files
barrnap_gff = sys.argv[1]  # Input GFF file from barrnap
mapping_file = sys.argv[2]  # Mapping file containing contig-to-locus tag mappings
output_gff = sys.argv[3]  # Output GFF file with updated IDs

# Load the mapping
id_mapping = {}  # Dictionary to store contig-to-locus tag mappings
with open(mapping_file, 'r') as mapping:
    for line in mapping:
        contig, acnh_id = line.strip().split()
        id_mapping[contig] = acnh_id

# Replace IDs in the barrnap GFF file
with open(barrnap_gff, 'r') as infile, open(output_gff, 'w') as outfile:
    for line in infile:
        if line.startswith('#'):
            outfile.write(line)
        else:
            fields = line.strip().split('\t')
            contig = fields[0]
            if contig in id_mapping:
                fields[0] = contig  # Keep the contig name in the first column
                fields[8] += f";Locus_tag={id_mapping[contig]}"  # Add locus tag to the last column
            outfile.write('\t'.join(fields) + '\n')
