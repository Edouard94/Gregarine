from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Input files
barrnap_gff_file = "barrnap_cleaned_no_fasta.gff3"
trnascan_gff_file = "trnascan_cleaned.gff3"
genbank_file = "Leidyana_gryllorum_copy.gbk"
output_file = "Leidyana_gryllorum_updated.gbk"

# Parse the existing GenBank file
with open(genbank_file, "r") as gbk_in:
    records = {record.id: record for record in SeqIO.parse(gbk_in, "genbank")}

# Function to parse GFF files and extract features
def parse_gff(gff_file, feature_type):
    features = {}
    with open(gff_file, "r") as gff_in:
        for line in gff_in:
            if line.startswith("#") or line.strip() == "":
                continue  # Skip comments and empty lines
            fields = line.strip().split("\t")
            contig, source, ftype, start, end, score, strand, phase, attributes = fields
            if ftype != feature_type:
                continue  # Only process the specified feature type
            qualifiers = {}
            for attribute in attributes.split(";"):
                if "=" in attribute:
                    key, value = attribute.split("=", 1)
                    qualifiers[key] = value
            strand = 1 if strand == "+" else -1
            feature = SeqFeature(
                location=FeatureLocation(int(start) - 1, int(end), strand=strand),
                type=ftype,
                qualifiers=qualifiers,
            )
            if contig not in features:
                features[contig] = []
            features[contig].append(feature)
    return features

# Parse the barrnap and tRNAscan GFF files
barrnap_features = parse_gff(barrnap_gff_file, "rRNA")
trnascan_features = parse_gff(trnascan_gff_file, "tRNA")

# Add features to the corresponding GenBank records
for contig, record in records.items():
    if contig in barrnap_features:
        record.features.extend(barrnap_features[contig])
    if contig in trnascan_features:
        record.features.extend(trnascan_features[contig])

# Write the updated GenBank file
with open(output_file, "w") as gbk_out:
    SeqIO.write(records.values(), gbk_out, "genbank")

print(f"Updated GenBank file written to {output_file}")
