# BLASTP Comparison between L. gryllorum and G. niphandrodes
#! /bin/bash
## Set paths to protein files
lgryllorum=~/GregarineStuff/Lgryllorum_Genome_analysis/Funannotate/Run_v2/species_train/annotate_results/Leidyana_gryllorum.proteins.fa
gniphandrodes=~/GregarineStuff/Lgryllorum_Genome_analysis/Funannotate/GCF_000223845.1/protein.faa
output_dir=Blast_Compare_Gni
nr_diamond=/data2/nr_diamond_DB/nr_diamond.dmnd

## Create output directory if it doesn't exist
mkdir -p $output_dir
cd $output_dir

## Create DIAMOND databases for both proteomes
diamond makedb --in $lgryllorum -d leidyana_gryllorum_dmd
diamond makedb --in $gniphandrodes -d gniphandrodes_dmd

## Run DIAMOND BLASTP in both directions
diamond blastp -d leidyana_gryllorum_dmd -q $gniphandrodes -o gregarine_vs_gniphandrodes.tab \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    --threads 20

diamond blastp -d gniphandrodes_dmd -q $lgryllorum -o gniphandrodes_vs_gregarine.tab \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    --threads 20

## Compare against NR database
diamond blastp -d $nr_diamond -q $lgryllorum -o lgryllorum_nr.tab \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    --threads 20

diamond blastp -d $nr_diamond -q $gniphandrodes -o gniphandrodes_nr.tab \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    --threads 20

## Extract top hits for comparison
cat > top_hit_extract.py << 'EOF'
import pandas as pd

# Define column names for the DIAMOND output
column_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle']

# Load and process L. gryllorum DIAMOND output
print("Processing L. gryllorum BLASTP results...")
lgryllorum_output = pd.read_csv('lgryllorum_nr.tab', sep='\t', header=None, names=column_names)
lgryllorum_top_hits = lgryllorum_output.sort_values(['qseqid', 'bitscore'], ascending=[True, False]).drop_duplicates(subset='qseqid', keep='first')
lgryllorum_top_hits.to_csv('top_hits_lgryllorum_nr.tab', sep='\t', index=False, header=False)

# Load and process G. niphandrodes DIAMOND output
print("Processing G. niphandrodes BLASTP results...")
gniphandrodes_output = pd.read_csv('gniphandrodes_nr.tab', sep='\t', header=None, names=column_names)
gniphandrodes_top_hits = gniphandrodes_output.sort_values(['qseqid', 'bitscore'], ascending=[True, False]).drop_duplicates(subset='qseqid', keep='first')
gniphandrodes_top_hits.to_csv('top_hits_gniphandrodes_nr.tab', sep='\t', index=False, header=False)

print("Top hits extraction completed for both species.")
EOF
python top_hit_extract.py

## Compare BLASTP results between both species
## compare_blast.py
cat > compare_blast.py << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import numpy as np

# Define column names for the BLASTp output
column_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle']

# Load BLASTp results with specified column names
print("Loading BLASTP results...")
lgryllorum_nr = pd.read_csv('top_hits_lgryllorum_nr.tab', sep='\t', header=None, names=column_names)
gniphandrodes_nr = pd.read_csv('top_hits_gniphandrodes_nr.tab', sep='\t', header=None, names=column_names)

# Filter results based on e-value
filtered_lgryllorum_nr = lgryllorum_nr[lgryllorum_nr['evalue'] < 1e-5]
filtered_gniphandrodes_nr = gniphandrodes_nr[gniphandrodes_nr['evalue'] < 1e-5]

# Calculate the number of shared and unique proteins based on sseqid
shared_proteins = set(filtered_lgryllorum_nr['sseqid']) & set(filtered_gniphandrodes_nr['sseqid'])
unique_to_lgryllorum = set(filtered_lgryllorum_nr['sseqid']) - shared_proteins
unique_to_gniphandrodes = set(filtered_gniphandrodes_nr['sseqid']) - shared_proteins

print(f"Shared Proteins: {len(shared_proteins)}")
print(f"Unique to L. gryllorum: {len(unique_to_lgryllorum)}")
print(f"Unique to G. niphandrodes: {len(unique_to_gniphandrodes)}")

# Calculate proteins without hits
proteins_without_hits_lgryllorum = len(lgryllorum_nr) - len(filtered_lgryllorum_nr)
proteins_without_hits_gniphandrodes = len(gniphandrodes_nr) - len(filtered_gniphandrodes_nr)

# Calculate percentages of proteins without hits
percentage_without_hits_lgryllorum = (proteins_without_hits_lgryllorum / len(lgryllorum_nr)) * 100
percentage_without_hits_gniphandrodes = (proteins_without_hits_gniphandrodes / len(gniphandrodes_nr)) * 100

# Create Venn diagram
plt.figure(figsize=(10, 6))
venn = venn2(
    subsets=(len(unique_to_lgryllorum), len(unique_to_gniphandrodes), len(shared_proteins)),
    set_labels=('L. gryllorum', 'G. niphandrodes')
)
plt.title('Shared and Unique Proteins Between Gregarine Species')
plt.savefig('gregarine_proteome_comparison.png', dpi=300)
print("Venn diagram saved as 'gregarine_proteome_comparison.png'")

# Generate summary report
with open('proteome_comparison_summary.txt', 'w') as report:
    report.write("PROTEOME COMPARISON SUMMARY\n")
    report.write("===========================\n\n")
    report.write(f"Total proteins analyzed from L. gryllorum: {len(lgryllorum_nr)}\n")
    report.write(f"Total proteins analyzed from G. niphandrodes: {len(gniphandrodes_nr)}\n\n")
    report.write(f"Shared proteins: {len(shared_proteins)} ({len(shared_proteins)/len(filtered_lgryllorum_nr)*100:.1f}% of L. gryllorum, {len(shared_proteins)/len(filtered_gniphandrodes_nr)*100:.1f}% of G. niphandrodes)\n")
    report.write(f"Unique to L. gryllorum: {len(unique_to_lgryllorum)} ({len(unique_to_lgryllorum)/len(filtered_lgryllorum_nr)*100:.1f}%)\n")
    report.write(f"Unique to G. niphandrodes: {len(unique_to_gniphandrodes)} ({len(unique_to_gniphandrodes)/len(filtered_gniphandrodes_nr)*100:.1f}%)\n\n")
    report.write(f"Proteins without hits in L. gryllorum: {proteins_without_hits_lgryllorum} ({percentage_without_hits_lgryllorum:.1f}%)\n")
    report.write(f"Proteins without hits in G. niphandrodes: {proteins_without_hits_gniphandrodes} ({percentage_without_hits_gniphandrodes:.1f}%)\n")

print("Summary report generated as 'proteome_comparison_summary.txt'")

# Fix e-values to avoid -inf or invalid values in log10
#filtered_lgryllorum_nr['evalue'] = filtered_lgryllorum_nr['evalue'].replace(0, 1e-300)
#filtered_gniphandrodes_nr['evalue'] = filtered_gniphandrodes_nr['evalue'].replace(0, 1e-300)

# Plot e-value distribution
#plt.figure(figsize=(12, 6))

# L. gryllorum e-value distribution
#plt.subplot(1, 2, 1)
#plt.hist(np.log10(filtered_lgryllorum_nr['evalue']), bins=20, alpha=0.7, color='blue')
#plt.title('L. gryllorum E-value Distribution (log10)')
#plt.xlabel('log10(E-value)')
#plt.ylabel('Count')

# G. niphandrodes e-value distribution
#plt.subplot(1, 2, 2)
#plt.hist(np.log10(filtered_gniphandrodes_nr['evalue']), bins=20, alpha=0.7, color='green')
#plt.title('G. niphandrodes E-value Distribution (log10)')
#plt.xlabel('log10(E-value)')
#plt.ylabel('Count')

#plt.tight_layout()
#plt.savefig('evalue_distribution.png', dpi=300)
#print("E-value distribution plot saved as 'evalue_distribution.png'")
EOF
python compare_blast.py

## Create a plot of the top protein frequencies for each species
## protein_frequency_plot.py

    cat > protein_frequency_plot.py << 'EOF'
    import pandas as pd
    import matplotlib.pyplot as plt
    import re
    import numpy as np
    from collections import Counter

    # Define column names for the BLAST output
    column_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle']

    # Load the top hits files
    print("Loading BLAST results...")
    lgryllorum_hits = pd.read_csv('top_hits_lgryllorum_nr.tab', sep='\t', header=None, names=column_names)
    gniphandrodes_hits = pd.read_csv('top_hits_gniphandrodes_nr.tab', sep='\t', header=None, names=column_names)

    # Filter by e-value significance
    lgryllorum_hits = lgryllorum_hits[lgryllorum_hits['evalue'] < 1e-5]
    gniphandrodes_hits = gniphandrodes_hits[gniphandrodes_hits['evalue'] < 1e-5]

    def clean_protein_name(name):
        """Extract the main protein name, removing species info and database identifiers"""
        # Remove anything in square brackets (typically species names)
        name = re.sub(r'\[.*?\]', '', name)
        
        # Remove database identifiers (like gi|, ref|, etc.)
        name = re.sub(r'\w+\|\S+\|', '', name)
        
        # Remove common prefixes
        name = re.sub(r'^(hypothetical protein|putative|predicted|probable) ', '', name)
        
        # Remove trailing spaces, periods, etc.
        name = name.strip(' .,;:')
        
        # Truncate very long names
        if len(name) > 50:
            name = name[:47] + "..."
            
        return name

    # Clean and count protein names
    lgryllorum_proteins = lgryllorum_hits['stitle'].apply(clean_protein_name)
    gniphandrodes_proteins = gniphandrodes_hits['stitle'].apply(clean_protein_name)

    lgryllorum_counts = Counter(lgryllorum_proteins)
    gniphandrodes_counts = Counter(gniphandrodes_proteins)

    # Get the top 15 most frequent proteins for each species
    top_lgryllorum = dict(lgryllorum_counts.most_common(15))
    top_gniphandrodes = dict(gniphandrodes_counts.most_common(15))

    # Create the plots
    plt.figure(figsize=(20, 10))

    # L. gryllorum plot
    plt.subplot(1, 2, 1)
    bars1 = plt.barh(range(len(top_lgryllorum)), list(top_lgryllorum.values()), color='skyblue')
    plt.yticks(range(len(top_lgryllorum)), list(top_lgryllorum.keys()), fontsize=9)
    plt.xlabel('Frequency')
    plt.title('Top 15 Most Frequent Proteins in L. gryllorum')
    plt.gca().invert_yaxis()  # Highest frequency at the top

    # Add count labels to bars
    for bar in bars1:
        width = bar.get_width()
        plt.text(width + 0.5, bar.get_y() + bar.get_height()/2, 
                f'{width}', ha='left', va='center')

    # G. niphandrodes plot
    plt.subplot(1, 2, 2)
    bars2 = plt.barh(range(len(top_gniphandrodes)), list(top_gniphandrodes.values()), color='lightgreen')
    plt.yticks(range(len(top_gniphandrodes)), list(top_gniphandrodes.keys()), fontsize=9)
    plt.xlabel('Frequency')
    plt.title('Top 15 Most Frequent Proteins in G. niphandrodes')
    plt.gca().invert_yaxis()  # Highest frequency at the top

    # Add count labels to bars
    for bar in bars2:
        width = bar.get_width()
        plt.text(width + 0.5, bar.get_y() + bar.get_height()/2, 
                f'{width}', ha='left', va='center')

    plt.tight_layout()
    plt.savefig('top_protein_frequencies.png', dpi=300, bbox_inches='tight')
    print("Protein frequency plot saved as 'top_protein_frequencies.png'")

    # Also create a table of all protein frequencies for reference
    all_proteins_df = pd.DataFrame({
        'Protein': list(lgryllorum_counts.keys()),
        'L_gryllorum_count': [lgryllorum_counts[p] for p in lgryllorum_counts.keys()]
    }).sort_values('L_gryllorum_count', ascending=False)

    all_proteins_df.to_csv('lgryllorum_protein_frequencies.tsv', sep='\t', index=False)

    all_proteins_df = pd.DataFrame({
        'Protein': list(gniphandrodes_counts.keys()),
        'G_niphandrodes_count': [gniphandrodes_counts[p] for p in gniphandrodes_counts.keys()]
    }).sort_values('G_niphandrodes_count', ascending=False)

    all_proteins_df.to_csv('gniphandrodes_protein_frequencies.tsv', sep='\t', index=False)
    print("Protein frequency tables saved as TSV files")
EOF
python protein_frequency_plot.py

## Create a plot of protein frequencies for each function
cat > function_bar_plot.py << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt
import re
import seaborn as sns
import numpy as np

# Define column names for the BLAST output
column_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle']

# Load the top hits files
print("Loading BLAST results...")
lgryllorum_hits = pd.read_csv('top_hits_lgryllorum_nr.tab', sep='\t', header=None, names=column_names)
gniphandrodes_hits = pd.read_csv('top_hits_gniphandrodes_nr.tab', sep='\t', header=None, names=column_names)

# Function to extract functional categories from titles
def extract_category(title):
    # Remove species and identifiers
    title = re.sub(r'\[.*?\]', '', title)
    title = re.sub(r'\w+\|\S+\|', '', title)

    # Define patterns for functional categories
    categories = {
        'metabolism': r'dehydrogenase|synthetase|kinase|metabolism|metabolic|reductase|transporter',
        'signaling': r'receptor|signal|transduction|kinase|phosphatase',
        'transcription': r'transcription|transcriptional|DNA-binding|zinc finger|helicase',
        'translation': r'translation|ribosom|elongation|initiation|tRNA',
        'cell cycle': r'cell cycle|mitosis|meiosis|division|proliferation',
        'structure': r'cytoskeleton|actin|tubulin|microtubule|flagell',
        'stress response': r'stress|heat shock|chaperone|antioxidant',
        'transport': r'transport|channel|pump|carrier|exchanger',
        'protein processing': r'protease|peptidase|proteasome|ubiquitin',
        'unknown': r'hypothetical|uncharacterized|unknown|predicted'
    }

    # Find the first matching category
    for category, pattern in categories.items():
        if re.search(pattern, title, re.IGNORECASE):
            return category

    # Default category if no match found
    return 'other'

# Handle missing values in the 'stitle' column
lgryllorum_hits['stitle'] = lgryllorum_hits['stitle'].fillna('unknown')
gniphandrodes_hits['stitle'] = gniphandrodes_hits['stitle'].fillna('unknown')

# Extract categories
lgryllorum_hits['category'] = lgryllorum_hits['stitle'].apply(extract_category)
gniphandrodes_hits['category'] = gniphandrodes_hits['stitle'].apply(extract_category)

# Count categories
lgry_category_counts = lgryllorum_hits.groupby('category').size().reset_index(name='count')
lgry_category_counts['species'] = 'L. gryllorum'

gnip_category_counts = gniphandrodes_hits.groupby('category').size().reset_index(name='count')
gnip_category_counts['species'] = 'G. niphandrodes'

# Combine data
combined = pd.concat([lgry_category_counts, gnip_category_counts])

# Create a pivot table to make sure values align correctly with the plot
pivot_data = combined.pivot(index='category', columns='species', values='count').reset_index()
pivot_data = pivot_data.fillna(0)  # Fill NaN values with 0

# Define a consistent order for categories
categories = sorted(combined['category'].unique())

# Create bar plot with a consistent order
plt.figure(figsize=(14, 10))
ax = plt.subplot(111)

# Get color map
colors = plt.cm.viridis(np.linspace(0, 1, len(categories)))

# Plot each category as a separate bar with explicit positions
x_pos = np.arange(2)  # Two species
bar_width = 0.8 / len(categories)
for i, category in enumerate(categories):
    category_data = combined[combined['category'] == category]
    positions = x_pos - 0.4 + (i + 0.5) * bar_width
    bars = ax.bar(positions, category_data['count'], width=bar_width, label=category, color=colors[i])
    
    # Add count labels directly on each bar
    for bar, value in zip(bars, category_data['count']):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                f"{int(value)}", ha='center', va='bottom', fontsize=10, fontweight='bold')

# Set x-axis labels
plt.xticks(x_pos, ['L. gryllorum', 'G. niphandrodes'])

# Customize plot
plt.title('Functional Categories of Proteins by Species', fontsize=16)
plt.xlabel('Species', fontsize=14)
plt.ylabel('Count', fontsize=14)
plt.legend(title='Functional Category', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig('functional_bar_plot.png', dpi=300, bbox_inches='tight')
print("Functional bar plot saved as 'functional_bar_plot.png'")

EOF
python function_bar_plot.py
