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

