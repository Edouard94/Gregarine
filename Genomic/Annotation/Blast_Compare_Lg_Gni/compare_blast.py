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
