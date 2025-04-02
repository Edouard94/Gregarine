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
