# %%
# Usage python merge_matrix.py infiles > outfile

import pandas as pd
import numpy as np
import sys

infiles = sys.argv[1:]

# Read the first file
df = pd.read_csv(infiles[0], sep='\t', index_col=None, header=None)
# Get the column names
chr = df.iloc[:,0]
pos = df.iloc[:,1]
strand = df.iloc[:,2]
# Create a new column with the combined index
df['index'] = chr.str.cat(pos.astype(str), sep=':').str.cat(strand, sep=':')
# Calculate the methylation ratio
df['methy_ratio'] = (df.iloc[:,4]/df.iloc[:,5])*100
# Filter out the columns we need
df = pd.DataFrame(df, columns=["index", "methy_ratio"])
# Rename the columns
df = df.rename(columns={'index': 'index', 'methy_ratio': infiles[0]})

# %%
# Iterate through the infiles list starting from the second element
for infile in infiles[1:]:
    # Read the CSV file with tab-delimited columns and no header
    df_tmp = pd.read_csv(infile, sep='\t', index_col=None, header=None)
    # Extract the chromosome, position, and strand columns
    chr = df_tmp.iloc[:,0]
    pos = df_tmp.iloc[:,1]
    strand = df_tmp.iloc[:,2]
    # Create a new column that combines the chromosome, position, and strand
    df_tmp['index'] = chr.str.cat(pos.astype(str), sep=':').str.cat(strand, sep=':')
    # Calculate the methylation ratio by dividing the methylated reads by the total reads
    df_tmp['methy_ratio'] = (df_tmp.iloc[:,4]/df_tmp.iloc[:,5])*100
    # Filter the DataFrame to contain only the required columns
    df_tmp = pd.DataFrame(df_tmp, columns=["index", "methy_ratio"])
    # Rename the columns
    df_tmp = df_tmp.rename(columns={'index': 'index', 'methy_ratio': infile})
    # Merge the current DataFrame with the previous ones on the index column
    df = pd.merge(df, df_tmp, on='index', how='outer')

# Write the merged DataFrame to stdout
df.to_csv(sys.stdout, sep='\t', index=False)


