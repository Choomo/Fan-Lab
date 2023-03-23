#Usage python calcu_methyl_level_NNCGNN_fa_bed.py inputfile outputfile

import pandas as pd
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(input_file, sep="\t", names=["chr", "start", "end", "index", "methy_level", "strand", "CG_type"])
df['CG_type'] = df["CG_type"].str.upper()
df = df.loc[:, ['CG_type', 'methy_level']]
df_mean = df.groupby("CG_type").mean()
CG_type_counts = df.groupby("CG_type").size()
df_n = pd.concat([df_mean, CG_type_counts], axis =1)
df_n.to_csv(output_file, sep="\t", header=None)