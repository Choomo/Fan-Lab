import os
import sys
import pandas as pd
import io

in_file = sys.argv[1]
out_file = sys.argv[2]

### Read file in and save as dict
with open(in_file, "r") as in_file:
    mt_dict = {}
    for line in in_file:
        line_list = line.split(sep="\t")
        if line[0].startswith('>'):
            seq_name = line_list[1]
            mt_dict[seq_name] = ''
        else:
            mt_dict[seq_name] += line


##########output transformed matrix
with open(out_file, "w") as w_file:
    for seq_name in mt_dict:
        data = pd.read_csv(io.StringIO(mt_dict[seq_name]), lineterminator='\n', sep='\t', header=None).T
        data.insert(loc=0, column="str", value=["A:", "C:", "G:", "T:"])
        data = data.to_string(header=None, index=None)
        data = data.strip()
        data = data.replace(" ", "\t")
        data = data.replace("\n\t", "\n")
        w_file.write(f'{seq_name}\n{data}\n')


