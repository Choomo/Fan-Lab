# %%
## Usage python script.py -i input.fa -o output.fa -k keyword

import argparse
import linecache

parser = argparse.ArgumentParser(description='''
Extract reads in fasta file with key words
''')
parser.add_argument('-i', type=str, help="input file (fasta format)")
parser.add_argument("-o", type=str, help="output file name (fasta format)")
parser.add_argument("-k", type=str, help="Key word in reads name")
agrs = parser.parse_args()

input_file = agrs.i # Input file path
output_file_name = agrs.o # Output file name
keyword = agrs.k # Keyword




# %%
fa_ID_list = []
with open(input_file, "r") as rfile:
    for lineID, line in enumerate(rfile):
        if line.startswith(keyword):
            fa_ID_list.append(lineID)
        else:
            continue
with open(output_file_name, "w") as wfile:
    for line_ID in fa_ID_list:
        line_ID += 1
        line = linecache.getline(input_file, line_ID)
        wfile.write(line)
        line_ID += 1
        line = linecache.getline(input_file, line_ID)
        while not line.startswith(">"):
            wfile.write(f"{line.rstrip()}\n")
            line_ID += 1
            line = linecache.getline(input_file, line_ID)


# %%



