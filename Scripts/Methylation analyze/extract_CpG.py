import re
import os
import sys
import linecache

fasta_file = sys.argv[1]

def wc_count(file_name):
    import subprocess
    out = subprocess.getoutput("wc -l %s" % file_name)
    return int(out.split()[0])

row_numbers = wc_count(fasta_file)
count = 0

for n in range(1,row_numbers):
    line = linecache.getline(fasta_file, n).strip()
    try:
        line_2 = linecache.getline(fasta_file, n+1).strip()
    except:
        None
    if line.startswith('>'):
        chr_name = line[1:].strip()
        count = 0
        continue
    else:
        line = line.strip().upper()
        seq = line
        it = re.finditer(r"CG", seq)
        match_results = []
        for i in it:
            match_results.append(i.start())
        if len(match_results) != 0:
            for a in match_results:
                start = int(a) + count
                end = start + 2
                print(f"{chr_name}\t{start}\t{end}")
        if line[-1] == "C" and line_2[0] == "G":
            start = 49 + count
            end = start + 2
            print(f"{chr_name}\t{start}\t{end}")
        count += 50