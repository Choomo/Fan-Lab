#Usage: python mlmf_to_bed.py input_file output_file
#Get NNCGNN bed files

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
cutoff = int(sys.argv[3])

strand_infor = "."

with open(output_file, "w") as wfile:
    with open(input_file, "r") as rfile:
        lines = rfile.readlines()
        lines = lines[1:]
        for line in lines:
            strs = line.split("\t")
            signal_counts = int(strs[4]) + int(strs[5]) + int(strs[6]) + int(strs[7]) + int(strs[8]) + int(strs[9]) + int(strs[10]) + int(strs[11])
            if signal_counts < cutoff:
                continue
            else:
                chr_infor = strs[0]
                start_infor = int(strs[1]) - 3
                end_infor = int(strs[1]) + 3
                CG_index = chr_infor + ":" + str(start_infor) + "-" + str(end_infor)
                methy_level = strs[2]
                wfile.write(f"{chr_infor}\t{start_infor}\t{end_infor}\t{CG_index}\t{methy_level}\t{strand_infor}\n")
