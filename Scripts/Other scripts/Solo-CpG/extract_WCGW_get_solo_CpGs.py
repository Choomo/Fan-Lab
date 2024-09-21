# Input bed file with fasta squences and output solo CpG sites in bed file form
# Usage python extract_WCGW_get_solo_CpGs.py input_bed_file_with_fasta output_file_prefix
# Input: bed file with fasta squences at the end of each line column1 (chr), column2 (CG_start), column3 (fasta_start), column4 (strand), column5 (fasta)
# Output: 1. All CpG sites in bed file form
#         2. Solo CpG sites in bed file form
#         3. Solo CpG sites in WCGW, WCGS, SCGW, SCGS in bed file form


import sys
import re
import math

input_file = sys.argv[1]
output_prefix = sys.argv[2]

overall_CG = open(f"{output_prefix}.bed", "w")
solo_CG = open(f"{output_prefix}.solo.bed", "w")
outfile_WCGW = open(f"{output_prefix}.solo.WCGW.bed", 'w')
outfile_WCGS = open(f"{output_prefix}.solo.WCGS.bed", 'w')
outfile_SCGW = open(f"{output_prefix}.solo.SCGW.bed", 'w')
outfile_SCGS = open(f"{output_prefix}.solo.SCGS.bed", 'w')

with open(input_file, "r") as rfile:
        for line in rfile:
            line = line.strip()
            lines = line.split("\t")
            fa = lines[len(lines)-1].upper()
            i = int(lines[2]) - int(lines[1])  # Calculate the start position of the CG site in this slop fasta sequence
            if fa[i:i+2] == "CG":
                number = len(re.findall(r'CG', fa))
                chr = lines[0]
                start = int(lines[2])
                end = start + 2
                CG_context = fa[i-1:i+3]  # Extact CG type
                if  CG_context == "ACGA" or CG_context == "ACGT" or CG_context == "TCGT"or CG_context == "TCGA":
                    CG_type = "WCGW"
                elif CG_context == "CCGC" or CG_context == "CCGG" or CG_context == "GCGG" or CG_context == "GCGC":
                    CG_type = "SCGS"
                elif CG_context == "ACGC" or CG_context == "ACGG" or CG_context == "TCGC"or CG_context == "TCGG":
                    CG_type = "WCGS"
                elif CG_context == "CCGA" or CG_context == "CCGT" or CG_context == "GCGA"or CG_context == "GCGT":
                    CG_type = "SCGW"
                overall_CG.write(f"{chr}\t{start}\t{end}\t{number-1}\t{CG_type}\n")
                if number - 1 == 0:
                    solo_CG.write(f"{chr}\t{start}\t{end}\t{number-1}\t{CG_type}\n")
                    if CG_type == "WCGW":
                        outfile_WCGW.write(f"{chr}\t{start}\t{end}\t{number-1}\t{CG_type}\n")
                    elif CG_type == "WCGS":
                        outfile_WCGS.write(f"{chr}\t{start}\t{end}\t{number-1}\t{CG_type}\n")
                    elif CG_type == "SCGW":
                        outfile_SCGW.write(f"{chr}\t{start}\t{end}\t{number-1}\t{CG_type}\n")
                    elif CG_type == "SCGS":
                        outfile_SCGS.write(f"{chr}\t{start}\t{end}\t{number-1}\t{CG_type}\n")
            else:
                continue