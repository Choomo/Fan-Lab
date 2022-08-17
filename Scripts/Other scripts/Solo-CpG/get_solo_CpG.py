# Input bed file with fasta squences and output solo CpG sites in bed file form
# Usage python get_solo_CpG.py input_bed_file output_file_name

import sys
import re



input_file = sys.argv[1]
output_file = sys.argv[2]

with open(output_file, "w") as wfile:
    with open(input_file, "r") as rfile:
        for line in rfile:
            line = line.strip()
            lines = line.split("\t")
            fa = lines[len(lines)-1].upper()
            if fa[35:37] == "CG":
                number = len(re.findall(r'CG', fa))
                chr = lines[0]
                start = int(lines[1]) + 36
                end = int(lines[2]) - 35
                methyl_call = lines[3].split(",")
                MM = int(methyl_call[0]) + int(methyl_call[4])
                MU = int(methyl_call[1]) + int(methyl_call[5])
                UU = int(methyl_call[2]) + int(methyl_call[6])
                UM = int(methyl_call[3]) + int(methyl_call[7])
                if fa[34:38] == "ACGA" or fa[34:38] == "ACGT" or fa[34:38] == "TCGT"or fa[34:38] == "TCGA":
                    CG_type = "WCGW"
                elif fa[34:38] == "CCGC" or fa[34:38] == "CCGG" or fa[34:38] == "GCGG" or fa[34:38] == "GCGC":
                    CG_type = "SCGS"
                elif fa[34:38] == "ACGC" or fa[34:38] == "ACGG" or fa[34:38] == "TCGC"or fa[34:38] == "TCGG":
                    CG_type = "WCGS"
                elif fa[34:38] == "CCGA" or fa[34:38] == "CCGT" or fa[34:38] == "GCGA"or fa[34:38] == "GCGT":
                    CG_type = "SCGW"
                wfile.write(f"{chr}\t{start}\t{end}\t{MM}\t{MU}\t{UU}\t{UM}\t{number-1}\t{CG_type}\n")
            else:
                continue