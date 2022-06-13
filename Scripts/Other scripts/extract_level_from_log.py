from operator import length_hint
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]


with open(output_file, "w") as wfile:
    with open(input_file, "r") as rfile:
        lines = rfile.readlines()
        length = range(int(len(lines)/3))
        for i in length:
            maintenance_ratio = lines[i*3+1].split("\t")[1].strip()
            de_novo_ratio = lines[i*3+2].split("\t")[1][:-2]
            wfile.write(f"{maintenance_ratio},{de_novo_ratio}\n")