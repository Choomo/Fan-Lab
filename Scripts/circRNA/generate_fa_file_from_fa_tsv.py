import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, "r") as rfile:
    with open(output_file, "w") as wfile:
        lines = rfile.readlines()
        for line in lines:
            line_list = line.split("\t")
            col_number = len(line_list) - 1
            circ_name = line_list[0].strip()
            fa_seq = line_list[col_number].strip()
            if fa_seq == "NA" or fa_seq == "partial":
                continue
            else:
                wfile.write(f">{circ_name}\n{fa_seq}\n")

