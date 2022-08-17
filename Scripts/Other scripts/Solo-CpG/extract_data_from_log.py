import sys
import linecache

input_file = sys.argv[1]
output_file = sys.argv[2]

def wc_count(file_name):
    import subprocess
    out = subprocess.getoutput("wc -l %s" % file_name)
    return int(out.split()[0])

row_numbers = wc_count(input_file)



with open(output_file, "w") as wfile:
    wfile.write("sample,ROI,maintennce_ratio,deno_vo_ratio\n")
    sample_lines_number = []
    for n in range(1,row_numbers):
        line = linecache.getline(input_file, n).strip()
        if line.startswith("Hairpin"):
            sample_lines_number.append(n)
    for i in sample_lines_number:
        sample = linecache.getline(input_file, i).strip()
        for a in range(i+1, row_numbers):
            if linecache.getline(input_file, a).strip().startswith("Hairpin"):
                break
            elif linecache.getline(input_file, a).strip().startswith("For"):
                ROI = linecache.getline(input_file, a-1).strip()
                maintennce_ratio = linecache.getline(input_file, a+1).strip().split("\t")[1]
                deno_vo_ratio = linecache.getline(input_file, a+2).strip().split("\t")[1][:-1]
                wfile.write(f"{sample},{ROI},{maintennce_ratio},{deno_vo_ratio}\n")
            else:
                continue
                