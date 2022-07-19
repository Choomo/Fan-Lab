#将bed with methyl_level & coverage & base in the position文件转换成如下格式，同时对正链和负链的CpG信号统一进行处理，使得一个CpG位点只有一个平均的甲基化信号, 基因组第一个碱基坐标为1
#chr    start   end index   methy_level coverage
#chr1   1000    1001    chr1_1000   40.2    20
#...
#Usage: python bed_to_metcov.py input_file output_file

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

print(f"Start processing {input_file}!")

with open(output_file, "w") as wfile:
    with open(input_file, "r") as rfile:
        lines = rfile.readlines()
        i = 0
        n = len(lines)
        while i <= n-1:
            if lines[i].split("\t")[5].strip() == "C" or lines[i].split("\t")[5].strip() == "c":  # 判断这个位点是否是C
                if int(lines[i].split("\t")[1]) + 1 == int(lines[i+1].split("\t")[1]) and lines[i].split("\t")[0] == lines[i+1].split("\t")[0]: # 是C并判断这一行是否与下一行连续
                    chromesome = lines[i].split("\t")[0]
                    start = lines[i].split("\t")[2]
                    end = int(start)+1
                    index = chromesome + "_" + start
                    coverage = int(lines[i].split("\t")[4])+int(lines[i+1].split("\t")[4])
                    methy_level = (int(float(lines[i].split("\t")[3])*float(lines[i].split("\t")[4]))+int(float(lines[i+1].split("\t")[3])*float(lines[i+1].split("\t")[4])))/(int(lines[i].split("\t")[4])+int(lines[i+1].split("\t")[4]))
                    wfile.write(f"{chromesome}\t{start}\t{end}\t{index}\t{methy_level}\t{coverage}\n")
                    i += 2
                else:  # 是C但是不连续
                    chromesome = lines[i].split("\t")[0]
                    start = lines[i].split("\t")[2]
                    end = int(start)+1
                    index = chromesome + "_" + start
                    coverage = lines[i].split("\t")[4]
                    methy_level = lines[i].split("\t")[3]
                    wfile.write(f"{chromesome}\t{start}\t{end}\t{index}\t{methy_level}\t{coverage}\n")
                    i += 1
            elif lines[i].split("\t")[5].strip() == "G" or lines[i].split("\t")[5].strip() == "g": # 是G的情况
                chromesome = lines[i].split("\t")[0]
                start = lines[i].split("\t")[1]
                end = int(start)+1
                index = chromesome + "_" + start
                coverage = lines[i].split("\t")[4]
                methy_level = lines[i].split("\t")[3]
                wfile.write(f"{chromesome}\t{start}\t{end}\t{index}\t{methy_level}\t{coverage}\n")
                i += 1
            if i == n-1:
                if lines[i].split("\t")[5].strip() == "C" or lines[i].split("\t")[5].strip() == "c":
                    chromesome = lines[i].split("\t")[0]
                    start = lines[i].split("\t")[2]
                    end = int(start)+1
                    index = chromesome + "_" + start
                    coverage = lines[i].split("\t")[4]
                    methy_level = lines[i].split("\t")[3]
                    wfile.write(f"{chromesome}\t{start}\t{end}\t{index}\t{methy_level}\t{coverage}\n")
                elif lines[i].split("\t")[5].strip() == "G" or lines[i].split("\t")[5].strip() == "g":
                    chromesome = lines[i].split("\t")[0]
                    start = lines[i].split("\t")[1]
                    end = int(start)+1
                    index = chromesome + "_" + start
                    coverage = lines[i].split("\t")[4]
                    methy_level = lines[i].split("\t")[3]
                    wfile.write(f"{chromesome}\t{start}\t{end}\t{index}\t{methy_level}\t{coverage}\n")
                break
            else:
                continue

print(f"{input_file} finished!\tOutput saved in {output_file}")