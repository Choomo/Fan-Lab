#将Bismark输出的coverage文件转换成如下格式，同时对正链和负链的CpG信号统一进行处理，使得一个CpG位点只有一个平均的甲基化信号
#chr    start   end index   methy_level coverage
#chr1   1000    1001    chr1_1000   40.2    20
#...
#Usage: python bismarkcov_to_methy_txt.py input_file output_file

import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(output_file, "w") as wfile:
    with open(input_file, "r") as rfile:
        lines = rfile.readlines()
        i = 0
        n = len(lines)
        while i <= n-1:
            if int(lines[i].split("\t")[1]) + 1 == int(lines[i+1].split("\t")[1]) and lines[i].split("\t")[0] == lines[i+1].split("\t")[0]:
                chromesome = lines[i].split("\t")[0]
                start = lines[i].split("\t")[1]
                end = int(start)+1
                index = chromesome + "_" + start
                coverage = int(lines[i].split("\t")[4])+int(lines[i].split("\t")[5])+int(lines[i+1].split("\t")[4])+int(lines[i+1].split("\t")[5])
                methy_level = ((int(lines[i].split("\t")[4])+int(lines[i+1].split("\t")[4]))/coverage)*100
                wfile.write(f"{chromesome}\t{start}\t{end}\t{index}\t{methy_level}\t{coverage}\n")
                i += 2    
            else:
                chromesome = lines[i].split("\t")[0]
                start = lines[i].split("\t")[1]
                end = int(start)+1
                index = chromesome + "_" + start
                coverage = int(lines[i].split("\t")[4])+int(lines[i].split("\t")[5])
                methy_level = lines[i].split("\t")[3]
                wfile.write(f"{chromesome}\t{start}\t{end}\t{index}\t{methy_level}\t{coverage}\n")
                i += 1
            if i == n-1:
                    chromesome = lines[i].split("\t")[0]
                    start = lines[i].split("\t")[1]
                    end = int(start)+1
                    index = chromesome + "_" + start
                    coverage = int(lines[i].split("\t")[4])+int(lines[i].split("\t")[5])
                    methy_level = lines[i].split("\t")[3]
                    wfile.write(f"{chromesome}\t{start}\t{end}\t{index}\t{methy_level}\t{coverage}\n")
                    break
            else:
                continue