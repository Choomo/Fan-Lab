import re
import os
import sys
import linecache

# 获取命令行参数
fasta_file = sys.argv[1]

# 定义一个函数，用于统计文件行数
def wc_count(file_name):
    import subprocess
    out = subprocess.getoutput("wc -l %s" % file_name)
    return int(out.split()[0])

# 获取fasta文件的总行数
row_numbers = wc_count(fasta_file)
count = 0

# 遍历每一行
for n in range(1,row_numbers):
    line = linecache.getline(fasta_file, n).strip()
    try:
        line_2 = linecache.getline(fasta_file, n+1).strip()
    except:
        None
    # 如果以">"开头，则获取染色体名称
    if line.startswith('>'):
        chr_name = line[1:].strip()
        count = 0
        continue
    else:
        # 将序列转换为大写
        line = line.strip().upper()
        seq = line
        # 查找序列中的"CG"
        it = re.finditer(r"CG", seq)
        match_results = []
        for i in it:
            match_results.append(i.start())
        # 如果找到"CG"，则输出染色体名称、起始位置和结束位置
        if len(match_results) != 0:
            for a in match_results:
                start = int(a) + count
                end = start + 2
                print(f"{chr_name}\t{start}\t{end}")
        # 如果当前行以"C"结尾，且下一行以"G"开头，则输出染色体名称、起始位置和结束位置
        if line[-1] == "C" and line_2[0] == "G":
            start = 49 + count
            end = start + 2
            print(f"{chr_name}\t{start}\t{end}")
        # 计数器加50
        count += 50