import os
import sys

#计算hammer-seq输出bed文件中每一个位点的甲基化水平
#使用方法 python calc_methyl_level_per_site.py input_bed_file output_file_name

input_file = sys.argv[1]
output_filename = sys.argv[2]

with open(output_filename, "w") as wfile:
    with open(input_file, "r") as rfile:
        for line in rfile:
            methyl_events = line.split("\t")[3] 
            sum = 0
            for i in methyl_events.split(","):
                sum += int(i)
            if sum == 0:
                continue
            else:
                events = methyl_events.split(",")
                methylation_level = (int(events[0])*2 + int(events[1]) + int(events[3]) + int(events[4]) + int(events[5]) + int(events[7]))/(sum*2)
                chro_info = line.split("\t")[0]
                pos_info = line.split("\t")[1]
                wfile.write(f"{chro_info}.{pos_info}\t{methylation_level}\n")
        
                    