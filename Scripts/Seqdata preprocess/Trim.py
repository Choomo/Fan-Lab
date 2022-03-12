import os

file_path_list=[]
file_out_list=[]
for home, dirnames, filenames in os.walk(os.getcwd()):
    for filename in filenames:
        if filename.endswith('.fastq.gz'):
            file_path_list.append(os.path.join(home, filename))
            file_out_list.append(filename[0:10]+'_out.fastq.gz')

for file_path in file_path_list:
    position=file_path_list.index(file_path)
    os.system('trimmomatic SE -threads 20 -phred33 -trimlog trim.logfile {a} /share/home/zhuxm/liaochh/data/scRNA_seq_retinal_organoid/GSE119343/clean_data/{b}  ILLUMINACLIP:/share/home/zhuxm/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:16 > trim.log'.format(a=file_path, b=file_out_list[position]))
