
#mm10.fa为基因组fasta文件，最后输出结果为第一列为染色体编号，第二列为length
cat mm10.fa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > mm10.genome
