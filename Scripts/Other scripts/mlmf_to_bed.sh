#!/bin/bash
# Inputfile_name = $1

sed "1d" $1 > $1.tmp
awk -v OFS="\t" '{print $1, $2-1, $2, $3, $5+$6+$7+$8}' $1.tmp > ${$1:0:33}.bed
wait
rm $1.tmp

for i in *.bed;
do
bedtools subtract -a $i -b mm10-blacklist.v2.bed > $i