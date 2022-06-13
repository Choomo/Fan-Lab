#!/bin/bash
ftp -v -n 210.22.109.162<<EOF
user huganlu huganlu
binary
cd /Project_huganlu_220311S
lcd /shareEMC/home/hugl/data/hairpin-5/
prompt
mget *
bye
EOF
echo "download from ftp successfully"
