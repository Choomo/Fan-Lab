import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='''
''')

parser.add_argument('-i', type= str, help="input file (DCC output)")
parser.add_argument("-o", type= str, help="output file name")

agrs = parser.parse_args()

input_file = agrs.i
output_file = agrs.o

DCC_data = pd.read_csv(input_file, sep="\t")
DCC_data.insert(loc=3, column="chr_index", value=DCC_data["Chr"].map(str)+":"+DCC_data["Start"].map(str)+"-"+DCC_data["End"].map(str))
DCC_duplicate = DCC_data[-DCC_data.duplicated(subset="chr_index", keep="first")]
DCC_duplicate.to_csv(output_file, sep="\t", index=False)