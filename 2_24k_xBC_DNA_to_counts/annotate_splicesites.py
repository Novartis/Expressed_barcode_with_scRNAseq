import sys

splice_sites = ["GGTAA", "GGTG"]

input_file = sys.argv[1]
barcode_col = 1

with open(input_file) as fh:
    print(fh.readline().strip()+"\tsplice_site")
    for line in fh:
        ls = line.strip().split("\t")
        barcode = ls[barcode_col]
        ss = any(item in barcode for item in splice_sites)
        print(line.strip()+"\t"+str(ss))
fh.close()

