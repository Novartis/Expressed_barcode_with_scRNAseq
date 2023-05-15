import sys

input_file = sys.argv[1]
col_names = sys.argv[2]

col_names_split = col_names.split(",")
sample_col = 0

use_neg_for_T = True

with open(input_file) as fh:
    print(fh.readline().strip()+"\t"+"\t".join(col_names_split))
    for line in fh:
        ls = line.strip().split("\t")
        annt = ls[sample_col].split("_")

        if use_neg_for_T:
            if annt[2] == "U":
                annt = [annt[0], "-2", annt[2], annt[3]]
            elif annt[2] == "T":
                annt = [annt[0], "0", annt[2], annt[3]]

        print(line.strip()+"\t"+"\t".join(annt))
