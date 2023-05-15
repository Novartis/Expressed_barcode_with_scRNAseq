import sys

input_file=sys.argv[1]
input_sample_col = 0
input_cpm_col = 10

ctg_info=sys.argv[2]
ctg_sample_col = 0
ctg_val_col = 3

ctg_dict = dict()

with open(ctg_info) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split(",")
        ctg_dict[ls[ctg_sample_col]] = float(ls[ctg_val_col])/100
fh.close()

with open(input_file) as fh:
    print(fh.readline().strip()+"\tctg_corrected_cpm")
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[input_sample_col]
        cpm = float(ls[input_cpm_col])
        print(line.strip() + "\t" + "%.3f" % (cpm * ctg_dict[sample]))
fh.close()
