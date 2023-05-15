import sys
import math

input_file = sys.argv[1]
d0_dict = dict()
group_dict = dict()

fraction_cutoff = 0.0001
cf_cutoff = 0.99
group_column = [8]
replicate_column = 0
barcode_column = 1
cell_line_column = 4
fraction_column = 3
cf = 0.0
prev_celline = "NA"
prev_sample = "NA"

print("cell_line\tbarcode\ttimepoint\tlog2fc\td0_cpm_avg\treplicate_cpm\tcum_fraction")


def print_dict_values(group_dict_param, d0_dict_param, prev_cell):
    for gp in group_dict_param:
        for bc in group_dict_param[gp]:
            for item in group_dict_param[gp][bc]:
                fr = item[0]
                cf = item[1]
                if bc in d0_dict_param:
                    d0_frac = sum(d0_dict_param[bc])/len(d0_dict_param[bc])
                    print(prev_cell + "\t" + bc + "\t" + gp.replace("_", "\t") + "\t" + "%.3f" % math.log(fr/d0_frac, 2) + "\t" + "%.3f" % (d0_frac * (10**6)) + "\t"+ "%.3f" % (fr * (10**6)) + "\t" + "%.3f" % cf)


with open(input_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[0]
        barcode = ls[barcode_column]
        group = ls[group_column[0]]# + "_" + ls[group_column[0]]
        replicate = ls[replicate_column]
        fraction = float(ls[fraction_column])
        cell_line = ls[cell_line_column]
        if sample != prev_sample:
            cf = 0.0
        if prev_celline != cell_line and prev_celline != "NA":
            print_dict_values(group_dict, d0_dict, prev_celline)
            group_dict = dict()
            d0_dict = dict()
        cf += fraction
        if group == "D0":
            if barcode not in d0_dict:
                d0_dict[barcode] = list()
            d0_dict[barcode].append(fraction)
        else:
            if group not in group_dict:
                group_dict[group] = dict()
            if barcode not in group_dict[group]:
                group_dict[group][barcode] = list()
            #if cf < cf_cutoff:
            group_dict[group][barcode].append([fraction, cf])
        prev_celline = cell_line
        prev_sample = sample

print_dict_values(group_dict, d0_dict, prev_celline)
fh.close()