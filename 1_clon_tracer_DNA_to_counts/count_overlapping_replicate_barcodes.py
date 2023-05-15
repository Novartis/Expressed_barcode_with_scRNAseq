import sys
import numpy as np

tall_skinny_file = sys.argv[1]

# fraction_cutoff = 0.0001
# cf_cutoff = 0.99
# group_column = [5, 8]
# replicate_column = 7
# barcode_column = 1
# cell_line_column = 4
# fraction_column = 3
#print("cell_line\tbarcode\ttimepoint\tdrug\treplicate\tquery_group\tcount")

# chemo
# fraction_cutoff = 250
# cf_cutoff = 950000
# group_column = [9, 10]
# replicate_column = 11
# sample_column = 0
# barcode_column = 16
# cell_line_column = 8
# fraction_column = 17
# exclusion_list = ["HCC4006_Aggr_EGF816-Doce_3", "HCC827_Aggr_EGF816-Doce_3"]
# print("cell_line\tbarcode\tmode\ttreatment\treplicate\tgroup\tquery_group\tcount")
# print("cell_line\tcloneid\tmode\ttreatment\tnumber_replicates\tgroup\tmean_cpm\tmedian_cpm", file=sys.stderr)

# ltr
# fraction_cutoff = 250
# cf_cutoff = 950000
# group_column = [6, 7]
# replicate_column = 8
# sample_column =0
# barcode_column = 9
# cell_line_column = 5
# fraction_column = 10
# exclusion_list = ["HCC4006_4_EGF816_1", "HCC4006_4_EGF816_4",
#                   "HCC4006_7_EGF816_1", "HCC4006_7_EGF816_4",
#                   "HCC4006_10_EGF816_1", "HCC4006_10_EGF816_2", "HCC4006_10_EGF816_3", "HCC4006_10_EGF816_4",
#                   "HCC4006_10_EGF816_5",
#                   "HCC827_4_EGF816_1", "HCC827_4_EGF816_4",
#                   "HCC827_5_EGF816_1", "HCC827_5_EGF816_4",
#                   "HCC827_7_EGF816_1", "HCC827_7_EGF816_4",
#                   "HCC827_6_EGF816-INC280_3", "HCC827_6_EGF816-INC280_4",
#                   "HCC827_7_EGF816-INC280_3", "HCC827_7_EGF816-INC280_4",
#                   "PC9_7_EGF816_1", "PC9_7_EGF816_4",
#                   "PC9_12_EGF816_2", "PC9_12_EGF816_3",
#                   "PC9_15_EGF816_5", "PC9_15_EGF816_2"]
# print("cell_line\tcloneid\ttreatment\ttimepoint\treplicate\tgroup\tquery_group\tcount")
# print("cell_line\tcloneid\ttreatment\ttimepoint\tnumber_replicates\tgroup\tmean_cpm\tmedian_cpm", file=sys.stderr)

# exclusion_list = ["HCC4006_10_EGF816_1", "HCC4006_10_EGF816_5",
#                   "HCC827_6_EGF816-INC280_3", "HCC827_6_EGF816-INC280_4",
#                   "PC9_12_EGF816_2", "PC9_15_EGF816_2"]

# CART CROPSeq
fraction_cutoff = 250
cf_cutoff = 800000
group_column = [7, 8, 9]
replicate_column = 10
sample_column = 12
barcode_column = 2
cell_line_column = "NAP"
fraction_column = 5
exclusion_list = ["Plasmid_control1", "Plasmid_control2", "Plasmid_control_R14B", "Plasmid_control_R13B", "", ""]
print("cell_line\tbarcode\tdonor\tsource\tclone\treplicate\tgroup\tquery_group\tcount")
print("cell_line\tcloneid\tdonor\tsource\tclone\tnumber_replicates\tgroup\tmean_cpm\tmedian_cpm", file=sys.stderr)

prev_celline = "NA"
prev_sample = "NA"
group_set = set()
barcode_group_replicate_dict = dict()
barcode_group_frac_dict = dict()
barcode_gp_replicate_set = set()


def print_dict_values(the_dict, frac_dict, prev_cell):
        for bc in the_dict:
            for gp in the_dict[bc]:
                for rp in the_dict[bc][gp]:
                    for super_group in group_set:
                        if super_group in the_dict[bc]:
                            print(prev_cell+"\t"+bc+"\t"+gp.replace("_", "\t")+"\t"+rp+"\t"+gp+"\t"+super_group+"\t"+str(len(the_dict[bc][super_group])))
                number_replicates = len(frac_dict[bc][gp])
                if number_replicates > 0:
                    mean_cpm = "%.3f" % np.mean([item for item in frac_dict[bc][gp].values()])
                    median_cpm = "%.3f" % np.median([item for item in frac_dict[bc][gp].values()])
                    print(prev_cell + "\t" + bc + "\t" + gp.replace("_", "\t") + "\t" + str(number_replicates) + "\t" + gp +
                          "\t" + mean_cpm + "\t" + median_cpm, file=sys.stderr)


with open(tall_skinny_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[sample_column]
        if sample in exclusion_list:
            continue
        barcode = ls[barcode_column]
        group = "_".join([ls[item] for item in group_column])
        replicate = ls[replicate_column]
        fraction = float(ls[fraction_column])
        cell_line = ls[cell_line_column] if cell_line_column != "NAP" else cell_line_column
        if sample != prev_sample:
            cf = 0.0
        if prev_celline != cell_line and prev_celline != "NA":
            print_dict_values(barcode_group_replicate_dict, barcode_group_frac_dict, prev_celline)
            barcode_group_replicate_dict = dict()
            barcode_group_frac_dict = dict()
        if barcode not in barcode_group_replicate_dict:
            barcode_group_replicate_dict[barcode] = dict()
            barcode_group_frac_dict[barcode] = dict()
        if group not in barcode_group_replicate_dict[barcode]:
            group_set.add(group)
            barcode_group_replicate_dict[barcode][group] = list()
            barcode_group_frac_dict[barcode][group] = dict()
        cf += fraction
        if cf <= cf_cutoff and (barcode, group, replicate) not in barcode_gp_replicate_set:
            barcode_group_replicate_dict[barcode][group].append(str(replicate))
            barcode_group_frac_dict[barcode][group][replicate] = fraction
            barcode_gp_replicate_set.add((barcode, group, replicate))
        prev_celline = cell_line
        prev_sample = sample
fh.close()

print_dict_values(barcode_group_replicate_dict, barcode_group_frac_dict, prev_celline)



