import sys

ts_file = sys.argv[1]

group_cols = []
# chemo
# fraction_cutoff = 250
# cf_cutoff = 900000
# group_column = [9, 10]
# replicate_column = 11
# barcode_column = 2
# cell_line_column = 8
# fraction_column = 4
# exclusion_list = ["HCC4006_Aggr_EGF816-Doce_3", "HCC827_Aggr_EGF816-Doce_3"]


# ltr
fraction_cutoff = 250
cf_cutoff = 900000
group_column = [6, 7]
replicate_column = 8
barcode_column = 9
cell_line_column = 5
fraction_column = 10
exclusion_list = ["HCC4006_10_EGF816_2", "HCC4006_10_EGF816_3", "HCC4006_10_EGF816_4", "HCC827_6_EGF816-INC280_3",
                  "HCC827_6_EGF816-INC280_4", "PC9_12_EGF816_2", "PC9_15_EGF816_2"]

with open(ts_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[0]
        if sample in exclusion_list:
            continue

