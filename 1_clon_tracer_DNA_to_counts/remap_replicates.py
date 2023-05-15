import sys

# chemo
# group_column = [9, 10]
# replicate_column = 11
# barcode_column = 16
# cell_line_column = 8
# print("cell_line\tbarcode\tmode\ttreatment\treplicate\tgroup\tquery_group\tcount")

# ltr
group_column = [2, 3]
replicate_column = 4
barcode_column = 1
cell_line_column = 0


tall_skinny_file = sys.argv[1]
replicate_remap_dict = dict()
remap_dict = dict()
prev_celline = "NA"

with open(tall_skinny_file) as fh:
    print(fh.readline().strip()+"\treplicate_remap")
    for line in fh:
        ls = line.strip().split("\t")
        barcode = ls[barcode_column]
        group = ls[group_column[1]] + "_" + ls[group_column[0]]
        replicate = ls[replicate_column]
        cell_line = ls[cell_line_column]
        if prev_celline != cell_line and prev_celline != "NA":
            replicate_remap_dict = dict()
            remap_dict = dict()
        if barcode not in replicate_remap_dict:
            replicate_remap_dict[barcode] = dict()
            remap_dict[barcode] = dict()
        if group not in replicate_remap_dict[barcode]:
            replicate_remap_dict[barcode][group] = set()
            remap_dict[barcode][group] = dict()
        if replicate not in replicate_remap_dict[barcode][group]:
            replicate_remap_dict[barcode][group].add(replicate)
            remap_dict[barcode][group][replicate] = len(replicate_remap_dict[barcode][group])
        print(line.strip() + "\t" + str(remap_dict[barcode][group][replicate]))
        prev_celline = cell_line
