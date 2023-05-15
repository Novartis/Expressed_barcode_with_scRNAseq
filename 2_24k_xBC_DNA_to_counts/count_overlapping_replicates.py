import sys


tall_skinny_file = sys.argv[1]

group_column = [1]
replicate_column = 2
barcode_column = 0
prev_sample = "NA"
group_set = set()

barcode_group_replicate_dict = dict()

print("barcode\tgroup\treplicate\tquery_group\tcount")


def print_dict_values(the_dict):
        for bc in the_dict:
            for gp in the_dict[bc]:
                for rp in the_dict[bc][gp]:
                    for super_group in group_set:
                        if super_group in the_dict[bc]:
                            print(bc+"\t"+gp+"\t"+rp+"\t"+super_group+"\t"+str(len(the_dict[bc][super_group])))


with open(tall_skinny_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        barcode = ls[barcode_column]
        group = ls[group_column[0]]
        replicate = ls[replicate_column]
        if barcode not in barcode_group_replicate_dict:
            barcode_group_replicate_dict[barcode] = dict()
        if group not in barcode_group_replicate_dict[barcode]:
            group_set.add(group)
            barcode_group_replicate_dict[barcode][group] = list()
        barcode_group_replicate_dict[barcode][group].append(replicate)
fh.close()

print_dict_values(barcode_group_replicate_dict)


