import sys

xbc_col = 1
cpm_col = 3
sub_col = 1
sample_col = 0

input_file = sys.argv[1]
clone_map = sys.argv[2]

bc_clonemap_dict = dict()
total = dict()
clone_cpm_dict = dict()

with open(clone_map) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        bc_clonemap_dict[ls[0]] = ls[1]
fh.close()


with open(input_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        xbc = ls[xbc_col]
        cpm = float(ls[cpm_col]) * 10**6
        sample = ls[sample_col]
        clone = "NA"

        if xbc in bc_clonemap_dict:
            clone = bc_clonemap_dict[xbc]
        else:
            clone = ls[sub_col]
            bc_clonemap_dict[xbc] = clone

        if sample not in clone_cpm_dict:
            clone_cpm_dict[sample] = dict()

        if clone not in clone_cpm_dict[sample]:
            clone_cpm_dict[sample][clone] = cpm
        else:
            clone_cpm_dict[sample][clone] = (clone_cpm_dict[sample][clone] + cpm) / 2
fh.close()


for sample in clone_cpm_dict:
    total[sample] = sum(clone_cpm_dict[sample].values())

with open(input_file) as fh:
    print(fh.readline().strip()+"\tcloneid\trevised_cellcount_cpm")
    for line in fh:
        ls = line.strip().split("\t")
        xbc = ls[xbc_col]
        sample = ls[sample_col]

        print(line.strip()+"\t"+bc_clonemap_dict[xbc]+"\t"+"%.3f" % ((clone_cpm_dict[sample][bc_clonemap_dict[xbc]] *
                                                                      10 ** 6) / total[sample]))
