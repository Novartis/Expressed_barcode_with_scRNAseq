import sys

cloneid_col = 9 # 2
ctg_cpm_col = 8 # 11
sub_col = 1
sample_col = 0

input_file = sys.argv[1]

bc_clonemap_dict = dict()
total = dict()
clone_cpm_dict = dict()


with open(input_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        cpm = float(ls[ctg_cpm_col])
        sample = ls[sample_col]
        clone = ls[cloneid_col]

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
    print(fh.readline().strip()+"\trevised_cellcount")
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[sample_col]
        cloneid = ls[cloneid_col]

        print(line.strip()+"\t"+"%.3f" % ((clone_cpm_dict[sample][cloneid] * 10**6) / total[sample]))
