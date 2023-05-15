import sys

tenx_tall_skinny = sys.argv[1]

tenx_dict = dict()
cw_cutoff = 0.78
tenx_reference = True

sample_col = 0
tenx_col = 1
exp_col = 2
tenx_ref_col = 3
count_col = 4
cw_col = 6


def print_dict(tenx_dict, sample, total):
    for tenx in tenx_dict:
        for exp in tenx_dict[tenx]:
            print("\t".join(map(str, [sample, tenx, exp, tenx_dict[tenx][exp], "%.2f" % ((tenx_dict[tenx][exp] * 10.0**6)/total)])))


with open(tenx_tall_skinny) as fh:
    fh.readline()
    sample = ""
    total = 0
    print("sample\ttenx_bc\texp_bc\tcount\tcpm")
    for line in fh:
        ls = line.strip().split("\t")
        if sample != ls[sample_col]:
            if sample != "":
                print_dict(tenx_dict, sample, total)
            tenx_dict = dict()
            total = 0
        if float(ls[cw_col]) > cw_cutoff or ls[tenx_ref_col] == "False":
            continue
        if ls[tenx_col] not in tenx_dict:
            tenx_dict[ls[tenx_col]] = dict()
        if ls[exp_col] not in tenx_dict[ls[tenx_col]]:
            tenx_dict[ls[tenx_col]][ls[exp_col]] = 0
        tenx_dict[ls[tenx_col]][ls[exp_col]] += int(ls[count_col])
        total += int(ls[count_col])
        sample = ls[sample_col]
    print_dict(tenx_dict, sample, total)
fh.close()
