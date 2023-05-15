import sys


# Step3: Input- Recounts file from merge similar. Keep only things that are true

tenx_recounts = sys.argv[1]
sample_tenx_exp_dict = dict()
total_counts_dict = dict()

s_col = 0
tenx_col = 1
exp_col = 2
counts_col = 3

ratio_cutoff = 0.35
count_cutoff = 1

with open(tenx_recounts) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[s_col]
        tenx = ls[tenx_col]
        exp = ls[exp_col]
        counts = float(ls[counts_col])
        if sample not in sample_tenx_exp_dict:
            sample_tenx_exp_dict[sample] = dict()
            total_counts_dict[sample] = 0

        if tenx not in sample_tenx_exp_dict[sample]:
            sample_tenx_exp_dict[sample][tenx] = list()
            sample_tenx_exp_dict[sample][tenx].append([exp, counts])
            total_counts_dict[sample] += counts
        else:
            if counts/sample_tenx_exp_dict[sample][tenx][0][1] >= ratio_cutoff:
                sample_tenx_exp_dict[sample][tenx].append([exp, counts])
                total_counts_dict[sample] += counts

print("sample\ttenx_bc\texp_bc\tcount\tcpm")
for sample in sample_tenx_exp_dict:
    for tenx in sample_tenx_exp_dict[sample]:
        for item in sample_tenx_exp_dict[sample][tenx]:
            print("\t".join(map(str, [sample, tenx, item[0], item[1], (item[1]*10**6)/total_counts_dict[sample]])))
