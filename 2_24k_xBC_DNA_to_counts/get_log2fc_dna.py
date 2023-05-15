import sys
import math


input_file = sys.argv[1]  # /da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20191219_4cellline_DNA/output/counts_summary_cloneid.txt
sample = 0
clone_id = 2
revised_cpm = 4

cpm_cutoff = 2000

sample_cpm_dict = dict()
ci_set = set()
tp_pairs_generic = [["14", "0"], ["16", "14"], ["1", "0"]]

with open(input_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        smp = ls[sample].split("_")
        ci = ls[clone_id]
        cp = float(ls[revised_cpm])
        cl = smp[0]
        tp = smp[1]

        ci_set.add(ci)

        if cl not in sample_cpm_dict:
            sample_cpm_dict[cl] = dict()
        if tp not in sample_cpm_dict[cl]:
            sample_cpm_dict[cl][tp] = dict()
        sample_cpm_dict[cl][tp][ci] = cp
fh.close()

for cl_l in sample_cpm_dict:
    with open(cl_l+"_log2fc.txt", 'w') as fhw:
        print("cellline\tclone_id\tcomparator_group\tnum_cpm\tden_cpm\tlog2fc", file=fhw)
        for pair in tp_pairs_generic:
            if pair[0] not in sample_cpm_dict[cl_l] or \
                    pair[1] not in sample_cpm_dict[cl_l]:
                continue

            numerator_dict = sample_cpm_dict[cl_l][pair[0]]
            denominator_dict = sample_cpm_dict[cl_l][pair[1]]

            for ci_l in ci_set:
                num_val = (numerator_dict[ci_l] + 1.0) if ci_l in numerator_dict else 1.0
                den_val = (denominator_dict[ci_l] + 1.0) if ci_l in denominator_dict else 1.0

                if max(num_val, den_val) <= cpm_cutoff:
                    continue

                comparator_grp = "_".join([pair[0], "vs", pair[1]])

                log2fc = math.log(num_val/den_val, 2)

                print("\t".join([cl_l, ci_l, comparator_grp, "%.3f" % num_val,
                                 "%.3f" % den_val, "%.3f" % log2fc]), file=fhw)
    fhw.close()
