import sys
import math


input_file = sys.argv[1]  # /da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20200410_holiday_retreat_DNA/output/counts_summary_annt_cloneid_ctg.txt
cell_line = 5
timepoint = 6
treatment = 7
replicate = 8
clone_id = 9
revised_cpm = 10
cpm_cutoff = 2000

sample_cpm_dict = dict()

cl_set = set()
tp_set = set()
tr_set = set()
rp_set = set()
ci_set = set()
tr_tp_dict = dict()

treatment_pairs_generic = [["T", "U"], ["HR", "U"], ["HR", "T"]]
treatment_pairs_tp_matched = [["HR", "H"]]


with open(input_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        cl = ls[cell_line]
        tp = ls[timepoint]
        tr = ls[treatment]
        rp = ls[replicate]
        ci = ls[clone_id]
        cp = float(ls[revised_cpm])

        cl_set.add(cl)
        tp_set.add(tp)
        tr_set.add(tr)
        rp_set.add(rp)
        ci_set.add(ci)

        if tr not in tr_tp_dict:
            tr_tp_dict[tr] = list()

        tr_tp_dict[tr].append(tp)

        if cl not in sample_cpm_dict:
            sample_cpm_dict[cl] = dict()
        if (tr, tp, rp) not in sample_cpm_dict[cl]:
            sample_cpm_dict[cl][(tr, tp, rp)] = dict()
        sample_cpm_dict[cl][(tr, tp, rp)][ci] = cp
fh.close()

for cl_l in sample_cpm_dict:
    with open(cl_l+"_log2fc.txt", 'w') as fhw:
        print("cellline\tclone_id\tcomparator_group\tcomparator_replicate\tnum_cpm\tden_cpm\tlog2fc", file=fhw)
        for pair in treatment_pairs_generic:
            for timepoint0 in tr_tp_dict[pair[0]]:
                for timepoint1 in tr_tp_dict[pair[1]]:
                    for rp_l in sorted(rp_set):
                        if (pair[0], timepoint0, rp_l) not in sample_cpm_dict[cl_l] or \
                                (pair[1], timepoint1, rp_l) not in sample_cpm_dict[cl_l]:
                            continue

                        numerator_dict = sample_cpm_dict[cl_l][(pair[0], timepoint0, rp_l)]
                        denominator_dict = sample_cpm_dict[cl_l][(pair[1], timepoint1, rp_l)]

                        for ci_l in ci_set:
                            num_val = numerator_dict[ci_l] if ci_l in numerator_dict else 1.0
                            den_val = denominator_dict[ci_l] if ci_l in denominator_dict else 1.0

                            if max(num_val, den_val) <= cpm_cutoff:
                                continue

                            comparator_grp = "_".join([pair[0], timepoint0, "vs", pair[1], timepoint1])
                            comparator_rep = "_".join([pair[0], timepoint0, rp_l, "vs", pair[1], timepoint1, rp_l])

                            if den_val == 0.0:
                                log2fc = 10
                            else:
                                log2fc = math.log(num_val/den_val, 2)

                            print("\t".join([cl_l, ci_l, comparator_grp, comparator_rep, "%.3f" % num_val,
                                             "%.3f" % den_val, "%.3f" % log2fc]), file=fhw)

        for pair in treatment_pairs_tp_matched:
            for timepoint0 in tr_tp_dict[pair[0]]:
                for timepoint1 in tr_tp_dict[pair[1]]:
                    if timepoint1 != timepoint0:
                        continue
                    for rp_l in sorted(rp_set):
                        if (pair[0], timepoint0, rp_l) not in sample_cpm_dict[cl_l] or \
                                (pair[1], timepoint1, rp_l) not in sample_cpm_dict[cl_l]:
                            continue

                        numerator_dict = sample_cpm_dict[cl_l][(pair[0], timepoint0, rp_l)]
                        denominator_dict = sample_cpm_dict[cl_l][(pair[1], timepoint1, rp_l)]

                        for ci_l in ci_set:
                            num_val = numerator_dict[ci_l] if ci_l in numerator_dict else 1.0
                            den_val = denominator_dict[ci_l] if ci_l in denominator_dict else 1.0

                            if max(num_val, den_val) <= cpm_cutoff:
                                continue

                            comparator_grp = "_".join([pair[0], timepoint0, "vs", pair[1], timepoint1])
                            comparator_rep = "_".join([pair[0], timepoint0, rp_l, "vs", pair[1], timepoint1, rp_l])

                            if den_val == 0.0:
                                log2fc = 10
                            else:
                                log2fc = math.log(num_val/den_val, 2)

                            print("\t".join([cl_l, ci_l, comparator_grp, comparator_rep, "%.3f" % num_val,
                                             "%.3f" % den_val, "%.3f" % log2fc]), file=fhw)
    fhw.close()
