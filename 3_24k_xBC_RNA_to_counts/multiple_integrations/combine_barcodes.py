import sys

input_file = "/da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20190805_4cellline/output_v3_DNA_info/cooccurence_bathroom_sink.txt"

pairs_sample_dict = dict()

sample_col = 0
xbc1_col = 1
xbc2_col = 2
xbc1_occ = 5
xbc2_occ = 6
co_occ = 11
rna_1_col = 7
rna_2_col = 8
dna_1_col = 9
dna_2_col = 10
xbc1_sc_col = 3
xbc2_sc_col = 4

dna_1_unt_cols = [16, 18, 20]
dna_2_unt_cols = [17, 19, 21]
dna_1_t_cols = [22, 24, 26]
dna_2_t_cols = [23, 25, 27]

ratio_cutoff = 0.25
min_cooccurence = 4
absolute_min_occurence = 2
min_frac_cooccurence = 0.2
max_cooccurence_override = 10

def get_cpm(line_split, col):
    if line_split[col] == "NA":
        return 0.0

    return float(line_split[col])


def get_ratios(val1, val2):
    if val1 == 0.0 and val2 == 0.0:
        return 1

    return abs(val1 - val2) / max(val1, val2)


def return_stats(line_split):

    stats = {"cooccur": float(line_split[co_occ]), "cooccur_frac_to_min": 0.0, "rna_diff": 0.0,
             "dna_diffs": [0.0, 0.0, 0.0], "min_rna_cpm": 0.0, "min_dna_cpm": 0.0}

    min_occurence = min(float(line_split[xbc1_occ]), float(line_split[xbc2_occ]))

    stats["cooccur_frac_to_min"] = abs(stats["cooccur"] - min_occurence) / min_occurence

    rna_1 = get_cpm(line_split, rna_1_col)
    rna_2 = get_cpm(line_split, rna_2_col)
    stats["min_rna_cpm"] = min(rna_1, rna_2)
    stats["rna_diff"] = get_ratios(rna_1, rna_2)

    dna_1 = get_cpm(line_split, dna_1_col)
    dna_2 = get_cpm(line_split, dna_2_col)
    stats["min_dna_cpm"] = min(dna_1, dna_2)
    stats["dna_diffs"][0] = get_ratios(dna_1, dna_2)

    dna_1_unt = [get_cpm(line_split, columns) for columns in dna_1_unt_cols]
    dna_2_unt = [get_cpm(line_split, columns) for columns in dna_2_unt_cols]
    dna_unt_ratios = [get_ratios(dna_1_unt[idx], dna_2_unt[idx]) for idx in range(0, len(dna_1_unt))]

    dna_1_t = [get_cpm(line_split, columns) for columns in dna_1_t_cols]
    dna_2_t = [get_cpm(line_split, columns) for columns in dna_2_t_cols]
    dna_t_ratios = [get_ratios(dna_1_t[idx], dna_2_t[idx]) for idx in range(0, len(dna_1_t))]

    stats["dna_diffs"][1] = min(dna_unt_ratios)
    stats["dna_diffs"][2] = min(dna_t_ratios)

    return stats


def get_best_max(val_list):
    return "%.3f" % max([float(item) for item in val_list])


def get_best_min(val_list):
    return "%.3f" % min([float(item) for item in val_list])


def check_pass(dict, bc1, bc2, sc_bc1, sc_bc2):
    total_passed = 0
    total_passed_lowoccurence = 0
    pair = (bc1, bc2, sc_bc1, sc_bc2)
    comment = "low_occurence_failure"

    co_occ_list = list()
    min_cooccurence_list = list()
    rna_ratio_list = list()
    dna_ratio_list = list()
    rna_cpm_list = list()
    dna_cpm_list = list()

    for samp in sorted(dict[pair]):
        co_occ_list.append("%.3f" % dict[pair][samp]["cooccur"])
        min_cooccurence_list.append("%.3f" % dict[pair][samp]["cooccur_frac_to_min"])
        rna_ratio_list.append("%.3f" % dict[pair][samp]["rna_diff"])
        dna_ratio_list.append("%.3f" % min(dict[pair][samp]["dna_diffs"]))
        rna_cpm_list.append(str(dict[pair][samp]["min_rna_cpm"]))
        dna_cpm_list.append(str(dict[pair][samp]["min_dna_cpm"]))

        if dict[pair][samp]["cooccur"] >= absolute_min_occurence:
            if dict[pair][samp]["cooccur_frac_to_min"] >= min_frac_cooccurence:
                if min(dict[pair][samp]["dna_diffs"]) < ratio_cutoff:
                    if dict[pair][samp]["cooccur"] >= min_cooccurence:
                        total_passed += 1
                    else:
                        total_passed_lowoccurence += 1
                else:
                    comment = "dna_ratio_cutoff_failed"
            else:
                comment = "min_frac_occurence_failed"
        else:
            comment = "occurence_failed"

    if total_passed >= 2:
        comment = "passed"
    elif total_passed_lowoccurence >= 2 or total_passed >= 1:
        comment = "passed_with_low_occurence"

    return "\t".join(map(str, [";".join(sorted(dict[pair])), len(dict[pair]), total_passed, comment,
                               ";".join(dna_ratio_list), get_best_min(dna_ratio_list), ";".join(rna_ratio_list),
                               get_best_min(rna_ratio_list), ";".join(co_occ_list), get_best_max(co_occ_list),
                               ";".join(min_cooccurence_list), get_best_max(min_cooccurence_list),
                               ";".join(dna_cpm_list), ";".join(rna_cpm_list)]))


with open(input_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")

        samp = ls[sample_col]
        bc1 = ls[xbc1_col]
        bc2 = ls[xbc2_col]
        sc_b1 = ls[xbc1_sc_col]
        sc_b2 = ls[xbc2_sc_col]

        dict_key = (bc1, bc2, sc_b1, sc_b2)

        if dict_key not in pairs_sample_dict:
            pairs_sample_dict[dict_key] = dict()

        pairs_sample_dict[dict_key][samp] = return_stats(ls)
fh.close()


print("xbc1\txbc2\tis_sc1\tis_sc2\tsamples\tnum_samp\tnum_samp_passed\tcomment\tdna_pct_diff\tbest_dna_pct_diff\trna_pct_diff\t"
      "best_rna_pct_diff\tcoccur\tbest_cooccur\tcoccur_pct_diff\tbest_coccur_pct_diff\tmin_dna_cpms\tmin_rna_cpms")

for (barcode_1, barcode_2, sc_bc1, sc_bc2) in pairs_sample_dict:

    passed_stats = check_pass(pairs_sample_dict, barcode_1, barcode_2, sc_bc1, sc_bc2)

    print("\t".join(map(str, [barcode_1, barcode_2, sc_bc1, sc_bc2, passed_stats])))





