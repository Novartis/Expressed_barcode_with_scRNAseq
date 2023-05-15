import sys
from itertools import combinations
from scipy.stats import pearsonr

input_file = sys.argv[1]
#input_file = "/da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20191219_4cellline_DNA/output/counts_summary.txt"
number_of_integrations = [2]

cell_line_condition_bc_dict = dict()
cell_line_barcode_dict = dict()
condition_set = set()
cell_lines_list = ["PC9", "HCC4006"]


def get_differences(wrk_dict, bc_m):
    cnd_pct_diff = dict()

    for cnd in condition_set:
        cpms = [wrk_dict[cnd][bc_m[0]] if bc_m[0] in wrk_dict[cnd] else 0.0,
                wrk_dict[cnd][bc_m[1]] if bc_m[1] in wrk_dict[cnd] else 0.0]
        if max(cpms) > 0 and min(cpms) > 0:
            pct_difference = abs(cpms[0]-cpms[1])/max(cpms) * 100
            cnd_pct_diff[cnd] = pct_difference, cpms[0], cpms[1]

    return cnd_pct_diff


with open(input_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[0]
        barcode = ls[2]
        cpm = float(ls[4]) * 10 ** 6

        # CPM should be at least 1/8th of expected
        if cpm < 250 or (len(cell_lines_list) > 0 and not any(item in sample for item in cell_lines_list)):
            continue

        sample_split = sample.split("_")
        cell_line = sample_split[0]
        condition = sample_split[1]

        condition_set.add(condition)

        if cell_line not in cell_line_condition_bc_dict:
            cell_line_condition_bc_dict[cell_line] = dict()
            cell_line_barcode_dict[cell_line] = set()

        if condition not in cell_line_condition_bc_dict[cell_line]:
            cell_line_condition_bc_dict[cell_line][condition] = dict()

        cell_line_condition_bc_dict[cell_line][condition][barcode] = cpm
        cell_line_barcode_dict[cell_line].add(barcode)
fh.close()

print("\t".join(["cell_line", "condition", "bc0", "bc1", "bc0_cpm", "bc1_cpm", "percent_difference_to_max"]))
for cl in cell_line_condition_bc_dict:
    working_dict = cell_line_condition_bc_dict[cl]
    working_barcodes = cell_line_barcode_dict[cl]
    for i in number_of_integrations:
        potential_barcode_combinations = combinations(list(working_barcodes), i)
        for bc_multi in potential_barcode_combinations:
            cnd_pct_diff = get_differences(working_dict, bc_multi)
            bc_0, bc_1 = sorted(bc_multi)
            for cnd in cnd_pct_diff:
                pct_diff = cnd_pct_diff[cnd][0]
                cpm_0 = cnd_pct_diff[cnd][1]
                cpm_1 = cnd_pct_diff[cnd][2]
                print("\t".join([cl, cnd, bc_0, bc_1, "%.3f" % cpm_0, "%.3f" % cpm_1, "%.3f" % pct_diff]))
