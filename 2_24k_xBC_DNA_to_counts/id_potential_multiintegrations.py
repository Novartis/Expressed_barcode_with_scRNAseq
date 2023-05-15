import sys
from itertools import combinations
from scipy.stats import pearsonr

#input_file = sys.argv[1]
input_file = "/da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20191219_4cellline_DNA/output/counts_summary.txt"
number_of_integrations = [4]#, 3, 2]

cell_line_condition_bc_dict = dict()
cell_line_barcode_dict = dict()
condition_set = set()


def get_correlations(wrk_dict, bc_m):
    cond_cpm_dict = dict()
    corr_list = list()

    for cnd in condition_set:
        cond_cpm_dict[cnd] = list()
        for bc in sorted(bc_m):
            cpm = wrk_dict[cnd][bc] if bc in wrk_dict[cnd] else 0.0
            cond_cpm_dict[cnd].append(cpm)

    for cnd1 in cond_cpm_dict:
        for cnd2 in cond_cpm_dict:
            if cnd1 == cnd2:
                continue
            r, _ = pearsonr(cond_cpm_dict[cnd1], cond_cpm_dict[cnd2])

            corr_list.append(r)

    return corr_list


with open(input_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[0]
        barcode = ls[2]
        cpm = float(ls[4]) * 10 ** 6

        # CPM should be at least 1/8th of expected
        if cpm < 250:
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

#print(cell_line_condition_bc_dict, file=sys.stderr)

for cl in cell_line_condition_bc_dict:
    working_dict = cell_line_condition_bc_dict[cl]
    working_barcodes = cell_line_barcode_dict[cl]
    for i in number_of_integrations:
        potential_barcode_combinations = combinations(list(working_barcodes), i)
        for bc_multi in potential_barcode_combinations:
            correlation_list = get_correlations(working_dict, bc_multi)