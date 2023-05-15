import sys

# Step2: Input- Counts file from summary. Merge similar looking expressed barcodes and merge them

input_file = sys.argv[1]

tenx_exp_dict = dict()
exp_set = set()
tenx_set = set()
prev_samp = "NA"

sample_col = 0
tenx_col = 1
exp_col = 2
count_col = 4
ref_tenx_col = 3

hd_cutoff = 3
cw_cutoff = 0.8


# Get basic hamming distance between 2 sequences
def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def merge_barcodes(ip_dict):
    cleaned_dict = dict()
    total = 0
    for exp_barcode in ip_dict:
        barcode_added = False
        for cleaned_expbc in cleaned_dict:
            if hamming_distance(cleaned_expbc, exp_barcode) <= hd_cutoff:
                cleaned_dict[cleaned_expbc] += ip_dict[exp_barcode]
                total += ip_dict[exp_barcode]
                barcode_added = True

        if not barcode_added:
            cleaned_dict[exp_barcode] = ip_dict[exp_barcode]
            total += ip_dict[exp_barcode]

    return cleaned_dict, total


def merge_and_print(ip_dict, sample):
    grand_total = 0
    cw = 0
    cleaned_dict = dict()

    for tenx_bc in ip_dict:
        out_dict, total = merge_barcodes(ip_dict[tenx_bc])
        grand_total += total
        for exp_bc in out_dict:
            cleaned_dict[(tenx_bc, exp_bc)] = out_dict[exp_bc]

    for (tenx_bc_e, exp_bc_e), counts_e in sorted(cleaned_dict.items(), key=lambda x: x[1], reverse=True):
        fraction = counts_e/grand_total
        cw += fraction
        if cw <= cw_cutoff:
            print("\t".join([sample, tenx_bc_e, exp_bc_e, "%.3f" % counts_e, "%.3f" % (fraction * 10 ** 6), "%.3f" % cw]))


print("\t".join(["sample", "tenx_bc", "exp_bc", "counts", "cpm", "cw"]))

with open(input_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        if prev_samp != ls[sample_col] and prev_samp != "NA":
            merge_and_print(tenx_exp_dict, prev_samp)
            tenx_exp_dict = dict()
            exp_set = set()
            tenx_set = set()
        tenx_bc = ls[tenx_col]
        exp_bc = ls[exp_col]
        counts = ls[count_col]
        prev_samp = ls[sample_col]
        ref_tenx = ls[ref_tenx_col]

        if ref_tenx == "False":
            continue

        if tenx_bc not in tenx_exp_dict:
            tenx_exp_dict[tenx_bc] = dict()

        tenx_exp_dict[tenx_bc][exp_bc] = int(counts)

    merge_and_print(tenx_exp_dict, prev_samp)
fh.close()

