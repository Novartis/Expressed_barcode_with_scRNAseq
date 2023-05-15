import sys


# Reads in the tenx_counts_file, tenx_recounts_file, CR_tp_file

#tenx_counts_file = sys.argv[1]
tenx_recounts_file = sys.argv[1]
CR_tp_file = sys.argv[2]

CR_10x_tp_barcodes = dict()
tenx_recounts_dict = dict()
tenx_counts_dict = dict()

exp_tenx_dict = dict()
tenx_exp_dict = dict()
tenx_cpm_dict = dict()
exp_cpm_dict = dict()

cpm_cutoff = 2.5
highfreq_cutoff = 20

count_cutoff = 2

# annotation_map = {"MGH707_tx_1": "MGH001_EGF816_14day_R1",
#                   "MGH707_tx_2": "MGH001_EGF816_14day_R2",
#                   "MGH707_ctrl_1": "MGH001_CTRL_14day_R1",
#                   "MGH707_ctrl_2": "MGH001_CTRL_14day_R2"}


annotation_map = dict()


def get_presence(barcode, dictionary):
    if barcode in dictionary:
        return "1"
    return "0"


def get_high_freq(sample, tenx):
    if tenx in tenx_cpm_dict[sample] and tenx_cpm_dict[sample][tenx] >= highfreq_cutoff:
        return "1"
    return "0"


def get_high_freq_exp(sample, exp):
    if exp in exp_cpm_dict[sample] and exp_cpm_dict[sample][exp] >= highfreq_cutoff:
        return "1"
    return "0"


with open(CR_tp_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split(",")
        sample_name = annotation_map[ls[0]] if len(annotation_map) > 0 else ls[0]
        cndn = "_".join(sample_name.split("_")[0:2])
        tenx_barcode = ls[1]
        if sample_name not in CR_10x_tp_barcodes:
            CR_10x_tp_barcodes[sample_name] = set()
            CR_10x_tp_barcodes[cndn] = set()
        CR_10x_tp_barcodes[sample_name].add(tenx_barcode)
        CR_10x_tp_barcodes[cndn].add(tenx_barcode)
fh.close()


with open(tenx_recounts_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[0]
        tenxb = ls[1]
        exp = ls[2]
        cpm = float(ls[4])
        count = float(ls[3])
        if sample not in tenx_recounts_dict:
            tenx_recounts_dict[sample] = set()

        if sample not in exp_tenx_dict:
            exp_tenx_dict[sample] = dict()
            tenx_exp_dict[sample] = dict()
            tenx_cpm_dict[sample] = dict()
            exp_cpm_dict[sample] = dict()

        if exp not in exp_tenx_dict[sample]:
            exp_tenx_dict[sample][exp] = set()
        if cpm >= cpm_cutoff and count >= count_cutoff:
            exp_tenx_dict[sample][exp].add(tenxb)

        if tenxb not in tenx_exp_dict[sample]:
            tenx_exp_dict[sample][tenxb] = set()
        if cpm >= cpm_cutoff and count >= count_cutoff:
            tenx_exp_dict[sample][tenxb].add(exp)
            tenx_recounts_dict[sample].add(tenxb)

        if tenxb not in tenx_cpm_dict[sample]:
            tenx_cpm_dict[sample][tenxb] = cpm

        if exp not in exp_cpm_dict[sample]:
            exp_cpm_dict[sample][exp] = cpm

fh.close()


# with open(tenx_counts_file) as fh:
#     fh.readline()
#     for line in fh:
#         ls = line.strip().split("\t")
#         sample = ls[0]
#         tenxb = ls[1]
#         if sample not in tenx_counts_dict:
#             tenx_counts_dict[sample] = set()
#         tenx_counts_dict[sample].add(tenxb)
# fh.close()

with open("tenx_exp_counts.txt", 'w') as fh:
    print("sample\ttenx\texp_counts\thigh_freq", file=fh)
    for sample in tenx_exp_dict:
        for tenx in tenx_exp_dict[sample]:
            print(sample+"\t"+tenx+"\t"+str(len(tenx_exp_dict[sample][tenx]))+"\t"+get_high_freq(sample, tenx), file=fh)
fh.close()

with open("exp_tenx_counts.txt", 'w') as fh:
    print("sample\texp\ttenx_counts\thigh_freq", file=fh)
    for sample in exp_tenx_dict:
        for exp in exp_tenx_dict[sample]:
            print(sample+"\t"+exp+"\t"+str(len(exp_tenx_dict[sample][exp]))+"\t"+get_high_freq_exp(sample, exp), file=fh)
fh.close()

# with open("tenx_recounts_recovery.txt", 'w') as fh:
#     print("sample\ttenx\tCR\txBC\thigh_freq", file=fh)
#     for sample in CR_10x_tp_barcodes:
#         if sample in tenx_recounts_dict:
#             barcode_set = set(CR_10x_tp_barcodes[sample]).union(set(tenx_recounts_dict[sample]))
#             for tenx in barcode_set:
#                 print(sample+"\t"+tenx+"\t"+get_presence(tenx, CR_10x_tp_barcodes[sample])+"\t"+get_presence(tenx, tenx_recounts_dict[sample])
#                       + "\t"+get_high_freq(sample, tenx), file=fh)
# fh.close()
#
# with open("tenx_counts_recovery.txt", 'w') as fh:
#     print("sample\ttenx\tCR\txBC", file=fh)
#     for sample in CR_10x_tp_barcodes:
#         if sample in tenx_counts_dict:
#             barcode_set = set(CR_10x_tp_barcodes[sample]).union(set(tenx_counts_dict[sample]))
#             for tenx in barcode_set:
#                 print(sample + "\t" + tenx + "\t" + get_presence(tenx, CR_10x_tp_barcodes[sample]) + "\t" + get_presence(tenx, tenx_counts_dict[sample]), file=fh)
# fh.close()