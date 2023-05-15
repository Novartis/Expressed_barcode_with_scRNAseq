import sys
import nwalign3 as nw
import os

# For reclustering, we are using 10x CR barcodes from Cell ranger assay
# We are also making sure that the 10x barcode in dominant in only that sample

tenx_tall_skinny = sys.argv[1]   # /da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20190421_Jen_Dave_Matt_RNA/output_v2/tenx_counts.txt
tenx_tp_barcodes = sys.argv[2]   # /da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20190421_Jen_Dave_Matt_RNA/output_v2/CR_10x_barcodes.txt
tenx_dna_barcodes = sys.argv[3] if len(sys.argv) == 4 else ""     # /da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20191219_4cellline_DNA/output/counts_summary.txt

#tenx_dominant_barcodes = sys.argv[3] if len(sys.argv) > 3 else "NA"  # /da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20190421_Jen_Dave_Matt_RNA/output_v2/tenx_dominant_samples.txt

# python $BITBUCKET/oncp-expressed/recluster_barcodes.py tenx_counts.txt ../CR_10x_barcodes ../../20191219_4cellline_DNA/output/counts_summary.txt

annotation_map = {"MGH707_tx_1": "MGH001_EGF816_14day_R1",
                  "MGH707_tx_2": "MGH001_EGF816_14day_R2",
                  "MGH707_ctrl_1": "MGH001_CTRL_14day_R1",
                  "MGH707_ctrl_2": "MGH001_CTRL_14day_R2"}

tenx_dict = dict()
cw_cutoff = 0.85
tenx_reference = True
dna_filter = True
dna_count_cutoff = 3

sample_col = 0
tenx_col = 1
exp_col = 2
tenx_ref_col = 3
count_col = 4
cw_col = 6
min_ham = 10  # 5

count_frac_cufoff = 0.01
count_cutoff = 1

tenx_barcode_dict = dict()
exp_set = set()

xBC_dna_barcode_dict = dict()
CR_10x_tp_barcodes = dict()


def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        print((s1, s2))
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def hamming_distance_with_indels(s1, s2, algorithm="global"):
    # This gives the best gapped alignment between 2 sequences. Match is scored 1 and
    # mis-match is scored -1. So the score should give the gapped hamming distance
    ham_tuple = nw.global_align(s1, s2, gap_open=-10, gap_extend=-5, match=1)
    return hamming_distance(ham_tuple[0], ham_tuple[1])


def merge_barcodes(tenx_dict):
    merged_tenx_dict = dict()

    for txbc in tenx_dict:
        sorted_exp_dict = [(k, tenx_dict[txbc][k]) for k in sorted(tenx_dict[txbc], key=tenx_dict[txbc].get, reverse=True)]

        original = ""
        merged_tenx_dict[txbc] = dict()
        for (exp_barcode, count) in sorted_exp_dict:
            merged = False
            if original == "":
                original = exp_barcode
                merged_tenx_dict[txbc][exp_barcode] = count
                original_count = count
            else:
                for exp_code in merged_tenx_dict[txbc]:
                    if hamming_distance(exp_code, exp_barcode) <= min_ham:
                        merged_tenx_dict[txbc][exp_code] += count
                        merged = True
                        break
                    elif hamming_distance_with_indels(exp_code, exp_barcode) <= min_ham:
                        merged_tenx_dict[txbc][exp_code] += count
                        merged = True
                        break
                if not merged:
                    #if count/original_count >= count_frac_cufoff:
                    merged_tenx_dict[txbc][exp_barcode] = count
    return merged_tenx_dict


def print_dict(tenx_dict, sample, total):
    for tenx in tenx_dict:
        for exp in tenx_dict[tenx]:
            print("\t".join(map(str, [sample, tenx, exp, tenx_dict[tenx][exp], "%.5f" % ((tenx_dict[tenx][exp] * 10.0**6)/total)])))


if os.path.isfile(tenx_dna_barcodes):
    with open(tenx_dna_barcodes) as fh:
        fh.readline()
        for line in fh:
            ls = line.strip().split("\t")
            if ls[0] not in xBC_dna_barcode_dict:
                xBC_dna_barcode_dict[ls[0]] = set()
            if float(ls[3]) >= dna_count_cutoff:
                xBC_dna_barcode_dict[ls[0]].add(ls[2])
    fh.close()


with open(tenx_tp_barcodes) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split(",")
        sample_name = annotation_map[ls[0]] if ls[0] in annotation_map else ls[0]
        tenx_barcode = ls[1]
        if sample_name not in CR_10x_tp_barcodes:
            CR_10x_tp_barcodes[sample_name] = set()
        CR_10x_tp_barcodes[sample_name].add(tenx_barcode.upper())
fh.close()


with open(tenx_tall_skinny) as fh:
    fh.readline()
    sample = ""
    total = 0
    print("sample\ttenx_bc\texp_bc\tcount\tcpm")
    for line in fh:
        ls = line.strip().split("\t")
        # Remove all combinations seen at a single count
        if float(ls[count_col]) <= count_cutoff:
            continue
        if ls[sample_col] not in tenx_barcode_dict:
            tenx_barcode_dict[ls[sample_col]] = dict()
            tenx_dict[ls[sample_col]] = dict()
        if sample != ls[sample_col]:
            if sample != "":
                tenx_merged_dict = merge_barcodes(tenx_dict[sample])
                for txb in tenx_merged_dict:
                    for exb in tenx_merged_dict[txb]:
                        if txb not in tenx_barcode_dict[sample]:
                            tenx_barcode_dict[sample][txb] = set()
                        tenx_barcode_dict[sample][txb].add(exb)
                print_dict(tenx_merged_dict, sample, total)
            total = 0
        if ls[tenx_col] not in CR_10x_tp_barcodes[ls[sample_col]]:
            sample = ls[sample_col]
            continue
        if ls[tenx_col] not in tenx_dict[ls[sample_col]]:
            tenx_dict[ls[sample_col]][ls[tenx_col]] = dict()
        # dna_sample = "_".join(ls[sample_col].split("_")[0:2])
        # if len(xBC_dna_barcode_dict) > 0 and ls[exp_col] not in xBC_dna_barcode_dict[dna_sample] or len(xBC_dna_barcode_dict) == 0:
        #     continue
        tenx_dict[ls[sample_col]][ls[tenx_col]][ls[exp_col]] = float(ls[count_col])
        total += float(ls[count_col])
        sample = ls[sample_col]
fh.close()


tenx_merged_dict = merge_barcodes(tenx_dict[sample])
for txb in tenx_merged_dict:
    for exb in tenx_merged_dict[txb]:
        if txb not in tenx_barcode_dict[sample]:
            tenx_barcode_dict[sample][txb] = set()
        tenx_barcode_dict[sample][txb].add(exb)
print_dict(tenx_merged_dict, sample, total)

# Get the cooccurence only within samples
printed_set = dict()

with open("co_occurence_summary.txt", 'w') as fhw:
    print("\t".join(["sample", "exp_pairs", "seen_together"]), file=fhw)
    for sample in tenx_barcode_dict:
        printed_set[sample] = set()
        co_occurence_pair = dict()
        for co_occurence_set in tenx_barcode_dict[sample].values():
            if len(co_occurence_set) > 0:
                for item1 in co_occurence_set:
                    exp_set.add(item1)
                    for item2 in co_occurence_set:
                        if (item1, item2) in co_occurence_pair:
                            co_occurence_pair[(item1, item2)] += 1
                        elif (item2, item1) in co_occurence_pair:
                            co_occurence_pair[(item2, item1)] += 1
                        else:
                            co_occurence_pair[(item1, item2)] = 1

        if not os.path.isdir("co_occurence"):
            os.mkdir("co_occurence")
        with open("co_occurence/co_occurence_"+sample+".txt", 'w') as fh:
            print_str = "barcode"
            for exp in sorted(exp_set):
                print_str += ("\t" + exp)
            print(print_str, file=fh)

            for exp1 in sorted(exp_set):
                print_str = exp1
                for exp2 in sorted(exp_set):
                    if (exp1, exp2) not in co_occurence_pair and (exp2, exp1) not in co_occurence_pair:
                        print_str += "\t0"
                    else:
                        counts = str((co_occurence_pair[(exp1, exp2)]) if (exp1, exp2) in co_occurence_pair else co_occurence_pair[(exp2, exp1)])
                        pairs = "_".join(sorted((exp1, exp2)))
                        if pairs not in printed_set[sample]:
                            print("\t".join([sample, pairs, counts]), file=fhw)
                        printed_set[sample].add(pairs)
                        print_str += "\t"+counts
                print(print_str, file=fh)
        fh.close()
fhw.close()
