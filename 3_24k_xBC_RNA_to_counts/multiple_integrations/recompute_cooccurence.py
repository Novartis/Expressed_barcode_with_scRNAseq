import sys


# Step 4: Recompute the number of time 2 expressed barcodes are seen associated with a single cell

tenx_tall_skinny = sys.argv[1]   # tenx_recounts2_sorted.txt
tenx_barcode_dict = dict()

sample_col = 0
tenx_col = 1
exp_col = 2

with open(tenx_tall_skinny) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[sample_col]
        txb = ls[tenx_col]
        exb = ls[exp_col]

        if sample not in tenx_barcode_dict:
            tenx_barcode_dict[sample] = dict()
        if txb not in tenx_barcode_dict[sample]:
            tenx_barcode_dict[sample][txb] = set()

        tenx_barcode_dict[sample][txb].add(exb)
fh.close()

printed_set = dict()

with open("co_occurence_summary2.txt", 'w') as fhw:
    print("\t".join(["sample", "exp_pairs", "seen_together"]), file=fhw)
    for sample in tenx_barcode_dict:
        printed_set[sample] = set()
        exp_set = set()
        co_occurence_pair = dict()
        for co_occurence_set in tenx_barcode_dict[sample].values():
            if len(co_occurence_set) > 0:
                for item1 in co_occurence_set:
                    exp_set.add(item1)
                    for item2 in co_occurence_set:
                        if item1 == item2:
                            continue
                        if (item1, item2) in co_occurence_pair:
                            co_occurence_pair[(item1, item2)] += 1
                        elif (item2, item1) in co_occurence_pair:
                            co_occurence_pair[(item2, item1)] += 1
                        else:
                            co_occurence_pair[(item1, item2)] = 1

        for exp1 in sorted(exp_set):
            print_str = exp1
            for exp2 in sorted(exp_set):
                if (exp1, exp2) not in co_occurence_pair and (exp2, exp1) not in co_occurence_pair:
                    print_str += "\t0"
                else:
                    counts = "%i" % ((co_occurence_pair[(exp1, exp2)]) if (exp1, exp2) in co_occurence_pair else co_occurence_pair[(exp2, exp1)] / 2)
                    pairs = "_".join(sorted((exp1, exp2)))
                    if pairs not in printed_set[sample]:
                        print("\t".join([sample, pairs, counts]), file=fhw)
                    printed_set[sample].add(pairs)
                    print_str += "\t" + counts
fhw.close()
