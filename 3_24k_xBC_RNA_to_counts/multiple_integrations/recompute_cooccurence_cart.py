import sys


# Step 4: Recompute the number of time 2 expressed barcodes are seen associated with a single cell
# python /da/onc/krishvi7/bitbucket/oncp-expressed/multiple_integrations/recompute_cooccurence_cart.py tenx_recounts2.annt.txt > occurence_ratio.txt
# After this. some work is done in R to separate singles, multiples and unsure ones

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

exp_barcode_associations = dict()

for sample in tenx_barcode_dict:
    for txb in tenx_barcode_dict[sample]:
        for barcode in tenx_barcode_dict[sample][txb]:
            if barcode not in exp_barcode_associations:
                exp_barcode_associations[barcode] = list()
            exp_barcode_associations[barcode].append(len(tenx_barcode_dict[sample][txb]) - 1)


print("exp_bc\tnum_occurence")
for barcode in exp_barcode_associations:
    print_list = [barcode+"\t"+str(num) for num in exp_barcode_associations[barcode]]
    print("\n".join(print_list))