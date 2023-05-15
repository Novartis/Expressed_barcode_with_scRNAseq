import sys

# Reads in the tenx_recounts_file, CR_tp_file
tenx_recounts_file = sys.argv[1]
CR_tp_file = sys.argv[2]
#tenx_distance_file = sys.argv[3]
#dna_data = sys.argv[4]

CR_10x_tp_barcodes = dict()

tenx_recounts_dict = dict()
final_associations = dict()

exp_tenx_dict = dict()
tenx_exp_dict = dict()
exp_count_dict = dict()

count_exp_barcode = dict()

cpm_cutoff = 2.5
highfreq_cutoff = 20
count_cutoff = 10

annotation_map = dict()


# Get the true positive 10x barcodes for the samples
with open(CR_tp_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split(",")
        name_condition = "_".join(ls[0].split("_")[0:2])
        tenx_barcode = ls[1]
        if name_condition not in CR_10x_tp_barcodes:
            CR_10x_tp_barcodes[name_condition] = set()
        CR_10x_tp_barcodes[name_condition].add(tenx_barcode)
fh.close()


# Open the recounts file and start munching
# The idea is that you store both cDNA and library information for a given sample
# We then try and use both to associate 10x and xBC

with open(tenx_recounts_file) as fh:
    fh.readline()
    for line in fh:

        ls = line.strip().split("\t")

        sample = ls[0].split("_")
        name_condition = "_".join(sample[0:2])
        source = sample[2]
        tenxb = ls[1]
        exp = ls[2]
        cpm = float(ls[4])
        count = int(ls[3])

        if tenxb in CR_10x_tp_barcodes[name_condition]:
            if name_condition not in exp_tenx_dict:
                exp_tenx_dict[name_condition] = dict()
                tenx_exp_dict[name_condition] = dict()
                exp_count_dict[name_condition] = dict()

            if source not in exp_tenx_dict[name_condition]:
                exp_tenx_dict[name_condition][source] = dict()
                tenx_exp_dict[name_condition][source] = dict()
                exp_count_dict[name_condition][source] = dict()

            if exp not in exp_tenx_dict[name_condition][source]:
                exp_tenx_dict[name_condition][source][exp] = set()

            if tenxb not in tenx_exp_dict[name_condition][source]:
                tenx_exp_dict[name_condition][source][tenxb] = set()

            # Only associate expressed barcodes and "True" 10x barcodes
            exp_tenx_dict[name_condition][source][exp].add(tenxb)
            tenx_exp_dict[name_condition][source][tenxb].add(exp)
            exp_count_dict[name_condition][source][(exp, tenxb)] = count
fh.close()


# Re-assigning 10x barcode to exp based on CR truth and uniqueness of 10x to that expressed
for cl_cnd in exp_tenx_dict:
    final_associations[cl_cnd] = dict()
    count_exp_barcode[cl_cnd] = 0

    # For a give cell line, condition, get all possible expressed barcodes
    exp_barcodes = set(list(exp_tenx_dict[cl_cnd]["cDNA"].keys()) + list(exp_tenx_dict[cl_cnd]["library"].keys()))

    for exp_barcode in exp_barcodes:

        final_associations[cl_cnd][exp_barcode] = list()

        tenx_barcodes_cDNA = list(exp_tenx_dict[cl_cnd]["cDNA"][exp_barcode]) if exp_barcode in exp_tenx_dict[cl_cnd]["cDNA"] else list()
        tenx_barcodes_library = list(exp_tenx_dict[cl_cnd]["library"][exp_barcode]) if exp_barcode in exp_tenx_dict[cl_cnd]["library"] else list()
        tenx_barcode_set = set(tenx_barcodes_cDNA + tenx_barcodes_library)

        for tenx_barcode in tenx_barcode_set:

            cnt_cDNA =  exp_count_dict[cl_cnd]["cDNA"][(exp_barcode, tenx_barcode)] if (exp_barcode, tenx_barcode) in exp_count_dict[cl_cnd]["cDNA"] else 0
            cnt_lib = exp_count_dict[cl_cnd]["library"][(exp_barcode, tenx_barcode)] if (exp_barcode, tenx_barcode) in exp_count_dict[cl_cnd]["library"] else 0

            if tenx_barcode in tenx_barcodes_cDNA and tenx_barcode in tenx_barcodes_library:
                if exp_barcode in tenx_exp_dict[cl_cnd]["cDNA"][tenx_barcode] and \
                        exp_barcode in tenx_exp_dict[cl_cnd]["library"][tenx_barcode]:
                    final_associations[cl_cnd][exp_barcode].append((tenx_barcode, cnt_cDNA + cnt_lib))
                    count_exp_barcode[cl_cnd] += cnt_cDNA + cnt_lib
                elif tenx_barcode in tenx_exp_dict[cl_cnd]["cDNA"] and \
                        len(tenx_exp_dict[cl_cnd]["cDNA"][tenx_barcode]) == 1 and \
                        exp_barcode in tenx_exp_dict[cl_cnd]["cDNA"][tenx_barcode]:
                    final_associations[cl_cnd][exp_barcode].append((tenx_barcode, cnt_cDNA + cnt_lib))
                    count_exp_barcode[cl_cnd] += cnt_cDNA + cnt_lib

            elif tenx_barcode in tenx_barcodes_cDNA:
                if cnt_cDNA >= count_cutoff:
                    final_associations[cl_cnd][exp_barcode].append((tenx_barcode, cnt_cDNA + cnt_lib))
                    count_exp_barcode[cl_cnd] += cnt_cDNA + cnt_lib
                elif tenx_barcode in tenx_exp_dict[cl_cnd]["cDNA"] and \
                        len(tenx_exp_dict[cl_cnd]["cDNA"][tenx_barcode]) == 1 and \
                        exp_barcode in tenx_exp_dict[cl_cnd]["cDNA"][tenx_barcode]:
                    final_associations[cl_cnd][exp_barcode].append((tenx_barcode, cnt_cDNA + cnt_lib))
                    count_exp_barcode[cl_cnd] += cnt_cDNA + cnt_lib

            elif tenx_barcode in tenx_barcodes_library:
                if tenx_barcode in tenx_exp_dict[cl_cnd]["library"] and \
                        len(tenx_exp_dict[cl_cnd]["library"][tenx_barcode]) == 1 and \
                        exp_barcode in tenx_exp_dict[cl_cnd]["library"][tenx_barcode] and \
                        tenx_barcode not in tenx_exp_dict[cl_cnd]["cDNA"]:
                    final_associations[cl_cnd][exp_barcode].append((tenx_barcode, cnt_cDNA + cnt_lib))
                    count_exp_barcode[cl_cnd] += cnt_cDNA + cnt_lib

with open("tenx_recounts2.txt", 'w') as fh:
    print("sample\ttenx_bc\texp_bc\tcount\tcpm", file=fh)
    for cl_cnd in final_associations:
        total = count_exp_barcode[cl_cnd]
        for exp_barcode in final_associations[cl_cnd]:
            for tenx_barcode, count in final_associations[cl_cnd][exp_barcode]:
                cpm = ((count*1.0/total) * 10**6) if total != 0 else 0
                print("\t".join(map(str, [cl_cnd, tenx_barcode, exp_barcode, count, "%.5f" % cpm])), file=fh)
fh.close()

