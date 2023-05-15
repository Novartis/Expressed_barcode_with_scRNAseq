import sys

cleaned_file = sys.argv[1]
combined_barcode_file = sys.argv[2]

# python /da/onc/krishvi7/bitbucket/oncp-expressed/multiple_integrations/create_clonal_map.py tenx_recounts2_sorted.txt cooccurence_bathroom_sink_combined_barcodes.txt > tenx_recounts2_sorted_clonemap.txt

bc_pair_dict = dict()
pairs = list()

clone_bc_dict = dict()
bc_clone_dict = dict()
barcode_set = set()
low_occurence_pairs = list()
clone_number = 0

best_cooccur_col = 13
ratio_best_cooccur = 0.2


def get_clone_id(bc):
    for clone in clone_bc_dict:
        if bc in clone_bc_dict[clone]:
            return clone
    return "NA"


with open(combined_barcode_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        bc1 = ls[0]
        bc2 = ls[1]
        result = ls[7]
        best_cooccur = float(ls[best_cooccur_col])
        if result == "passed":

            if bc1 not in bc_pair_dict:
                bc_pair_dict[bc1] = dict()
            if bc2 not in bc_pair_dict:
                bc_pair_dict[bc2] = dict()

            bc_pair_dict[bc1][bc2] = best_cooccur
            bc_pair_dict[bc2][bc1] = best_cooccur

        if result == "passed_with_low_occurence":
            low_occurence_pairs.append(sorted((bc1, bc2)))
fh.close()

bc_seen = set()

for bc in bc_pair_dict:
    good_matches = list()
    best_value = 10000
    for key, value in sorted(bc_pair_dict[bc].items(), key=lambda item: item[1], reverse=True):
        if len(good_matches) == 0:
            best_value = value
            good_matches.append((key, value))
        elif value / best_value >= ratio_best_cooccur:
            if sorted((bc, key)) not in low_occurence_pairs:
                good_matches.append((key, value))

    for bc2 in good_matches:
        pair = sorted([bc, bc2[0]])
        if pair not in pairs:
            pairs.append(pair)

for bc1, bc2 in pairs:
    bc1_prev = True
    bc2_prev = True

    if bc1 in bc_pair_dict:
        if bc2 not in bc_pair_dict:
            continue
    elif bc2 in bc_pair_dict:
        if bc1 not in bc_pair_dict:
            continue

    if bc1 not in barcode_set:
        barcode_set.add(bc1)
        bc1_prev = False

    if bc2 not in barcode_set:
        barcode_set.add(bc2)
        bc2_prev = False

    if bc1_prev:
        clone = get_clone_id(bc1)
    elif bc2_prev:
        clone = get_clone_id(bc2)
    else:
        clone_number += 1
        clone = "clone" + str(clone_number).rjust(5, "0")
        clone_bc_dict[clone] = set()

    clone_bc_dict[clone].add(bc1)
    clone_bc_dict[clone].add(bc2)


with open(cleaned_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        bc = ls[2]

        clone = "NA"
        if bc in barcode_set:
            clone = get_clone_id(bc)
        else:
            clone_number += 1
            clone = "clone" + str(clone_number).rjust(5, "0")
            clone_bc_dict[clone] = set()
            barcode_set.add(bc)

        clone_bc_dict[clone].add(bc)
fh.close()


for clone in clone_bc_dict:
    for bc in clone_bc_dict[clone]:
        bc_clone_dict[bc] = clone

with open(cleaned_file) as fh:
    print(fh.readline().strip()+"\tcloneid")
    for line in fh:
        ls = line.strip().split("\t")
        bc = ls[2]
        print(line.strip()+"\t"+bc_clone_dict[bc])
fh.close()

with open("clone_map.txt", 'w') as fhw:
    print("xbc\tcloneid", file=fhw)
    for bc in bc_clone_dict:
        print(bc + "\t" + bc_clone_dict[bc], file=fhw)
fhw.close()