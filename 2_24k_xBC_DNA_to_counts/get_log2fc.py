import sys
import math
import operator
import numpy as np

input_file = sys.argv[1] # /da/onc/BFx/research/krishvi7/barcoding/expressed_barcode/20190923_HCC827_drug_DNA/output/counts_mat_cloneid.txt
group_col_dict = dict()

name_sub_dict = {
    "EGF816_Single": "ES",
    "EGF816_INC280": "EI",
    "EGF816_Trametinib": "ET",
    "EGF816_BGJ398": "EB",
    "EGF816_Carboplatin-Pemetrexed": "EC",
    "Carboplatin-Pemetrexed_Single": "CS",
    "Day6_Untreated": "U6",
    "Day0_Untreated": "U0"
}

EGF_single = "ES"
Chemo_single = "CS"
INC_combo = "EI"
Tram_combo = "ET"
BGJ_combo = "EB"
Chemo_combo = "EC"
Day6_group = "U6"
Day0_group = "U0"

min_clone_cutoff = 100
clone_cutoff = 1000


def sort_representation_dict(rep_dict):
    return [k for k, v in sorted(rep_dict.items(), key=operator.itemgetter(1))]


def get_stats(arr, group):

    lt = [arr[f] for f in group_col_dict[group]]
    return max(np.median(lt), 0.1), max(min(lt), 0.1), lt


def get_significant_groups(arr, groups):
    representation_dict = dict()

    for grp in groups:
        if "Day" in grp:
            continue
        representation = get_stats(arr, grp)
        representation_dict[grp] = representation

    return sort_representation_dict(representation_dict)


with open("log2fc_cloneid.txt", 'w') as fhw:
    print("clonal_barcode\tgroup\treplicate\tlog2fc", file=fhw)
    with open(input_file) as fh:
        header = fh.readline().strip().split("\t")
        bc = header[0]
        for index in range(1, len(header)):
            group = name_sub_dict["_".join(header[index].split("_")[1:3])]
            if group not in group_col_dict:
                group_col_dict[group] = list()
            group_col_dict[group].append(index-1)
        for line in fh:
            ls = line.strip().split("\t")

            try:
                float_array = list(map(float, ls[1:]))

                # Check if it has decent representation in at least one of the samples
                if max(float_array) >= min_clone_cutoff:
                    # Has good representation to begin with
                    clone_annt = list()

                    key_groups = get_significant_groups(float_array, group_col_dict.keys())

                    median_Day0, min_Day0, Day0_replicate_cpm = get_stats(float_array, Day0_group)
                    median_EGF_single, min_EGF_single, EGF_replicate_cpm = get_stats(float_array, EGF_single)
                    median_Day6, min_Day6, Day6_replicate_cpm = get_stats(float_array, Day6_group)

                    for key_group_index in range(0, len(key_groups)):
                        median_key_grp, min_key_grp, key_grp_rep_cpm = get_stats(float_array, key_groups[key_group_index])

                        for cpm_index in range(0, len(key_grp_rep_cpm)):
                            rep_cpm = max(key_grp_rep_cpm[cpm_index], 0.5)
                            rep_fc = math.log(rep_cpm/min_Day6, 2)
                            print(ls[0]+"\t"+key_groups[key_group_index]+"\t"+str(cpm_index)+"\t"+"%.3f" % rep_fc, file=fhw)

            except:
                print(line.strip(), file=sys.stderr)
                raise
    fh.close()
fhw.close()