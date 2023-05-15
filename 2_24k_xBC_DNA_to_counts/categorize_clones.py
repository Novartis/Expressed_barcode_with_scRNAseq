import sys
import math
import operator
import numpy as np

input_file = sys.argv[1]
group_col_dict = dict()

# EGF_single = "EGF816_Single"
# Chemo_single = "Carboplatin_Single"
# INC_combo = "EGF816_INC280"
# Tram_combo = "EGF816_Trametinib"
# BGJ_combo = "EGF816_BGJ398"
# Chemo_combo = "EGF816_Carboplatin"
# Day6_group = "Day6_Untreated"
# Day0_group = "Day0_Untreated"



name_sub_dict = {
    "EGF816_Single": "ES",
    "EGF816_INC280": "EI",
    "EGF816_Trametinib": "ET",
    "EGF816_BGJ398": "EB",
    "EGF816_Carboplatin": "EC",
    "Carboplatin_Single": "CS",
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

log2fc_cutoff_high = -6.5
log2fc_cutoff = -5
log2fc_cutoff_mild = -3

min_clone_cutoff = 100
clone_cutoff = 1000

base_category = 127

# weight_dict = {
#     "EGF816_Single": 2,
#     "EGF816_INC280": 4,
#     "EGF816_Trametinib": 8,
#     "EGF816_BGJ398": 16,
#     "EGF816_Carboplatin": 32,
#     "Carboplatin_Single": 64
# }

weight_dict = {
    "ES": 2,
    "EI": 4,
    "ET": 8,
    "EB": 16,
    "EC": 32,
    "CS": 64
}


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


with open("replicate_sensitivity.txt", 'w') as fhw:
    print("clonal_barcode\tgroup\treplicate\tlog2fc", file=fhw)
    print("clonal_barcode\tbinary_category\tcategory")
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
                clone_category = base_category

                float_array = list(map(float, ls[1:]))
                min_representation = 0
                insignificant = False
                print_all = False

                if ls[0] == "GTGTACACCATGCACAACCATG":
                    print_all = False

                # Check if it has decent representation in at least one of the samples
                if max(float_array) >= min_clone_cutoff:
                    # Has good representation to begin with
                    clone_annt = list()
                    if get_stats(float_array, Day6_group)[1] >= clone_cutoff:
                        clone_category = 129
                        #clone_annt.append("Veh")
                    else:
                        clone_category = 128
                        #clone_annt.append("nVeh")

                    key_groups = get_significant_groups(float_array, group_col_dict.keys())

                    median_EGF_single, min_EGF_single, EGF_replicate_cpm = get_stats(float_array, EGF_single)
                    median_Day6, min_Day6, Day6_replicate_cpm = get_stats(float_array, Day6_group)

                    max_Day6 = max(max(Day6_replicate_cpm), 0.5)

                    if math.log(median_EGF_single/median_Day6, 2) <= log2fc_cutoff_high and \
                            min_Day6 >= clone_cutoff:
                        clone_category |= 2
                        clone_annt.append("ES")

                    if print_all:
                        print(key_groups, file=sys.stderr)

                    for key_group_index in range(0, len(key_groups)):
                        median_key_grp, min_key_grp, key_grp_rep_cpm = get_stats(float_array, key_groups[key_group_index])
                        med_fc = math.log(median_key_grp / median_Day6, 2)
                        max_fc = math.log(min_key_grp / median_Day6, 2)

                        for cpm_index in range(0, len(key_grp_rep_cpm)):
                            rep_cpm = max(key_grp_rep_cpm[cpm_index], 0.5)
                            rep_fc = math.log(rep_cpm/min_Day6, 2)
                            if rep_fc <= log2fc_cutoff_high:
                                print(ls[0]+"\t"+key_groups[key_group_index]+"\t"+str(cpm_index)+"\t"+"%.3f" % rep_fc,
                                      file=fhw)

                        if key_group_index == 0:
                            if math.log(median_key_grp / median_Day6, 2) <= log2fc_cutoff and \
                                    min_Day6 >= clone_cutoff:
                                min_log2fc = med_fc
                                clone_category |= weight_dict[key_groups[key_group_index]]
                                clone_annt.append(key_groups[key_group_index])
                            else:
                                insignificant = True
                        else:
                            if print_all:
                                print(key_groups[key_group_index], file=sys.stderr)
                                print(weight_dict[key_groups[key_group_index]], file=sys.stderr)
                                print(med_fc, min_log2fc, file=sys.stderr)
                            if med_fc - min_log2fc <= log2fc_cutoff_mild or max_fc <= log2fc_cutoff_high:
                                clone_category |= weight_dict[key_groups[key_group_index]]
                                clone_annt.append(key_groups[key_group_index])

                        if insignificant:
                            break

                    if clone_category > 129:
                        print(ls[0]+"\t"+"{0:b}".format(clone_category)+"\t"+";".join(sorted(clone_annt)))

            except:
                print(line.strip(), file=sys.stderr)
                raise
    fh.close()
fhw.close()