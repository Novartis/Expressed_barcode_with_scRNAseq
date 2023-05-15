import sys
import numpy as np
from scipy import stats

cw_file = sys.argv[1]
sample_col = 0
cw_col = 1
unique_barcodes_col = 2

# group_col_name = ["group"]
# exclude_pattern = ["exvivo"]
# reference_group_pattern = "cell_pellet"
# annotation_col_name = ["model", "source", "tissue", "timepoint", "condition"]
# aggregator_col_name = ["model"]

#group_col_name = ["source", "collection_day", "tumor_growth", "condition", "mouse", "donor"]
group_col_name = ["GROUP"]
exclude_pattern = ["Plasmid"]
exclude_sample = ["HA-91-GK71", "ZA-90-UN79"]
reference_group_pattern = "Pellet"
annotation_col_name = ["CONDITION", "SOURCE", "GLUCOSE", "Tumor_weight/Cell_number"]
#annotation_col_name = ["Sample", "source", "spl.source", "cellline", "collection_day", "tumor_growth", "treatment_days",
#                       "condition", "group", "mouse", "donor"]
aggregator_col_name = ["CELL_LINE"]
#
# group_col_name = ["GROUP1"]
# exclude_pattern = []
# reference_group_pattern = "invitro.Baseline"
# annotation_col_name = ["CELL_LINE", "POOL_NAME", "SAMPLE_TYPE", "TREATMENT", "GLUCOSE"]
# aggregator_col_name = ["CELL_LINE"]
cw_range = [90, 95]


group_barcodenum_dict = dict()
reference_barcodenum_dict = dict()

with open(cw_file) as fh:
    header = fh.readline().strip().split("\t")
    group_cols = [header.index(gp_col_name) for gp_col_name in group_col_name]
    aggregator_cols = [header.index(agg_col_name) for agg_col_name in aggregator_col_name]
    annotations_cols = [header.index(ann_col_name) for ann_col_name in annotation_col_name]
    for line in fh:
        ls = line.strip().split("\t")
        cw = int(ls[cw_col])

        if ls[sample_col] in exclude_sample:
            continue

        if cw_range[0] <= cw <= cw_range[1]:
            unique_barcodes = ls[unique_barcodes_col]
            annotations = "\t".join([ls[annotations_col] for annotations_col in annotations_cols])
            agg = ".".join([ls[aggregator_col] for aggregator_col in aggregator_cols])
            group = ".".join([ls[gp_col] for gp_col in group_cols])
            if any([exclude in group for exclude in exclude_pattern]):
                continue

            if agg not in group_barcodenum_dict:
                group_barcodenum_dict[agg] = dict()
                reference_barcodenum_dict[agg] = dict()

            if group not in group_barcodenum_dict[agg]:
                group_barcodenum_dict[agg][group] = [dict(), ""]

            if cw not in group_barcodenum_dict[agg][group][0]:
                group_barcodenum_dict[agg][group][0][cw] = list()

            if cw not in reference_barcodenum_dict[agg]:
                reference_barcodenum_dict[agg][cw] = list()

            if reference_group_pattern in group:
                reference_barcodenum_dict[agg][cw].append(int(unique_barcodes))
            else:
                group_barcodenum_dict[agg][group][0][cw].append(int(unique_barcodes))
                group_barcodenum_dict[agg][group][1] = annotations

print(".".join(aggregator_col_name) + "\ttakerate_mean\ttakerate_stdmin\ttakerate_stdmean\t" +
      "\t".join(annotation_col_name))

for agg in group_barcodenum_dict:
    for group in group_barcodenum_dict[agg]:
        if agg in reference_barcodenum_dict and len(reference_barcodenum_dict[agg]) > 0:
            reference_barcodes_dict = reference_barcodenum_dict[agg]

            try:
                for cw in group_barcodenum_dict[agg][group]:
                    group_barcode_dict = group_barcodenum_dict[agg][group][0]

                    ratio_list = list()

                    for cw in group_barcode_dict:
                        if cw in reference_barcodes_dict:
                            for ref_item in reference_barcodes_dict[cw]:
                                for group_item in group_barcode_dict[cw]:
                                    ratio_list.append((group_item*100.0)/ref_item)

                    annt = group_barcodenum_dict[agg][group][1]

                    if "M9" in annt:
                        print(agg, file=sys.stderr)
                        print(group, file=sys.stderr)
                        print(group_barcode_dict, file=sys.stderr)
                        print(reference_barcodes_dict, file=sys.stderr)

                    if len(ratio_list) >= 2:
                        res_mean, res_var, res_std = stats.bayes_mvs(ratio_list, alpha=0.95)
                        mean_val = float(res_mean.statistic)
                        std_min = float(res_std.minmax[0])
                        std_max = float(res_std.minmax[1])

                        print(agg + "\t" + "%.5f" % mean_val + "\t" + "%.5f" % std_min + "\t" + "%.5f" % std_max + "\t"
                              + annt)
            except:
                print(type(str(res_mean.statistic)), file=sys.stderr)
                print(type(str(res_std.minmax[1])), file=sys.stderr)
                print(type(str(res_std.minmax[0])), file=sys.stderr)
                raise

