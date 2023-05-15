import sys

kitchen_sink_file = sys.argv[1]
subsidiary_file = sys.argv[2]

subsidiary_dict = dict()
timepoint_set = set()

# python $BITBUCKET/oncp-expressed/append_bathroom_sink.py co_occurence_summary2.txt exp_tenx_counts.txt tenx_recounts2.txt

with open(subsidiary_file) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[0]
        bc = ls[2]
        cpm = float(ls[4]) * 10 ** 6

        if not "00_U" in sample and not "02_T" in sample:
            continue

        timepoint_orig = "0" if "00_U" in sample else "14"
        sample_orig = sample.split("_")[0] + "_" + timepoint_orig

        timepoint_replicate = "_".join(sample.split("_")[1:4])
        timepoint_set.add(timepoint_replicate)

        if sample_orig not in subsidiary_dict:
            subsidiary_dict[sample_orig] = dict()

        if timepoint_replicate not in subsidiary_dict[sample_orig]:
            subsidiary_dict[sample_orig][timepoint_replicate] = dict()

        if cpm > 0:
            subsidiary_dict[sample_orig][timepoint_replicate][bc] = "%.3f" % cpm
fh.close()

with open(kitchen_sink_file) as fh:
    print_header = fh.readline().strip()
    for tp in sorted(timepoint_set):
        print_header += "\txbc1_" + tp + "_HR_cpm\txbc2_" + tp + "_HR_cpm"

    print(print_header)

    for line in fh:
        ls = line.strip().split("\t")
        samp_orig = ls[0]
        bc1 = ls[1]
        bc2 = ls[2]

        if samp_orig in subsidiary_dict:
            print_line = line.strip()

            for tp in sorted(timepoint_set):
                if tp in subsidiary_dict[samp_orig]:
                    print_line += "\t" + ((subsidiary_dict[samp_orig][tp][bc1]) if bc1 in subsidiary_dict[samp_orig][tp] else "0.0")
                    print_line += "\t" + ((subsidiary_dict[samp_orig][tp][bc2]) if bc2 in subsidiary_dict[samp_orig][tp] else "0.0")
                else:
                    print_line += "\tNA\tNA"

            print(print_line)
        else:
            nas = ["NA\tNA" for tp in sorted(timepoint_set)]
            print(line.strip()+"\t"+"\t".join(nas))
