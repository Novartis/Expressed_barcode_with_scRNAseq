import sys

tall_skinny = sys.argv[1]
waterfall_file = sys.argv[2]
sample_col = 0
barcode_col = 1
freq_col = 2

barcode_cl_dict = dict()


with open(tall_skinny) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[sample_col]
        barcode = ls[barcode_col]
        fraction = float(ls[freq_col])
        #group = "_".join([ls[item] for item in [5, 8]])
        #sample = "_".join(ls[4].split("_")[0])
        #group = "_".join(ls[4].split("_")[0:2])


        if barcode not in barcode_cl_dict:
            #barcode_cl_dict[barcode] = [sample, group, "%.3f" % (fraction * 10**6)]
            barcode_cl_dict[barcode] = [sample, "%.3f" % (fraction * 10 ** 6)]
        else:
            #if fraction > float(barcode_cl_dict[barcode][2])/10**6:
                #barcode_cl_dict[barcode] = [sample, group, "%.3f" % (fraction * 10**6)]
            if fraction > float(barcode_cl_dict[barcode][1]) / 10 ** 6:
                barcode_cl_dict[barcode] = [sample, "%.3f" % (fraction * 10 ** 6)]
fh.close()

with open(waterfall_file) as fh:
    #print(fh.readline().strip()+"\tdominant_cellline\tdominant_condition\tmax_cpm")
    print(fh.readline().strip() + "\tdominant_cellline\tmax_cpm")
    for line in fh:
        ls = line.strip().split("\t")
        barcode = ls[barcode_col]
        print(line.strip()+"\t"+"\t".join(barcode_cl_dict[barcode]))
fh.close()