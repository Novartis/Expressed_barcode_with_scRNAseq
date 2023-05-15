import sys

sensitivity_info = sys.argv[1] #log2fc_sensitivities_HCC827_filt.txt
tall_skinny = sys.argv[2]
cell_line = sys.argv[3]

cloneid_sensitivity = dict()

with open(sensitivity_info) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        if ls[3] == cell_line and ls[7] == "t0":
            cloneid_sensitivity[ls[0]] = ls[8]
fh.close()

prev_sample = ""
cw = 0
cw_cutoff = 10000000
tp_clone_dict = dict()

with open(tall_skinny) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split()
        sample = ls[0]
        if not sample.startswith(cell_line):
            continue
        cpm = float(ls[10])
        tp_treat = ls[6] + "_" + ls[7]
        if tp_treat == "-2_U":
            tp_treat = "0_U"
        if tp_treat not in tp_clone_dict:
            tp_clone_dict[tp_treat] = dict()
        cloneid = ls[9]
        if cpm == 0.0:
            continue
        if prev_sample != ls[0]:
            cw = 0
        prev_sample = ls[0]
        cw += cpm
        if cw <= cw_cutoff:
            if cloneid not in tp_clone_dict[tp_treat]:
                tp_clone_dict[tp_treat][cloneid] = cw
fh.close()


with open(tall_skinny) as fh:
    print(fh.readline().strip()+"\tcw\tsensitivity\tcleaned_tp_tr")
    for line in fh:
        ls = line.strip().split()
        sample = ls[0]
        if not sample.startswith(cell_line):
            continue
        cpm = float(ls[10])
        cloneid = ls[9]
        tp_treat = ls[6] + "_" + ls[7]
        if tp_treat == "-2_U":
            tp_treat = "0_U"
        sensitivity = "NA" if cloneid not in cloneid_sensitivity else cloneid_sensitivity[cloneid]
        if cloneid in tp_clone_dict[tp_treat]:
            print(line.strip() + "\t" + "%.3f" % tp_clone_dict[tp_treat][cloneid] + "\t" + sensitivity + "\t" + tp_treat)
fh.close()
