import sys

co_occurence_summary = sys.argv[1]
single_clones = sys.argv[2]
exp_tenx_associations = sys.argv[3]
tenx_recounts2 = sys.argv[4]
dna_counts = sys.argv[5]

# python $BITBUCKET/oncp-expressed/annotate_cooccurence.py co_occurence_summary2.txt ../../20200415_long_term_resistance_invivo_DNA/output/counts_summary_6_clones.txt exp_tenx_counts.txt tenx_recounts2_sorted.txt ../../20191219_4cellline_DNA/output/counts_summary.txt

sc_clone_dict = dict()
sc_barcode_set = set()
samp_exp_counts = dict()
rna_exp_counts = dict()
dna_exp_counts = dict()

with open(single_clones) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        if ls[0] not in sc_clone_dict:
            sc_clone_dict[ls[2]] = list()
        sc_clone_dict[ls[2]].append([ls[0], "%.2f" % float(ls[4])])
        sc_barcode_set.add(ls[2])
fh.close()

with open(exp_tenx_associations) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        if ls[0] not in samp_exp_counts:
            samp_exp_counts[ls[0]] = dict()
        samp_exp_counts[ls[0]][ls[1]] = ls[2]
fh.close()


with open(tenx_recounts2) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[0]
        exp_bc = ls[2]
        if sample not in rna_exp_counts:
            rna_exp_counts[sample] = dict()
        if exp_bc not in rna_exp_counts[sample]:
            rna_exp_counts[sample][exp_bc] = 0
        rna_exp_counts[sample][exp_bc] += float(ls[3])
fh.close()


with open(dna_counts) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[0]
        exp_bc = ls[2]
        if sample not in dna_exp_counts:
            dna_exp_counts[sample] = dict()
        if exp_bc not in dna_exp_counts[sample]:
            dna_exp_counts[sample][exp_bc] = 0
        dna_exp_counts[sample][exp_bc] += float(ls[4])*(10**6)
fh.close()


print("sample\txbc1\txbc2\txbc1_pc9_single\txbc2_pc9_single\txbc1_occurence_cleaned\txbc2_occurence_cleaned"
      "\txbc1_rna_cpm\txbc2_rna_cpm\txbc1_dna_cpm\txbc2_dna_cpm\tco_occurence\tpc9_single_1\tdna_frac_in_pc9_s1"
      "\tpc9_single_2\tdna_frac_in_pc9_s2")

with open(co_occurence_summary) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[0]
        bcs = sorted(ls[1].split("_"))

        if bcs[0] == bcs[1]:
            continue

        bc_presence = [barcode in sc_barcode_set for barcode in bcs]

        if any(bc_presence) and "PC9" in sample:
            clone_name = list()
            for bc in bcs:
                clone_name.append(("\t".join(sc_clone_dict[bc][0])) if bc in sc_clone_dict else "NA\t0.0")
        else:
            clone_name = ["NA\t0.0", "NA\t0.0"]

        samp_exp = [(samp_exp_counts[sample][bc] if bc in samp_exp_counts[sample] else "0") for bc in bcs]
        rna_repr = [(str(rna_exp_counts[sample][bc]) if bc in rna_exp_counts[sample] else "1") for bc in bcs]
        dna_repr = [("%.3f" % (dna_exp_counts[sample][bc]) if bc in dna_exp_counts[sample] else "1") for bc in bcs]

        print("\t".join([sample, "\t".join(bcs), "\t".join(map(str, bc_presence)), "\t".join(samp_exp),
                         "\t".join(rna_repr), "\t".join(dna_repr), ls[2], "\t".join(clone_name)]))
fh.close()