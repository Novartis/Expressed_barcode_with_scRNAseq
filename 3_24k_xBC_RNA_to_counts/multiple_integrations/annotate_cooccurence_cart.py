import sys

co_occurence_summary = sys.argv[1]
exp_tenx_associations = sys.argv[2]
tenx_recounts2 = sys.argv[3]

# python $BITBUCKET/oncp-expressed/annotate_cooccurence_cart.py co_occurence_summary2.txt exp_tenx_counts.txt tenx_recounts2.txt > co

samp_exp_counts = dict()
rna_exp_counts = dict()

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


print("sample\txbc1\txbc2\txbc1_occurence_cleaned\txbc2_occurence_cleaned"
      "\txbc1_rna_cpm\txbc2_rna_cpm\tco_occurence")

with open(co_occurence_summary) as fh:
    fh.readline()
    for line in fh:
        ls = line.strip().split("\t")
        sample = ls[0]
        bcs = sorted(ls[1].split("_"))

        if bcs[0] == bcs[1]:
            continue

        samp_exp = [(samp_exp_counts[sample][bc] if bc in samp_exp_counts[sample] else "0") for bc in bcs]
        rna_repr = [(str(rna_exp_counts[sample][bc]) if bc in rna_exp_counts[sample] else "1") for bc in bcs]

        print("\t".join([sample, "\t".join(bcs), "\t".join(samp_exp),
                         "\t".join(rna_repr), ls[2]]))
fh.close()