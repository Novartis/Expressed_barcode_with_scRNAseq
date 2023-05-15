__author__ = 'krishvi7'
import os.path
import glob

## Input file format (for individual sample)
# Barcode                         Count   Fraction
# ACACTCTGTGTGAGAGAGACTGTGTCACAC  14      0.00106772422209
# TGTGAGTGAGTCACTGACTCAGACAGAGAG  8       0.000610128126907
# ACTGTCACAGACTGTGTGTGTCTCAGTGAC  8       0.000610128126907
# TCTGACACAGAGACAGTCAGTCAGTGACAG  7       0.000533862111043

## output file format- matrix fraction (similar for counts)
# Barcode                         Sample1 Sample2                 Sample3
# ACACACACAGACTGTCAGAGAGTCTGAGTG  0.0     0.0                     0.000147754137116
# ACACACACAGAGAGTGTGTGAGTGTGTGTG  0.0     0.00035253472467        0.0
# ACACACACAGTGAGAGACACAGTCAGTCTG  0.0     0.000282027779736       0.0
# ACACACACTGTGTGACAGAGACTCACTGTC  0.0     0.0                     0.000147754137116


def matricize(output_directory, matrix_count, matrix_fraction):

    # Initialize variables
    bcsamp2count = {}
    bcsamp2frac = {}
    samplist = []
    barcodelist = set()

    try:
        for filename in glob.glob(output_directory+"/*barcode.txt"):
            samp = os.path.basename(filename).replace(".barcode.txt", "")
            samplist.append(samp)
            fh = open(filename)
            fh.readline()
            rank = 0
            for line in fh:
                rank += + 1
                vals = line.strip().split()
                barcodelist.add(vals[0])
                keyval = (vals[0], samp)
                bcsamp2count[keyval] = vals[1]
                bcsamp2frac[keyval] = vals[2]

        with open(output_directory+"/"+matrix_count, 'w') as f:
            f.write("Barcode\t" + "\t".join(samplist))
            f.write("\n")
            for barcode in sorted(barcodelist):
                f.write(barcode + "\t" + "\t".join([bcsamp2count[(barcode, samp)] if (barcode, samp) in bcsamp2count
                                                        else "0" for samp in samplist]))
                f.write("\n")

        with open(output_directory+"/"+matrix_fraction, 'w') as f:
            f.write("Barcode\t" + "\t".join(samplist))
            f.write("\n")
            for barcode in sorted(barcodelist):
                f.write(barcode + "\t" + "\t".join(["%.3f" % (float(bcsamp2frac[(barcode, samp)]) * 10**6) if (barcode, samp) in bcsamp2frac
                                                       else "0.0" for samp in samplist]))
                f.write("\n")

    except:
        raise
