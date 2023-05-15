__author__ = 'krishvi7'
import glob
import os.path

## Input file format (for individual sample)
# Barcode                         Count   Fraction
# ACACTCTGTGTGAGAGAGACTGTGTCACAC  14      0.00106772422209
# TGTGAGTGAGTCACTGACTCAGACAGAGAG  8       0.000610128126907
# ACTGTCACAGACTGTGTGTGTCTCAGTGAC  8       0.000610128126907
# TCTGACACAGAGACAGTCAGTCAGTGACAG  7       0.000533862111043

## output file format
# Sample  Barcode                         Count   Fraction
# Sample1 ACACTCTGTGTGAGAGAGACTGTGTCACAC  14      0.00106772422209
# Sample1 TGTGAGTGAGTCACTGACTCAGACAGAGAG  8       0.000610128126907
# Sample1 ACTGTCACAGACTGTGTGTGTCTCAGTGAC  8       0.000610128126907
# Sample1 TCTGACACAGAGACAGTCAGTCAGTGACAG  7       0.000533862111043
# Sample1 AGAGAGTGAGAGTCAGAGTGAGAGTGTCTG  5       0.000381330079317
# Sample1 ACAGTGTCACACTGTGACTGTCTCTGAGAC  5       0.000381330079317
# Sample1 AGACAGTGTCAGACAGTGAGTGTCACTGAG  5       0.000381330079317


def tall_skinny(output_directory, annotation_dict, annotation_header, tall_skinny_output):
    try:
        with open(output_directory+"/"+tall_skinny_output, 'w') as f:
            f.write("Sample\tBarcode\tCount\tFraction\t"+annotation_header)
            f.write("\n")
            for filename in glob.glob(output_directory+"/*.barcode.txt"):
                samp = os.path.basename(filename).replace(".barcode.txt", "")
                fh = open(filename)
                fh.readline()
                for line in fh:
                    sample_annt = annotation_dict[samp] if samp in annotation_dict else ""
                    f.write(samp+"\t"+line.strip()+"\t"+sample_annt)
                    f.write("\n")
                fh.close()
    except:
        raise

