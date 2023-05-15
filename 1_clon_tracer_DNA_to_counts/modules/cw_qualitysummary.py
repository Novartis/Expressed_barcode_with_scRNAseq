__author__ = 'KRISHVI7'
import glob
import os.path
import sys


def cw_qualitysummary(output_directory, annotation_dict, annotation_header, cw_summary, quality_summary):

    try:

        with open(output_directory+"/"+quality_summary, 'w') as f:
            f.write("Sample\tCategory\tRead_Count\t"+annotation_header)
            f.write("\n")
            for filename in glob.glob(output_directory+"/*barcode_verbose.txt"):
                samp = os.path.basename(filename).replace(".barcode_verbose.txt", "")
                with open(filename) as fh:
                    for line in fh:
                        if "#" in line and ":" in line:
                            line_split = line.strip().split(":")
                            sample_annt = annotation_dict[samp] if samp in annotation_dict else ""
                            f.write(samp+"\t"+line_split[0][1:]+"\t"+line_split[1].strip()+"\t"+sample_annt)
                            f.write("\n")
                        else:
                            break
                fh.close()

        with open(output_directory+"/"+cw_summary, 'w') as f:
            f.write("Sample\tCumulativeWealth\tUniqueBarcodeCount\t"+annotation_header)
            f.write("\n")
            for filename in glob.glob(output_directory+"/*cumulative_wealth.txt"):
                samp = os.path.basename(filename).replace(".cumulative_wealth.txt", "")
                prev_unique_barcode_count = "0"
                fh = open(filename)
                fh.readline()
                for line in fh:
                    ls = line.strip().split("\t")
                    if len(ls) == 2:
                        prev_unique_barcode_count = ls[1]
                    sample_annt = annotation_dict[samp] if samp in annotation_dict else ""
                    f.write(samp+"\t"+ls[0]+"\t"+prev_unique_barcode_count+"\t"+sample_annt)
                    f.write("\n")
    except:
        raise
