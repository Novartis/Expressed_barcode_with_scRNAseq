author__ = 'krishvi7'
import math


def simplify_barcodes(barcode_file, barcode_output, cumulative_output):
    wealthresults = [0] * 101

    fh = open(barcode_file)
    curwealth = 0
    cumbc = 0
    linenum = 0
    total = 0

    try:
        with open(barcode_output, 'w') as f:
            f.write("Barcode\tCount\tFraction")
            f.write("\n")
            for line in fh:
                if "Total_good_reads" in line:
                    total = int(line.strip().split("\t")[1])
                elif "Total_less_than_two" in line:
                    total += int(line.strip().split("\t")[1])

                if "#" in line:
                    continue

                linenum += 1
                vals = line.strip().split()
                f.write(vals[0] + "\t" + vals[3] + "\t" + vals[5])
                f.write("\n")
                cumbc += int(vals[3])
                if 100 * cumbc / total >= curwealth:
                    for i in range(curwealth, int(math.floor(100 * cumbc / total)) + 1):
                        wealthresults[i] = linenum
                    curwealth = int(math.floor(100 * cumbc / total)) + 1

            for i in range(curwealth, 101):
                wealthresults[i] = ""

        with open(cumulative_output, 'w') as f:
            f.write("PercentageOfBarcodesRead\tUniqueBarcodesRequiredToReachPercentage")
            f.write("\n")
            for idx, val in enumerate(wealthresults):
                f.write(str(idx) + "\t" + str(val))
                f.write("\n")

    except:
        raise
