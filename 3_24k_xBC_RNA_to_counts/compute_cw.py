import math
import sys


def simplify_barcodes(tenx_file, stats_file):
    wealthresults = [0] * 101
    curwealth = 0
    cumbc = 0
    linenum = 0
    total = 0

    try:
        with open(stats_file) as fh:
            for line in fh:
                if "good" in line:
                    total += int(line.strip().split(":")[1])
        fh.close()

        with open(tenx_file) as fh:
            fh.readline()
            for line in fh:
                linenum += 1
                vals = line.strip().split("\t")
                cumbc += int(vals[2])
                if 100 * cumbc / total >= curwealth:
                    for i in range(curwealth, int(math.floor(100 * cumbc / total)) + 1):
                        wealthresults[i] = linenum
                    curwealth = int(math.floor(100 * cumbc / total)) + 1
        fh.close()

        sample = stats_file.replace("_stats.txt", "")

        for i in range(curwealth, 101):
            wealthresults[i] = ""

        print("Sample\tPercentageOfBarcodesRead\tUniqueBarcodesRequiredToReachPercentage")
        for idx, val in enumerate(wealthresults):
            print(sample + "\t" + str(idx) + "\t" + str(val))


    except:
        raise

simplify_barcodes(sys.argv[1], sys.argv[2])
