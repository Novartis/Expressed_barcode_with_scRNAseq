import subprocess
from pickle import FALSE

import classes

charvals = list()
charvals.append({'A': 0, 'T': 1})
charvals.append({'G': 0, 'C': 1})

intvals = list()
intvals.append({0: 'A', 1: 'T'})
intvals.append({0: 'G', 1: 'C'})

global totalcount


class zippedFile:

    def unzip(self, fnames):

        try:
            if fnames[0][-4:] == ".bz2":
                f = subprocess.Popen(["bunzip2", "-c"] + fnames, stdout=subprocess.PIPE, close_fds=True)
            elif fnames[0][-3:] == ".gz":
                f = subprocess.Popen(["zcat", "-f"] + fnames, stdout=subprocess.PIPE, close_fds=True)
            elif fnames[0][-6:] == ".fastq" or fnames[0][-3:] == ".fq":
                f = subprocess.Popen(["cat"] + fnames, stdout=subprocess.PIPE, close_fds=True)
            else:
                raise classes.Common.print_error("Unrecognized file format")
            try:
                for line in f.stdout:
                    yield line
            finally:
                f.terminate()
        except:
            raise

    def __init__(self, fname):
        self.fname = fname
        self.fnames = fname.split(",")
        self.generator = self.unzip(self.fnames)

    def close(self):
        self.generator.close()


class Barcode:
    global totalcount

    def __init__(self, barcode, count, nearestHam, nearestHamDist):
        self.barcodes = [barcode]
        self.counts = [count]
        self.nearestHam = nearestHam
        self.nearestHamDist = nearestHamDist
        self.totalcount = count

    def add(self, barcode, count):
        self.barcodes.append(barcode)
        self.counts.append(count)
        self.totalcount += count

    def __str__(self):
        return int2string(self.barcodes[0])+"\t"+str(self.counts[0])+"\tTotal\t" + str(self.totalcount) + \
               "\tFraction_of_all_barcodes\t" + str(self.totalcount/totalcount) + "\t"+str(self.nearestHamDist) + \
               "_From_" + self.nearestHam + "\t" + "\t".join(int2string(b) + '\t'+str(c) for b, c
                                                             in zip(self.barcodes[1:], self.counts[1:]))

    def hamdist(self, otherint, numbits=30):
        xor = self.barcodes[0] ^ otherint
        return sum([(xor >> i) & 1 for i in xrange(numbits)])

    def __lt__(self, other):
        return self.totalcount > other.totalcount

    def __eq__(self, other):
        return self.barcodes[0]

    def __hash__(self):
        return self.barcodes[0]

powersoftwo = [2**i for i in range(31)]


def similarBarcodes(barcode, numbits=30, hamDist=2, allowIndels=False):
    if allowIndels:
        for bc in mutate_with_indels(barcode, 0, numbits, hamDist):
            yield bc
    else:
        for bc in mutate(barcode, 0, numbits, hamDist):
            yield bc


def mutate(barcode, bitstart, bitend, hamDist):
    if hamDist == 1:
        for i in xrange(bitstart, bitend):
            yield barcode ^ powersoftwo[i]
    else:
        for i in xrange(bitstart, bitend):
            temp = barcode ^ powersoftwo[i]
            for result in mutate(temp, i+1, bitend, hamDist - 1):
                yield result


def mutate_with_indels(barcode, bitstart, bitend, hamDist):
    if hamDist == 1:
        for i in xrange(bitstart, bitend):
            yield barcode ^ powersoftwo[i]
        tested = set()
        for bmi in xrange(bitend-2, 1, -1):
            temp = (barcode >> bmi) % 4
            temp2 = (barcode >> bmi-2) % 4
            if temp == temp2:
                #2bp repeat
                bc = (barcode >> bmi << bmi) + ((barcode % (1<<bmi-2)) << 2)
                if not bc in tested:
                    tested.add(bc)
                    yield bc
                #2bp insertion
                bc = (barcode >> bmi << bmi) + ((barcode % (1<<bmi)) >> 2) + (temp << bmi-2)
                if not bc in tested:
                    tested.add(bc)
                    yield bc
                #4bp insertion
                if bmi>=4:
                    bc = (barcode >> bmi << bmi) + ((barcode % (1<<bmi)) >> 4) + (temp << bmi-2) + (temp << bmi-4)
                    if not bc in tested:
                        tested.add(bc)
                        yield bc
        for bmi in xrange(bitend-4, 3, -1):
            temp = (barcode >> bmi) % 16
            temp2 = (barcode >> bmi-4) % 16
            if temp==temp2:
                #4bp repeat
                bc = (barcode >> bmi << bmi) + ((barcode % (1<<bmi)) >> 4) + (temp << bmi-4)
                if not bc in tested:
                    tested.add(bc)
                    yield bc
                #8bp insertion
                if bmi >= 8:
                    bc = (barcode >> bmi << bmi) + ((barcode % (1<<bmi)) >> 8) + (temp << bmi-4) + (temp << bmi-8)
                    if not bc in tested:
                        tested.add(bc)
                        yield bc
    else:
        for i in xrange(bitstart, bitend):
            temp = barcode ^ powersoftwo[i]
            for result in mutate_with_indels(temp, i+1, bitend, hamDist - 1):
                yield result


def string2int(line):
    cur = 0
    toreturn = 0
    for char in line:
        if not char in charvals[cur]:
            return None
        toreturn = (toreturn << 1) + charvals[cur][char]
        if cur == 0:
            cur = 1
        else:
            cur = 0
    return toreturn


def int2string(input, mylen=30):
    cur = 0
    if mylen % 2 == 0:
        cur += 1
    toreturn = ''
    for i in xrange(mylen):
        toreturn += intvals[cur][input & 1]
        input = input >> 1
        if cur == 0:
            cur = 1
        else:
            cur = 0
    return toreturn[::-1]


def barcode_count(input_file_list, output_file, barcode_index, known_contaminant="", phred_format=33,
                  library_mode=False, barcode_length=30, bad_seq_reads=True):

    global totalcount
    has_error = False

    try:

        linenum = 0
        known_contaminant_set = set()
        counts = {}
        totalcount = 0
        totaltries = 0
        numbadend = 0
        numbadseq = 0
        numbadqualbase = 0
        numbadqualavg = 0
        numexceptions = 0
        known_contaminant_count = 0

        base_cutoff = 10
        average_cutoff = 30
        barcode_index_length = len(barcode_index)

        myval=""

        if known_contaminant != "":
            fh = open(known_contaminant)
            for line in fh:
                known_contaminant_set.add(line.strip())
            fh.close()


        fh = zippedFile(input_file_list)#open(sys.argv[1])

        for line in fh.generator:
            linenum += 1

            if linenum >= 4:
                linenum = 0

            if linenum != 2 and linenum != 0:
                continue

            if linenum == 2:
                strippedline = line.strip()
                consider_qual = True
                totaltries += 1

                if known_contaminant != "" and (strippedline[:barcode_length] in known_contaminant_set):
                    known_contaminant_count += 1
                    consider_qual = False
                    continue

                if barcode_index_length > 0 and len(strippedline) >= barcode_length+barcode_index_length and strippedline[barcode_length:barcode_length +
                        barcode_index_length] != barcode_index:
                    numbadend += 1
                    consider_qual = False
                    continue

                myval = string2int(strippedline[:barcode_length])
                if bad_seq_reads and myval is None:
                        numbadseq += 1
                        consider_qual = False
                        continue

            if linenum == 0 and consider_qual:

                qualline = line.strip()
                bad_base = False
                # Initialize
                total_score = 0
                count_chars = 0

                # Get the average quality of all the nucleotides
                for char_index in range(0, barcode_length):
                    phred_score = ord(qualline[char_index]) - phred_format
                    if phred_score < base_cutoff:
                        bad_base = True
                        break
                    total_score += phred_score
                    count_chars += 1

                if bad_base:
                    numbadqualbase += 1
                    continue

                try:

                    if (total_score * 1.0)/barcode_length < average_cutoff:
                        numbadqualavg += 1
                        continue
                except ZeroDivisionError:
                    numexceptions += 1
                    continue

                totalcount += 1
                if myval in counts:
                    counts[myval] += 1
                else:
                    counts[myval] = 1

        fh.close()

        totalcount = float(totalcount)
        barcodes = sorted(counts.items(), key=lambda x: x[1], reverse=True)
        toprint = []
        addedBarcodes = {}
        similarBarcodeDict1off = {}
        similarBarcodeDict2off = {}
        lessThanTwo = 0
        lessThanTwoList = []

        if library_mode:
            for item in barcodes:
                toprint.append(Barcode(item[0], item[1], "None.0", 0))
        else:
            for i in xrange(len(barcodes)):

                if barcodes[i][1] < 2:
                    lessThanTwo += 1
                    lessThanTwoList.append(barcodes[i][0])
                    continue
                if barcodes[i][0] in similarBarcodeDict1off:
                    similarBarcodeDict1off[barcodes[i][0]].add(barcodes[i][0], barcodes[i][1])
                    continue
                if barcodes[i][0] in similarBarcodeDict2off:
                    similarBarcodeDict2off[barcodes[i][0]].add(barcodes[i][0], barcodes[i][1])
                    continue

                nearestHam = ""
                nearestHamDist = 30

                for barcode in toprint:
                    if barcode.counts[0] < 50*barcodes[i][1]:
                        break
                    myHamDist = barcode.hamdist(barcodes[i][0])
                    if myHamDist < nearestHamDist:
                        nearestHam = barcode
                        nearestHamDist = myHamDist
                if nearestHamDist <= 3:
                    nearestHam.add(barcodes[i][0], barcodes[i][1])
                    continue

                if nearestHam == "":
                    mybarcode = Barcode(barcodes[i][0], barcodes[i][1], "None.0", nearestHamDist)
                else:
                    mybarcode = Barcode(barcodes[i][0], barcodes[i][1], int2string(nearestHam.barcodes[0])+"."+str(nearestHam.totalcount), nearestHamDist)
                toprint.append(mybarcode)
                addedBarcodes[barcodes[i][0]] = mybarcode
                if barcodes[i][1] >= 8: #need to see main barcode at least 8 times to expect to see a 1-off at least twice
                    for bc in similarBarcodes(barcodes[i][0], hamDist=1):
                        if not bc in similarBarcodeDict1off:
                            similarBarcodeDict1off[bc] = mybarcode
                if barcodes[i][1] >= 40: #need to see main barcode at least 40 times to expect to see a 2-off at least twice
                    for bc in similarBarcodes(barcodes[i][0], hamDist=2):
                        if not bc in similarBarcodeDict2off:
                            similarBarcodeDict2off[bc] = mybarcode

        with open(output_file, 'w') as f:

            f.write("#Total_sequences:\t" + str(int(totaltries)))
            f.write("\n")
            f.write("#Total_good_reads:\t" + str(int(totalcount) - int(lessThanTwo)))
            f.write("\n")
            f.write("#Total_bad_end_reads:\t" + str(int(numbadend)))
            f.write("\n")
            f.write("#Total_bad_seq_reads:\t" + str(int(numbadseq)))
            f.write("\n")
            f.write("#Total_bad_base_quals:\t" + str(int(numbadqualbase)))
            f.write("\n")
            f.write("#Total_bad_avg_qual:\t" + str(int(numbadqualavg)))
            f.write("\n")
            f.write("#Total_exceptions_qual:\t" + str(int(numexceptions)))
            f.write("\n")
            f.write("#Total_less_than_two:\t" + str(int(lessThanTwo)))
            f.write("\n")
            f.write("#Total_known_contaminants:\t" + str(int(known_contaminant_count)))
            f.write("\n")

            for barcode in toprint:
                f.write(str(barcode))
                f.write("\n")

            f.write("## less than two")
            f.write("\n")

            for barcode in lessThanTwoList:
                f.write("# "+int2string(barcode))
                f.write("\n")
    except:
        has_error = True

    return has_error
