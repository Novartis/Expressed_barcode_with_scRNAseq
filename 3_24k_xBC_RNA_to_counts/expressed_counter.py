import sys
import getopt
import datetime
import os
import multiprocessing
from operator import itemgetter
import nwalign as nw
from Bio import pairwise2
import subprocess
import warnings
import gzip

warnings.filterwarnings("ignore")

__author__ = 'krishvi7'


# To run express counter
# python express_only.py -i input.fastq.gz -t tenx.fastq.gz -s sample_name -o op

# Takes in 1 and 2 read and name of the exon and returns the barcode tag
def get_pos_base_counts(output, input_files, sample_name, pipeline, tenx_dict, possible_tenx_reference_set):
    stats = {"total_reads": 0,
             "good.exact_match": 0,
             "good.hd_lt_" + str(Default.pipeline_config[pipeline]["max_merging_hd"]): 0,
             "bad.missing.start": 0,
             "bad.missing.end": 0,
             "bad.missing.last3.ecor1": 0,
             "bad.missing."+Default.pipeline_config[pipeline]["end"][0]+".end_or_missing.AC.at.4": 0,
             "bad.hd_indel_gt_"+str(Default.pipeline_config[pipeline]["max_merging_hd"]): 0,
             "bad." + str(Default.pipeline_config[pipeline]["start_seed"]) + "_start_seed_match": 0,
             "bad.hd_tie": 0}

    tag_dict = dict()
    read_tracker = list()
    return_message = {"exp": False, "mess": ""}
    reference_tag_name_dict, reference_seed_tag_dict = Default.get_reference_tag_names(pipeline)
    reference_tag_set = reference_tag_name_dict.keys()
    tenx_exp_dict = dict()
    tenx_cr_list = list()

    try:
        for tag in reference_tag_set:
            tag_dict[tag] = 0

        # Open the fastq.gz file directly
        counter = 0
        with gzip.open(input_files[0]) as fh:
            for line in fh:
                counter += 1
                # Ignore blank lines
                if not line:
                    continue

                if counter % 4 == 1:
                    header = line.strip().split(" ")[0]
                elif counter % 4 == 2:
                    ac_missing = False
                    end_missing = False
                    read = line.strip()
                    # Add this to total reads
                    stats["total_reads"] += 1

                    if stats["total_reads"] % 25000 == 0:
                        print >> sys.stderr, str(stats["total_reads"]*1.0/1000)+"K reads processed"

                    if "ecor1" in Default.pipeline_config[pipeline]:
                        ecor1_site = Default.pipeline_config[pipeline]["ecor1"][0]
                        ecor1_start = Default.pipeline_config[pipeline]["ecor1"][1]
                        ecor1_len = len(ecor1_site)

                        # Check for ecor1 site
                        if ecor1_site != read[ecor1_start:ecor1_start + ecor1_len]:
                            stats["bad.missing.last3.ecor1"] += 1
                            continue

                    if "AC" in Default.pipeline_config[pipeline]:
                        ac_site = Default.pipeline_config[pipeline]["AC"][0]
                        ac_start = Default.pipeline_config[pipeline]["AC"][1]
                        ac_len = len(ac_site)

                        ac_seq = read[ecor1_start+ecor1_len+ac_start:ecor1_start+ecor1_len+ac_start+ac_len]

                        # Check if the read ends with expected sequence we have AC present at
                        if Default.pipeline_config[pipeline]["AC"][0] != ac_seq:
                            ac_missing = True

                    if "start" in Default.pipeline_config[pipeline]:
                        start_site = Default.pipeline_config[pipeline]["start"][0]
                        start_site_len = len(start_site)
                        start_start = Default.pipeline_config[pipeline]["start"][1]

                        start_seq = read[start_start:start_start+start_site_len]

                        # Check if the read ends with expected sequence we have AC present at
                        if start_site != start_seq:
                            cat = "bad.missing.start"
                            if Default.v:
                                read_tracker.append("\t".join(map(str, [header, cat, "NA", "NA"])))
                            stats[cat] += 1
                            continue

                    end_sequence = read[Default.pipeline_config[pipeline]["end"][1]:
                                        Default.pipeline_config[pipeline]["end"][1]+len(Default.pipeline_config[pipeline]["end"][0])]

                    if Default.end_check and Default.pipeline_config[pipeline]["end"][0] != end_sequence:
                        end_missing = True

                    if ac_missing or end_missing:
                        cat = "bad.missing."+Default.pipeline_config[pipeline]["end"][0]+".end_or_missing.AC.at.4"
                        if "AC" not in Default.pipeline_config[pipeline]:
                            cat = "bad.missing.end"
                        if Default.v:
                            read_tracker.append("\t".join(map(str, [header, cat, "NA", "NA"])))
                        stats[cat] += 1
                        continue

                    if "ecor1" in Default.pipeline_config[pipeline]:
                        tag_start = ecor1_start + len(Default.pipeline_config[pipeline]["ecor1"][0])
                    else:
                        tag_start = start_start + start_site_len

                    tag = read[tag_start:tag_start+Default.pipeline_config[pipeline]["barcode_length"]]

                    # Check for exact match
                    if tag in tag_dict:
                        tag_dict[tag] += 1
                        stats["good.exact_match"] += 1
                        tenx_bc = tenx_dict[header]
                        if Default.v:
                            read_tracker.append("\t".join(map(str, [header, "good_exact", "NA", "NA"])))
                        if (tenx_bc, tag) not in tenx_exp_dict:
                            tenx_exp_dict[(tenx_bc, tag)] = 0
                        tenx_exp_dict[(tenx_bc, tag)] += 1

                    elif Default.pipeline_config[pipeline]["merge_tags"]:
                        min_ham = 300
                        hd_dict = dict()

                        # Check for exact seed match
                        if tag[0:Default.pipeline_config[pipeline]["start_seed"]] in reference_seed_tag_dict:
                            subset_tag_dict = reference_seed_tag_dict[
                                tag[0:Default.pipeline_config[pipeline]["start_seed"]]]

                            # Compare tags with min_ham distance
                            for ref_tag in subset_tag_dict:
                                ham_distance = hamming_distance(tag, ref_tag)

                                # Changes the reference only if a lower hamming distance is found
                                if min_ham >= ham_distance:
                                    preferred = ref_tag
                                    if ham_distance not in hd_dict:
                                        hd_dict[ham_distance] = list()
                                    hd_dict[ham_distance].append(preferred)
                                    min_ham = ham_distance

                            # If there is a good match with a unique cluster member, assign the count
                            if len(hd_dict[min_ham]) == 1 and min_ham <= Default.pipeline_config[pipeline][
                                "max_merging_hd"]:
                                stats["good.hd_lt_" + str(Default.pipeline_config[pipeline]["max_merging_hd"])] += 1
                                tenx_bc = tenx_dict[header]
                                if (tenx_bc, tag) not in tenx_exp_dict:
                                    tenx_exp_dict[(tenx_bc, tag)] = 0
                                tenx_exp_dict[(tenx_bc, tag)] += 1

                                tag_dict[hd_dict[min_ham][0]] += 1
                                if Default.v:
                                    read_tracker.append("\t".join(map(str, [header, "bad.hd_tie",
                                                                        ";".join(hd_dict[min_ham]),
                                                                        "min_ham:" + str(min_ham)])))
                            # If there is a tie, discard (indel correction won't help tag assignment)
                            elif len(hd_dict[min_ham]) > 1:
                                stats["bad.hd_tie"] += 1
                                if Default.v:
                                    read_tracker.append("\t".join(map(str, [header, "bad.hd_tie",
                                                                        ";".join(hd_dict[min_ham]),
                                                                        "min_ham:" + str(min_ham)])))
                        else:
                            stats["bad." + str(Default.pipeline_config[pipeline]["start_seed"]) + "_start_seed_match"] += 1
                            if Default.v:
                                read_tracker.append(
                                    "\t".join(
                                        map(str, [header, "bad." + str(Default.pipeline_config[pipeline]["start_seed"])
                                                  + "_start_seed_match", "NA", "NA"])))
        fh.close()

        sorted_tags = sorted(tag_dict, key=tag_dict.get, reverse=True)
        sum_total = sum(tag_dict.values())

        if sum_total == 0:
            with open(output + "/" + sample_name + "_stats.txt", 'w') as fh:
                for key in sorted(stats.keys()):
                    print >> fh, "#" + key + ":" + str(stats[key])
                for item in read_tracker:
                    print >> fh, item
            fh.close()
            print >> sys.stderr, "No reads found"
            sys.exit(1)

        sorted_10x_tags = sorted(tenx_exp_dict, key=tenx_exp_dict.get, reverse=True)
        tenx_sum_total = sum(tenx_exp_dict.values())

        cw = 0.0
        with open(output + "/" + sample_name + "_counts.txt", 'w') as fh:
            print >> fh, "name\ttag\tcount\tfraction\tcw"
            for tag in sorted_tags:
                frac = tag_dict[tag] * 1.0 / sum_total
                cw += frac
                print >> fh, "\t".join(map(str, [reference_tag_name_dict[tag], tag, tag_dict[tag],
                                                 '%.3f' % frac, '%.3f' % cw]))
        fh.close()

        with open(output + "/" + sample_name + "_stats.txt", 'w') as fh:
            for key in sorted(stats.keys()):
                print >> fh, "#" + key + ":" + str(stats[key])
            for item in read_tracker:
                print >> fh, item
        fh.close()

        tenx_list = list()

        cw = 0.0
        with open(output + "/" + sample_name + "_tenx_associations.txt", 'w') as fh:
            print >> fh, "sample\t10x_barcode\texp_barcode\treference_10x\tcount\tfraction\tcw\tsc_10x\tsc_10x_fraction\tcr_found"
            for ten_x, exp in sorted_10x_tags:
                frac = tenx_exp_dict[(ten_x, exp)] * 1.0 / tenx_sum_total
                cw += frac
                tenx_list.append(ten_x)
                print >> fh, sample_name+"\t"+ten_x+"\t"+exp+"\t"+str(ten_x in possible_tenx_reference_set)+"\t"+str(
                    tenx_exp_dict[(ten_x, exp)])+"\t"+'%.8f' % frac+"\t"+'%.3f' % cw +"\t" + "\t" + ("True" if ten_x in tenx_cr_list else "False")
            for tenx_cr_bc in tenx_cr_list:
                if tenx_cr_bc not in tenx_list:
                    print >> fh, sample_name+"\t"+tenx_cr_bc + "\tNA\t" + str(ten_x in possible_tenx_reference_set) + \
                    "\t0\t0.0\t0.0\t" + "\tTrue"

        fh.close()

    except ValueError as exp:
        return_message["exp"] = True
        return_message["mess"] = exp

    return return_message


# Get basic hamming distance between 2 sequences
def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


# Perform pairwise alignment between 2 sequences and run in global or local mode
#
# Algorithm:
# ----------
# 1. Global alignment is optimized for hamming distance calculation with indels.
#    Best for alignment with reference sequences.
# 2. Local alignment is gap opening penalized.
#    Best for creating a single tag from paired end reads. Use score carefully
#
def hamming_distance_with_indels(s1, s2, algorithm="global"):
    # This gives the best gapped alignment between 2 sequences. Match is scored 1 and
    # mis-match is scored -1. So the score should give the gapped hamming distance
    if algorithm == "global":
        # Perform global alignment
        return nw.global_align(s1, s2, gap_open=-10, gap_extend=-5, match=1)
    else:
        # Perform gapped alignment. parameters (seq1, seq2, match_score, mismatch_score,
        # gap_opening_score, gap_extension_score)
        return pairwise2.align.localms(s1, s2, 2, 0, -4, -2)[0]


class Common:
    def __init__(self):
        self.dataType_dict = {"input": str,
                              "output": str,
                              "sample_name": str,
                              "pipeline": str,
                              "tenx_fastq": str,
                              "tenx_cr": str}
        self.range_dict = {"pipeline": ["expressed_barcode", "cropseq"]}

    @staticmethod
    def reverse_compliment(sequence):
        return "".join(Default.compliment.get(base, base) for base in reversed(sequence))

    @staticmethod
    def just_compliment(sequence):
        return "".join(Default.compliment.get(base, base) for base in sequence)

    @staticmethod
    def print_error(message):
        print >> sys.stderr, "ERROR: " + message
        sys.exit(1)

    @staticmethod
    def get_timestamp():
        return datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')

    @staticmethod
    def super_usage():
        print("\nExpress counter pipeline- Express counter pipeline is used to map expressed barcodes to references")


# Contains all the default parameters
class Default:

    o = "output_" + Common().get_timestamp()
    v = False
    end_check = False
    paired_end = False
    cr_cellular_barcode = True

    # relative_length # Length of the expected reference in the sequence read
    # max_merging_hd  # Max hamming distance between 2 seqs to be merged after indel correction
    # merge_tags = True # Enable tag merging
    # start_seed  # Search for exact seed match

    flash = "/usr/prog/onc/seqtools/flash/1.2.11/flash"
    pipeline_config = {"expressed_barcode": {"reference": "/da/onc/krishvi7/bitbucket/oncp-expressed/references/expressed_bc_24K.txt",
                                             "barcode_length": 22, "max_merging_hd": 3, "merge_tags": True,
                                             "ecor1": ["TTC", 70], "end": ["GTT", -30], "AC": ["AC", 4], "start_seed": 3,
                                             "10x_start": 0, "10x_length": 16,
                                             "10x_references": "/da/onc/krishvi7/bitbucket/oncp-expressed/references/tenx_barcodes.txt"},
                       "cropseq": {
                           "reference": "/da/onc/krishvi7/bitbucket/oncp-expressed/references/cropseq.txt",
                           "barcode_length": 20, "max_merging_hd": 3, "merge_tags": True,
                           "start": ["CACCG", 15], "end": ["GTTT", -48], "start_seed": 3,
                           "10x_start": 0, "10x_length": 16,
                           "10x_references": "/da/onc/krishvi7/bitbucket/oncp-expressed/references/tenx_barcodes_v3.txt.gz"}
                       }

    @staticmethod
    def get_reference_file(pipeline):
        # dataframe with pos:base as index, samples as col names and fraction
        # of each samples for given pos:base as data
        return Default.pipeline_config[pipeline]["reference"]

    @staticmethod
    def get_reference_tag_names(pipeline):
        reference_tag_name = dict()
        reference_tag_seed_dict = dict()
        # Load expected to a dictionary
        with open(Default.get_reference_file(pipeline)) as fh:
            fh.readline()
            for line in fh:
                ls = line.strip().split("\t")

                name = ls[0]
                seq = ls[1]

                seq_seed = seq[0:Default.pipeline_config[pipeline]["start_seed"]]

                if seq_seed not in reference_tag_seed_dict:
                    reference_tag_seed_dict[seq_seed] = list()

                reference_tag_seed_dict[seq_seed].append(seq)
                reference_tag_name[seq] = name
        fh.close()

        return reference_tag_name, reference_tag_seed_dict

    # Don't change this
    compliment = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}





# Contains code to validate a list of input files
class ExpressCounter:
    def __init__(self, input_file="", sample_name="", output=Default.o, pipeline="", tenx_fastq="", tenx_cr=""):
        self.input = input_file
        self.sample_name = sample_name
        self.output = output
        self.pipeline = pipeline
        self.tenx_fastq = tenx_fastq
        self.tenx_cr = tenx_cr
        self.file_extension_list = [".fastq.gz"]
        self.args_with_params = ["output", "sample_name", "input", "tenx_fastq", "pipeline", "tenx_cr"]
        self.args_without_params = []
        self.args_dict = {"output": "o", "input": "i", "tenx_fastq": "t", "sample_name": "s", "pipeline": "p", "tenx_cr": "c"}
        self.reverse_args_dict = {value: key for key, value in self.args_dict.iteritems()}
        self.args_description = {"i": "Comma separated list of input fastq files Supports " +
                                      ",".join(self.file_extension_list) + ".",
                                 "s": "Sample name",
                                 "o": "Output directory to store counts",
                                 "t": "10X fastq file",
                                 "c": "10X barcodes detected from single cell RNAseq 10x run",
                                 "p": "Pipeline to run"
                                 }
        self.options = ["s", "o", "p", "c", "t"]

    def usage(self):
        Common().super_usage()
        print("\nUsage:")
        print("-------")
        print("\tpython express_counter.py [options]* -i||--input <Single end fq file or Paired end fq files separated by comma> -p|--pipeline <Pipeline to be run>")
        print("\noptions:")
        print("--------")
        for option in self.options:
            print(" \t-" + option + "|--" + self.reverse_args_dict[option] + "\t: " + self.args_description[option] +
                  (". Default- " + str(getattr(Default, option) if option in Default.__dict__.keys() else "")))
        print("")
        sys.exit(1)

    # Check if file exists and is of compatible extension
    def check_input(self):

        # Check if input files are provided
        if self.pipeline == "":
            Common().print_error("Pipeline for reference based counting is not supplied")

        # Check if input files are provided
        if self.input == "":
            Common().print_error("Input file(s) for reference based counting is not supplied")

        # Check if input files are provided
        if self.tenx_fastq == "":
                Common().print_error("10X fastq file for reference based counting is not supplied")

        input_file_split = self.input.split(",")

        if len(input_file_split) == 1:
            print >> sys.stderr, "Running in single end mode"

        for input_file in input_file_split:
            if not os.path.isfile(input_file):
                Common().print_error("File " + input_file + " does not exist")

            if not any(input_file.endswith(x) for x in self.file_extension_list):
                Common().print_error(
                    "Invalid file: " + input_file + ", Must be a *" + ' ,*'.join(self.file_extension_list))

        if not os.path.isfile(self.tenx_fastq):
            Common().print_error("File " + self.tenx_fastq + " does not exist")

    # Check output
    def check_output(self):
        if not self.output:
            print("No output")
            sys.exit(1)

        if not os.path.isdir(self.output):
            os.makedirs(self.output)

    # Appends tag files with sample information
    def process(self):
        input_files = self.input.split(",")
        possible_tenx_reference_set = set()
        tenx_dict = dict()
        header = None

        if self.sample_name == "":
            self.sample_name = os.path.basename(input_files[0].split(".fastq")[0])

        smfuid = os.path.dirname(input_files[0].split(".fastq")[0])

        # Disabled for single thread running
        # pool = multiprocessing.Pool(2)

        counter = 0
        with gzip.open(self.tenx_fastq) as fh:
            for line in fh:
                counter += 1
                if not line:
                    continue
                if counter % 4 == 1:
                    header = line.strip().split(" ")[0]
                elif counter % 4 == 2:
                    tenx_dict[header] = line.strip()[Default.pipeline_config[self.pipeline]["10x_start"]: \
                                                     Default.pipeline_config[self.pipeline]["10x_start"] + \
                                                     Default.pipeline_config[self.pipeline]["10x_length"]]
        fh.close()

        if "gz" in Default.pipeline_config[self.pipeline]["10x_references"]:
            with gzip.open(Default.pipeline_config[self.pipeline]["10x_references"]) as fh:
                for line in fh:
                    possible_tenx_reference_set.add(line.strip())
            fh.close()

        else:
            with open(Default.pipeline_config[self.pipeline]["10x_references"]) as fh:
                for line in fh:
                    possible_tenx_reference_set.add(line.strip())
            fh.close()

        # Create a results list and send them for execution
        results = list()
        results.append(get_pos_base_counts(self.output, input_files, self.sample_name, self.pipeline, tenx_dict,
                                           possible_tenx_reference_set))
        # results.append(pool.apply_async(get_pos_base_counts, (self.output, input_files, self.sample_name)))

        # Get the results for both exons
        for result in results:
            # return_message = result.get()
            return_message = result
            if return_message["exp"]:
                print >> sys.stderr, "***********************************"
                print >> sys.stderr, "Exception occurred while processing"
                print >> sys.stderr, return_message["mess"]
                print >> sys.stderr, "***********************************"
                sys.exit(1)


def main(argv):
    express_counter = ExpressCounter()

    try:
        opts, args = getopt.getopt(argv, ":".join(express_counter.args_dict.values()) + ":", map(lambda x: x + "=",
                                                                                         express_counter.args_dict.keys()))
    except getopt.GetoptError as exp:
        Common().print_error("Invalid input argument: " + exp.msg)

    if len(opts) == 0:
        express_counter.usage()

    for opt, arg in opts:
        if opt in map(lambda x: "-" + x, express_counter.args_dict.values()) or opt in map(lambda x: "--" + x, express_counter.
                args_dict.keys()):
            if opt in map(lambda x: "-" + x, express_counter.args_dict.values()):
                clte_member = express_counter.args_dict.keys()[express_counter.args_dict.values().index(opt[1:])]
            else:
                clte_member = opt[2:]

            if clte_member in express_counter.args_without_params:
                setattr(express_counter, clte_member, not getattr(Default(), express_counter.args_dict[clte_member]))
            else:
                tryout = None
                try:
                    tryout = Common().dataType_dict[clte_member](arg)
                except ValueError:
                    Common().print_error("Invalid argument " + arg + " for parameter " + opt)
                if clte_member in Common().range_dict.keys():
                    if tryout not in Common().range_dict[clte_member]:
                        raise Common().print_error("Argument for " + opt + " must be " +
                                                   ",".join(map(str, Common().range_dict[clte_member])))
                setattr(express_counter, clte_member, tryout)
        else:
            Common().print_error("Invalid input argument " + opt)

    print >> sys.stderr, "Evaluating arguments."
    express_counter.check_input()
    express_counter.check_output()

    print >> sys.stderr, "Processing "+express_counter.sample_name
    express_counter.process()
    print >> sys.stderr, express_counter.sample_name + " processing complete"


if __name__ == "__main__":
    main(sys.argv[1:])
