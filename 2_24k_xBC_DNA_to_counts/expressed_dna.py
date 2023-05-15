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


# To run leuk
# python expressed_dns.py -i input.fastq.gz -s sample_name -o op


# Takes in 1 and 2 read and name of the exon and returns the barcode tag
def get_pos_base_counts(output, input_files, sample_name, pipeline):
    stats = {"total_reads": 0,
             "good.exact_match": 0,
             "good.hd_lt_" + str(Default.pipeline_config[pipeline]["max_merging_hd"]): 0,
             "good.hd_indel_lt_" + str(Default.pipeline_config[pipeline]["max_merging_hd"]): 0,
             "bad." + str(Default.pipeline_config[pipeline]["start_seed"]) + "_start_seed_match": 0,
             "bad.hd_tie": 0,
             "bad.hd_indel_gt_"+str(Default.pipeline_config[pipeline]["max_merging_hd"]): 0}

    tag_dict = dict()
    read_tracker = list()
    return_message = {"exp": False, "mess": ""}
    reference_tag_name_dict, reference_seed_tag_dict = Default.get_reference_tag_names(pipeline)
    reference_tag_set = reference_tag_name_dict.keys()

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
                    read = line.strip()
                    tag = read[0:Default.pipeline_config[pipeline]["relative_length"]]

                    # Add this to total reads
                    stats["total_reads"] += 1

                    if stats["total_reads"] % 25000 == 0:
                        print >> sys.stderr, str(stats["total_reads"]*1.0/1000)+"K reads processed"

                    # Check for exact match
                    if tag in tag_dict:
                        tag_dict[tag] += 1
                        # read_tracker.append("\t".join(map(str, [header, "good.exact_match", tag, "NA"])))
                        stats["good.exact_match"] += 1
                        continue
                    elif Default.pipeline_config[pipeline]["merge_tags"]:
                        min_ham = 300
                        hd_dict = dict()

                        # Check for exact seed match
                        if tag[0:Default.pipeline_config[pipeline]["start_seed"]] in reference_seed_tag_dict:
                            subset_tag_dict = reference_seed_tag_dict[tag[0:Default.pipeline_config[pipeline]["start_seed"] ]]

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
                            if len(hd_dict[min_ham]) == 1 and min_ham <= Default.pipeline_config[pipeline]["max_merging_hd"]:
                                stats["good.hd_lt_"+str(Default.pipeline_config[pipeline]["max_merging_hd"])] += 1
                                tag_dict[hd_dict[min_ham][0]] += 1
                                read_tracker.append("\t".join(map(str, [header, "bad.hd_tie",
                                                                        ";".join(hd_dict[min_ham]),
                                                                        "min_ham:" + str(min_ham)])))
                            # If there is a tie, discard (indel correction won't help tag assignment)
                            elif len(hd_dict[min_ham]) > 1:
                                stats["bad.hd_tie"] += 1
                                read_tracker.append("\t".join(map(str, [header, "bad.hd_tie",
                                                                        ";".join(hd_dict[min_ham]),
                                                                        "min_ham:"+str(min_ham)])))
                            # If the ham distance is more than allowed, pick the potential reference and
                            # perform indel correction
                            else:
                                pot_reference = hd_dict[min_ham][0]
                                hd_tuple = hamming_distance_with_indels(tag, pot_reference)
                                ham_distance = hamming_distance(hd_tuple[0], hd_tuple[1])

                                if ham_distance > Default.pipeline_config[pipeline]["max_merging_hd"]:
                                    stats["bad.hd_indel_gt_"+str(Default.pipeline_config[pipeline]["max_merging_hd"])] += 1
                                    read_tracker.append("\t".join(map(str, [header, read,
                                                                            "bad.hd_indel_gt_"
                                                                            +str(Default.pipeline_config[pipeline]["max_merging_hd"]),
                                                                            pot_reference,
                                                                            "min_ham:" + str(min_ham)])))
                                else:
                                    # Merge the tag with the winner and note the hamming distance as well
                                    tag_dict[hd_dict[min_ham][0]] += 1
                                    stats["good.hd_indel_lt_" + str(Default.pipeline_config[pipeline]["max_merging_hd"])] += 1
                                    if Default.pipeline_config[pipeline]["indel_correction"]:
                                        read_tracker.append("\t".join(map(str, [header, "good.hd_indel_lt_" +
                                                                                str(Default.pipeline_config[pipeline]["max_merging_hd"]),
                                                                                pot_reference,
                                                                                ";".join(["min_ham:" + str(min_ham),
                                                                                "corrected_ham:"+str(ham_distance)])])))
                        else:
                            stats["bad."+str(Default.pipeline_config[pipeline]["start_seed"])+"_start_seed_match"] += 1
                            if Default.pipeline_config[pipeline]["indel_correction"]:
                                read_tracker.append(
                                    "\t".join(map(str, [header, "bad."+str(Default.pipeline_config[pipeline]["start_seed"])
                                                        +"_start_seed_match", "NA", "NA"])))

        fh.close()

        sorted_tags = sorted(tag_dict, key=tag_dict.get, reverse=True)
        sum_total = sum(tag_dict.values())

        with open(output + "/" + sample_name + "_counts.txt", 'w') as fh:
            print >> fh, "name\ttag\tcount\tfraction"
            for tag in sorted_tags:
                print >> fh, "\t".join(map(str, [reference_tag_name_dict[tag], tag, tag_dict[tag],
                                                 tag_dict[tag] * 1.0 / sum_total]))
        fh.close()

        with open(output + "/" + sample_name + "_stats.txt", 'w') as fh:
            for key in sorted(stats.keys()):
                print >> fh, "#" + key + ":" + str(stats[key])
            if Default.v:
                for item in read_tracker:
                    print >> fh, item
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
                              "pipeline": str}
        self.range_dict = {"pipeline": ["expressed_barcode"]}

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
        print("\nExpressedDNA pipeline- Expressed DNA pipeline is used to map tags to variants of a reference exon")


# Contains all the default parameters
class Default:

    o = "output_" + Common().get_timestamp()
    v = False
    paired_end = False

    # relative_length # Length of the expected reference in the sequence read
    # max_merging_hd  # Max hamming distance between 2 seqs to be merged after indel correction
    # merge_tags = True # Enable tag merging
    # indel_correction # Enable indel based merging
    # start_seed  # Search for exact seed match

    flash = "/usr/prog/onc/seqtools/flash/1.2.11/flash"
    pipeline_config = {"expressed_barcode": {"reference": os.environ["EXPRESSED_DNA_HOME"]+"/references/expressed_bc_24K.txt",
                                             "relative_length": 22, "max_merging_hd": 3, "merge_tags": True,
                                             "indel_correction": True, "start_seed": 3}}

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

                if len(ls) == 4:
                    name_orig = ls[0]
                    name_opt = ls[0]+"_opt"
                    seq_orig = ls[1]
                    seq_opt = ls[3]

                    seq_seed_opt = seq_opt[0:Default.pipeline_config[pipeline]["start_seed"]]
                    seq_seed_orig = seq_orig[0:Default.pipeline_config[pipeline]["start_seed"]]

                    if seq_seed_opt not in reference_tag_seed_dict:
                        reference_tag_seed_dict[seq_seed_opt] = list()
                    if seq_seed_orig not in reference_tag_seed_dict:
                        reference_tag_seed_dict[seq_seed_orig] = list()

                    if seq_opt == seq_orig:
                        reference_tag_name[seq_orig] = name_opt+";"+name_orig
                        reference_tag_seed_dict[seq_seed_orig].append(seq_orig)
                    else:
                        reference_tag_name[seq_orig] = name_orig
                        reference_tag_name[seq_opt] = name_opt
                        reference_tag_seed_dict[seq_seed_orig].append(seq_orig)
                        reference_tag_seed_dict[seq_seed_opt].append(seq_opt)
                elif len(ls) == 2:
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
class ExpressedDNA:
    def __init__(self, input_file="", sample_name="", output=Default.o, pipeline=""):
        self.input = input_file
        self.sample_name = sample_name
        self.output = output
        self.pipeline = pipeline
        self.file_extension_list = [".fastq.gz"]
        self.args_with_params = ["output", "sample_name", "input", "pipeline"]
        self.args_without_params = []
        self.args_dict = {"output": "o", "input": "i", "sample_name": "s", "pipeline": "p"}
        self.reverse_args_dict = {value: key for key, value in self.args_dict.iteritems()}
        self.args_description = {"i": "Comma separated list of input fastq files Supports " +
                                      ",".join(self.file_extension_list) + ".",
                                 "s": "Sample name",
                                 "o": "Output directory to store counts",
                                 "p": "Pipeline to run"
                                 }
        self.options = ["s", "o", "p"]

    def usage(self):
        Common().super_usage()
        print("\nUsage:")
        print("-------")
        print("\tpython expressed_dna.py [options]* -i||--input <Single end fq file or Paired end fq files separated by comma> -p|--pipeline <Pipeline to be run>")
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

        input_file_split = self.input.split(",")

        if len(input_file_split) == 1:
            print >> sys.stderr, "Running in single end mode"

        for input_file in input_file_split:
            if not os.path.isfile(input_file):
                Common().print_error("File " + input_file + " does not exist")

            if not any(input_file.endswith(x) for x in self.file_extension_list):
                Common().print_error(
                    "Invalid file: " + input_file + ", Must be a *" + ' ,*'.join(self.file_extension_list))



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

        if self.sample_name == "":
            self.sample_name = os.path.basename(input_files[0].split(".fastq")[0])

        # Disabled for single thread running
        # pool = multiprocessing.Pool(2)

        # Create a results list and send them for execution
        results = list()
        results.append(get_pos_base_counts(self.output, input_files, self.sample_name, self.pipeline))
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
    expressed_dna = ExpressedDNA()

    try:
        opts, args = getopt.getopt(argv, ":".join(expressed_dna.args_dict.values()) + ":", map(lambda x: x + "=",
                                                                                         expressed_dna.args_dict.keys()))
    except getopt.GetoptError as exp:
        Common().print_error("Invalid input argument: " + exp.msg)

    if len(opts) == 0:
        expressed_dna.usage()

    for opt, arg in opts:
        if opt in map(lambda x: "-" + x, expressed_dna.args_dict.values()) or opt in map(lambda x: "--" + x, expressed_dna.
                args_dict.keys()):
            if opt in map(lambda x: "-" + x, expressed_dna.args_dict.values()):
                clte_member = expressed_dna.args_dict.keys()[expressed_dna.args_dict.values().index(opt[1:])]
            else:
                clte_member = opt[2:]

            if clte_member in expressed_dna.args_without_params:
                setattr(expressed_dna, clte_member, not getattr(Default(), expressed_dna.args_dict[clte_member]))
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
                setattr(expressed_dna, clte_member, tryout)
        else:
            Common().print_error("Invalid input argument " + opt)

    print >> sys.stderr, "Evaluating arguments."
    expressed_dna.check_input()
    expressed_dna.check_output()

    print >> sys.stderr, "Processing "+expressed_dna.sample_name
    expressed_dna.process()
    print >> sys.stderr, expressed_dna.sample_name + " processing complete"


if __name__ == "__main__":
    main(sys.argv[1:])
