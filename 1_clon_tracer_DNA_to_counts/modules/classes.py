__author__ = 'krishvi7'
import os
import sys
import os.path
import re
import datetime
import barcode_count
import simplify
import matricize
import tall_skinny
import cw_qualitysummary
import multiprocessing
import math
import glob

# Contains functions used across the package
class Common():

    # Default values for paramters
    output_prefix = "output"
    barcode_index = "AGCAGAG"
    phred_format = 33
    barcode_length = 30

    def __init__(self):
        self.dataType_dict = {"input_files": str,
                              "output_prefix": str,
                              "metadata_file": str,
                              "output_directory": str,
                              "phred_format": int,
                              "barcode_index": str,
                              "sample_pattern_regex": str,
                              "known_contaminant": str,
                              "suppress_matrix": bool,
                              "suppress_quality_cw_summary": bool,
                              "suppress_tall_skinny": bool,
                              "library_mode": bool,
                              "barcode_length": int}
        self.range_dict = {"phred_format": [33, 64]}
        self.paramless_dict = {"suppress_quality_cw_summary": True,
                               "suppress_matrix": True,
                               "suppress_tall_skinny": True,
                               "library_mode": True}
        self.metadata_dict = {"SMFID": "sample_name",
                              "Fastq": "fastq_files"}
        return

    @staticmethod
    def print_error(message):
        print >> sys.stderr, "ERROR: "+message
        sys.exit(1)

    @staticmethod
    def get_timestamp():
        return datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')

    @staticmethod
    def super_usage():
        print("\nClonTracer package V1.0- ClonTracer package is used to analyze barcode data obtained from ClonTracer"
              " library")


    @staticmethod
    def get_file_by_pattern(input):
        return glob.glob(input)[0]


# Contains code to validate a list of input files
class ClonTracerCountBarcodes():

    def __init__(self, delimiter=".", input_files="", input_file_list=list(), output_prefix="output",
                 barcode_index=Common.barcode_index, known_contaminant="", phred_format=Common.phred_format, library_mode=False,
                 barcode_length=Common.barcode_length):
        self.delimiter = delimiter
        self.input_file_list = input_file_list
        self.input_files = input_files
        self.output_prefix = output_prefix
        self.barcode_index = barcode_index
        self.phred_format = phred_format
        self.module = "clonTracer_countBarcodes"
        self.barcode_file = ""
        self.test_mode = False
        self.library_mode = library_mode
        self.known_contaminant = known_contaminant
        self.barcode_length = barcode_length
        self.file_extension_list = [".fastq", ".fq", ".fastq.gz", ".fastq.bz2"]
        self.args_with_params = {"input_files", "output_prefix", "barcode_index", "known_contaminant"}
        self.args_dict = {"input_files": "i", "output_prefix": "o", "phred_format": "q", "barcode_index": "b", "barcode_length": "x",
                          "known_contaminant": "k"}
        self.args_dict_paramless = {"library_mode": "l"}
        self.args_dict_all = self.args_dict.copy()
        self.args_dict_all.update(self.args_dict_paramless)
        self.reverse_args_dict_all = {value: key for key, value in self.args_dict_all.iteritems()}
        self.args_description = {"i": "Comma separated list of fastq files. Supports "+",".join(self.file_extension_list),
                                 "o": "Output prefix",
                                 "q": "phred quality format. must be 33 or 64",
                                 "b": "primer index following WS barcode",
                                 "k": "Text file containing known contaminant barcodes",
                                 "x": "Length of the barcode",
                                 "l": "Library mode excludes barcode merging step"}
        self.options = ["o", "q", "b", "k", "l"]
        self.args_default = {"o": "output", "b": "AGCAGAG", "q": "33"}

    def usage(self):
        Common.super_usage()
        print("\nUsage:")
        print("-------")
        print("  python clonTracer_countBarcodes.py [options]* -i|--input_files <ip_files>")
        print("\n  -i|--input_files "+self.args_description["i"])
        print("\noptions:")
        print("--------")
        for option in self.options:
            print ("  -"+option+"|--"+self.reverse_args_dict_all[option]+"\t: "+self.args_description[option] +
                   (". Default- "+self.args_default[option] if option in self.args_default.keys() else ""))
        print ("")
        sys.exit(1)

    def check_input_files(self):
        if self.input_files == "":
            Common.print_error("Input file(s) for barcode counting not supplied")

        for input_file in self.input_files.split(","):
            if not os.path.isfile(input_file):
                Common.print_error("Input file "+input_file+" does not exist")

            if not any(input_file.endswith(x) for x in self.file_extension_list):
                Common.print_error("Invalid input "+input_file+", Must be a *"+' ,*'.join(self.file_extension_list))

    def count_barcodes(self):
        self.barcode_file = self.output_prefix+".barcode_verbose.txt"
        print("\nCounting barcodes on sample "+self.input_files+"....")
        barcode_count.barcode_count(self.input_files, self.barcode_file, self.barcode_index, self.known_contaminant,
                                    self.phred_format, self.library_mode)

    def simplify_barcodes(self):
        print("Simplifying barcode from sample "+self.output_prefix+"....")
        simplify.simplify_barcodes(self.barcode_file, self.output_prefix+".barcode.txt", self.output_prefix +
                                   ".cumulative_wealth.txt")


# Contains code to validate an input directory and submits individual samples
class ClonTracerCountBarcodesExperiment():

    def __init__(self):
        self.metadata_file = ""
        self.output_directory = Common.output_prefix+"_"+Common.get_timestamp()
        self.smfuid_col = -1
        self.fastq_col = -1
        self.suppress_quality_cw_summary = False
        self.library_mode = False
        self.suppress_matrix = False
        self.suppress_tall_skinny = False
        self.barcode_index = Common.barcode_index
        self.known_contaminant = ""
        self.phred_format = Common.phred_format
        self.sample_pattern_regex = ""
        self.test_mode = False
        self.annotation_dict = dict()
        self.annotation_header = ""
        self.barcode_length = Common.barcode_length
        self.module = "clonTracer_countBarcodes_experiment"
        self.args_dict = {"metadata_file": "i", "output_directory": "o", "phred_format": "q",
                          "barcode_index": "b",  "known_contaminant": "k", "sample_pattern_regex": "p",
                          "barcode_length": "x"}
        self.args_dict_paramless = {"suppress_quality_cw_summary": "c", "suppress_matrix": "m",
                                    "suppress_tall_skinny": "s", "library_mode": "l"}
        self.args_dict_all = self.args_dict.copy()
        self.args_dict_all.update(self.args_dict_paramless)
        self.reverse_args_dict_all = {value: key for key, value in self.args_dict_all.iteritems()}
        self.args_description = {"i": "sample key containing information of sample name and input fastq files",
                                 "o": "Output directory",
                                 "c": "Do not produce the quality and cumulative summary file",
                                 "m": "Do not produce barcodes x samples matrix of barcode fraction",
                                 "s": "Do not produce tall-skinny version of barcode-sample matrix",
                                 "q": "phred quality format. must be 33 or 64",
                                 "b": "primer index following WS barcode",
                                 "k": "Text file containing known contaminant barcodes",
                                 "p": "Sample pattern regex",
                                 "x": "Length of the barcode",
                                 "l": "Library mode excludes barcode merging step"}
        self.options = ["o", "c", "m", "s", "q", "b", "k", "p", "z", "l"]
        self.args_default = {"o": "./output_<datetimevalue>", "b": "AGCAGAG", "q": "33"}

    def usage(self):
        Common.super_usage()
        print("\nUsage:")
        print("-------")
        print("  python clonTracer_countBarcodes_experiment.py [options]* -i|--input_directory <ip_dir>")
        print("\n  -i|--input_directory: "+self.args_description["i"])
        print("\noptions:")
        print("--------")
        for option in self.options:
            print ("  -"+option+"|--"+self.reverse_args_dict_all[option]+"\t: "+self.args_description[option] +
                   (". Default- "+self.args_default[option] if option in self.args_default.keys() else ""))

        print("")
        sys.exit(1)

    # Checks the existence of input directory and prints appropriate messages
    def check_metadata_file(self):

        fastq_col_found = False
        smfuid_col_found = False

        print("Checking metadata file")
        if self.metadata_file == "":
            Common.print_error("Metadata file for barcode counting not supplied")

        if not os.path.isfile(self.metadata_file):
            Common.print_error("Metadata file "+self.metadata_file + " does not exist")

        with open(self.metadata_file) as fh:
            header = fh.readline().strip().split("\t")

            for col in header:
                if col == Common().metadata_dict["Fastq"]:
                    fastq_col_found = True
                    self.fastq_col = header.index(col)
                elif col == Common().metadata_dict["SMFID"]:
                    smfuid_col_found = True
                    self.smfuid_col = header.index(col)
        fh.close()

        if not fastq_col_found or not smfuid_col_found:
            Common.print_error("Metadata file doesn't have expected columns")

        print("Metadata file "+self.metadata_file+" looks good")

    # Checks the existence of the output directory
    def check_output_directory(self):
        if not os.path.exists(self.output_directory):
            if not self.test_mode:
                os.makedirs(self.output_directory)
        print("Output would be stored in "+self.output_directory)

    # Checks the existence of the output directory
    def check_known_contaminant_file(self):
        if self.known_contaminant != "":
            if not os.path.exists(self.known_contaminant):
                Common.print_error("Contaminant file " + self.known_contaminant+" does not exist")

    # Calls method to generate cw and quality summary
    def cw_quality_summary(self):
        print("Collating quality and cumulative summary for the experiment....")
        cw_qualitysummary.cw_qualitysummary(self.output_directory, self.annotation_dict, self.annotation_header,
                                            "cw_summary.txt", "quality_summary.txt")

    # Calls method to generate tall skinny file
    def tall_skinny(self):

        if not self.test_mode:
            print("Generating tall-skinny barcode fractions....")
            tall_skinny.tall_skinny(self.output_directory, self.annotation_dict, self.annotation_header,
                                    "tall_skinny.txt")
        else:
            print("Count and Fraction matrix would be generated")

    # Call method to create barcode matrix
    def matricize(self):

        if not self.test_mode:
            print("Generating count and fraction matrix....")
            matricize.matricize(self.output_directory, "matrix_count.txt", "matrix_cpm.txt")
        else:
            print("Count and Fraction matrix would be generated")

    # Populate the sample_key dict
    def populate_sample_key(self):
        if os.path.exists(self.metadata_file):
            fh = open(self.metadata_file)
            first_line = fh.readline().strip().split("\t")
            self.annotation_header = "\t".join(first_line[3:])
            for line in fh:
                line_split = line.strip().split("\t")
                self.annotation_dict[line_split[1]] = "\t".join(line_split[3:])
            fh.close()
        else:
            Common.print_error("Sample key file "+self.metadata_file+" not found")

    # Go over all the fastq files in the directory and count the barcodes
    def count_barcodes_on_samples(self):

        clt_obj_list = list()

        with open(self.metadata_file) as fh:
            fh.readline()
            for line in fh:

                ls = line.strip().split("\t")

                sample = ls[self.smfuid_col]
                fastq_files = ls[self.fastq_col].split(",")

                # Check for regex if any
                if self.sample_pattern_regex != "":
                    if not re.match(self.sample_pattern_regex, sample):
                        continue

                sample_name = sample
                clt_object = ClonTracerCountBarcodes(input_file_list=list(),
                                                     output_prefix=self.output_directory+"/"+sample_name,
                                                     barcode_index=self.barcode_index,
                                                     phred_format=self.phred_format,
                                                     known_contaminant=self.known_contaminant,
                                                     library_mode=self.library_mode,
                                                     barcode_length=self.barcode_length)

                for fastq_file in fastq_files:
                    clt_object.input_file_list.append(Common.get_file_by_pattern(fastq_file.replace("\"", "")))

                clt_object.barcode_file = clt_object.output_prefix + ".barcode_verbose.txt"
                clt_object.input_files = ",".join(clt_object.input_file_list)
                clt_obj_list.append(clt_object)

        fh.close()

        pools = multiprocessing.Pool(max(int(math.ceil(len(clt_obj_list)/2)), 6))
        results = list()

        for obj in clt_obj_list:
            print("\nCounting barcodes on sample " + obj.input_files + "....")
            # results.append(pools.apply_async(barcode_count.barcode_count, (obj.input_files, obj.barcode_file,
            #                                                                obj.barcode_index, obj.known_contaminant,
            #                                                                obj.phred_format, obj.library_mode)))

            barcode_count.barcode_count(obj.input_files, obj.barcode_file, obj.barcode_index, obj.known_contaminant,
                                        obj.phred_format, obj.library_mode)

        # for result in results:
        #     has_error = result.get()
        #     if has_error:
        #         print >> sys.stderr, "Exception while counting barcodes"
        #         sys.exit(1)

        # results = list()

        for obj in clt_obj_list:
            print("\nSimplifying barcodes from sample " + obj.output_prefix + "....")
            # results.append(pools.apply_async(simplify.simplify_barcodes, (obj.barcode_file,
            #                                                               obj.output_prefix + ".barcode.txt",
            #                                                               obj.output_prefix + ".cumulative_wealth.txt"))
            #                )
            simplify.simplify_barcodes(obj.barcode_file, obj.output_prefix + ".barcode.txt",
                                         obj.output_prefix + ".cumulative_wealth.txt")
        # for result in results:
        #     has_error = result.get()
        #     if has_error:
        #         print >> sys.stderr, "Exception while counting barcodes"
        #         sys.exit(1)


