clonTracer_v0.0a

Installation:
	Requirements
		- Mac or Unix OS
		- Python 2 (does not work with Python 3)
			Modules:
				- subprocess
				- re	
				- datetime
				- glob
				- numpy
				- math
				- getopt

	Download the clonTracer_analyze zip file and extract the content


Running the scripts:

	Experiment level:

		Metadata file with sample_name and fastq_files columns

	    python clonTracer_countBarcodes_experiment.py [options]* -i|--input_directory <ip_dir>

		 -i|--metadata_file: Metadata file containing columns as described in the documentation
 		
		options:
		--------
		 -o|--output_directory 		      : Output directory. Default- ./output_<datetimevalue>
		 -c|--suppress_quality_cw_summary : Do not produce the quality and cumulative summary file
		 -m|--suppress_matrix_fractions   : Do not produce barcodes x samples matrix of barcode fraction
		 -q|--phred_format		          : phred quality format. must be 33 or 64. Default- 33
		 -s|--suppress_tall_skinny        : Do not produce tall-skinny version of barcode-sample matrix
		 -b|--barcode_index 		      : primer index following WS barcode. Default- AGCAGAG
		 -k|--known_contaminant 	      : Text file containing known contaminant barcodes
		 -p|--sample_pattern_regex 	      : Sample pattern regex. Python style regex


	Sample level:

	    python clonTracer_countBarcodes.py [options]* -i|--input_file_list <ip_file_list>
	
       	 -i|--input_file_list Comma separated list of fastq files
 
		options:
	 	--------
		 -o|--output_prefix		: Output prefix
		 -q|--phred_format		: phred quality format. must be 33 or 64. Default- 33
		 -b|--barcode_index		: primer index following  barcode. Default- AGCAGAG
		 -k|--known_contaminant	: Text file containing known contaminant barcodes


Sample Run:
	
	Experiment level
	
		python <clonTracer_home>/clonTracer_countBarcodes_experiment.py -i test
	
	Sample level

		python <clonTracer_home>/clonTracer_countBarcodes.py -i <clonTracer_home>/test/Sample1/Sample1.fastq.gz

	where <clonTracer_home> is the downloaded clonTracer directory
	

Output files:

	Output files by default will have the same prefix as the input file. For a given sample, the following files are generated
	
	Sample level:
	
		a. *.barcode_verbose.txt	: Contains barcode, count, merged_count (with 1 or 2 hamming distance), fraction and details of barcodes merged in to parent. It also contains the quality summary header
		b. *.barcode.txt         	: Summary of *.barcode_verbose.txt. Contains barcode, merged_count and fraction
		c. *.cumulative_wealth.txt  : Contains the information of number of unique barcodes needed to reach a cumulative wealth of 'X'%
	
	Experiment level:
	
		a. quality_summary.txt  : Summary of quality header for all the samples
		b. cw_summary.txt       : Summary of cumulative wealth for all the samples
		c. matrix_count.txt     : Barcode count matrix of barcode x sample
		d. matrix_fraction.txt  : Barcode fraction matrix of barcode x sample
		e. tall_skinny.txt      : Tall skinny version of count and fraction matrices