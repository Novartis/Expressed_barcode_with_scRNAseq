import sys
import multiprocessing
import os
import subprocess
import glob

# to run wrapper
# python expressed_dna_wrapper.py <metadata_file> <Dave_ruddy's_sequencing_directory> <output_directory>
# Sample_Key,Name
# Sample_Key- Custom name given to the sample

config = sys.argv[1]  # Provide the full path
sequence_dir = sys.argv[2]  # Provide the full path
output_dir = sys.argv[3]  # Provide the full path
append = False  # Whether to  Sample_ /Project_ to the beginning of every sequence
query_dict = dict()
sep = "\t"
csv = True
qsub_mode = True
pipeline = "expressed_barcode"

with open(config, "r") as fh:
    fh.readline()
    for line in fh:
        if csv:
            sep = ','
        line_strip = line.strip().split(sep)
        sample = line_strip[0]
        fastq_sets = glob.glob(os.path.join(sequence_dir, sample+"*fastq.gz"))
        #print(fastq_sets)
        if len(fastq_sets) == 1:
            fastq = fastq_sets
            query_dict[sample] = fastq
        #fastq = line_strip[2]
        # if len(line_strip) >= 3:
        #     fastq = line_strip[2]
        #     query_dict[sample] = (name, fastq)
        # else:
        #query_dict[sample] = fastq

    # Check if a job needs to be submitted to the cluster
    if qsub_mode:
        input_prefix = sequence_dir
        for sample_key in query_dict:
            # if len(query_dict[sample_key]) == 1:
            sample_name = sample_key
            input_files = query_dict[sample_key]
            #input_files = [input_prefix + "/" + query_dict[sample_key]]
            # elif len(query_dict[sample_key]) == 2:
            #     sample_name = query_dict[sample_key][0]
            #     input_files = [input_prefix + "/" + query_dict[sample_key][1]]
            # else:
            #     input_files = [input_prefix + "/" + f for f in os.listdir(input_prefix) if ".fastq.gz" in f]
            #     sample_name = query_dict[sample_key]

            # print(" ".join(["qsub", "-o", output_dir+"/"+sample_name+"_job_op.txt",
            #                           "-e", output_dir+"/"+sample_name+"_job_err.txt",
            #                           os.environ["EXPRESSED_DNA_HOME"]+"/expressed_dna_jobsubmitter.sh",
            #                           sample_name,
            #                           ",".join(sorted(input_files)),
            #                           output_dir,
            #                           pipeline]))
            cmd = " ".join(["qsub", "-o", output_dir+"/"+sample_name+"_job_op.txt",
                                      "-e", output_dir+"/"+sample_name+"_job_err.txt",
                                      os.environ["EXPRESSED_DNA_HOME"]+"/expressed_dna_jobsubmitter.sh",
                                      sample_name,
                                      ",".join(sorted(input_files)),
                                      output_dir+"/"+sample_name,
                                      pipeline])
            #print(cmd)
            subprocess.call(cmd, shell=True)
    else:
        # Run locally by creating pool
        pool_query = multiprocessing.Pool(2)
        for sample_key in query_dict.keys():
            results_query = list()
            input_prefix = sequence_dir
            if len(query_dict[sample_key]) == 2:
                sample_name = query_dict[sample_key][0]
                input_files = [input_prefix + "/" + query_dict[sample_key][1]]
            else:
                input_files = [input_prefix + "/" + f for f in os.listdir(input_prefix) if ".fastq.gz" in f]
                sample_name = query_dict[sample_key]

            results_query.append(pool_query.apply_async(subprocess.call, ([["python", os.environ["EXPRESSED_DNA_HOME"]
                                                                            + "/expressed_dna.py",
                                                                            "--sample_name", sample_name,
                                                                            "--pipeline", pipeline,
                                                                            "--input", ",".join(sorted(input_files)),
                                                                            "--output", output_dir]])))
        for result in results_query:
            # TODO- Handle return messages
            return_message = result.get()






