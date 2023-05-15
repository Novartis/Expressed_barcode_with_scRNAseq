import sys
import multiprocessing
import os
import subprocess

metadata = sys.argv[1]
input_directory = sys.argv[2]
output_dir = sys.argv[3]
sep = ","
csv = False
qsub_mode = False
#pipeline = "expressed_barcode"
pipeline = "cropseq"

query_dict = dict()

with open(metadata, "r") as fh:
    fh.readline()
    for line in fh:
        if csv:
            sep = ','
        line_strip = line.strip().split(sep)
        sample = line_strip[0]
        name = line_strip[1]
        input_prefix = input_directory
        fastq = [input_prefix + "/" + f for f in os.listdir(input_prefix) if ".fastq.gz" in f and "R2_001" in f and sample in f][0]
        query_dict[sample] = (name, fastq)

#pool_query = multiprocessing.Pool(2)
#results_query = list()

for sample_key in query_dict:
    output_dir_sample = output_dir+"/"+sample_key
    if not os.path.isdir(output_dir_sample):
        os.mkdir(output_dir_sample)

    cmd = ["python", os.environ["EXPRESSED_HOME"]+"/expressed_counter.py",  "--sample_name", query_dict[sample_key][0],
           "--input", query_dict[sample_key][1], "--tenx_fastq", query_dict[sample_key][1].replace("_R2_001", "_R1_001"),
           "--pipeline", pipeline,  "--output", output_dir_sample]

    logfile = open(output_dir_sample+ "/" + "command_log.txt", 'w')
    run_command = " ".join(["qsub", "-o", output_dir_sample + "/job_op.txt",
                            "-e", output_dir_sample + "/job_err.txt",
                            os.environ["EXPRESSED_HOME"] + "/expressed_barcode_jobsubmitter.sh", query_dict[sample_key][0],
                            query_dict[sample_key][1], query_dict[sample_key][1].replace("_R2_001", "_R1_001"),
                            output_dir_sample, pipeline])
    print(" ".join(cmd), file=logfile)
    print(run_command, file=logfile)
    logfile.close()
    subprocess.call(run_command, shell=True)
    #results_query.append(pool_query.apply_async(subprocess.call, ([cmd])))

# for result in results_query:
#     # TODO- Handle return messages
#     return_message = result.get()