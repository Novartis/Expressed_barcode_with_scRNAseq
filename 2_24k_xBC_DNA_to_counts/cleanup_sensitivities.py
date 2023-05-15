import sys

log2_sensitivity = sys.argv[1]
counts_summary = sys.argv[2]

sensitivity_dict = dict()

with open(counts_summary) as fh:
    for line in fh:
        ls = line.strip().split()
        treatment = ls[7]
        cloneid = ls[9]
        cleaned_tp = ls[14]
        sensitivity = ls[13]
        if sensitivity != "NA":
            if cloneid not in sensitivity_dict:
                sensitivity_dict[cloneid] = dict()
            if treatment not in sensitivity_dict[cloneid]:
                sensitivity_dict[cloneid][treatment] = dict()
            if cleaned_tp not in sensitivity_dict[cloneid][treatment]:
                sensitivity_dict[cloneid][treatment][cleaned_tp] = sensitivity
fh.close()

with open(log2_sensitivity) as fh:
    print(fh.readline().strip()+"\tsensitivity\tt20h\tt10h\tt20hr\tt10hr")
    for line in fh:
        ls = line.strip().split()
        cloneid = ls[0]
        cleaned_tp = ls[7]

        sensitivity = "NA"
        t20h = "False"
        t20hr = "False"
        t10h = "False"
        t10hr = "False"

        if cloneid in sensitivity_dict:
            # Check if the clone has the holiday timepoint
            if "T" in sensitivity_dict[cloneid]:
                sensitivity = sensitivity_dict[cloneid]["T"]["0t"]
            if "H" in sensitivity_dict[cloneid] and "20" in sensitivity_dict[cloneid]["H"]:
                t20h = "True"
            if "H" in sensitivity_dict[cloneid] and "10" in sensitivity_dict[cloneid]["H"]:
                t10h = "True"
            if "HR" in sensitivity_dict[cloneid] and "20" in sensitivity_dict[cloneid]["HR"]:
                t20hr = "True"
            if "HR" in sensitivity_dict[cloneid] and "10" in sensitivity_dict[cloneid]["HR"]:
                t10hr = "True"

            treatments = ";".join([treatment for treatment in sensitivity_dict[cloneid]])

        print(line.strip()+"\t"+sensitivity+"\t"+t20h+"\t"+t10h+"\t"+t20hr+"\t"+t10hr)
fh.close()