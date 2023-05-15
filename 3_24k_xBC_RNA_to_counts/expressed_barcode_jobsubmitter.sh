#!/bin/bash

#$ -M viveksagar.kr@novartis.com
#$ -l h_rt=360000
#$ -pe smp 2
#$ -l m_mem_free=40G

module load pulsar

EXPRESSED_HOME=/da/onc/krishvi7/bitbucket/oncp-expressed

/usr/prog/onc/pulsar/venv_py2/bin/python $EXPRESSED_HOME/expressed_counter.py --sample_name $1 --input $2 --tenx_fastq $3 --pipeline $5 --output $4