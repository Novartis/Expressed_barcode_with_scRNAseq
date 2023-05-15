#!/bin/bash

#$ -M viveksagar.kr@novartis.com
#$ -m a
#$ -l h_rt=36000
#$ -pe smp 2
#$ -l m_mem_free=20G

module load pulsar

/usr/prog/onc/pulsar/venv_py2/bin/python $EXPRESSED_DNA_HOME/expressed_dna.py --sample_name $1 --input $2 --output $3 --pipeline $4