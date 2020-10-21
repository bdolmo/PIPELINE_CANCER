#!/usr/bin/env python3

import os
import sys
import re
import logging
import gzip
import shutil
import csv
from collections import defaultdict
from pathlib import Path
from scipy import stats
import numpy as np
import subprocess
from modules import utils as u
from modules import params as p
from modules import trimming as t

def do_var_call():

    if p.analysis_env['VARIANT_CLASS'] == 'somatic':

        do_mutect2()

        do_manta()

def do_manta():

    for sample in p.sample_env:

        # Create BAM_FOLDER if not present
        p.sample_env[sample]['VCF_FOLDER'] = \
         p.sample_env[sample]['SAMPLE_FOLDER'] + "/" + "VCF_FOLDER"

        vcf_folder_path = Path(p.sample_env[sample]['VCF_FOLDER'])
        if not vcf_folder_path.is_dir():
            os.mkdir(vcf_folder_path)

        # This is the predefined output from manta
        p.sample_env[sample]['MANTA_TUMOR'] = \
            p.sample_env[sample]['VCF_FOLDER'] + "/results/variants/tumorSV.vcf.gz"
        p.sample_env[sample]['MANTA_TUMOR_NAME'] = "tumorSV.vcf.gz"

        # Renamed manta output for our sample
        p.sample_env[sample]['MANTA_VCF'] = \
            p.sample_env[sample]['VCF_FOLDER'] + "/" + sample + ".manta.vcf"
        p.sample_env[sample]['MANTA_VCF_NAME'] = sample + ".manta.vcf"

        # Ready SV file
        p.sample_env[sample]['READY_SV_VCF'] = \
            p.sample_env[sample]['VCF_FOLDER'] + "/" + sample + ".manta.vcf"
        p.sample_env[sample]['READY_SV_VCF_NAME'] = sample + ".manta.vcf"

        if not os.path.isfile(p.sample_env[sample]['MANTA_VCF']):
            
            bashCommand = ('{} --tumorBam {} --referenceFasta {} --exome --runDir {}').format(
             p.system_env['MANTA_CONFIG'], p.sample_env[sample]['READY_BAM'], \
             p.aux_env['GENOME_FASTA'], p.sample_env[sample]['VCF_FOLDER'] )
            print (bashCommand)
            logging.info(bashCommand)
            process = subprocess.Popen(bashCommand,#.split(),
                shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            output, error = process.communicate()

            if not error.decode('UTF-8'):
                msg = " INFO: Manta config was created successfully for sample " + sample
                print (msg)
                logging.info(msg)
            else:
                msg = " ERROR: Manta config could not be created for sample " + sample
                print (msg)
                logging.error(msg)

            run_manta_script = p.sample_env[sample]['VCF_FOLDER'] + "/" + "runWorkflow.py"

            bashCommand = ('{} --j {}').format(run_manta_script, p.analysis_env['THREADS'])
            msg = " INFO: Calling SVs with Manta on sample " + sample
            print(msg)
            logging.info(msg)
            logging.info(bashCommand)
            print (bashCommand)

            process = subprocess.Popen(bashCommand,#.split(),
                shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            output, error = process.communicate()

            if not error.decode('UTF-8'):
                msg = " INFO: Manta was executed successfully for sample " + sample
                print (msg)
                logging.info(msg)
            else:
                msg = " ERROR: Manta could not be executed for sample " + sample
                print (msg)
                logging.error(msg)

        if not os.path.isfile(p.sample_env[sample]['MANTA_VCF']):
            print(p.sample_env[sample]['MANTA_TUMOR'])
            input = gzip.GzipFile(p.sample_env[sample]['MANTA_TUMOR'], 'rb')
            s = input.read()
            input.close()

            output = open(p.sample_env[sample]['MANTA_VCF'], 'wb')
            output.write(s)
            output.close()
        else:
            msg = " INFO: Skipping Manta analysis for "+ sample
            print(msg)
            logging.info(msg)

def do_mutect2():

    for sample in p.sample_env:

        # Create BAM_FOLDER if not present
        p.sample_env[sample]['VCF_FOLDER'] = \
         p.sample_env[sample]['SAMPLE_FOLDER'] + "/" + "VCF_FOLDER"

        vcf_folder_path = Path(p.sample_env[sample]['VCF_FOLDER'])
        if not vcf_folder_path.is_dir():
            os.mkdir(vcf_folder_path)

        p.sample_env[sample]['MUTECT2_VCF'] = p.sample_env[sample]['VCF_FOLDER'] + \
            "/" + sample + ".mutect2.vcf"
        p.sample_env[sample]['MUTECT2_VCF_NAME'] =  sample +  ".mutect2.vcf"

        p.sample_env[sample]['READY_SNV_VCF'] = p.sample_env[sample]['MUTECT2_VCF']
        p.sample_env[sample]['READY_SNV_VCF_NAME'] = p.sample_env[sample]['MUTECT2_VCF_NAME']

        bashCommand = ('{} run -v {}:/bam_data/ -v {}:/vcf_data/ -v {}:/bundle/ '
            '-v {}:/panel_data/ -it {} gatk Mutect2 -I /bam_data/{} '
            '-O /vcf_data/{} -L /panel_data/{} -R /bundle/{} --native-pair-hmm-threads {}'.format(p.system_env['DOCKER'],
            p.sample_env[sample]['BAM_FOLDER'],  p.sample_env[sample]['VCF_FOLDER'], p.defaults['BUNDLE_FOLDER'], p.defaults['PANEL_FOLDER'],
            p.docker_env['GATK'], p.sample_env[sample]['READY_BAM_NAME'], p.sample_env[sample]['MUTECT2_VCF_NAME'], \
            p.analysis_env['PANEL_LIST_NAME'], p.aux_env['GENOME_NAME'], p.analysis_env['THREADS']
            ))
        
        if not os.path.isfile(p.sample_env[sample]['MUTECT2_VCF']):
            msg = " INFO: Running Mutect2 for sample "+ sample
            print(msg)
            logging.info(msg)
            logging.info(bashCommand)

            process = subprocess.Popen(bashCommand,#.split(),
                shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = process.communicate()
            if not error.decode('UTF-8'):
                if re.search(r'Mutect2 done.', output.decode('UTF-8')):
                    msg = " INFO: Mutect2 varcall ended OK for sample " + sample
                    print (msg)
                    logging.info(msg)
                    ok = True
                else:
                    msg = " ERROR: Something went wrong with Mutect2 varcall for sample " + sample
                    print (msg)
                    logging.error(msg)
            else:
                msg = " ERROR: Something went wrong with Mutect2 varcall for sample " + sample
                print (msg)
                logging.error(msg)
        else:
            msg = " INFO: Skipping Mutect2 analysis for "+ sample
            print(msg)
            logging.info(msg)