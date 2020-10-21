#!/usr/bin/env python3

import os
import sys
import shutil
import re
import logging
import gzip
import csv
from collections import defaultdict
from pathlib import Path
from scipy import stats
import numpy as np
import subprocess
from civicpy import civic, exports
import requests

from modules import utils as u
from modules import params as p
from modules import trimming as t
from modules import sqlite as s

def do_annotation():

    do_vep()

    #do_civic()


def do_civic():

    for sample in p.sample_env:

        p.sample_env[sample]['CIVIC_VCF'] = \
            p.sample_env[sample]['READY_SNV_VCF'].replace(".vcf", ".civic.vcf")
        p.sample_env[sample]['CIVIC_VCF_NAME'] = \
            os.path.basename(p.sample_env[sample]['READY_SNV_VCF_NAME']).replace(".vcf", ".vep.vcf")

        p.sample_env[sample]['READY_SNV_VCF'] = p.sample_env[sample]['CIVIC_VCF']
        p.sample_env[sample]['READY_SNV_VCF_NAME'] = p.sample_env[sample]['CIVIC_VCF_NAME']        

def do_vep():
    '''
    VEP annotation of vcf files
    '''

    for sample in p.sample_env:

        p.sample_env[sample]['VEP_VCF'] = \
            p.sample_env[sample]['READY_SNV_VCF'].replace(".vcf", ".vep.vcf")

        p.sample_env[sample]['VEP_VCF_NAME'] = \
            os.path.basename(p.sample_env[sample]['READY_SNV_VCF_NAME']).replace(".vcf", ".vep.vcf")

        vep_output_vcf = p.defaults['VEP_DATA_OUTPUT'] + "/" + p.sample_env[sample]['VEP_VCF_NAME']

        bashCommand = ('{} run -t -i -v {}:/opt/vep/.vep ensemblorg/ensembl-vep'
        ' perl vep --cache --offline --dir_cache /opt/vep/.vep/ --dir_plugins /opt/vep/.vep/Plugins/'
        ' --input_file /opt/vep/.vep/input/{} --output_file /opt/vep/.vep/output/{} '
        ' --af_1kg --af_gnomad '
        ' --format vcf --vcf --hgvs --hgvsg --max_af --pubmed --gene_phenotype --ccds --sift b --polyphen b'
        ' --symbol --force_overwrite --fork {} --canonical '
        .format(p.system_env['DOCKER'], p.defaults['VEP_DATA'], p.sample_env[sample]['READY_SNV_VCF_NAME'],\
        p.sample_env[sample]['VEP_VCF_NAME'], p.analysis_env['THREADS']))

        if not os.path.isfile(p.sample_env[sample]['VEP_VCF']):

            # First, copy vcf to vep input dir
            shutil.copy2(p.sample_env[sample]['READY_SNV_VCF'], p.defaults['VEP_DATA_INPUT']+"/"+p.sample_env[sample]['READY_SNV_VCF_NAME'])

            msg = " INFO: Annotating sample " + sample + " with VEP"
            print (msg)
            print (bashCommand)
            logging.info(msg)
            logging.info(bashCommand)

            process = subprocess.Popen(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = process.communicate()
            if not error.decode('UTF-8'):
                pass
            else:
                msg = " ERROR: VEP annotation failed for sample" + sample
                print(msg)
                logging.error(msg)

            shutil.copy2(vep_output_vcf, p.sample_env[sample]['VEP_VCF'])
        else:
            msg = " INFO: Skipping VEP annotation for "+ sample
            print(msg)
            logging.info(msg)

        p.sample_env[sample]['READY_SNV_VCF'] = p.sample_env[sample]['VEP_VCF']
        p.sample_env[sample]['READY_SNV_VCF_NAME'] = p.sample_env[sample]['VEP_VCF_NAME']        
