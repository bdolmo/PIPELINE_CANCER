#!/usr/bin/env python3

import os
import sys
import re
import logging
from collections import defaultdict
from pathlib import Path
import subprocess
from civicpy import civic, exports
from modules import utils as u
from modules import sqlite as s

global sample_env 
sample_env = defaultdict(dict)

def update_variant_databases():

    if not os.path.isfile(analysis_env['CIVIC_VCF']):
        with open(analysis_env['CIVIC_VCF'], 'w', newline='') as file:
            w = exports.VCFWriter(file)
            all_variants = civic.get_all_variants()
            w.addrecords(all_variants)
            w.writerecords()


def set_analysis_env(args):


    global analysis_env 
    #analysis_env = defaultdict(dict)

    analysis_env = {
        'INPUT_DIR'      : os.path.abspath(args.input_dir),
        'OUTPUT_DIR'     : os.path.abspath(args.output_dir),
        'OUTPUT_NAME'    : os.path.basename(args.output_dir),
        'PANEL'          : os.path.abspath(args.panel),
        'PANEL_NAME'     : os.path.basename(args.panel),
        'PANEL_LIST'     : defaults['PANEL_FOLDER'] + "/" \
            + os.path.basename(args.panel).replace(".bed", ".list"),
        'PANEL_LIST_NAME': os.path.basename(os.path.abspath(args.panel).replace(".bed", ".list")),
        'ROI_NUMBER'     : '.',
        'GENOME_VERSION' : args.reference,
        'VARIANT_CLASS'  : args.var_class,
        'CIVIC_VCF'      : defaults['CIVIC_FOLDER'] + "/" + "civic.variants.vcf",
        'THREADS'        : args.threads,
    }
    # Creating output directory
    output_path = Path(analysis_env['OUTPUT_DIR'])
    if not output_path.is_dir():
        os.mkdir(analysis_env['OUTPUT_DIR'])

    global logging
    log_file = analysis_env['OUTPUT_DIR'] + "/" + analysis_env['OUTPUT_NAME'] + ".pipeline.log"
    logging.basicConfig(filename=log_file, filemode='w', 
        format='PID:%(process)d\t%(asctime)s\t%(message)s')
    logging.getLogger().setLevel(logging.INFO)

    # Checking existance of input
    input_path = Path(analysis_env['INPUT_DIR'])
    if not input_path.is_dir():
        msg = "ERROR: Input folder " +  analysis_env['INPUT_DIR'] + " does not exist"
        print (msg)
        logging.error(msg)
        sys.exit()

def set_defaults(main_dir):
    ''' 
        Set default parameters
    '''
    global defaults 
    defaults = {
        # Folder with all binaries needed
        'BIN_FOLDER'    : main_dir +  "/BIN_FOLDER",
        # Folder with ancillary files
        'BUNDLE_FOLDER' : main_dir + "/BUNDLE_FOLDER",
        # Folder with ancillary files
        'PANEL_FOLDER'  : main_dir + "/PANEL_FOLDER",
        # Folder with annotation files
        'ANNOTATION_FOLDER' : main_dir + "/ANNOTATION_FOLDER",

        'CIVIC_FOLDER' : main_dir + "/ANNOTATION_FOLDER/CIVIC",

        # Setting VEP directories
        'SQLITE_DB_FOLDER' : main_dir + "/SQLITE_DB_FOLDER",
        'VEP_DATA'        :  main_dir + "/VEP_DATA",
        'VEP_DATA_INPUT'  :  main_dir + "/VEP_DATA/input",
        'VEP_DATA_OUTPUT' :  main_dir + "/VEP_DATA/output"
    }
    # Check all folders are present
    for folder in defaults:
        folder_path = Path(defaults[folder])
        if folder_path.is_dir():
            pass
        else:
            msg = " ERROR: folder " +  defaults[folder] + " is not available"
            print (msg)
            logging.error(msg)
            sys.exit()

def set_panel_configuration(main_dir):
    ''' 
        Set ancillary files (gnome fasta , etc)
    '''    
    global roi_env 
    roi_env = defaultdict(dict)

    roi_env = s.load_panel_transcripts(analysis_env['PANEL_NAME'])


def set_auxfiles_env():
    ''' 
        Set ancillary files (gnome fasta , etc)
    '''
    global aux_env 
    aux_env = defaultdict(dict)
    if analysis_env['GENOME_VERSION'] == 'hg19':
        aux_env['GENOME_FASTA']=  defaults['BUNDLE_FOLDER'] + "/ucsc.hg19.fasta"
        aux_env['GENOME_NAME'] = "ucsc.hg19.fasta"
        aux_env['GENOME_DICT'] = defaults['BUNDLE_FOLDER'] + "/ucsc.hg19.dict"
        aux_env['GENOME_DICT_NAME'] = "ucsc.hg19.dict"


def set_system_env():
    ''' 
        Set system binary parameters
    '''
    global system_env 
    system_env = defaultdict(dict)
    system_env = {
        'DOCKER' : '/usr/bin/docker',
        'SAMTOOLS' : defaults['BIN_FOLDER'] +  "/samtools",
        'BEDTOOLS' : defaults['BIN_FOLDER'] +  "/bedtools",
        'FREEBAYES': defaults['BIN_FOLDER'] +  "/freebayes",
        'FASTP'    : defaults['BIN_FOLDER'] +  "/fastp",
        'BWA'      : defaults['BIN_FOLDER'] +  "/bwa",
        'MOSDEPTH' : defaults['BIN_FOLDER'] +  "/mosdepth",
        'MANTA_CONFIG' : defaults['BIN_FOLDER'] + '/manta-1.6.0.centos6_x86_64/bin/configManta.py',
        'BGZIP'    : defaults['BIN_FOLDER'] +  "/bgzip",
        'TABIX'    : defaults['BIN_FOLDER'] +  "/tabix",
        'GUNZIP'   : defaults['BIN_FOLDER'] +  "/gunzip"
    }

    for binary in system_env:
        bin_path = Path(system_env[binary])
        if bin_path.is_file():
            msg = " INFO:" +  system_env[binary] + " is available"
            print (msg)
            logging.info(msg)
            pass
        else:
            msg = " ERROR: binary " +  system_env[binary] + " is not available"
            print (msg)
            logging.error(msg)
            sys.exit()

    global docker_env 
    docker_env = defaultdict(dict)

    gatkimage   = 'broadinstitute/gatk:4.1.3.0'
    picardimage = 'broadinstitute/picard'
    #vepimage    = 'ensembl/vep:latest

    docker_env = {
        'GATK'   : gatkimage,
        'PICARD' : picardimage
    }

    for image in docker_env:
        u.check_docker_images(system_env['DOCKER'], docker_env[image])

    # Now check all images are available
    # Now check that every binary can be executed