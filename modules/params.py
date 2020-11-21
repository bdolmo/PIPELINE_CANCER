#!/usr/bin/env python3

import os
import sys
import re
import logging
from collections import defaultdict
from pathlib import Path
import subprocess
import docx
import csv
from civicpy import civic, exports
from modules import utils as u
from modules import sqlite as s
from datetime import datetime
import pandas as pd

global sample_env 
sample_env = defaultdict(dict)


def set_analysis_env(args):


    global analysis_env 
    #analysis_env = defaultdict(dict)
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    analysis_env = {
        'INPUT_DIR'      : os.path.abspath(args.input_dir),
        'OUTPUT_DIR'     : os.path.abspath(args.output_dir),
        'OUTPUT_NAME'    : os.path.basename(args.output_dir),
        'PANEL'          : os.path.abspath(args.panel),
        'PANEL_NAME'     : os.path.basename(args.panel),
        'PANEL_LIST'     : defaults['PANEL_FOLDER'] + "/" \
            +  os.path.basename(args.panel) + "/" + os.path.basename(args.panel).replace(".bed", ".list"),
        'PANEL_LIST_NAME': os.path.basename(os.path.abspath(args.panel).replace(".bed", ".list")),
        'ROI_NUMBER'     : '.',
        'GENOME_VERSION' : args.reference,
        'VARIANT_CLASS'  : args.var_class,
        'THREADS'        : args.threads,
        'ANALYSIS_DATE'  : dt_string,
        'LANGUAGE'       : args.language,
        'MIN_FUSION_SIZE': args.min_fusion_size
    }
 
    if os.path.isfile(args.sample_data) :
        analysis_env['SAMPLE_DATA'] = os.path.abspath(args.sample_data)
    else:
        analysis_env['SAMPLE_DATA'] = "."
    if os.path.isfile(args.lab_data):
        analysis_env['LAB_DATA'] = os.path.abspath(args.lab_data)

    # Default annotation files for each genome version
    if analysis_env['VARIANT_CLASS'] == "somatic":
        if analysis_env['GENOME_VERSION'] == "hg19":
            analysis_env['CHIMERKB_BED'] = defaults['CHIMERKB_FOLDER'] + "/" \
                + "chimerKB_hg19_fusions.bed"
            analysis_env['CHIMERKB_BED_NAME'] = "chimerKB_hg19_fusions.bed"
            analysis_env['CGI_BIOMARKERS'] = defaults['CGI_FOLDER'] + "/" \
                + "cgi_biomarkers_latest" + "/" + "cgi_biomarkers_per_variant.tsv"
            analysis_env['CGI_BIOMARKERS_NAME'] = "cgi_biomarkers_per_variant.tsv"
            analysis_env['GNOMAD_AF_VCF'] = defaults['GNOMAD_FOLDER'] + "/" \
                + "somatic-b37_af-only-gnomad.raw.sites.vcf"
            analysis_env['GNOMAD_AF_VCF_NAME'] = "somatic-b37_af-only-gnomad.raw.sites.chr.vcf.gz"
    if analysis_env['LANGUAGE'] == "cat":
        defaults['JASPERREPORT_FOLDER'] = defaults['JASPERREPORT_FOLDER'] + "/cat"

    # Creating output directory
    output_path = Path(analysis_env['OUTPUT_DIR'])
    if not output_path.is_dir():
        os.mkdir(analysis_env['OUTPUT_DIR'])

    global logging
    log_file = os.path.abspath(analysis_env['OUTPUT_DIR']) + "/" + analysis_env['OUTPUT_NAME'] + ".pipeline.log"
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

    global sample_data
    sample_data = defaultdict(dict)

    # Load sample data; id's, tumor purity, etc
    if os.path.isfile(analysis_env['SAMPLE_DATA']):

        if not analysis_env['SAMPLE_DATA'].endswith(".docx"):
            msg = " ERROR: --sample_data " + analysis_env['SAMPLE_DATA'] + " is not a docx file"
            print(msg)
            logging.error(msg)
            sys.exit()

        doc = docx.Document(analysis_env['SAMPLE_DATA'])
        all_t = doc.tables
        is_date = False
        is_sample = False
        petition_date = ""
        for t in all_t:
            row_idx = 0
            for row in t.rows:
                cell_idx = 0
                for cell in row.cells:
                    cell.text = cell.text.rstrip("\n")
                    if cell.text == "":
                        continue                    
                    if cell.text.startswith('DATA'):
                        is_date = True
                       # break
                        continue
                    if cell.text.startswith('NOM'):
                        is_date = False
                        is_sample= True
                        break
                        continue
                    if is_date:
                        if row_idx == 0:
                            petition_date = cell.text
                            is_date = False
                    if is_sample:
                        if cell_idx == 0:
                            sample = cell.text
                            sample_data[sample]['PETITION_DATE'] = petition_date
                        if cell_idx == 1:
                            purity = re.search('^\d+%', cell.text)
                            if purity:
                                purity = purity.group(0)
                            else:
                                purity = 0
                            sample_data[sample]['PURITY'] = purity
                    cell_idx+=1
            row_idx += 1
    else:
        msg = " WARNING: sample data (.docx) file was not found"
        print (msg)
        logging.error(msg)

    global lab_data
    lab_data = defaultdict(dict)
    if os.path.isfile(analysis_env['LAB_DATA']):
        lab_df = pd.read_excel(analysis_env['LAB_DATA'], skiprows=10)
        for index,row in lab_df.iterrows():
            lab_id  = str(row['#SAMPLE IDIBGI ID'])
            if lab_id == "nan":
                continue
            ext1_id = str(row['SAMPLE FIC ID'])
            ext2_id = str(row['HC PACIENT'])

           # if lab_id in sample_data:             
            if ext1_id in sample_data:
               # sample_data[lab_id]= defaultdict(dict)
                lab_data[lab_id]['AP_CODE'] = ext1_id
                lab_data[lab_id]['HC_CODE'] = ext2_id
                lab_data[lab_id]['PURITY']  = sample_data[ext1_id]['PURITY']
                lab_data[lab_id]['PETITION_DATE'] = sample_data[ext1_id]['PETITION_DATE']
            else:
                lab_data[lab_id]['AP_CODE'] = ext1_id
                lab_data[lab_id]['HC_CODE'] = ext2_id
                lab_data[lab_id]['PURITY']  = '.'
                lab_data[lab_id]['PETITION_DATE'] = '.'
    for sample in lab_data:
        print ("mostra " + sample + " " + lab_data[sample]['AP_CODE'] + " " + lab_data[sample]['PURITY'])

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
        # chimerKB folder
        'CHIMERKB_FOLDER' : main_dir + "/ANNOTATION_FOLDER/chimerKB",
        # gnomAD folder
        'GNOMAD_FOLDER' : main_dir + "/ANNOTATION_FOLDER/gnomAD",
        # CGI folder
        'CGI_FOLDER' : main_dir + "/ANNOTATION_FOLDER/Cancer_Genome_Interpreter",
        # Setting JASPER directories
        'JASPERREPORT_FOLDER' : main_dir +  "/BIN_FOLDER/JASPERREPORTS/MyReports",
        'JDBC_FOLDER' : main_dir +  "/BIN_FOLDER/JASPERREPORTS/JDBC/",

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

def get_panel_configuration(main_dir):
    ''' 
        Set ancillary files (gnome fasta , etc)
    '''    
    global roi_env 
    roi_env = defaultdict(dict)
    roi_env = s.load_panel_transcripts(analysis_env['PANEL_NAME'])

    global biomarker_env
    biomarker_env = defaultdict(dict)
    biomarker_env = s.load_panel_biomarkers(analysis_env['PANEL_NAME'])

    global disclaimers_env
    disclaimers_env = defaultdict(dict)
    disclaimers_env = s.load_panel_disclaimers(analysis_env['PANEL_NAME'], analysis_env['LANGUAGE'])

    global cna_plot_env
    cna_plot_env = defaultdict(dict)
    cna_plot_env = s.load_cna(analysis_env['PANEL_NAME'])
    print(str(cna_plot_env))

def set_auxfiles_env():
    ''' 
        Set ancillary files (gnome fasta , etc)
    '''
    global aux_env 
    aux_env = defaultdict(dict)
    if analysis_env['GENOME_VERSION'] == 'hg19':
        defaults['BUNDLE_FOLDER'] = defaults['BUNDLE_FOLDER'] + "/hg19/"
        aux_env['GENOME_FASTA']     =  defaults['BUNDLE_FOLDER'] + "/ucsc.hg19.fasta"
        aux_env['GENOME_NAME']      = "ucsc.hg19.fasta"
        aux_env['GENOME_DICT']      = defaults['BUNDLE_FOLDER'] + "/ucsc.hg19.dict"
        aux_env['GENOME_DICT_NAME'] = "ucsc.hg19.dict"
        aux_env['GENE_LIST']        = defaults['BUNDLE_FOLDER'] + "/genelist.hg19.bed.gz"
        aux_env['GENE_LIST_NAME']   = "genelist.hg19.bed.gz"
        aux_env['NORMALS_REF_CNA']  = defaults['PANEL_FOLDER'] + "/" +  analysis_env['PANEL_NAME'].replace(".bed", "")\
          + "/" + analysis_env['PANEL_NAME'].replace(".bed", ".cnn")
        aux_env['NORMALS_REF_CNA_NAME']  = analysis_env['PANEL_NAME'].replace(".bed", "cnn")
    if analysis_env['LANGUAGE'] == "cat": 
        aux_env['REPORT_JRXML'] = defaults['BIN_FOLDER'] \
            +"/JASPERREPORTS/MyReports/cat/LungCancer_Report_v1_cat.jrxml"
    if analysis_env['LANGUAGE'] == "en": 
        aux_env['REPORT_JRXML'] = defaults['BIN_FOLDER'] \
            +"/JASPERREPORTS/MyReports/en/LungCancer_Report_v1_en.jrxml"

def set_system_env():
    ''' 
        Set system binary parameters
    '''
    global system_env 
    system_env = defaultdict(dict)
    system_env = {
        'DOCKER' : '/usr/bin/docker',
        'BCFTOOLS' : defaults['BIN_FOLDER'] +  "/bcftools",
        'SAMTOOLS' : defaults['BIN_FOLDER'] +  "/samtools",
        'BEDTOOLS' : defaults['BIN_FOLDER'] +  "/bedtools",
        'FREEBAYES': defaults['BIN_FOLDER'] +  "/freebayes",
        'FASTP'    : defaults['BIN_FOLDER'] +  "/fastp",
        'BWA'      : defaults['BIN_FOLDER'] +  "/bwa",
        'MOSDEPTH' : defaults['BIN_FOLDER'] +  "/mosdepth",
        'MANTA_CONFIG' : defaults['BIN_FOLDER'] + '/manta-1.6.0.centos6_x86_64/bin/configManta.py',
        'PURECN'   : defaults['BIN_FOLDER'] + "/PureCN/extdata/PureCN.R",
        'CNVKIT'   : defaults['BIN_FOLDER'] + "/cnvkit/cnvkit.py",
        'BGZIP'    : defaults['BIN_FOLDER'] +  "/bgzip",
        'TABIX'    : defaults['BIN_FOLDER'] +  "/tabix",
        'GUNZIP'   : defaults['BIN_FOLDER'] +  "/gunzip",
        'JASPERSTARTER' : defaults['BIN_FOLDER'] +"/JASPERREPORTS/jasperstarter/bin/jasperstarter"
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
    cnvkitimage = 'etal/cnvkit'
    #vepimage    = 'ensembl/vep:latest

    docker_env = {
        'GATK'   : gatkimage,
        'PICARD' : picardimage,
        'CNVKIT' : cnvkitimage
    }

    for image in docker_env:
        u.check_docker_images(system_env['DOCKER'], docker_env[image])

    # Now check all images are available
    # Now check that every binary can be executed