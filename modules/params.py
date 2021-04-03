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
from modules import utils as u
from modules import sqlite as s
from datetime import datetime
import pandas as pd
import openpyxl

global sample_env 
sample_env = defaultdict(dict)

def set_labdata_env():
    '''
    TEMPORARY, we need a petition manager to integrate all data from lab and external
    Lab data (docx) and sample_data (xlsx) loading.
    These two files are required to get code relationship (Lab, AP, HC codes)
    and others  like tumour purity estimates,
    '''
    # Sample data from petitions (either from DB or docx)  
    global sample_data
    sample_data = defaultdict(dict)

    # Load lab data, internal ids and relationship with external ids
    global lab_data
    lab_data = defaultdict(dict)

    if analysis_env['DB_DIR']:

        all_petitions = s.load_petitions()
        if all_petitions:
            msg = " INFO: Loading petitions from database"
            logging.info(msg)
            print(msg)
            for petition in all_petitions:
                print("AP CODE: " + petition.AP_code)
                sample_data[petition.AP_code] = defaultdict(dict)
                sample_data[petition.AP_code]['PETITION_ID'] = petition.Petition_id
                sample_data[petition.AP_code]['AP_CODE'] = petition.AP_code
                sample_data[petition.AP_code]['HC_CODE'] = petition.HC_code
                sample_data[petition.AP_code]['PURITY']  = petition.Tumour_pct
                sample_data[petition.AP_code]['PETITION_DATE'] = petition.Date
                sample_data[petition.AP_code]['VOLUME'] = petition.Volume
                sample_data[petition.AP_code]['CONC_NANODROP'] = petition.Conc_nanodrop
                sample_data[petition.AP_code]['RATIO_NANODROP'] = petition.Ratio_nanodrop
                sample_data[petition.AP_code]['MEDICAL_DOCTOR'] = petition.Medical_doctor
                sample_data[petition.AP_code]['BILLING_UNIT'] = petition.Billing_unit
        else:
            msg = " INFO: No petitions were found"
            logging.info(msg)
            print(msg)

    elif 'SAMPLE_DATA' in analysis_env:

        doc = docx.Document(analysis_env['SAMPLE_DATA'])

        # Getting all tables  
        all_t = doc.tables
        is_date   = False
        is_sample = False
        
        petition_date  = ""
        petition_id    = "" 
        ap_code        = ""
        purity         = ""
        residual_volume= ""
        nanodrop_conc  = ""
        nanodrop_ratio = ""
        hc_number      = ""
        medical_doctor = ""
        billing_unit   = ""

        abs_idx = 0
        for t in all_t:
            row_idx = 0
            for row in t.rows:
                abs_idx+=1
                cell_idx = 0
                if is_sample:
                    if len(row.cells) > 9:
                        msg = " ERROR: row " + str(abs_idx) + " has more thant 8 cells"
                        print (msg)
                        logging.error(msg)
                        break
                for cell in row.cells:
                    cell.text = cell.text.rstrip("\n")
                    if cell.text == "":
                        continue                    
                    if cell.text.startswith('DATA'):
                        is_date = True
                        continue
                    if cell.text.startswith('CODI'):
                        is_date = False
                        is_sample= True
                    if is_date:
                        if row_idx == 0:
                            petition_date = cell.text
                            tmp_date = petition_date.split("/")
                            if len(tmp_date) != 3:
                                msg = " ERROR: Date format is incorrect. It should follow dd/mm/yyyy format"
                                print (msg)
                                logging.error(msg)
                            if len(tmp_date) == 3:
                                days   = int(tmp_date[0])
                                month  = int(tmp_date[1])
                                year   = int(tmp_date[2])
                                if days > 31 or month > 12:
                                    msg = " ERROR: Date format is incorrect. It should follow dd/mm/yyyy format"
                                    print (msg)
                                    logging.error(msg)
                            is_date = False
                    if is_sample:
                        if row_idx > 1:
                            # AP code, must be alphanumeric
                            if cell_idx == 0:
                                ap_code = cell.text
                                
                                # Ensuring AP specific format  
                                if re.search('\d+[A-Z]+\d+?\s+[A-Z]\d+', ap_code):
                                    pass
                                else:
                                    msg = " ERROR: Incorrect AP code at row " + str(abs_idx)
                                    print (msg)
                                    logging.error(msg)
                                if not ap_code:
                                    msg = " ERROR: Missing AP code at row " + str(abs_idx)
                                    print (msg)
                                    logging.error(msg)                                
                            # Tumour purity (%) 
                            if cell_idx == 1:
                                purity = re.search('^\d+%', cell.text)
                                if purity:
                                    purity = purity.group(0).replace("%", "")
                                    if re.search('\d+', purity):
                                        if int(purity) >= 0 and int(purity) <= 100:
                                            pass
                                        else:
                                            msg = " ERROR: Incorrect %Tumoral value for sample " + ap_code
                                            print (msg)
                                            logging.error(msg)                                
                                    else:
                                        msg = " ERROR: Non-numeric %Tumoral value was found for sample " + ap_code
                                        print (msg)
                                        logging.error(msg)
                                else:
                                    msg = " ERROR: Missing %Tumoral value for sample " + ap_code
                                    print (msg)
                                    logging.error(msg)                                    

                            # Residual volume 
                            if cell_idx == 2:
                                residual_volume = cell.text
                                if re.search('\d+', purity):
                                    pass
                                else:
                                    msg = " ERROR: Non-numeric residual value for was found for sample " + ap_code
                                    print (msg)
                                    logging.error(msg)
                                if not residual_volume:
                                    msg = " ERROR: Missing residual volume for sample " + ap_code
                                    print (msg)
                                    logging.error(msg)

                            # Nanodrop concentration 
                            if cell_idx == 3:
                                nanodrop_conc = cell.text.replace(",", ".")
                                if re.search('\d+?\.?\d+', nanodrop_conc):
                                    pass
                                else:
                                    msg = " ERROR: Invalid Nanodrop concentration value for sample " + ap_code
                                    print(msg)
                                    logging.error(msg)
                                if not nanodrop_conc:
                                    msg = " ERROR: Missing Nanodrop concentration for sample " + ap_code
                                    print (msg)
                                    logging.error(msg)

                            # Nanodrop ratio 
                            if cell_idx == 4:
                                nanodrop_ratio = cell.text.replace(",", ".")
                                if re.search('\d+?\.?\d+', nanodrop_ratio):
                                    pass
                                else:
                                    msg = " ERROR: Invalid value for Nanodrop ratio 260/280 for sample " + ap_code
                                    print(msg)
                                    logging.error(msg)
                                if not nanodrop_ratio:
                                    msg = " ERROR: Missing Nanodrop ratio 260/280 value for sample " + ap_code
                                    print(msg)
                                    logging.error(msg)                           
                            # HC number
                            if cell_idx == 5:
                                hc_number = cell.text
                                if not hc_number:
                                    msg = " ERROR: Missing HC (clinic history) number for sample " + ap_code
                                    print(msg)
                                    logging.error(msg) 
                            # Medical doctor solicitant 
                            if cell_idx == 6:
                                medical_doctor = cell.text
                                if not medical_doctor:
                                    msg = " ERROR: Missing Physician solicitant for sample " + ap_code
                                    print(msg)
                                    logging.error(msg)                                     
                            # Billing unit 
                            if cell_idx == 7:
                                billing_unit = cell.text
                                if not billing_unit:
                                    msg = " ERROR: Missing billing unit for sample " + ap_code
                                    print(msg)
                                    logging.error(msg)                                    
                            sample_dict[ap_code] = defaultdict(dict)
                            sample_dict[ap_code]['Date']           = petition_date
                            sample_dict[ap_code]['Petition_id']    = '.'  
                            sample_dict[ap_code]['AP_code']        = ap_code  
                            sample_dict[ap_code]['HC_number']      = hc_number  
                            sample_dict[ap_code]['Tumour_purity']  = purity  
                            sample_dict[ap_code]['Residual_volume']= residual_volume  
                            sample_dict[ap_code]['Nanodrop_conc']  = nanodrop_conc  
                            sample_dict[ap_code]['Nanodrop_ratio'] = nanodrop_ratio  
                            sample_dict[ap_code]['Medical_doctor'] = medical_doctor  
                            sample_dict[ap_code]['Billing_unit']   = billing_unit  

                    cell_idx+=1
                row_idx += 1

    if 'LAB_DATA' in analysis_env:

        msg = " INFO: Loading lab data from " + analysis_env['LAB_DATA']
        logging.info(msg)
        print(msg)

        lab_df = pd.read_excel(analysis_env['LAB_DATA'],engine='openpyxl',skiprows=10)
        for index,row in lab_df.iterrows():
            lab_id  = str(row['#SAMPLE IDIBGI ID'])
            if lab_id == "nan":
                continue
            ext1_id = str(row['SAMPLE FIC ID'])
            ext2_id = str(row['HC PACIENT'])
            if ext1_id in sample_data:
                lab_data[lab_id]['AP_CODE'] = ext1_id

                lab_data[lab_id]['HC_CODE'] = ext2_id
                lab_data[lab_id]['PURITY']  = sample_data[ext1_id]['PURITY']
                lab_data[lab_id]['PETITION_DATE'] = sample_data[ext1_id]['PETITION_DATE']
            else:
                lab_data[lab_id]['AP_CODE'] = ext1_id
                lab_data[lab_id]['HC_CODE'] = ext2_id
                lab_data[lab_id]['PURITY']  = '.'
                lab_data[lab_id]['PETITION_DATE'] = '.'
            if lab_data[lab_id]['AP_CODE'] == lab_id:
                lab_data[lab_id]['AP_CODE'] = '.'
    else:
        msg = " ERROR: lab data (.xlsx) file was not found"
        print (msg)
        logging.error(msg)

def set_analysis_env(args, mode):

    global analysis_env 
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    analysis_env = {
        'ANALYSIS_MODE'  : mode,
        'INPUT_DIR'      : os.path.abspath(args.input_dir),
        'OUTPUT_DIR'     : os.path.abspath(args.output_dir),
        'ROI_NUMBER'     : '.',
        'SEQUENCING'     : args.sequencing,
        'GENOME_VERSION' : args.reference,
        'VARIANT_CLASS'  : args.var_class,
        'THREADS'        : args.threads,
        'ANALYSIS_DATE'  : dt_string,
        'DB_DIR'         : os.path.abspath(args.db_dir),
        'ANN_DIR'        : os.path.abspath(args.ann_dir),
        'REF_DIR'        : os.path.abspath(args.ref_dir),
        'LANGUAGE'       : args.language,
        'MIN_FUSION_SIZE': args.min_fusion_size
    }

    # Creating output directory
    output_path = Path(analysis_env['OUTPUT_DIR'])
    if not output_path.is_dir():
        os.mkdir(analysis_env['OUTPUT_DIR'])

    # Setting output name (just the basename of the output dir) 
    analysis_env['OUTPUT_NAME'] = os.path.basename(analysis_env['OUTPUT_DIR'])

    # Define log file for the analysis
    global logging
    log_file = os.path.abspath(analysis_env['OUTPUT_DIR']) + "/" \
        + analysis_env['OUTPUT_NAME'] + ".pipeline.log"
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

    # Setting targeted (panel, exome) analysis params     
    if args.sequencing == "targeted":
        analysis_env['SEQ_APPLICATION'] = "targeted"
        if args.panel:
            tmp_panel = os.path.basename(args.panel).split("\.")
            panel_name = tmp_panel[0]
            analysis_env['PANEL'] = args.panel
            analysis_env['PANEL_NAME'] = panel_name
            analysis_env['PANEL_WORKDIR'] =  defaults['PANEL_FOLDER'] +  "/" \
                + os.path.basename(args.panel).replace(".bed", "")
            if not os.path.isdir(analysis_env['PANEL_WORKDIR']):
                os.mkdir(analysis_env['PANEL_WORKDIR'])
            analysis_env['PANEL_LIST'] = analysis_env['PANEL_WORKDIR'] + "/" \
                + os.path.basename(args.panel).replace(".bed", ".list")
            analysis_env['PANEL_LIST_NAME']  = \
                os.path.basename(analysis_env['PANEL_LIST'].replace(".bed", ".list"))
        else:
            msg = "ERROR: It is required a gene panel (--panel) for targeted sequencing analysis"
            print (msg)
            logging.error(msg)
            sys.exit()
    else:
        analysis_env['SEQ_APPLICATION'] = "wgs"
        
    # Define user_id if you want to link with database information 
    if args.user_id:
        analysis_env['USER_ID'] = args.user_id
    else:
        analysis_env['USER_ID'] = "."

    # Load sample data (docx), just for somatic! 
    if args.sample_data:
        if os.path.isfile(args.sample_data) :
            analysis_env['SAMPLE_DATA'] = os.path.abspath(args.sample_data)
        else:
            analysis_env['SAMPLE_DATA'] = "."

    # Load lab sample data (from xlsx file) 
    if args.lab_data:
        if os.path.isfile(args.lab_data):
            analysis_env['LAB_DATA'] = os.path.abspath(args.lab_data)

    # Check existance of the database directory 
    if not os.path.isdir(analysis_env['DB_DIR'] ):
        msg = " ERROR: Missing input database directory (--db_dir)"
        print (msg)
        logging.info(msg)
        sys.exit()

    # Check existance of the annotation directory 
    if not os.path.isdir(analysis_env['ANN_DIR'] ):
        msg = " ERROR: Missing input annotation directory (--ann_dir)"
        print (msg)
        logging.info(msg)
        sys.exit()
        
    # Check existance of the reference directory 
    if not os.path.isdir(analysis_env['REF_DIR'] ):
        msg = " ERROR: Missing input reference directory (--ref_dir)"
        print (msg)
        logging.info(msg)
        sys.exit()

    # Check out input missense predictors
    try:
        getattr(args, missense_predictors)
    except:
        default_predictors = ['sift','polyphen2','mutationtaster2','provean','fathmm','revel','mutpred']
        analysis_env['MISSENSE_PREDICTORS'] = ','.join(default_predictors)
    else:
        valid_predictors = ['sift','polyphen2','mutationtaster2','provean','fathmm','revel','mutpred']  
        input_list = args.missense_predictors.split(",")
        for predictor in input_list:
            if predictor not in valid_predictors:
                msg = " ERROR: Input missense predictor " + predictor + " is not valid"
                print (msg)
                logging.info(msg)
                sys.exit()
        analysis_env['MISSENSE_PREDICTORS'] = args.missense_predictors

    # Check out input missense predictors
    try:
        getattr(args, genomewide_predictors)
    except:
        default_predictors = ['cadd','ncer']  
        analysis_env['GENOMEWIDE_PREDICTORS'] = ','.join(default_predictors)
    else:
        valid_predictors = ['cadd','ncer']  
        input_list = args.genomewide_predictors.split(",")
        for predictor in input_list:
            if predictor not in valid_predictors:
                msg = " ERROR: Input genomewide predictor " + predictor + " is not valid"
                print (msg)
                logging.info(msg)
                sys.exit()
        analysis_env['GENOMEWIDE_PREDICTORS'] = args.genomewide_predictors

    # Check out input missense predictors
    try:
        getattr(args, splicing_predictors)
    except:
        default_predictors = ['spliceai','maxentscan','dbscsnv']  
        analysis_env['SPLICING_PREDICTORS'] = ','.join(default_predictors)
    else:
        valid_predictors = ['spliceai','maxentscan','dbscsnv']  
        input_list = args.splicing_predictors.split(",")
        for predictor in input_list:
            if predictor not in valid_predictors:
                msg = " ERROR: Input splicing predictor " + predictor + " is not valid"
                print (msg)
                logging.info(msg)
                sys.exit()
        analysis_env['SPLICING_PREDICTORS'] = args.splicing_predictors

    # Check out conservation scores
    try:
        getattr(args, conservation)
    except:
        default_conservation = ['phastcons', 'phylop', 'gerp'] 
        analysis_env['CONSERVATION_SCORES'] = ','.join(default_conservation)
    else:
        valid_conservation = ['phastcons', 'phylop', 'gerp'] 
        input_list = args.conservation.split(",")
        for predictor in input_list:
            if predictor not in valid_conservation:
                msg = " ERROR: Input conservation score " + predictor + " is not valid"
                print (msg)
                logging.info(msg)
                sys.exit()
        analysis_env['CONSERVATION_SCORES'] = args.conservation

    # Check out gnomAD annotation
    try:
        getattr(args, gnomad)
    except:
        # Default 
        analysis_env['GNOMAD'] = True
    else:
        if args.gnomad == True:
            analysis_env['GNOMAD'] = True
        else:
            analysis_env['GNOMAD'] = False

    # Check out 1KG Annotation
    try:
        getattr(args, thousand_genomes)
    except:
        analysis_env['1KG'] = True
    else:
        if args.thousand_genomes == True:
            analysis_env['1KG'] = True
        else:
            analysis_env['1KG'] = False

def set_defaults(main_dir):
    ''' 
        Set default parameters
    '''
    global defaults 
    defaults = {
        # Folder with all binaries needed
        'BIN_FOLDER'    : main_dir +  "/BIN_FOLDER",
        # Folder with ancillary files
        'PANEL_FOLDER'  : main_dir + "/PANEL_FOLDER", 
        # Setting JASPER directories
        'JASPERREPORT_FOLDER' : main_dir +  "/BIN_FOLDER/JASPERREPORTS/MyReports",
        'JDBC_FOLDER' : main_dir +  "/BIN_FOLDER/JASPERREPORTS/JDBC/",
        #'SQLITE_DB_FOLDER' : main_dir + "/SQLITE_DB_FOLDER",
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
        Set panel configuration if available
    '''

    # Loading panel defined transcript id's    
    global roi_env 
    roi_env = defaultdict(dict)
    roi_env = s.load_panel_transcripts(analysis_env['PANEL_NAME'])
    if not roi_env:
        msg = " WARNING: Unable to load transcript info for panel " + analysis_env['PANEL_NAME']
        print(msg)
        logging.error(msg)
    else:
        msg = " INFO: loading transcript info for panel " + analysis_env['PANEL_NAME'] + " ended successfully"
        print(msg)
        logging.info(msg)

    # Panel version 
    panel_info = s.Roi.query.filter_by(panel=analysis_env['PANEL_NAME']).first()
    if panel_info:
        analysis_env['PANEL_VERSION'] = panel_info.panel_version
    else: 
        analysis_env['PANEL_VERSION'] = '.'

    # Loading panel defined biomarkers
    global biomarker_env
    biomarker_env = defaultdict(dict)
    biomarker_env = s.load_panel_biomarkers(analysis_env['PANEL_NAME'])
    if not roi_env:
        msg = " WARNING: Unable to load biomarker info for panel " + analysis_env['PANEL_NAME']
        print(msg)
        logging.error(msg)
    else:
        msg = " INFO: loading biomarker info for panel " + analysis_env['PANEL_NAME'] + " ended successfully"
        print(msg)
        logging.info(msg)

    # Loading panel disclaimers
    global disclaimers_env
    disclaimers_env = defaultdict(dict)
    disclaimers_env = s.load_panel_disclaimers(analysis_env['PANEL_NAME'], analysis_env['LANGUAGE'])
    if not roi_env:
        msg = " WARNING: Unable to load disclaimers info for panel " + analysis_env['PANEL_NAME']
        print(msg)
        logging.error(msg)

    # Loading panel defined genes for CNA analysis
    global cna_env
    cna_env = defaultdict(dict)
    cna_env = s.load_cna(analysis_env['PANEL_NAME'])
    if not roi_env:
        msg = " WARNING: Unable to load CNA genes for panel " + analysis_env['PANEL_NAME']
        print(msg)
        logging.error(msg)

def set_auxfiles_env():
    ''' 
        Set ancillary files (genome fasta , etc). Better with yaml file => Coming soon!
    '''
    global aux_env 
    aux_env = defaultdict(dict)
    if analysis_env['GENOME_VERSION'] == 'hg19':

        # Setting genome folder and checking that it exists 
        aux_env['GENOME_FOLDER']    = analysis_env['REF_DIR'] + "/hg19"
        if not os.path.isdir(aux_env['GENOME_FOLDER'] ):
            msg = " ERROR: Missing genome folder " + aux_env['GENOME_FOLDER']
            print (msg)
            logging.error(msg)
            sys.exit()

        # Setting genome fasta file and checking that it exists 
        aux_env['GENOME_FASTA']     = aux_env['GENOME_FOLDER'] + "/ucsc.hg19.fasta"
        if not os.path.isfile(aux_env['GENOME_FASTA'] ):
            msg = " ERROR: Missing FASTA genome " + aux_env['GENOME_FASTA']
            print (msg)
            logging.error(msg)
            sys.exit()
        aux_env['GENOME_NAME']      = os.path.basename(aux_env['GENOME_FASTA'])

        # Setting genome fasta dictionary and checking that it exists 
        aux_env['GENOME_DICT']      = aux_env['GENOME_FOLDER'] + "/ucsc.hg19.dict"
        if not os.path.isfile(aux_env['GENOME_DICT'] ):
            msg = " ERROR: Missing FASTA dict " + aux_env['GENOME_DICT']
            print (msg)
            logging.error(msg)
            sys.exit()
        aux_env['GENOME_DICT_NAME'] = os.path.basename(aux_env['GENOME_DICT'])

        # Setting gene list and checking that it exists 
        aux_env['GENE_LIST']        = aux_env['GENOME_FOLDER'] + "/genelist.hg19.bed.gz"
        if not os.path.isfile(aux_env['GENE_LIST'] ):
            msg = " ERROR: Missing gene list bed file " + aux_env['GENE_LIST']
            print (msg)
            logging.error(msg)
            sys.exit()
        aux_env['GENE_LIST_NAME']   = "genelist.hg19.bed.gz"

        if analysis_env['SEQ_APPLICATION'] == "targeted": 
            # Setting normals reference and checking that it exists 
            aux_env['NORMALS_REF_CNA']  = defaults['PANEL_FOLDER'] + "/" \
                +  analysis_env['PANEL_NAME'].replace(".bed", "")\
                + "/" + analysis_env['PANEL_NAME'].replace(".bed", ".cnn")
            if not os.path.isfile(aux_env['NORMALS_REF_CNA'] ):
                msg = " ERROR: Missing normals reference (.cnn) " + aux_env['NORMALS_REF_CNA']
                print (msg)
                logging.error(msg)
                sys.exit()
            aux_env['NORMALS_REF_CNA_NAME']  = os.path.basename(aux_env['NORMALS_REF_CNA'])

        # Annotation files
        if analysis_env['VARIANT_CLASS'] == "somatic":

            # Setting VEP and checking that it exists 
            aux_env['VEP_FOLDER'] = analysis_env['ANN_DIR'] + "/VEP"
            if not os.path.isdir(aux_env['VEP_FOLDER']):
                msg = " ERROR: Missing VEP directory " + aux_env['VEP_FOLDER']
                print (msg)
                logging.error(msg)
                sys.exit()
            
            # Where input vcfs will be placed
            aux_env['VEP_FOLDER_INPUT'] = aux_env['VEP_FOLDER'] + "/input"
            if not os.path.isdir(aux_env['VEP_FOLDER_INPUT']):
                msg = " ERROR: Missing VEP input directory " + aux_env['VEP_FOLDER_INPUT']
                print (msg)
                logging.error(msg)
                sys.exit()

            # Where output (annotated) vcfs will be placed
            aux_env['VEP_FOLDER_OUTPUT'] = aux_env['VEP_FOLDER'] + "/output"
            if not os.path.isdir(aux_env['VEP_FOLDER_OUTPUT']):
                msg = " ERROR: Missing VEP output directory " + aux_env['VEP_FOLDER_OUTPUT']
                print (msg)
                logging.error(msg)
                sys.exit()

            # Setting chimerKB and checking that everything exists 
            aux_env['CHIMERKB_FOLDER'] = analysis_env['ANN_DIR'] + "/chimerKB"
            if not os.path.isdir(aux_env['CHIMERKB_FOLDER']):
                msg = " ERROR: Missing chimerKB directory " + aux_env['CHIMERKB_FOLDER']
                print (msg)
                logging.error(msg)
                sys.exit()
            aux_env['CHIMERKB_BED'] = aux_env['CHIMERKB_FOLDER'] + "/" \
                + "chimerKB_hg19_fusions.bed"
            if not os.path.isfile(aux_env['CHIMERKB_BED']):
                msg = " ERROR: Missing chimerKB directory " + aux_env['CHIMERKB_FOLDER']
                print (msg)
                logging.error(msg)
                sys.exit()
            aux_env['CHIMERKB_BED_NAME'] = os.path.basename(aux_env['CHIMERKB_BED'])

            # Setting CGI and checking that everything exists 
            aux_env['CGI_FOLDER'] = analysis_env['ANN_DIR'] + "/CGI"
            if not os.path.isdir(aux_env['CGI_FOLDER']):
                msg = " ERROR: Missing CGI directory " + aux_env['CGI_FOLDER']
                print (msg)
                logging.error(msg)
                sys.exit()
            aux_env['CGI_BIOMARKERS'] = aux_env['CGI_FOLDER'] + "/" \
                + "cgi_biomarkers_latest" + "/" + "cgi_biomarkers_per_variant.tsv"
            if not os.path.isfile(aux_env['CGI_BIOMARKERS']):
                msg = " ERROR: Missing CGI biomarker file " + aux_env['CGI_BIOMARKERS']
                print (msg)
                logging.error(msg)
                sys.exit()                
            aux_env['CGI_BIOMARKERS_NAME'] = os.path.basename(aux_env['CGI_BIOMARKERS'])

            # Setting CGI and checking that everything exists 
            aux_env['GNOMAD_FOLDER'] = analysis_env['ANN_DIR'] + "/gnomAD"
            if not os.path.isdir(aux_env['GNOMAD_FOLDER']):
                msg = " ERROR: Missing gnomAD directory: " + aux_env['GNOMAD_FOLDER']
                print (msg)
                logging.error(msg)
                sys.exit()

            # Setting gnomAD and checking that everything exists 
            aux_env['GNOMAD_AF_VCF'] = aux_env['GNOMAD_FOLDER'] + "/" \
                + "somatic-b37_af-only-gnomad.raw.sites.chr.vcf.gz"
            if not os.path.isfile(aux_env['GNOMAD_AF_VCF']):
                msg = " ERROR: Missing gnomAD VCF: " + aux_env['GNOMAD_AF_VCF']
                print (msg)
                logging.error(msg)
                sys.exit()            
            aux_env['GNOMAD_AF_VCF_NAME'] = os.path.basename(aux_env['GNOMAD_AF_VCF'])

            # Setting gnomAD and checking that everything exists 
            aux_env['GNOMAD'] = aux_env['GNOMAD_FOLDER'] + "/gnomad.genomes.r2.1.1.sites.only_af.vcf.gz"
            if not os.path.isfile(aux_env['GNOMAD']):
                msg = " ERROR: Missing gnomAD: " + aux_env['GNOMAD']
                print (msg)
                logging.error(msg)
                sys.exit()            
            aux_env['GNOMAD_FILENAME'] = os.path.basename(aux_env['GNOMAD'])
            # Setting 1KG
            aux_env['1KG_FOLDER'] =  analysis_env['ANN_DIR'] + "/1000Genomes"
            if not os.path.isdir(aux_env['1KG_FOLDER']):
                msg = " ERROR: Missing 1000Genomes folder at " + analysis_env['ANN_DIR']
                print (msg)
                logging.error(msg)
                sys.exit()
            aux_env['1KG'] = aux_env['1KG_FOLDER'] + "/" \
                + "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz"
            aux_env['1KG_FILENAME'] = os.path.basename(aux_env['1KG'])

            # Setting dbNSFP
            aux_env['DBNSFP_FOLDER'] = analysis_env['ANN_DIR'] + "/dbNSFP"
            aux_env['DBNSFP'] = aux_env['DBNSFP_FOLDER'] + "/dbNSFP4.1a_grch37.gz"
            aux_env['DBNSFP_FILENAME'] = os.path.basename(aux_env['DBNSFP'])
            if not os.path.isfile(aux_env['DBNSFP']):
                msg = " ERROR: Missing " + aux_env['DBNSFP']
                print (msg)
                logging.error(msg)
                sys.exit()  

            # Setting SpliceAI
            aux_env['SPLICEAI_FOLDER'] =  analysis_env['ANN_DIR'] + "/spliceAI"
            if not os.path.isdir(aux_env['SPLICEAI_FOLDER']):
                msg = " ERROR: Missing SpliceAI folder at " + analysis_env['ANN_DIR']
                print (msg)
                logging.error(msg)
                sys.exit()
            aux_env['SPLICEAI_SNV']  =  aux_env['SPLICEAI_FOLDER'] + "/spliceai_scores.raw.snv.hg19.vcf.gz"
            if not os.path.isfile(aux_env['SPLICEAI_SNV']):
                msg = " ERROR: Missing " + aux_env['SPLICEAI_SNV']
                print (msg)
                logging.error(msg)
                sys.exit()  
            aux_env['SPLICEAI_SNV_FILENAME'] = os.path.basename(aux_env['SPLICEAI_SNV'])
        
            aux_env['SPLICEAI_INDEL']=  aux_env['SPLICEAI_FOLDER'] + "/spliceai_scores.raw.indel.hg19.vcf.gz"
            if not os.path.isfile(aux_env['SPLICEAI_INDEL']):
                msg = " ERROR: Missing " + aux_env['SPLICEAI_INDEL']
                print (msg)
                logging.error(msg)
                sys.exit()
            aux_env['SPLICEAI_INDEL_FILENAME'] = os.path.basename(aux_env['SPLICEAI_INDEL'])

            # Setting MaxEntScan
            aux_env['MAXENT_FOLDER'] =  analysis_env['ANN_DIR'] + "/MaxEntScan"
            if not os.path.isdir(aux_env['MAXENT_FOLDER']):
                msg = " ERROR: Missing MaxEntScan folder at " + analysis_env['ANN_DIR']
                print (msg)
                logging.error(msg)
                sys.exit()
            # Setting Conservation
            aux_env['CONSERVATION_FOLDER'] =  analysis_env['ANN_DIR'] + "/Conservation"
            if not os.path.isdir(aux_env['CONSERVATION_FOLDER']):
                msg = " ERROR: Missing Conservation folder at " + analysis_env['ANN_DIR']
                print (msg)
                logging.error(msg)
                sys.exit()
            aux_env['PHASTCONS'] = aux_env['CONSERVATION_FOLDER'] + "/" + "hg19.100way.phastCons.bw"
            aux_env['PHASTCONS_FILENAME'] = os.path.basename(aux_env['PHASTCONS'])
            aux_env['PHYLOP'] = aux_env['CONSERVATION_FOLDER'] + "/" + "hg19.100way.phyloP100way.bw"
            aux_env['PHYLOP_FILENAME'] = os.path.basename(aux_env['PHYLOP'])

    # Setting jrxlm files depending on language
    if analysis_env['LANGUAGE'] == "cat":
        defaults['JASPERREPORT_FOLDER'] = defaults['JASPERREPORT_FOLDER'] + "/cat/"
        aux_env['REPORT_JRXML'] = defaults['BIN_FOLDER'] \
            +"/JASPERREPORTS/MyReports/cat/LungCancer_Report_v1_cat.jrxml"
    if analysis_env['LANGUAGE'] == "en":
        defaults['JASPERREPORT_FOLDER'] = defaults['JASPERREPORT_FOLDER'] + "/en/"
        aux_env['REPORT_JRXML'] = defaults['BIN_FOLDER'] \
            +"/JASPERREPORTS/MyReports/en/LungCancer_Report_v1_en.jrxml"
    if not os.path.isfile(aux_env['REPORT_JRXML'] ):
        msg = " WARNING: Missing jrxml file " + aux_env['REPORT_JRXML']
        print (msg)
        logging.error(msg)

def set_system_env():
    ''' 
        Set system binaries and other NGS utils. 
        would be better to be renamed ngs_utils_env?
    '''
    global system_env 
    system_env = defaultdict(dict)
    system_env = {
        'DOCKER'   : u.get_bin_path('docker'),
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

    # Now setting docker environment 
    global docker_env 
    docker_env = defaultdict(dict)

    # Images 
    gatkimage   = 'broadinstitute/gatk:4.1.3.0'
    picardimage = 'broadinstitute/picard'
    cnvkitimage = 'etal/cnvkit'
    vepimage    = 'ensemblorg/ensembl-vep'

    docker_env = {
        'GATK'   : gatkimage,
        'PICARD' : picardimage,
        'CNVKIT' : cnvkitimage,
        'VEP'    : vepimage
    }
    # Checking images are available 
    for image in docker_env:
        u.check_docker_images(system_env['DOCKER'], docker_env[image])