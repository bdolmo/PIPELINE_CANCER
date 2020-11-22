#!/usr/bin/env python3

import os
import sys
import re
import logging
import gzip
import csv
from collections import defaultdict
from pathlib import Path
from scipy import stats
import numpy as np
import subprocess
import shutil
from modules import utils as u
from modules import params as p
from modules import trimming as t

def do_all():

    do_generate_list_file()

    # Mapping phase
    map_fastq()

    # Removing duplicates
    remove_duplicates()

    # Get coverage metrics
    extract_coverage_metrics()

    # Get mapping metrics
    extract_mapping_metrics()

    # Gather all metrics and create a summary qc
    create_summary_qc()


def create_summary_qc():

    msg = " INFO: Creating summary QC"
    print (msg)
    logging.info(msg)

    p.analysis_env['SUMMARY_QC'] = p.analysis_env['OUTPUT_DIR'] + \
        "/" + "summary_qc.tsv"

    if os.path.isfile(p.analysis_env['SUMMARY_QC']):
        msg = " INFO: Skipping summary QC joining"
        print (msg)
        logging.info(msg)
        return

    header = ['Fic ID','Lab ID','Read L','Panel','N Ex','Total Reads','20X','30X','100X','20X','30X','100X','Mean_cov','%ROI','%Rmdup']

    with open(p.analysis_env['SUMMARY_QC'], 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        tsv_writer.writerow(header)
        for sample in p.sample_env:
            info = []
            info.append("9999999")
            info.append(sample)
            info.append("76")
            info.append(p.analysis_env['PANEL_NAME'])
            info.append(str(p.analysis_env['ROI_NUMBER']))
            info.append(str(u.num_to_human(p.sample_env[sample]['TOTAL_READS'])))
            for threshold in p.sample_env[sample]['CALL_RATE']:
                if threshold in header:
                    info.append(str(p.sample_env[sample]['CALL_RATE'][threshold]))
            for threshold in p.sample_env[sample]['LOST_EXONS']:
                if threshold in header:
                    info.append(str(p.sample_env[sample]['LOST_EXONS'][threshold]))
            info.append(str(p.sample_env[sample]['MEAN_COVERAGE']))
            info.append(str(p.sample_env[sample]['ROI_PERCENTAGE']))
            info.append(str(p.sample_env[sample]['PCR_DUPLICATES_PERCENTAGE']))
           
            tsv_writer.writerow(info)

def extract_mapping_metrics():

    p.analysis_env['SUMMARY_QC'] = p.analysis_env['OUTPUT_DIR'] + \
        "/" + "summary_qc.tsv"

    if os.path.isfile(p.analysis_env['SUMMARY_QC']):
        msg = " INFO: Skipping mapping metric extraction"
        print (msg)
        logging.info(msg)
        with open(p.analysis_env['SUMMARY_QC']) as f:
            for line in f:
                line = line.rstrip('\n')
                tmp = line.split('\t')
                sample = tmp[1] 
                if 'Lab ID' in line:
                    continue
                else:
                    p.sample_env[sample]['TOTAL_READS'] = tmp[5] 
                    p.sample_env[sample]['MEAN_COVERAGE'] = tmp[12] 
                    p.sample_env[sample]['ROI_PERCENTAGE'] = tmp[13] 
                    p.sample_env[sample]['PCR_DUPLICATES_PERCENTAGE'] = tmp[14] 
        return
        
    for sample in p.sample_env:

        p.sample_env[sample]['TOTAL_READS']        = "."
        p.sample_env[sample]['ON_TARGET_READS']    = "."
        p.sample_env[sample]['ROI_PERCENTAGE']     = "."
        p.sample_env[sample]['MEAN_INSERT_SIZE']   = "."
        p.sample_env[sample]['SD_INSERT_SIZE']     = "."
  
        # Getting total number of reads
        bashCommand = ('{} view -c {}').format(p.system_env['SAMTOOLS'], \
        p.sample_env[sample]['READY_BAM'])
        msg = " INFO: Extracting mapping metrics for sample "+ sample
        print (msg)
        logging.info(bashCommand)

        process = subprocess.Popen(bashCommand,#.split(),
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()

        if not error.decode('UTF-8'):
            p.sample_env[sample]['TOTAL_READS'] = output.decode('UTF-8')
        else:
            msg = " ERROR: Could not extract the total number read of sample "+ sample
            logging.error(msg)

        # Getting total number of on-target reads
        bashCommand = ('{} view -c {} -L {}').format(p.system_env['SAMTOOLS'], \
        p.sample_env[sample]['READY_BAM'], p.analysis_env['PANEL'])
        msg = " INFO: Getting total number of reads for sample " + sample
        logging.info(msg)
        logging.info(bashCommand)
        print (msg)

        process = subprocess.Popen(bashCommand,#.split(),
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()

        if not error.decode('UTF-8'):
            p.sample_env[sample]['ON_TARGET_READS'] = output.decode('UTF-8')
        else:
            msg = " ERROR: Could not extract the total on-target reads of sample "+ sample
            logging.error(msg)

        if p.sample_env[sample]['TOTAL_READS'] != '.' \
            and p.sample_env[sample]['TOTAL_READS'] != '.':
            p.sample_env[sample]['ROI_PERCENTAGE'] = \
               round(100*(int(p.sample_env[sample]['ON_TARGET_READS'])/int(p.sample_env[sample]['TOTAL_READS'])),3)

        # Getting insert sizes from the first N aligments
        N = 5000
        bashCommand = ('{} view {} | head -{}').format(p.system_env['SAMTOOLS'], \
        p.sample_env[sample]['READY_BAM'], N)

        msg = " INFO: Getting insert size stats for " + sample
        logging.info(msg)
        logging.info(bashCommand)
        print (msg)

        process = subprocess.Popen(bashCommand,#.split(),
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()

        alignments = []
        i_sizes    = []
        if not error.decode('UTF-8'):
            alignments = output.decode('UTF-8').split('\n')
            for aln in alignments:
                tmp = aln.split('\t')
                if len(tmp) < 10:
                    continue
                isize = abs(int(tmp[8]))
                i_sizes.append(isize)
            array = np.array(i_sizes)
            mean_isize = stats.trim_mean(array, 0.1)
            sd_isize   = stats.tstd(array)
            #print(mean_isize+ " "+ sd_isize)
            p.sample_env[sample]['MEAN_INSERT_SIZE'] = str(mean_isize)
            p.sample_env[sample]['SD_INSERT_SIZE']   = str(sd_isize)

        else:
            msg = " ERROR: Could not extract insert sizes from sample "+ sample
            logging.error(msg)


# def calculate_hts_metrics():

#     do_generate_list_file()

#     for sample in p.sample_env:

#         # Create QC directory
#         p.sample_env[sample]['QC_FOLDER'] = \
#             p.sample_env[sample]['SAMPLE_FOLDER'] + "/" + "QC_FOLDER"

#         qc_path = Path(p.sample_env[sample]['QC_FOLDER'])
#         if not qc_path.is_dir():
#             os.mkdir(qc_path)

#         p.sample_env[sample]['SUMMARY_QC'] = \
#             p.sample_env[sample]['QC_FOLDER'] + "/" + sample + "_summary_qc.csv"      

#         p.sample_env[sample]['HS_METRICS'] = \
#             p.sample_env[sample]['QC_FOLDER'] + "/" + sample + ".hs.metrics.txt"  

#         p.sample_env[sample]['HS_METRICS_NAME'] = sample + ".hs.metrics.txt"  

#         if not os.path.isfile(p.sample_env[sample]['SUMMARY_QC']):

#             bashCommand = ('{} run -v {}:/bam_data/ -v {}:/panel_data/ -v {}:/qc_data/' \
#             ' -it {} gatk CollectHsMetrics -BI /panel_data/{} -I /bam_data/{} -O /qc_data/{} -TI /panel_data/{}') \
#                 .format(p.system_env['DOCKER'], p.sample_env[sample]['BAM_FOLDER'], \
#                  p.defaults['PANEL_FOLDER'], p.sample_env[sample]['QC_FOLDER'], \
#                  p.docker_env['GATK'], p.analysis_env['PANEL_LIST_NAME'], \
#                  p.sample_env[sample]['READY_BAM_NAME'], p.sample_env[sample]['HS_METRICS_NAME'], \
#                  p.analysis_env['PANEL_LIST_NAME'])

#             msg = " INFO: Extracting HS metrics for sample " + sample
#             logging.info(msg)
#             logging.info(bashCommand)

#             process = subprocess.Popen(bashCommand,#.split(),
#                 shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#             output, error = process.communicate()
            
#             if not error.decode('UTF-8'):
#                 msg = " INFO: HS metrics extraction ended OK"
#                 logging.error(msg)
#             else:
#                 msg = " ERROR: Something went wrong with HS metrics extraction"
#                 logging.info(msg)


# Required Arguments:

# --BAIT_INTERVALS,-BI <File>   An interval list file that contains the locations of the baits used.  This argument must
#                               be specified at least once. Required. 

# --INPUT,-I <File>             An aligned SAM or BAM file.  Required. 

# --OUTPUT,-O <File>            The output file to write the metrics to.  Required. 

# --TARGET_INTERVALS,-TI <File> An interval list file that contains the locations of the targets.  This argument must be
#                               specified at least once. Required. 

def extract_coverage_metrics():

    p.analysis_env['SUMMARY_QC'] = p.analysis_env['OUTPUT_DIR'] + \
        "/" + "summary_qc.tsv"

    if os.path.isfile(p.analysis_env['SUMMARY_QC']):
        msg = " INFO: Skipping coverage metric extraction"
        print (msg)
        logging.info(msg)
       # return

    for sample in p.sample_env:

        # Create QC directory
        p.sample_env[sample]['QC_FOLDER'] = \
            p.sample_env[sample]['SAMPLE_FOLDER'] + "/" + "QC_FOLDER"

        qc_path = Path(p.sample_env[sample]['QC_FOLDER'])
        if not qc_path.is_dir():
            os.mkdir(qc_path)

        p.sample_env[sample]['SUMMARY_QC'] = \
             p.sample_env[sample]['QC_FOLDER'] + "/" + sample + "_summary_qc.csv"

        p.sample_env[sample]['COV_THRESHOLDS_FILE'] = \
            p.sample_env[sample]['QC_FOLDER'] + "/" + sample + ".thresholds.bed.gz"

        # Getting cov metrics with Mosdepth
        if not os.path.isfile(p.sample_env[sample]['COV_THRESHOLDS_FILE']):
            bashCommand = ('{} --thresholds 1,10,20,30,100,200 --fast-mode --by {} {} {}') \
                .format(p.system_env['MOSDEPTH'], p.analysis_env['PANEL'],  \
                p.sample_env[sample]['QC_FOLDER'] + "/" + sample, p.sample_env[sample]['READY_BAM'])

            msg = " INFO: Extracting coverage metrics of sample " + sample
            print (msg)
            logging.info(msg)
            logging.info(bashCommand)

            process = subprocess.Popen(bashCommand,#.split(),
                shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = process.communicate()
            
            if not error.decode('UTF-8'):
                msg = " INFO: Coverage metrics extraction ended OK"
                print (msg)
                logging.error(msg)
            else:
                msg = " ERROR: Something went wrong with coverage metrics extraction"
                print (msg)
                logging.info(msg)
        else:
            msg = " INFO: Skipping coverage metrics analysis"
            print (msg)
            logging.info(msg)

        # Now getting call_rate and lost_exons from thresholds file
        call_rate_dict = defaultdict(dict)
        lost_exons_dict= defaultdict(dict)
        total_bases = 0
        total_rois  = 0
        header = []
        with gzip.open(p.sample_env[sample]['COV_THRESHOLDS_FILE'],'rt') as fin:        
            for line in fin:
                line = line.rstrip('\n')
                if line.startswith('#'):
                    tmp = line.split('\t')
                    for idx in range (4,len(tmp)-1):
                        header.append(tmp[idx])
                else:
                    tmp = line.split('\t')
                    length = int(tmp[2])-int(tmp[1])
                    total_bases+=length
                    i = 0
                    for x in range (4,len(tmp)-1):
                        field = header[i]
                        # Checking number of bases at the threshold
                        if field not in call_rate_dict:
                            call_rate_dict[field] = 0
                        call_rate_dict[field] += int(tmp[x])
                        if field not in lost_exons_dict:
                            lost_exons_dict[field] = 0
                        # Now checking lost exons at the threshold
                        if int(tmp[x]) < length:
                            lost_exons_dict[field]+=1
                        i+=1
                total_rois+=1
        p.analysis_env['ROI_NUMBER'] = str(total_rois)
        if not 'CALL_RATE' in p.sample_env[sample]:
            p.sample_env[sample]['CALL_RATE'] = {}
        if not 'LOST_EXONS' in p.sample_env[sample]:
            p.sample_env[sample]['LOST_EXONS'] = {}                 
        for field in call_rate_dict:
            p.sample_env[sample]['CALL_RATE'][field] = round(100*(call_rate_dict[field]/total_bases),3)

            if not field in p.sample_env[sample]['LOST_EXONS']:
                p.sample_env[sample]['LOST_EXONS'][field] = 0
            p.sample_env[sample]['LOST_EXONS'][field] = lost_exons_dict[field]
        # Now get the mean coverage
        p.sample_env[sample]['MOSDEPTH_SUMMARY'] =  \
            p.sample_env[sample]['QC_FOLDER'] + "/" + sample + ".mosdepth.summary.txt"
        
        with open(p.sample_env[sample]['MOSDEPTH_SUMMARY']) as f:
            for line in f:
                line = line.rstrip('\n')
                if line.startswith('total_region'):
                    tmp = line.split('\t')
                    p.sample_env[sample]['MEAN_COVERAGE'] = tmp[3]
        f.close()

def remove_duplicates():

  for sample in p.sample_env:

    # Setting rmdup_bam
    p.sample_env[sample]['RMDUP_BAM'] = \
      p.sample_env[sample]['RAW_BAM'].replace(".bam", ".rmdup.bam")

    p.sample_env[sample]['RMDUP_BAM_NAME'] = \
      os.path.basename(p.sample_env[sample]['RMDUP_BAM'])

    # Now ready bam is rmdup_bam
    p.sample_env[sample]['READY_BAM'] = p.sample_env[sample]['RMDUP_BAM']

    p.sample_env[sample]['READY_BAM_NAME'] = \
      os.path.basename(p.sample_env[sample]['READY_BAM'])

    # Create QC directory
    p.sample_env[sample]['QC_FOLDER'] = \
      p.sample_env[sample]['SAMPLE_FOLDER'] + "/" + "QC_FOLDER"

    qc_path = Path(p.sample_env[sample]['QC_FOLDER'])
    if not qc_path.is_dir():
        os.mkdir(qc_path)

    # Set picard metrics file
    p.sample_env[sample]['PICARD_METRICS'] = \
      p.sample_env[sample]['QC_FOLDER'] + "/" + sample + ".picard.metrics.txt"

    p.sample_env[sample]['PICARD_METRICS_NAME'] = sample + ".picard.metrics.txt"

    # Now map fastq files with bwa mem
    bashCommand = ('{} run -v {}:/work_data/ -it {} gatk MarkDuplicates '
        '-I /work_data/{} -O /work_data/{} -M /work_data/{} '
        '-REMOVE_DUPLICATES true -ASSUME_SORTED true'.format(
        p.system_env['DOCKER'], p.sample_env[sample]['BAM_FOLDER'],
        p.docker_env['GATK'], p.sample_env[sample]['RAW_BAM_NAME'],
        p.sample_env[sample]['RMDUP_BAM_NAME'], p.sample_env[sample]['PICARD_METRICS_NAME']))

    if not os.path.isfile(p.sample_env[sample]['RMDUP_BAM']) \
        or not os.path.isfile(p.sample_env[sample]['PICARD_METRICS']):

        msg = " INFO: Removing duplicates of sample " + sample
        print (msg)
        logging.info(msg)
        logging.info(bashCommand)

        process = subprocess.Popen(bashCommand,#.split(),
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()
        #check if everything is ok
        if not error.decode('UTF-8'):
            if re.search(r'MarkDuplicates done.', output.decode('UTF-8')):
                #rm dup ok
                msg = " INFO: RM duplicated process ended OK for sample " + sample
                print (msg)
                logging.info(msg)

                bam_summary_picard = p.sample_env[sample]['BAM_FOLDER'] + "/" + sample + ".picard.metrics.txt"
                qc_summary_picard  = p.sample_env[sample]['QC_FOLDER'] + "/" + sample + ".picard.metrics.txt"

                shutil.copy2(bam_summary_picard, qc_summary_picard)
            else:
                msg = " ERROR: Something went wrong with rm duplicates " + sample
                print (msg)
                logging.error(msg)
        else:
            msg = " ERROR: Something went wrong with rm duplicates " + sample
            print (msg)
            logging.error(msg)
    else:
        msg = " INFO: Skipping duplicate removal for sample " + sample
        print (msg)
        logging.error(msg)

    index_bam(p.sample_env[sample]['READY_BAM'])

    # now get the number of pcr duplicates
    with open (p.sample_env[sample]['PICARD_METRICS']) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith(sample) or line.startswith("Unknown"):
                tmp = line.split('\t')
                perc_duplicates = round(100*float(tmp[8]),3)
                p.sample_env[sample]['PCR_DUPLICATES_PERCENTAGE'] = str(perc_duplicates)
    f.close()

def do_generate_list_file():
    '''generate list file from bed file
    '''

    bashCommand = ('{} run -v {}:/bundle/ -v {}:/panels/ '
        '-it {} gatk BedToIntervalList -I /panels/{} -O /panels/{} '
        '-SD /bundle/{} --SORT --UNIQUE'.format(p.system_env['DOCKER'],
        p.defaults['BUNDLE_FOLDER'], p.analysis_env['PANEL_WORKDIR'], p.docker_env['GATK'], \
        p.analysis_env['PANEL_NAME'], p.analysis_env['PANEL_LIST_NAME'], p.aux_env['GENOME_DICT_NAME']))


    if not os.path.isfile(p.analysis_env['PANEL_LIST']):

        msg = " INFO: Generating list for the gene panel "+ p.analysis_env['PANEL_NAME']
        print (msg)
        logging.info(msg)
        logging.info(bashCommand)

        process = subprocess.Popen(bashCommand,#.split(),
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()

        if not error.decode('UTF-8'):
            if re.search(r'BedToIntervalList done.', output.decode('UTF-8')):
                msg = " INFO: list file generation for " + p.analysis_env['PANEL_NAME'] + "ended OK."
                print (msg)
                logging.info(msg)
            else:
                msg = " ERROR: Something went wrong with list file generation for " + p.analysis_env['PANEL_NAME'] 
                print (msg)
                logging.error(msg)
        else:
            msg = " ERROR: Something went wrong with list file generation for " + p.analysis_env['PANEL_NAME'] 
            print (msg)
            logging.error(msg)
    else:
        msg = " INFO: Skipping list file generation for panel " + p.analysis_env['PANEL_NAME'] + " . File is already available"
        print (msg)
        logging.info(msg)

def map_fastq():

  for sample in p.sample_env:

    # Create BAM_FOLDER if not present
    p.sample_env[sample]['BAM_FOLDER'] = \
      p.sample_env[sample]['SAMPLE_FOLDER'] + "/" + "BAM_FOLDER"

    bam_folder_path = Path(p.sample_env[sample]['BAM_FOLDER'])
    if not bam_folder_path.is_dir():
      os.mkdir(bam_folder_path)

    # Defining raw bam
    p.sample_env[sample]['RAW_BAM'] = \
      p.sample_env[sample]['BAM_FOLDER'] + "/" + sample + ".bam"
    p.sample_env[sample]['RAW_BAM_NAME'] = sample + ".bam"

    # READY_BAM will be overwritten after each post-processing step (e.g remove dups, etc) 
    p.sample_env[sample]['READY_BAM'] = p.sample_env[sample]['RAW_BAM']
    p.sample_env[sample]['READY_BAM_NAME'] = sample + ".bam"

    # Now map fastq files with bwa mem 
    # Plus samtools sort and conversion (avoid threading in sam->bam if you have less than 16Gb of RAM)
    bashCommand = ('{} mem {} -R \'@RG\\tID:{}\\tSM:{}\' -M -t {} {} {} | {} view -Shu - |' 
      '{} sort -T {} -o {}').format(
      p.system_env['BWA'], 
      p.aux_env['GENOME_FASTA'],
      sample,
      sample,
      p.analysis_env['THREADS'], 
      p.sample_env[sample]['READY_FQ1'], 
      p.sample_env[sample]['READY_FQ2'],
      p.system_env['SAMTOOLS'], 
      p.system_env['SAMTOOLS'], 
      'TMP',
      p.sample_env[sample]['RAW_BAM']
      )

    if not os.path.isfile(p.sample_env[sample]['RAW_BAM']):

      msg = " INFO: Mapping sample " +  sample
      print (msg)
      p.logging.info(msg)
      p.logging.info(bashCommand)
    
      process = subprocess.Popen(bashCommand, shell=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
      
      output, error = process.communicate()
    else:
        msg = " INFO: Skipping mapping for sample " + sample + " BAM file is already available"
        print (msg)
        p.logging.info(msg)
    index_bam(p.sample_env[sample]['RAW_BAM'])

def index_bam(bam):

    bai = bam + ".bai"
    #generate index
    bashCommand = '{} index {}'.format(p.system_env['SAMTOOLS'], bam)
    
    if not os.path.isfile(bai):

        msg = " INFO: Indexing bam " + bam
        print (msg)
        logging.info(bashCommand)

        process = subprocess.Popen(bashCommand,#.split(),
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()
        if not error.decode('UTF-8') and not output.decode('UTF-8'):
            msg = " INFO: Bam indexing for " + bam + " ended OK"
            logging.info(msg)
        else:
            msg = " ERROR: Something went wrong with Indexing process for " + bam
            logging.error(msg)
