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
import json
from modules import utils as u
from modules import params as p
from modules import trimming as t
from modules import sqlite as s

def do_report():

    clinical_vcf()

    clinical_report_csv()

def clinical_report_csv():

    for sample in p.sample_env:

        p.sample_env[sample]['REPORT_FOLDER'] = p.sample_env[sample]['SAMPLE_FOLDER'] + "/" + "REPORT_FOLDER"

        report_folder_path = Path(p.sample_env[sample]['REPORT_FOLDER'])
        if not report_folder_path.is_dir():
            os.mkdir(report_folder_path)

        p.sample_env[sample]['CLINICAL_REPORT_CSV'] = \
            p.sample_env[sample]['REPORT_FOLDER'] + "/" + sample + ".clinical.csv"
        o =  open (p.sample_env[sample]['CLINICAL_REPORT_CSV'], "w")

        report_header = [] 
        report_header.append("Lab_ID")
        report_header.append("Ext_ID")
        report_header.append("Gene")
        report_header.append("P_CODE")
        report_header.append("G_CODE")
        report_header.append("C_CODE")
        report_header.append("Exon")
        report_header.append("ENST")
        report_header.append("VAF")
        report_header.append("Identifiers")
        report_header.append("MAX_AF")
        report_header.append("MAX_AF_POP")
        report_header.append("Drugs")
        report_header.append("Clinical_trials")
        o.write('\t'.join(report_header)+"\n")
        var_dict = defaultdict(dict)
        with open(p.sample_env[sample]['READY_MERGED_JSON']) as f:
            var_dict = json.load(f)
        
        for variant in var_dict['variants']:
          results_list = [] 
          lab_id = sample
          results_list.append(lab_id)
          ext_id = "999999"
          results_list.append(ext_id)

          if 'SVTYPE' in var_dict['variants'][variant]['INFO']:

            if 'PR' in var_dict['variants'][variant]:
              pe_info = var_dict['variants'][variant]['PR'].split(",")
              pe_ref = int(pe_info[0])
              pe_alt = int(pe_info[1])
            else:
              pe_ref = 0
              pe_alt = 0
            if 'SR' in var_dict['variants'][variant]:
              sr_info = var_dict['variants'][variant]['SR'].split(",")
              sr_ref = int(sr_info[0])
              sr_alt = int(sr_info[1])
            else:
              sr_ref = 0
              sr_alt = 0
            fusion_out = False
            if (pe_alt > 0 and sr_alt > 0):
              fusion_out = True
            if sr_alt+sr_ref > 0:
              VAF = str(round((sr_alt/(sr_alt+sr_ref)),3))
            elif pe_ref+pe_alt >0:
              VAF = str(round((pe_alt/(pe_alt+pe_ref)),3))
            else:
              VAF = "0"
            fusion = '.'
            if var_dict['variants'][variant]['INFO']['FUSION']['PARTNERS'] != '.' and \
              var_dict['variants'][variant]['INFO']['FUSION']['PARTNERS'] != "" :
              fusion = var_dict['variants'][variant]['INFO']['FUSION']['PARTNERS'] + " FUSION"
            if fusion != '.' or fusion_out == True:
              
              results_list.append(fusion)
              p_code = '.'
              results_list.append(p_code)
              g_code = '.'
              results_list.append(g_code)
              c_code = '.'
              results_list.append(c_code)
              exon = '.'
              results_list.append(exon)
              enst_id = '.'
              results_list.append(enst_id)
              results_list.append(VAF)
              clin_sig = '.'
              results_list.append(clin_sig)
              rs_id = '.'
              results_list.append(rs_id)
              max_af = '.'
              results_list.append(max_af)
              max_af_pop = '.'
              results_list.append(max_af_pop)
              drugs = '.'
              results_list.append(drugs)
              clin_trials = '.'
              results_list.append(clin_trials)
              o.write('\t'.join(results_list)+ "\n")
          else:
            gene   = var_dict['variants'][variant]['INFO']['CSQ']['SYMBOL']
            results_list.append(gene)
            p_code = var_dict['variants'][variant]['INFO']['CSQ']['HGVSp']
            if p_code:
              tmp_p_code = p_code.split(":")
              if len(tmp_p_code)>1:
                p_code = tmp_p_code[1] 
            results_list.append(p_code)
            g_code = var_dict['variants'][variant]['INFO']['CSQ']['HGVSg'] 
            results_list.append(g_code)
            c_code = var_dict['variants'][variant]['INFO']['CSQ']['HGVSc']
            if c_code:
              tmp_c_code = c_code.split(":")
              if len(tmp_c_code)> 1:
                c_code = tmp_c_code[1]
            results_list.append(c_code)
            exon   = var_dict['variants'][variant]['INFO']['CSQ']['EXON']
            results_list.append(exon)
            enst_id= var_dict['variants'][variant]['INFO']['CSQ']['Feature']
            results_list.append(enst_id)
            VAF    = var_dict['variants'][variant]['AF']
            results_list.append(VAF)
            clin_sig = var_dict['variants'][variant]['INFO']['CSQ']['CLIN_SIG']
          # results_list.append(clin_sig)
            rs_id  = var_dict['variants'][variant]['INFO']['CSQ']['Existing_variation']
            rs_id = rs_id.replace('&', ',')
            results_list.append(rs_id)
            max_af = var_dict['variants'][variant]['INFO']['CSQ']['MAX_AF']
            results_list.append(max_af)
            max_af_pop = var_dict['variants'][variant]['INFO']['CSQ']['MAX_AF_POPS']
            max_af_pop = max_af_pop.replace('&', ',')
            if max_af != '.':
              if float(max_af) == 0:
                max_af_pop = '.'
            results_list.append(max_af_pop)
            drugs_list = []
            clin_trials_list = []  
            if 'CIVIC' in var_dict['variants'][variant]['INFO']:
              for ev_id in var_dict['variants'][variant]['INFO']['CIVIC']:
                if 'EV_DRUGS' in var_dict['variants'][variant]['INFO']['CIVIC']:
                  drugs = var_dict['variants'][variant]['INFO']['CIVIC']['EV_DRUGS']  
                  drugs_list.append(drugs)
                if 'EV_CLINICAL_TRIALS' in var_dict['variants'][variant]['INFO']['CIVIC']:
                  clin_trials = var_dict['variants'][variant]['INFO']['CIVIC']['EV_CLINICAL_TRIALS']
                  clin_trials_list.append(clin_trials)

            if drugs_list :
              drugs_str = ','.join(list(set(drugs_list)))
            else:
              drugs_str = '.'
            if clin_trials_list:
              clintrials_str = ','.join(list(set(clin_trials_list)))
            else:
              clintrials_str = '.'
            results_list.append(drugs_str)
            results_list.append(clintrials_str)

            if 'pathogenic' in clin_sig:
              o.write('\t'.join(results_list)+ "\n")

        o.close()





def clinical_vcf():

    for sample in p.sample_env:

        p.sample_env[sample]['CLINICAL_SNV_VCF'] = p.sample_env[sample]['READY_SNV_VCF'].replace(".vcf", ".clinical.vcf")
        o = open(p.sample_env[sample]['CLINICAL_SNV_VCF'], 'w')

        vep_dict = defaultdict(dict)
        vep_list = []

        civic_dict = defaultdict(dict)
        civic_list = []

        with open (p.sample_env[sample]['READY_SNV_VCF']) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    o.write(line+"\n")
                    if re.search("ID=CSQ", line):
                        vep_dict, vep_list = u.create_vep_dict(line)
                    if re.search("ID=CIVIC", line):
                        civic_dict, civic_list = u.create_vep_dict(line)
                else:
                    tmp = line.split('\t')
                    chr = tmp[0]
                    pos = tmp[1]
                    id  = tmp[2]
                    ref = tmp[3]
                    alt = tmp[4]
                    qual= tmp[5]
                    filter = tmp[6]
                    info = tmp[7]
                    format_tag = tmp[8]
                    format = tmp[9]

                    info_list = info.split(';')
                    idx = 0
                    for item in info_list:
                        if item.startswith('CSQ'):
                            tmp_transcript = item.split(",")
                            for transcript_info in tmp_transcript:
                                transcript_info = transcript_info.replace("CSQ=", "")
                                transcript_list = transcript_info.split("|")
                                gene = transcript_list[vep_dict['SYMBOL']]
                                ensg_id = transcript_list[vep_dict['Gene']]
                                enst_id = transcript_list[vep_dict['Feature']]
                                max_af = 0.00
                                if transcript_list[vep_dict['MAX_AF']] != "":
                                    max_af = float(transcript_list[vep_dict['MAX_AF']])

                                if gene in p.roi_env or ensg_id in p.roi_env:
                                    # Select the annotations from the wanted transcript

                                    if enst_id == p.roi_env[gene] or enst_id == p.roi_env[ensg_id]:
                                        # Now, dump the variant if present at civic

                                        if transcript_list[civic_dict['EV_DIRECTION']] == "Supports" \
                                            or 'pathogenic' in transcript_list[vep_dict['CLIN_SIG']]:
                                            line = '\t'.join(tmp)
                                            o.write(line+"\n")
                                            info_list[idx] = transcript_info
                                            break
                                    else:
                                        continue
                                else:
                                    continue
                            #item = tmp_transcript[tidx]
                        idx+=1

                    info = ';'.join(info_list)
                    tmp[7] = info
                    line = '\t'.join(tmp)
                    #o.write(line+"\n")
        f.close()
        o.close()

