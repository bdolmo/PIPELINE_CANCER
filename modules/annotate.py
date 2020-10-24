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
import pprint
import json

from modules import utils as u
from modules import params as p
from modules import trimming as t
from modules import sqlite as s

def do_annotation():

    do_vep()

    do_civic()

class Civic:
  '''civic class
  '''

  def __init__(self):
    # Defining attributes
    self.variants_dict = defaultdict(dict)
    self.evidence_dict = defaultdict(dict)

  def queryCivic(self, gene, variant, vartype, exon, conseq):
    '''
       return list of evidences given a gene,variant,exon and conseq
    '''
    ev_list = []
    candidate_var = variant

    candidate_var_list = []
    candidate_var_list.append(variant)
    candidate_var_list.append(variant[:-1])
    candidate_var_list.append(variant[1:])
    candidate_var_list.append(variant[1:-1])    
    candidate_var_list.append(variant[:-1]+"X")

    if conseq == 'intron_variant' or conseq == 'splice_region_variant':
      v = "EXON " +  str(exon)
      candidate_var_list.append(v)
    if vartype == 'DELETION':
      v = "EXON " + str(exon) + " DELETION"
      candidate_var_list.append(v)
    if vartype == 'INSERTION':
      v ="EXON " + str(exon) + " INSERTION"
      candidate_var_list.append(v)

    chosen_var_id = '.'
    if gene in self.variants_dict:
      # Now look for all variants
      
      for var in candidate_var_list:
        for mut in self.variants_dict[gene]:
          mut_uc = mut.upper()
          if var in mut_uc:
            chosen_var_id = self.variants_dict[gene][mut]
            break
        if chosen_var_id != '.':
          break

      ev_list = []
      # Get all evidences
      for ev_id in self.evidence_dict[chosen_var_id]:
        ev_direction   =  self.evidence_dict[chosen_var_id][ev_id]['evidence_direction']
        ev_level       =  self.evidence_dict[chosen_var_id][ev_id]['evidence_level']
        ev_significance=  self.evidence_dict[chosen_var_id][ev_id]['clinical_significance']
        ev_drugs       =  self.evidence_dict[chosen_var_id][ev_id]['drugs']
        ev_disease     =  self.evidence_dict[chosen_var_id][ev_id]['disease']
        ev_pmid        =  self.evidence_dict[chosen_var_id][ev_id]['pmid']
        ev_clintrials  =  self.evidence_dict[chosen_var_id][ev_id]['clinical_trials']
        ev_info = []
        ev_info.append(ev_direction)
        ev_info.append(ev_level)
        ev_info.append(ev_significance)
        ev_info.append(ev_drugs)
        ev_info.append(ev_disease)
        ev_info.append(ev_pmid)
        ev_info.append(ev_clintrials)
        ev_str = '|'.join(ev_info)
        ev_list.append(ev_str)
    return ev_list

  def loadCivic(self):
    '''load civic data in memory
    '''
    URL = "https://civicdb.org/api/variants?count=10000"
    r = requests.get(url = URL) 

    if not r.ok:
      r.raise_for_status()
      sys.exit()
    decoded = r.json()

    # Load all variants
    for variant in decoded['records']:
      gene_name = variant['entrez_name']
      variant_name = variant['name']
      self.variants_dict[gene_name][variant_name] = variant['id']

    URL = "https://civicdb.org/api/evidence_items?count=10000"
    r = requests.get(url = URL) 
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    decoded = r.json()

    # Load all evidences
    for evidence in decoded['records']:
      variant_id = evidence['variant_id']
      if not variant_id in self.evidence_dict:
        self.evidence_dict[variant_id] = defaultdict(dict)

      ev_id = evidence['name']
      self.evidence_dict[variant_id][ev_id]['evidence_direction'] = evidence['evidence_direction']
      self.evidence_dict[variant_id][ev_id]['evidence_level'] = evidence['evidence_level']
      self.evidence_dict[variant_id][ev_id]['clinical_significance'] = evidence['clinical_significance']
      drugs_list = []
      for drug in evidence['drugs']:
        drugs_list.append(drug['name'])
      drug_set = set(drugs_list)
      drugs_list = list(drug_set)
      self.evidence_dict[variant_id][ev_id]['drugs']  = ','.join(drugs_list) 
      self.evidence_dict[variant_id][ev_id]['disease'] = evidence['disease']['name'].replace(" ", "_")
      self.evidence_dict[variant_id][ev_id]['pmid'] = evidence['source']['citation_id']
      clintrials_list = []
      for trial in evidence['source']['clinical_trials']:
        clintrials_list.append(trial['nct_id'])
      self.evidence_dict[variant_id][ev_id]['clinical_trials'] = ','.join(clintrials_list)

def do_civic():

    msg = " INFO: Loading CIViC database"
    print(msg)
    logging.info(msg)

    #Instantiante CIViC class
    civic = Civic()

    # In-memory loading of CIViC database
    civic.loadCivic()

    for sample in p.sample_env:

        p.sample_env[sample]['CIVIC_VCF'] = \
            p.sample_env[sample]['READY_SNV_VCF'].replace(".vcf", ".civic.vcf")
        p.sample_env[sample]['CIVIC_VCF_NAME'] = \
            os.path.basename(p.sample_env[sample]['READY_SNV_VCF_NAME']).replace(".vcf", ".vep.vcf")

        print(p.sample_env[sample]['READY_SNV_VCF'])
        sys.exit()

        vep_dict = defaultdict(dict)
        vep_list = []
        o = open(p.sample_env[sample]['CIVIC_VCF'], 'w')

        civic_fields = []
        civic_fields.append("EV_ID")
        civic_fields.append("EV_DIRECTION")
        civic_fields.append("EV_LEVEL")
        civic_fields.append("EV_SIGNIFICANCE")
        civic_fields.append("EV_DRUGS")
        civic_fields.append("EV_DISEASE")
        civic_fields.append("EV_PMID")
        civic_fields.append("EV_CLINICAL_TRIALS")
        civic_info_header = "##INFO=<ID=CIVIC,Number=.,Type=String,Description=\"Civic evidence. Format: " + '|'.join(civic_fields) + "\">"
        with open (p.sample_env[sample]['READY_SNV_VCF']) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    if re.search("ID=CSQ", line):
                        vep_dict, vep_list = u.create_vep_dict(line)
                        o.write(line+"\n")
                        o.write(civic_info_header+"\n")
                        continue
                    o.write(line+"\n")
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

                    vartype = '.'
                    if len(ref)==1 and len(alt) ==1:
                      vartype = "SNV"
                    elif len(ref) >1 and len(alt) > 1 \
                     and len(ref) == len(alt):
                      vartype = "MNV"
                    elif len(ref)>len(alt):
                      vartype = "DELETION"
                    elif len(ref)<len(alt):
                      vartype = "INSERTION"

                    info_list = info.split(';')
                    idx = 0
                    civic_annotation = '.'
                    for item in info_list:
                        if item.startswith('CSQ'):
                            tmp_transcript = item.split(",")
                            for transcript_info in tmp_transcript:
                                transcript_info = transcript_info.replace("CSQ=", "")
                                transcript_list = transcript_info.split("|")
                                gene = transcript_list[vep_dict['SYMBOL']]
                                ensg_id = transcript_list[vep_dict['Gene']]
                                enst_id = transcript_list[vep_dict['Feature']]

                                exon = re.search("\d+", transcript_list[vep_dict['EXON']])
                                consequence = transcript_list[vep_dict['Consequence']]
                                aminoacids = transcript_list[vep_dict['Amino_acids']]
                                protein_position = transcript_list[vep_dict['Protein_position']]
                                variant = "."
                                if not exon:
                                    exon = '.'
                                if not consequence:
                                    consequence = '.'
                                if aminoacids != "" and protein_position != "":
                                  if "/" in aminoacids:
                                    variant = aminoacids.replace("/", protein_position)
                                  else:
                                    variant = aminoacids + protein_position
                                  civic_ev_list = civic.queryCivic(gene, variant, vartype, exon, consequence)
                                  civic_annotation = ','.join(civic_ev_list)
                            #item = tmp_transcript[tidx]
                        idx+=1

                    info = ';'.join(info_list)
                    tmp[7] = info + ";CIVIC="+civic_annotation
                    line = '\t'.join(tmp)
                    o.write(line+"\n")
        f.close()
        o.close()
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

        bashCommand = ('{} run -t -i -v {}:/civic_folder/ -v {}:/opt/vep/.vep ensemblorg/ensembl-vep'
        ' perl vep --cache --offline --dir_cache /opt/vep/.vep/ --dir_plugins /opt/vep/.vep/Plugins/'
        ' --input_file /opt/vep/.vep/input/{} --output_file /opt/vep/.vep/output/{} '
        ' --af_1kg --af_gnomad '
        ' --format vcf --vcf --hgvs --hgvsg --max_af --pubmed --gene_phenotype --ccds --sift b --polyphen b'
        ' --symbol --force_overwrite --fork {} --canonical '
        .format(p.system_env['DOCKER'], p.defaults['CIVIC_FOLDER'], p.defaults['VEP_DATA'], p.sample_env[sample]['READY_SNV_VCF_NAME'],\
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
