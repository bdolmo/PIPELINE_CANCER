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

    # Variant effect Predict annotation
    do_vep()

    # Select VEP-matching transcripts with our panel-defined transcripts
    filter_transcripts()

    if p.analysis_env['VARIANT_CLASS'] == "somatic":

        # CIViC annotation for SNV/Indels
        do_civic()

        # Cancer Genome Interpreter annotation
        do_cgi()

        # Fusion annotation with chimerKB  
        annotate_known_fusions()

        # Add flanking genes to SVs
        add_edge_genes()

        # CIViC annotation for CNAs
        annotate_cnas()
        
        # Merge snv and fusion vcf  
        merge_vcfs()
        # sys.exit()

def annotate_cnas():
  
  msg = " INFO: Annotating CNAs with CIViC"
  logging.info(msg)
  print(msg)

  do_ann = True
  for sample in p.sample_env:
    p.sample_env[sample]['READY_CNA_VCF'] = \
      p.sample_env[sample]['CNV_VCF'].replace(".vcf", ".annotated.vcf")
    if os.path.isfile(p.sample_env[sample]['READY_CNA_VCF']):
      do_ann = False
    else:
      do_ann = True
      break
  
  if do_ann == False:
    msg = " INFO: Skipping gene annotation to CNAs"
    logging.info(msg)
    print(msg)
    return

  # Instantiante CIViC class
  civic = Civic()

  # In-memory loading of CIViC database
  civic.loadCivic()

  for sample in p.sample_env:

    p.sample_env[sample]['READY_CNA_VCF'] = \
      p.sample_env[sample]['CNV_VCF'].replace(".vcf", ".annotated.vcf")

    if not os.path.isfile(p.sample_env[sample]['READY_CNA_VCF']):

      tmp_intersect =  p.sample_env[sample]['READY_CNA_VCF'].replace(".vcf", ".bed")

      bashCommand = '{} intersect -a {} -b {} -header -wa -wb | cut -f 1,2,3,4,5,6,7,8,9,10,14 | uniq > {}' \
        .format(p.system_env['BEDTOOLS'], p.sample_env[sample]['CNV_VCF'], p.analysis_env['PANEL'], tmp_intersect)

      logging.info(bashCommand)
      p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      output = p1.stdout.decode('UTF-8')
      error  = p1.stderr.decode('UTF-8')
      if not error and not output:
        msg = " INFO: CNA gene annotation for sample " + sample + " ended successfully"
        print(msg)
        logging.info(msg)
      else:
        msg = " ERROR: Something went wrong with CNA gene annotation for sample " + sample
        print(msg)
        logging.error(msg)

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
      genes_info_header = "##INFO=<ID=CNA_GENES,Number=1,Type=String,Description=\"Genes affected by a CNA\">" 
      seen = defaultdict(dict)
      o = open(p.sample_env[sample]['READY_CNA_VCF'], "w")
      with open(tmp_intersect) as f:
        for line in f:
          line = line.rstrip("\n")
          if not line in seen:
            seen[line] = 0
          else:
            seen[line]+=1
            if seen[line] > 1:
              continue   
          tmp = line.split("\t")
          if line.startswith("#"):
            if line.startswith("#CHROM"):
              o.write(civic_info_header + "\n")
              o.write(genes_info_header + "\n")
              o.write(line+"\n")
            else:
              o.write(line+"\n")
          else:
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
            gene = tmp[-1] 
            vartype = "CNA"
            variant  = ""
            tmp_info = info.split(";")
            for field in tmp_info:
              if 'SVTYPE' in field:
                variant = field.replace("SVTYPE=","")
            exon = '.'
            consequence = '.'
            civic_ev_list = civic.queryCivic(gene, variant, vartype, exon, consequence)
            civic_annotation = ','.join(civic_ev_list)
            
            if civic_annotation == "":
              civic_annotation = "."
            tmp[7] = tmp[7] + ";CNA_GENES=" + gene + ";CIVIC="+civic_annotation
            del tmp[-1]

            line = '\t'.join(tmp)
            o.write(line+"\n")
      o.close()
      if not os.path.isfile(p.sample_env[sample]['READY_CNA_VCF']):
        msg = " ERROR: " + p.sample_env[sample]['READY_CNA_VCF'] + " was not created"
        print(msg)
        logging.error(msg)

      u.convert_vcf_2_json(p.sample_env[sample]['READY_CNA_VCF'])
      os.remove(tmp_intersect)
    else:
        msg = " INFO: Skipping CNA gene annotation for sample " + sample
        print(msg)
        logging.info(msg)     


def add_edge_genes():

  msg = " INFO: Adding flanking genes to SVs"
  logging.info(msg)
  print(msg)

  for sample in p.sample_env:

    fusions_bed = u.vcf_2_bed(p.sample_env[sample]['READY_SV_VCF'])
    intersect_file  = p.sample_env[sample]['VCF_FOLDER'] + "/" + "intersect.bed"
    fusions_vcf = p.sample_env[sample]['READY_SV_VCF'].replace(".vcf", ".tmp.vcf")
    annotated_vcf = p.sample_env[sample]['READY_SV_VCF'].replace(".vcf", ".annotated.vcf")

    if os.path.isfile(annotated_vcf):
      bashCommand = 'grep "EDGE_GENES" {}'.format(annotated_vcf)
      p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      output = p1.stdout.decode('UTF-8')
      error  = p1.stderr.decode('UTF-8')
      if re.search("EDGE_GENES", output):
        continue

    o = open(fusions_vcf, 'w')

    bashCommand = '{} pairtopair -a {} -b {} -type either | sort -V | uniq > {}' \
      .format(p.system_env['BEDTOOLS'], fusions_bed, p.aux_env['GENE_LIST'], intersect_file)
      
    p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p1.stdout.decode('UTF-8')
    error  = p1.stderr.decode('UTF-8')
    if not error and not output:
      msg = " INFO: SV gene annotation for sample " + sample + " ended successfully"
      print(msg)
      logging.info(msg)
    else:
      msg = " ERROR: Something went wrong with SV gene annotation for sample " + sample
      print(msg)
      logging.error(msg)

    fusions_dict = defaultdict(dict)
    with open (intersect_file, "r") as f:
      for line in f:
        line = line.rstrip("\n")
        tmp = line.split("\t")
        coordinate = tmp[0]+"\t"+tmp[1]

        if not coordinate in fusions_dict:
          fusions_dict[coordinate] = defaultdict(dict)
          fusions_dict[coordinate]['LEFT_FLANK'] = []
          fusions_dict[coordinate]['RIGHT_FLANK'] = []

        chrA = tmp[0]
        posA = tmp[1]
        endA = tmp[2]

        chrB = tmp[3]
        posB = tmp[4]      
        endB = tmp[5]

        chr = tmp[7]
        pos = tmp[8]
        end = tmp[9]
        gene = tmp[10]
        if chr == chrA and chr != chrB:
          if posA >= pos and endA <= end:
            fusions_dict[coordinate]['LEFT_FLANK'].append(gene)
        elif chr == chrB and chr != chrA:
          if posB >= pos and endB <= end:
            fusions_dict[coordinate]['RIGHT_FLANK'].append(gene)
        else:
          if posA >= pos and posA <= end:
            fusions_dict[coordinate]['LEFT_FLANK'].append(gene)
          if posB >= pos and posB <= end:
            fusions_dict[coordinate]['RIGHT_FLANK'].append(gene)                     
    f.close()

    fusion_fields = []
    fusion_info_header = "##INFO=<ID=EDGE_GENES,Number=.,Type=String,Description=\"Genes at SV edges\">"

    with open (p.sample_env[sample]['READY_SV_VCF']) as f:
      for line in f:
        line = line.rstrip("\n")
        tmp = line.split("\t")
        if line.startswith("#"):
          if re.search("cmdline", line):
            o.write(line+"\n")
            o.write(fusion_info_header+"\n")
          else:
            o.write(line+"\n")
          continue
        coordinate = tmp[0]+"\t"+tmp[1]
        if coordinate in fusions_dict:
            fusion_ann_list = []
            left_flank = ','.join(list(set( fusions_dict[coordinate]['LEFT_FLANK'])))
            right_flank = ','.join(list(set( fusions_dict[coordinate]['RIGHT_FLANK'])))
            annot = left_flank+"-"+right_flank
          # fusion_ann_list.append(genes)
            tmp[7]+=";EDGE_GENES=" + annot
        else:
          tmp[7]+=";EDGE_GENES=."
        line  = '\t'.join(tmp)
        o.write(line+"\n")
    f.close()

    os.remove(p.sample_env[sample]['READY_SV_VCF'])
    os.rename(fusions_vcf, p.sample_env[sample]['READY_SV_VCF'])
    u.convert_vcf_2_json(p.sample_env[sample]['READY_SV_VCF'])

def filter_transcripts():
    '''
        Select VEP-matching transcripts with our panel
    '''

    for sample in p.sample_env:

        p.sample_env[sample]['TMP_SNV_VCF'] = p.sample_env[sample]['READY_SNV_VCF'].replace(".vcf", ".tmp.vcf")
        o = open(p.sample_env[sample]['TMP_SNV_VCF'], 'w')

        print (p.sample_env[sample]['TMP_SNV_VCF'])

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
                                
                               # if not gene in p.roi_env: 
                                if gene in p.roi_env or ensg_id in p.roi_env:
                                    # Select the annotations from the desired transcript
                                    if enst_id in p.roi_env[gene] or enst_id in p.roi_env[ensg_id]:
                                        info_list[idx] ="CSQ="+transcript_info
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
                    o.write(line+"\n")
        f.close()
        o.close()
        os.remove(p.sample_env[sample]['READY_SNV_VCF'])
        os.rename(p.sample_env[sample]['TMP_SNV_VCF'], p.sample_env[sample]['READY_SNV_VCF'])
          # Create JSON file
       # u.convert_vcf_2_json(p.sample_env[sample]['CLINICAL_SNV_VCF'])

def merge_vcfs():
    '''
        merge VCFs with bcftools
    '''
    for sample in p.sample_env:

        vcf1 = p.sample_env[sample]['READY_SV_VCF']
        vcf2 = p.sample_env[sample]['READY_SNV_VCF']
        vcf3 = p.sample_env[sample]['READY_CNA_VCF']  
        msg = " INFO: Merging " + sample + " vcfs"
        print(msg)
        logging.info(msg)

        merged_vcf = p.sample_env[sample]['VCF_FOLDER'] + "/" + sample + ".merged.variants.vcf"
        # Compress vcf and index
        if not '.gz' in vcf1:
            vcf1 = u.compress_vcf(vcf1)
            u.index_vcf(vcf1)
        else:
            vcf1 = vcf1 + ".gz"
        if not '.gz' in vcf2:
            vcf2 = u.compress_vcf(vcf2)
            u.index_vcf(vcf2)
        else:
            vcf2 = vcf2 + ".gz"
        if not '.gz' in vcf3:
            vcf3 = u.compress_vcf(vcf3)
            u.index_vcf(vcf3)
        else:
            vcf3 = vcf3 + ".gz"
        bashCommand = ('{} merge {} {} {} --force-samples > {}').format(p.system_env['BCFTOOLS'], vcf1, vcf2, vcf3, merged_vcf)

        if not os.path.isfile(merged_vcf):
            msg = " INFO: " + bashCommand
            logging.info(msg)
            p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output = p1.stdout.decode('UTF-8')
            error  = p1.stderr.decode('UTF-8')
            if not error:
                pass
            else:
                if re.search("Failed", error):
                  print(error)
                  msg = " ERROR: Could not merge " + vcf1 + " and" + vcf2 + " and" + vcf3
                  print(msg)
                  logging.error(msg)
                  sys.exit()
        p.sample_env[sample]['READY_MERGED_VCF'] = merged_vcf
        p.sample_env[sample]['READY_MERGED_JSON'] = merged_vcf.replace(".vcf", ".json")
        p.sample_env[sample]['READY_MERGED_VCF_NAME'] = sample + ".merged.variants.vcf"
        u.convert_vcf_2_json(p.sample_env[sample]['READY_MERGED_VCF'])

class Civic:
  '''civic class
  '''

  def __init__(self):
    # Defining attributes
    self.variants_dict = defaultdict(dict)
    self.evidence_dict = defaultdict(dict)

  def queryCivic(self, gene, variant, vartype, exon, conseq):
    '''
      return list of evidences given a gene, variant, exon and conseq
    '''
    ev_list = []
    candidate_var = variant
    candidate_var_list = []

    if vartype == "CNA":
      if variant == "DEL":
        candidate_var_list.append("DELETION")
        candidate_var_list.append("LOSS")
      if variant == "DUP":
        candidate_var_list.append("AMPLIFICATION")
        candidate_var_list.append("GAIN")
    else:
      candidate_var_list.append(variant)
      if not 'synonymous_variant' in conseq:
        if 'splice' in conseq:
          var = "EXON "  + str(exon) + " SKIPPING MUTATION"
          candidate_var_list.append(var)
        if vartype == "SNV" or vartype == "MNV":
          candidate_var_list.append(variant[:-1])
          candidate_var_list.append(variant[1:])
          candidate_var_list.append(variant[1:-1])    
          candidate_var_list.append(variant[:-1]+"X")

    tmp_var = variant.split("_")
    if vartype == 'DELETION':
      v = "EXON " + str(exon) + " DELETION"
      candidate_var_list.append(v)
      v = variant.replace("del", "")
      candidate_var_list.append(v)
      if len(tmp_var) > 1:          
        candidate_var_list.append(tmp_var[0]+"DEL")
      v = "EXON " + str(exon) + " MUTATION"
      candidate_var_list.append(v)        
    if vartype == 'INSERTION':
      v ="EXON " + str(exon) + " INSERTION"
      candidate_var_list.append(v)
      if len(tmp_var) > 1:
        v = tmp_var[0]+"INS"    
        candidate_var_list.append(v)
        v = variant.replace("ins", "")  
        candidate_var_list.append(v)
        v = variant.replace("dup", "")  
        candidate_var_list.append(v)  
      v = "EXON " + str(exon) + " MUTATION"
      candidate_var_list.append(v)

    chosen_var_id = '.'
    if gene in self.variants_dict:
      # Now look for all variants
      for var in candidate_var_list:
        for mut in self.variants_dict[gene]:
          mut_uc = mut.upper()
          var = var.upper()
          if vartype == "CNA":
            if var == mut_uc:
              chosen_var_id = self.variants_dict[gene][mut]
              break
          else:
            if var == mut_uc:
              chosen_var_id = self.variants_dict[gene][mut]
              print("FOUND " + var + " " + mut)
              break
        if chosen_var_id != '.':
          break
      if chosen_var_id == '.':
        return ev_list

      # Get all evidences
      num = 0
      for ev_id in self.evidence_dict[chosen_var_id]:
        ev_direction   =  self.evidence_dict[chosen_var_id][ev_id]['evidence_direction']
        ev_level       =  self.evidence_dict[chosen_var_id][ev_id]['evidence_level']
        ev_significance=  self.evidence_dict[chosen_var_id][ev_id]['clinical_significance']
        ev_drugs       =  self.evidence_dict[chosen_var_id][ev_id]['drugs']
        ev_disease     =  self.evidence_dict[chosen_var_id][ev_id]['disease']
        ev_pmid        =  self.evidence_dict[chosen_var_id][ev_id]['pmid']
        ev_clintrials  =  self.evidence_dict[chosen_var_id][ev_id]['clinical_trials']

        ev_info = []
        if ev_id is not None:
          ev_info.append(ev_id)
        if ev_direction is not None:
          ev_info.append(ev_direction)
        if ev_level is not None:
          ev_info.append(ev_level)
        if ev_significance is not None:
          ev_info.append(ev_significance)
        if ev_drugs is not None:
          ev_info.append(ev_drugs)
        if ev_disease is not None:
          ev_info.append(ev_disease)
        if ev_pmid is not None:
          ev_info.append(ev_pmid)
        if ev_clintrials is not None:
          ev_info.append(ev_clintrials)
        ev_str = '|'.join(ev_info)
        if num <= 20:
          ev_list.append(ev_str)
        num+=1
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
      if evidence['evidence_direction'] is None:
        continue
      if not variant_id in self.evidence_dict:
        self.evidence_dict[variant_id] = defaultdict(dict)

      ev_id = evidence['name']
      self.evidence_dict[variant_id][ev_id]['ev_id'] = evidence['name']
      self.evidence_dict[variant_id][ev_id]['evidence_direction'] = evidence['evidence_direction']
      self.evidence_dict[variant_id][ev_id]['evidence_level'] = evidence['evidence_level']
      self.evidence_dict[variant_id][ev_id]['clinical_significance'] = evidence['clinical_significance']
      drugs_list = []
      for drug in evidence['drugs']:
        drugs_list.append(drug['name'])
      drug_set = set(drugs_list)
      drugs_list = list(drug_set)
      self.evidence_dict[variant_id][ev_id]['drugs']  = '&'.join(drugs_list) 
      self.evidence_dict[variant_id][ev_id]['disease'] = evidence['disease']['name'].replace(" ", "_")
      self.evidence_dict[variant_id][ev_id]['pmid'] = evidence['source']['citation_id']
      clintrials_list = []
      for trial in evidence['source']['clinical_trials']:
        clintrials_list.append(trial['nct_id'])
      self.evidence_dict[variant_id][ev_id]['clinical_trials'] = '&'.join(clintrials_list)

def do_civic():

    proceed = True

    for sample in p.sample_env:
        p.sample_env[sample]['CIVIC_VCF'] = \
            p.sample_env[sample]['READY_SNV_VCF'].replace(".vcf", ".civic.vcf")
        p.sample_env[sample]['CIVIC_VCF_NAME'] = \
            os.path.basename(p.sample_env[sample]['READY_SNV_VCF_NAME']).replace(".vcf", ".civic.vcf")
        if os.path.isfile(p.sample_env[sample]['CIVIC_VCF']):
            p.sample_env[sample]['READY_SNV_VCF'] = p.sample_env[sample]['CIVIC_VCF']
            p.sample_env[sample]['READY_SNV_JSON'] = p.sample_env[sample]['CIVIC_VCF'].replace(".vcf", ".json")
            p.sample_env[sample]['READY_SNV_VCF_NAME'] = p.sample_env[sample]['CIVIC_VCF_NAME']  
            proceed = False
        else:
          proceed = True
          break

    if proceed == False:
        msg = " INFO: Skipping CIViC annotation"
        print(msg)
        logging.info(msg)
        return

    msg = " INFO: Loading CIViC database"
    print(msg)
    logging.info(msg)

    # Instantiante CIViC class
    civic = Civic()

    # In-memory loading of CIViC database
    civic.loadCivic()

    for sample in p.sample_env:

        msg = " INFO: Annotating sample "+ sample + " with CIViC database"
        print(msg)
        logging.info(msg)

        p.sample_env[sample]['CIVIC_VCF'] = \
            p.sample_env[sample]['READY_SNV_VCF'].replace(".vcf", ".civic.vcf")
        p.sample_env[sample]['CIVIC_VCF_NAME'] = \
            os.path.basename(p.sample_env[sample]['READY_SNV_VCF_NAME']).replace(".vcf", ".civic.vcf")

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
                    civic_ann_list = []
                    civic_annotation = '.'
                    for item in info_list:
                        if item.startswith('CSQ'):
                            tmp_transcript = item.split(",")
                            for transcript_info in tmp_transcript:
                                transcript_info = transcript_info.replace("CSQ=", "")
                                transcript_list = transcript_info.split("|")

                                # Get features from VEP annotation
                                gene = transcript_list[vep_dict['SYMBOL']]
                                ensg_id = transcript_list[vep_dict['Gene']]
                                enst_id = transcript_list[vep_dict['Feature']]
                                hgvs_p = transcript_list[vep_dict['HGVSp']]
                                exon = re.search("\d+", transcript_list[vep_dict['EXON']])
                                if not exon:
                                  exon = re.search("\d+", transcript_list[vep_dict['INTRON']])                                
                                if exon is not None:
                                  exon = exon.group(0)
                                else:
                                  exon= '.'

                                consequence = transcript_list[vep_dict['Consequence']]
                                aminoacids = transcript_list[vep_dict['Amino_acids']]
                                protein_position = transcript_list[vep_dict['Protein_position']]
                                variant = "."

                                # Do civic query
                                if aminoacids != "" and protein_position != "":
                                  if "/" in aminoacids:
                                    variant = aminoacids.replace("/", protein_position)
                                  else:
                                    variant = aminoacids + protein_position
                                  if vartype == "INSERTION" or vartype == "DELETION" and hgvs_p != '.':

                                   # ENSP00000275493.2:p.Glu746_Ala750del
                                    tmp_p_code = hgvs_p.split(":")
                                    if len(tmp_p_code) > 1:
                                      p_code = tmp_p_code[1].replace("p.", "")
                                      for aa in u.aa_dict:
                                        if aa in p_code:
                                          p_code = p_code.replace(aa, u.aa_dict[aa])
                                      variant = p_code
                                else:
                                  if 'splice' in consequence:
                                    variant = "SKIPPING MUTATION"
                                if variant != '.':

                                  civic_ev_list = civic.queryCivic(gene, variant, vartype, exon, consequence)
                                  civic_annotation = ','.join(civic_ev_list)
                                  if len(civic_ev_list)>0:
                                    civic_ann_list.append(civic_annotation)
                            #item = tmp_transcript[tidx]
                        idx+=1

                    info = ';'.join(info_list)
                    civic_out = "."
                    if len(civic_ann_list) > 0:
                      civic_out = ','.join(civic_ann_list)
                    else:
                      civic_out = civic_annotation
                    tmp[7] = info + ";CIVIC="+civic_out

                    line = '\t'.join(tmp)
                    o.write(line+"\n")
        f.close()
        o.close()
        p.sample_env[sample]['READY_SNV_VCF'] = p.sample_env[sample]['CIVIC_VCF']
        p.sample_env[sample]['READY_SNV_JSON'] = p.sample_env[sample]['CIVIC_VCF'].replace(".vcf", ".json")
        p.sample_env[sample]['READY_SNV_VCF_NAME'] = p.sample_env[sample]['CIVIC_VCF_NAME']     
        
        # Create JSON file
        u.convert_vcf_2_json(p.sample_env[sample]['READY_SNV_VCF'])

def do_vep():
    '''
    VEP annotation of vcf files
    '''

    dbnsfp_dict = {
      "sift" : "SIFT_score,SIFT_pred",
      "polyphen2" : "Polyphen2_HDIV_score,Polyphen2_HDIV_pred",
      "mutationtaster2" : "MutationTaster_score,MutationTaster_pred",
      "provean" : "PROVEAN_score,PROVEAN_pred",
      "fathmm"  : "FATHMM_score,FATHMM_pred",
      "revel"   : "REVEL_score",
      "mutpred" : "MutPred_score",
    } 

    for sample in p.sample_env:

        p.sample_env[sample]['VEP_VCF'] = \
            p.sample_env[sample]['READY_SNV_VCF'].replace(".vcf", ".vep.vcf")

        p.sample_env[sample]['VEP_VCF_NAME'] = \
            os.path.basename(p.sample_env[sample]['READY_SNV_VCF_NAME']).replace(".vcf", ".vep.vcf")

        dbnsfp_ann    = ''
        gnomad_ann    = ''
        ncer_ann      = ''
        cadd_ann      = ''
        spliceai_ann  = ''
        maxent_ann    = ''
        phylop_ann    = ''
        phastcons_ann = ''

        # Missense predictors from dbNSFP 
        if p.analysis_env['MISSENSE_PREDICTORS']:
          missense_list = p.analysis_env['MISSENSE_PREDICTORS'].split(',')
          tmp = [] 
          for predictor in missense_list:
            tmp.append(dbnsfp_dict[predictor])
          dbnsfp_ann = " --plugin dbNSFP,/anndir/dbNSFP/" + \
            p.aux_env['DBNSFP_FILENAME'] + "," + ','.join(tmp)
        # Setting splicing predictors 
        if p.analysis_env['SPLICING_PREDICTORS']:
          splicing_list = p.analysis_env['SPLICING_PREDICTORS'].split(',')
          for predictor in splicing_list:
            #  Including MaxEntScan
            if predictor == "maxentscan":
              maxent_ann = " --plugin MaxEntScan,/anndir/MaxEntScan"
            #  Including SpliceAI for SNV and Indels
            if predictor == "spliceai":
              spliceai_ann = " --plugin SpliceAI,snv=/anndir/spliceAI/" + \
                p.aux_env['SPLICEAI_SNV_FILENAME'] + ",indel=/anndir/spliceAI/" + \
                p.aux_env['SPLICEAI_INDEL_FILENAME']
        # Setting genomewide predictors 
        if p.analysis_env['GENOMEWIDE_PREDICTORS']:
          pass
        if p.analysis_env['GNOMAD'] == True:
          gnomad_fields = "AC,AN,AF,rf_tp_probability,AC_afr,AN_afr,AF_afr,AC_eas,AF_eas,AC_nfe,AN_nfe,AF_nfe,AC_fin,AN_fin,AF_fin,AC_asj,AN_asj,AF_asj,AC_oth,AN_oth,AF_oth,popmax,AC_popmax,AN_popmax,AF_popmax"
          gnomad_ann = ' --custom /anndir/gnomAD/' + p.aux_env['GNOMAD_FILENAME'] + ",GNOMAD,vcf,exact,0" + gnomad_fields

        # Setting conservation scores
        if p.analysis_env['CONSERVATION_SCORES']:
          cons_list = p.analysis_env['CONSERVATION_SCORES'].split(',')
          for conservation in cons_list:
            if conservation == 'phastcons':
              phastcons_ann = ' --custom /anndir/Conservation/' + p.aux_env['PHASTCONS_FILENAME'] + ",Phastcons" + ',bigwig,overlap'
            if conservation == 'phylop':
              phylop_ann = ' --custom /anndir/Conservation/' + p.aux_env['PHYLOP_FILENAME'] + ",PhyloP" + ',bigwig,overlap'

       # --plugin Conservation,/path/to/bigwigfile.bw
       # --plugin SpliceAI,snv=/path/to/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/path/to/spliceai_scores.raw.indel.hg38.vcf.gz
       # --plugin MaxEntScan,/path/to/maxentscan/fordownload,SWA,NCSS
        vep_output_vcf = p.aux_env['VEP_FOLDER_OUTPUT'] + "/" + p.sample_env[sample]['VEP_VCF_NAME']
        dbnsfp_fields = ',Ensembl_transcriptid,gnomAD_genomes_POPMAX_AC,gnomAD_genomes_POPMAX_AN,gnomAD_genomes_flag,gnomAD_genomes_POPMAX_AF,GERP++_RS,phyloP100way_vertebrate'
        bashCommand = ('{} run -t -i -v {}:/genomedir/ -v {}:/anndir/ -v {}:/opt/vep/.vep {}'
        ' perl vep --cache --offline --dir_cache /opt/vep/.vep/ --dir_plugins /opt/vep/.vep/Plugins/'
        ' --input_file /opt/vep/.vep/input/{} --output_file /opt/vep/.vep/output/{} '
        ' --af_1kg --af_gnomad --cache_version 101  --canonical '
        ' --format vcf --vcf --hgvs --hgvsg --max_af --pubmed --gene_phenotype --ccds --sift b --polyphen b --symbol --force_overwrite --fork {}'
        ' {} {} {} {} {} {} --fasta /genomedir/{} '
        .format(p.system_env['DOCKER'], p.aux_env['GENOME_FOLDER'], p.analysis_env['ANN_DIR'], p.aux_env['VEP_FOLDER'], p.docker_env['VEP'],\
        p.sample_env[sample]['READY_SNV_VCF_NAME'], p.sample_env[sample]['VEP_VCF_NAME'], p.analysis_env['THREADS'],\
        dbnsfp_ann, spliceai_ann, maxent_ann, phastcons_ann, phylop_ann, gnomad_ann, p.aux_env['GENOME_NAME'] ))

        if not os.path.isfile(p.sample_env[sample]['VEP_VCF']):

            # First, copy vcf to vep input dir
            shutil.copy2(p.sample_env[sample]['READY_SNV_VCF'], p.aux_env['VEP_FOLDER_INPUT']+"/"+p.sample_env[sample]['READY_SNV_VCF_NAME'])

            msg = " INFO: Annotating sample " + sample + " with VEP"
            print (msg)
            logging.info(msg)
            logging.info(bashCommand)

            p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output = p1.stdout.decode('UTF-8')
            error  = p1.stderr.decode('UTF-8')
            if not error:
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
        p.sample_env[sample]['READY_SNV_JSON'] = p.sample_env[sample]['VEP_VCF'].replace(".vcf", ".json")
        p.sample_env[sample]['READY_SNV_VCF_NAME'] = p.sample_env[sample]['VEP_VCF_NAME']     

        # Create JSON file
        u.convert_vcf_2_json(p.sample_env[sample]['READY_SNV_VCF'])

def load_cgi_data():
  '''
  Load Cancer Genome Interpreter
  '''
  cgi_dict = defaultdict(dict)

  with open (p.aux_env['CGI_BIOMARKERS'], 'r') as f:
    for line in f:
      line = line.rstrip("\n")
      tmp = line.split("\t")
      if line.startswith("Alteration"):
        continue
      else:
        tmp_alteration = tmp[0].split(":")
        if len(tmp_alteration) < 2:
          continue
        gene = tmp_alteration[0]
        biomarkers = tmp_alteration[1].split(",")

        if gene not in cgi_dict:
          cgi_dict[gene] = defaultdict(dict)

        for biomarker in biomarkers:
          cgi_dict[gene][biomarker]['association']= tmp[3]
          cgi_dict[gene][biomarker]['drug']   = tmp[8]
          cgi_dict[gene][biomarker]['p_code'] = biomarker
          cgi_dict[gene][biomarker]['c_code'] = tmp[21]
          cgi_dict[gene][biomarker]['g_code'] = tmp[22]
          cgi_dict[gene][biomarker]['pmid']   = tmp[17].replace("PMID:", "").replace(";", "&")
          cgi_dict[gene][biomarker]['tumor_type'] = tmp[27]
  return cgi_dict

def query_cgi_data(cgi_dict, gene, p_code):

  cgi_list = []
  if gene in cgi_dict:

    for biomarker in cgi_dict[gene]:
      if biomarker == p_code.upper():
        cgi_info = []
        cgi_association = cgi_dict[gene][biomarker]['association']
        cgi_drug  = cgi_dict[gene][biomarker]['drug']
        cgi_pcode = cgi_dict[gene][biomarker]['p_code']
        cgi_ccode = cgi_dict[gene][biomarker]['c_code']
        cgi_gcode = cgi_dict[gene][biomarker]['g_code']
        cgi_pmid  = cgi_dict[gene][biomarker]['pmid']
        cgi_tumor_type= cgi_dict[gene][biomarker]['tumor_type']
        cgi_info.append(cgi_association)
        cgi_info.append(cgi_drug)
        cgi_info.append(cgi_pcode)
        cgi_info.append(cgi_ccode)
        cgi_info.append(cgi_gcode)
        cgi_info.append(cgi_pmid)
        cgi_info.append(cgi_tumor_type)
        cgi_str = '|'.join(cgi_info)
        cgi_list.append(cgi_str)
  return cgi_list

def do_cgi():

    msg = " INFO: Loading Cancer Genome Interpreter database"
    logging.info(msg)
    print(msg)

    cgi_dict = load_cgi_data()
    if not cgi_dict:
        msg = " ERROR: CGI data could not be loaded"
        logging.error(msg)
        print(msg)

    for sample in p.sample_env:

        msg = " INFO: Annotating sample "+ sample + " with CGI database"
        print(msg)
        logging.info(msg)

        p.sample_env[sample]['CGI_VCF'] = \
            p.sample_env[sample]['READY_SNV_VCF'].replace(".vcf", ".cgi.vcf")
        p.sample_env[sample]['CGI_VCF_NAME'] = \
            os.path.basename(p.sample_env[sample]['READY_SNV_VCF_NAME']).replace(".vcf", ".cgi.vcf")

        vep_dict = defaultdict(dict)
        vep_list = []
        o = open(p.sample_env[sample]['CGI_VCF'], 'w')

        cgi_fields = []
        cgi_fields.append("ASSOCIATION")
        cgi_fields.append("DRUG")
        cgi_fields.append("PCODE")
        cgi_fields.append("CCODE")
        cgi_fields.append("GCODE")
        cgi_fields.append("PMID")
        cgi_fields.append("TUMOR_TYPE")
        cgi_info_header = "##INFO=<ID=CGI,Number=.,Type=String,Description=\"CGI evidence. Format: " + '|'.join(cgi_fields) + "\">"
        with open (p.sample_env[sample]['READY_SNV_VCF']) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    if re.search("ID=CSQ", line):
                        vep_dict, vep_list = u.create_vep_dict(line)
                        o.write(line+"\n")
                        o.write(cgi_info_header+"\n")
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
                    cgi_ann_list = []
                    cgi_annotation = '.'
                    for item in info_list:
                        if item.startswith('CSQ'):
                            tmp_transcript = item.split(",")
                            for transcript_info in tmp_transcript:
                                transcript_info = transcript_info.replace("CSQ=", "")
                                transcript_list = transcript_info.split("|")

                                # Get features from VEP annotation
                                gene = transcript_list[vep_dict['SYMBOL']]
                                ensg_id = transcript_list[vep_dict['Gene']]
                                enst_id = transcript_list[vep_dict['Feature']]
                                exon = re.search("\d+", transcript_list[vep_dict['EXON']])
                                consequence = transcript_list[vep_dict['Consequence']]
                                aminoacids = transcript_list[vep_dict['Amino_acids']]
                                protein_position = transcript_list[vep_dict['Protein_position']]
                                variant = "."

                                # Do cgi query
                                if aminoacids != "" and protein_position != "":
                                  if "/" in aminoacids:
                                    variant = aminoacids.replace("/", protein_position)
                                  else:
                                    variant = aminoacids + protein_position
                                  cgi_ev_list = query_cgi_data(cgi_dict, gene, variant)
                                  cgi_annotation = ','.join(cgi_ev_list)
                                  if len(cgi_ev_list)>0:
                                    cgi_ann_list.append(cgi_annotation)

                        idx+=1

                    info = ';'.join(info_list)
                    cgi_ann_list = list(set(cgi_ann_list))
                    cgi_annotation = ','.join(cgi_ann_list)

                    if not cgi_ann_list:
                      cgi_annotation = "."

                    tmp[7] = info + ";CGI="+cgi_annotation
                    line = '\t'.join(tmp)
                    o.write(line+"\n")
        f.close()
        o.close()  
        p.sample_env[sample]['READY_SNV_VCF'] = p.sample_env[sample]['CGI_VCF']
        p.sample_env[sample]['READY_SNV_JSON'] = p.sample_env[sample]['CIVIC_VCF'].replace(".vcf", ".json")
        p.sample_env[sample]['READY_SNV_VCF_NAME'] = p.sample_env[sample]['CGI_VCF_NAME']

        # Create JSON file
        u.convert_vcf_2_json(p.sample_env[sample]['READY_SNV_VCF'])

def annotate_known_fusions():

    for sample in p.sample_env:

        fusions_bed = u.vcf_2_bed(p.sample_env[sample]['READY_SV_VCF'])
        intersect_file  = p.sample_env[sample]['VCF_FOLDER'] + "/" + "intersect.bed"
        fusions_vcf = p.sample_env[sample]['READY_SV_VCF'].replace(".vcf", ".fusions.vcf")
        
        if os.path.isfile(fusions_vcf):
          msg = " INFO: Skipping fusion annotation for sample "+ sample
          print(msg)
          logging.info(msg)
          p.sample_env[sample]['READY_SV_VCF'] = fusions_vcf
          u.convert_vcf_2_json(fusions_vcf)
          continue

        o = open(fusions_vcf, 'w')

        bashCommand = ('{} pairtopair -a {} -b {} | sort -V | uniq > {}')\
        .format(p.system_env['BEDTOOLS'], fusions_bed, p.aux_env['CHIMERKB_BED'], intersect_file)
            
        p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = p1.stdout.decode('UTF-8')
        error  = p1.stderr.decode('UTF-8')
        if not error and not output:
            msg = " INFO: Fusion annotation for sample "+ sample +" ended successfully"
            print(msg)
            logging.info(msg)
        else:
            msg = " ERROR: Something went wrong with Fusion annotation for sample " + sample 
            print(msg)
            logging.error(msg)

        fusions_dict = defaultdict(dict)
        with open (intersect_file, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                tmp = line.split("\t")
                info = tmp[13]
                tmp_info = info.split(";")
                
                coordinate = tmp[0]+"\t"+tmp[1]
                fusion_pair = "."
                disease     = "."
                source      = "."
                if len(tmp_info) > 1:
                    fusion_pair = tmp_info[0].replace(" ", "")
                    source      = tmp_info[1]
                    disease     = tmp_info[2]
                    if not coordinate in fusions_dict:
                        fusions_dict[coordinate] = defaultdict(dict)
                        fusions_dict[coordinate]['PARTNERS'] = []
                        fusions_dict[coordinate]['SOURCES']  = []
                        fusions_dict[coordinate]['DISEASES'] = []
                        fusions_dict[coordinate]['PARTNERS'].append(fusion_pair)
                        fusions_dict[coordinate]['SOURCES'].append(source)
                        fusions_dict[coordinate]['DISEASES'].append(disease)
        f.close()

        fusion_fields = []
        fusion_fields.append("PARTNERS")
        fusion_fields.append("SOURCES")
        fusion_fields.append("DISEASES")
        fusion_info_header = "##INFO=<ID=FUSION,Number=.,Type=String,Description=\"Fusion evidence. Format:" + '|'.join(fusion_fields) + "\">"

        with open (p.sample_env[sample]['READY_SV_VCF']) as f:
            for line in f:
                line = line.rstrip("\n")
                tmp = line.split("\t")
                if line.startswith("#"):
                    if re.search("cmdline", line):
                        o.write(line+"\n")
                        o.write(fusion_info_header+"\n")
                    else:
                        o.write(line+"\n")
                    continue
                coordinate = tmp[0]+"\t"+tmp[1]
                if coordinate in fusions_dict:
                    fusion_ann_list = []
                    unique_partners = '&'.join(list(set(fusions_dict[coordinate]['PARTNERS'])))
                    unique_sources  = '&'.join(list(set(fusions_dict[coordinate]['SOURCES'])))
                    unique_diseases = '&'.join(list(set(fusions_dict[coordinate]['DISEASES'])))
                    fusion_ann_list.append(unique_partners)
                    fusion_ann_list.append(unique_sources)
                    fusion_ann_list.append(unique_diseases)
                    tmp[7]+=";FUSION=" + '|'.join(fusion_ann_list)
                else:
                    tmp[7]+=";FUSION=||"
                line  = '\t'.join(tmp)
                o.write(line+"\n")
        f.close()

        p.sample_env[sample]['READY_SV_VCF'] = fusions_vcf
        p.sample_env[sample]['READY_SV_JSON'] = fusions_vcf.replace(".vcf", "json")
        p.sample_env[sample]['READY_SV_VCF_NAME'] = os.path.basename(fusions_vcf)
        # Create JSON file
        u.convert_vcf_2_json(p.sample_env[sample]['READY_SV_VCF'])
