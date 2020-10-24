#!/usr/bin/env python3

import sys
import re
import gzip
import binascii
import os.path
import glob
from os import path
from collections import defaultdict
from pprint import pprint
import json

from pathlib import Path
import subprocess
from modules import params as p


def vcf_2_bed(vcf):

  bed = vcf.replace(".vcf", ".bed")
  o.open(bed, 'w')
  with open (vcf, "r") as f:
    for line in f:
      line = line.rstrip("\n")
      if line.startswith("#"):
        continue
      else:
        tmp = line.split("\t")
        #chr6	128129334	MantaBND:84:0:1:0:0:0:0	G	[chr6:117660945[G
        info = tmp[7]
        chr_A = tmp[0]
        pos_A = tmp[1]
        end_A = int(pos_A)+1

        chr_B = "."
        pos_B = "."

        info_list = info.split(";")
        for field in info_list:
          if field.startswith("END="):
            chr_B = chr_A
            pos_B = field.replace("END=", "")
        if pos_B == ".":
          m = re.search('chr(\d+)+:(\d+)', tmp[4])
          if m:
            tmp_posB = m.split(":")
            chr_B = tmp_posB[0]
            pos_B = tmp_posB[1]
        out_list = []
        out_list.append(str(chr_A))
        out_list.append(str(pos_A))
        out_list.append(str(end_A))
        out_list.append(str(chr_B))
        out_list.append(str(pos_B))
        out_list.append(str(end_B))
        out_list.append(str(info))
        o.write("\t".join(out_list))
  o.close()


def create_vep_dict(info):

    vep_dict = dict()
    vep_list = list()
    if info:
        tmp = info.split('Format:')

        # remove unwanted characters
        rawfields = re.sub('(">)', "", tmp[1])
        fields = rawfields.split('|')
        i = 0
        for field in fields:
            vep_list.append(field)
            vep_dict[field] = i
            i += 1

    return vep_dict, vep_list

def convert_vcf_2_json(vcf):

    out_json = vcf.replace(".vcf", ".json")

    vcf_dict =  defaultdict(dict)
    vcf_dict['header'] =  defaultdict(list)
    vcf_dict['variants'] =  defaultdict(list)

    ##fileformat=VCFv4.2
    ##fileDate=20200418
    ##source=Mutect2,freeBayes
    ##reference=/home/bdelolmo/Desktop/PIPELINE/BUNDLE/ucsc.hg19.fasta
    n_var = 0
    with open (vcf) as f:
        for line in f:
            line = line.rstrip('\n')
            # Parsing header
            if line.startswith('#'):
               # print (line)
                tmp = line.split('=')
                if re.search("##fileformat", line):
                    vcf_dict['header']['fileFormat'] = tmp[1]
                if re.search("##fileDate", line):
                    vcf_dict['header']['fileDate'] = tmp[1]
                if re.search("##source", line):
                    vcf_dict['header']['source'] = tmp[1]
                if re.search("##reference", line):
                    vcf_dict['header']['reference'] = tmp[1]
                if re.search("##contig", line):
                    ###contig=<ID=chr10,length=135534747>
                    line = line.replace("##contig=<ID=", "")
                    line = line.replace("length=", "")
                    line = line.replace(">","")
                    tmp = line.split(",")
                    contig_dict = defaultdict(dict)
                    contig_dict['ID'] = tmp[0]
                    contig_dict['length'] = tmp[1]
                    vcf_dict['header']['contig'].append(contig_dict)

                if re.search("ID=CSQ", line):
                    vep_dict, vep_list = create_vep_dict(line)
                if re.search("ID=CIVIC", line):
                    civic_dict, civic_list = create_vep_dict(line)
                if re.search('##FORMAT', line):
                    line = line.replace("##FORMAT=<ID=", "")
                    line = line.replace(">","")
                    tmp = line.split(',')
                    format_dict = defaultdict(dict)
                    format_id = ''
                    n = 0
                    for item in tmp:
                        tmp2= item.split('=')
                        if n == 0:
                            format_id = item
                        if len(tmp2) >1:
                            format_dict[format_id][tmp2[0]] = tmp2[1]
                        n=n+1
                    vcf_dict['header']['FORMAT'].append(format_dict)
                continue
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

                var_name = "var_" + str(n_var)
                vcf_dict['variants'][var_name] = defaultdict(dict)
                vcf_dict['variants'][var_name]['CHROM'] = chr
                vcf_dict['variants'][var_name]['POS'] = pos
                vcf_dict['variants'][var_name]['ID'] = id
                vcf_dict['variants'][var_name]['REF'] = ref
                vcf_dict['variants'][var_name]['ALT'] = alt
                vcf_dict['variants'][var_name]['QUAL'] = qual
                vcf_dict['variants'][var_name]['FILTER'] = filter

                info_list = info.split(';')
                vcf_dict['variants'][var_name]['INFO'] = defaultdict(dict)
                for item in info_list:
                    tmp_item = item.split('=')
                    if item.startswith('CSQ'):
                        vcf_dict['variants'][var_name]['INFO']['CSQ'] = defaultdict(dict)
                        tmp_multidim = item.split(',')
                        n = 0
                        for subitem in tmp_multidim:
                            subitem = subitem.replace("CSQ=", "")
                            conseq_name = "consequence_"+str(n)
                            tmp_subfield = subitem.split('|')
                            num_field = 0
                            for subfield in tmp_subfield:
                                if subfield == "":
                                    subfield = '.'
                                subitem_name = vep_list[num_field]
                                #print (subfield + " " + subitem_name + " " + str(num_field))
                                vcf_dict['variants'][var_name]['INFO']['CSQ'][subitem_name] = subfield
                                num_field=num_field+1
                            n =n+1
                    elif item.startswith('CIVIC'):
                        tmp_multidim = item.split(',')
                        for subitem in tmp_multidim:
                            civic_evidence = "Evidence_"+str(n)
                            tmp_subfield = subitem.split('|')
                            num_field = 0
                            for subfield in tmp_subfield:
                                if subfield == "":
                                    subfield = '.'
                                subitem_name = civic_list[num_field]
                                vcf_dict['variants'][var_name]['INFO']['CIVIC_EVIDENCE'][subitem_name] = subfield
                                num_field=num_field+1
                            n =n+1
                    else:
                        vcf_dict['variants'][var_name]['INFO'][tmp_item[0]] = tmp_item[1]
 
                tmp_format_tag = format_tag.split(':')
                tmp_format = format.split(':')
                j = 0
                for val in tmp_format:
                    vcf_dict['variants'][var_name][tmp_format_tag[j]] = val
                    j=j+1

                n_var = n_var+1
                #chr1	65301110	.	C	T	.	.

    with open(out_json, 'w') as fp:
        json.dump(vcf_dict, fp)



def decompress_vcf(input_vcf):
  '''
    Decompress a bgzipped vcf file
  '''
  output_vcf = input_vcf.replace(".gz.vcf", ".vcf")
  bashCommand = ('{} -c {} > {}').format(p.system_env['GUNZIP'], input_vcf, output_vcf)
  if not os.path.isfile(output_vcf):

    process = subprocess.Popen(bashCommand, shell=True, stdout=subprocess.PIPE,
      stderr=subprocess.PIPE)
    output, error = process.communicate()
    if not error.decode('UTF-8'):
      pass
    else:
      msg = " ERROR: Could not create a gzip for vcf " + input_vcf
      print(msg)
      logging.error(msg)

def index_vcf(input_vcf):
  '''
    Index a bgzipped vcf
  '''
  output_vcf_tbi = input_vcf + ".tbi"
  bashCommand = ('){} -p vcf {}').format(p.system_env['TABIX'], input_vcf)
  if not os.path.isfile(output_vcf_tbi):

    process = subprocess.Popen(bashCommand, shell=True, stdout=subprocess.PIPE,
      stderr=subprocess.PIPE)
    output, error = process.communicate()
    if not error.decode('UTF-8'):
      pass
    else:
      msg = " ERROR: Could not create an index file for  for vcf " + input_vcf
      print(msg)
      logging.error(msg)

def compress_vcf(input_vcf):
  '''
    Compress a plain text vcf
  '''
  output_vcf = input_vcf.replace(".vcf", ".gz.vcf")
  bashCommand = ('{} -c {} > {}').format(p.system_env['BGZIP'], input_vcf, output_vcf)
  if not os.path.isfile(output_vcf):

    process = subprocess.Popen(bashCommand, shell=True, stdout=subprocess.PIPE,
      stderr=subprocess.PIPE)
    output, error = process.communicate()
    if not error.decode('UTF-8'):
      pass
    else:
      msg = " ERROR: Could not create a gzip for vcf " + input_vcf
      print(msg)
      logging.error(msg)



def get_fastq_files():
  '''
    Get fastq files from input directory
  '''

  fastq_list = []
  fastq_list = glob.glob(p.analysis_env['INPUT_DIR'] +  "/*.fastq.gz")
  if not fastq_list:
    fastq_list = glob.glob(p.analysis_env['INPUT_DIR'] +  "/*.fa.gz")
  if not fastq_list:
    msg = " ERROR: No input fastq files were detected"
    print (msg)
    p.logging.error(msg)
    sys.exit()

  #out_fq_list = []
  #for fq in fastq_file:

  return fastq_list

def check_docker_images(command, image):
    '''simple check for docker images installed
    '''
    bashCommand = '{} image ls {}'.format(command, image)

    p.logging.info(bashCommand)
    process = subprocess.Popen(bashCommand,#.split(),
      shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    output, error = process.communicate()

    if error.decode('UTF-8') == '':
        output = output.decode('UTF-8')
        if output:
            c_lines = 0
            for line in output.split('\n'):
                c_lines += 1
            if c_lines !=3:
                msg = " ERROR: docker image for " + image +  " was not found" 
                print (msg)
                p.logging.error(msg)
            else:
                msg = " INFO: docker image for " + image +  " was found" 
                print (msg)
                p.logging.info(msg)
    else:
        msg = " ERROR: docker image for " + image +  " was not found" 
        print (msg)
        p.logging.error(msg)
