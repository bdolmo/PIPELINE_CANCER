#!/usr/bin/env python3

import os
import sys
import re
import logging
from collections import defaultdict
from pathlib import Path
import subprocess
from modules import params as p
from modules import utils as u
from modules import fastq as f

def trim_fastqs():
  '''Trim fastq files
  '''

  # load fastq files
  fastq_files = u.get_fastq_files()

  # Check fastq files are well formated
  seen_fq = defaultdict(dict)
  for fq in fastq_files:

    fq1 = ""
    fq2 = ""
    if '_R1_' in fq:
      fq1 = fq
      fq2 = fq1.replace("R1", "R2")
    else :
      fq2 = fq
      fq1 = fq2.replace("R2", "R1")
      
    fq_tmp = os.path.basename(fq1).split("_")
    sample_name  = fq_tmp[0]

    if not sample_name in seen_fq:
      seen_fq[sample_name] = 0
    else:
      seen_fq[sample_name] += 1
    if seen_fq[sample_name] > 1:
      continue

    f_pair = f.Fastq(fq1, fq2, True)

    if f_pair.isValid():
      pass
    else:
      msg = " ERROR: Fastq pairs " +  fq1 + " " + fq2 + " are invalid"
      logging.error(msg)
      sys.exit()
  
    fq_tmp = os.path.basename(fq1).split("_")
    sample_name  = fq_tmp[0]

    if not sample_name in seen_fq:
      seen_fq[sample_name] = 1
    else:
      seen_fq[sample_name] += 1
    if seen_fq[sample_name] > 1:
      continue

    p.sample_env[sample_name]['SAMPLE_FOLDER'] = \
      p.analysis_env['OUTPUT_DIR'] + "/" + sample_name

    sample_path = Path(p.sample_env[sample_name]['SAMPLE_FOLDER'])
    if not sample_path.is_dir():
      os.mkdir(sample_path)

    # Creating FASTQ directory
    p.sample_env[sample_name]['FASTQ_FOLDER'] \
      = p.sample_env[sample_name]['SAMPLE_FOLDER'] + "/" + "FASTQ_FOLDER"
    sample_fastq_path = Path(p.sample_env[sample_name]['FASTQ_FOLDER'])
    if not sample_fastq_path.is_dir():
      os.mkdir(sample_fastq_path)

    trimmed_fq1 = os.path.basename(fq1).replace(".fastq.gz", ".trimmed.fastq.gz")
    trimmed_fq2 = os.path.basename(fq2).replace(".fastq.gz", ".trimmed.fastq.gz")

    p.sample_env[sample_name]['READY_FQ1'] = \
      p.sample_env[sample_name]['FASTQ_FOLDER'] + "/" + trimmed_fq1
    p.sample_env[sample_name]['READY_FQ2'] = \
      p.sample_env[sample_name]['FASTQ_FOLDER'] + "/" + trimmed_fq2

    # Now trimming with fastp
    bashCommand = ('{} -i  {} -I {} -o {} -O {} -w {}') \
      .format(p.system_env['FASTP'], fq1, fq2, p.sample_env[sample_name]['READY_FQ1'], 
      p.sample_env[sample_name]['READY_FQ2'], p.analysis_env['THREADS'] )

    if not os.path.isfile(p.sample_env[sample_name]['READY_FQ1']) \
      and not os.path.isfile(p.sample_env[sample_name]['READY_FQ2']):

      msg = " INFO: Trimming sample " +  sample_name 
      print(msg)
      p.logging.info(msg)
      p.logging.info(bashCommand)

      process = subprocess.Popen(bashCommand, shell=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
      
      output, error = process.communicate()
    else:
      msg = " INFO: Skipping trimming of sample " +  sample_name 
      print(msg)
      p.logging.info(msg)