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

def do_report():

    select_clinical_variants()

def select_clinical_variants():

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
                        civic_dict, civic_list = create_vep_dict(line)
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

                                        if transcript_list[vep_dict['EV_DIRECTION']] == "Supports" \
                                            or max_af < 0.01:
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

