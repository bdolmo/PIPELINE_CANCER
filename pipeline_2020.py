#!/usr/bin/env python3

import os
import sys
import re
import argparse
import logging
from pathlib import Path

from modules import params as p
from modules import trimming as t
from modules import map as m
from modules import var_call as v
from modules import annotate as a
from modules import sqlite as s
from modules import report as r

main_dir = os.path.dirname(os.path.abspath(__file__))

def main(args):

    # Setting default directories
    p.set_defaults(main_dir)

    # Setting analysis variables
    p.set_analysis_env(args)

    # Settingg panel configuration
    p.set_panel_configuration(main_dir)

    #p.update_variant_databases()

    # Setting system binary variables
    p.set_system_env()

    # Setting aux files
    p.set_auxfiles_env()

    # Trim fastq files
    t.trim_fastqs()

    # Map fastq files, rmdup, bam qc
    m.do_all()

    # Perform var calling
    v.do_var_call()

    # Perform vcf annotation
    a.do_annotation()

    r.do_report()


def parse_arguments():
    '''parsing input arguments
    '''
    parser = argparse.ArgumentParser(description="Pipeline for NGS analysis v.1.0")
    parser.add_argument("--panel", type=str, required=True,
        help="panel to analyze options today [genes]", dest='panel')
    parser.add_argument("-r", "--reference", required=True, type=str, choices=['hg19', 'hg38'], 
        default='hg19', help="Genome reference to do mapping var calling and annotation",  dest='reference')
    parser.add_argument("-t", "--threads", type=int, default=4,
        help="Number of CPU threads to operate", dest='threads')
    parser.add_argument("--var_class", required=True, type=str, choices=['somatic', 'germline'],
        default='somatic', help="Perform somatic or germline variant analysis", dest='var_class')
    parser.add_argument("-o", "--output_dir", required=True, type=str,
        help="Output directory", dest='output_dir')
    parser.add_argument("-i", "--input_dir", required=True, type=str,
        help="Input directory", dest='input_dir')

    # Now subparsers
    subparsers = parser.add_subparsers(help='sub-command help', dest='subparser_name')
    parser_mapping = subparsers.add_parser('mapping', help='mapping pipeline')
    parser_mapping.add_argument("--rm_dup", action='store_true',
        help="Remove duplicates in bam. Default = True")
    parser_mapping.add_argument('--mapper', choices=['bwa'], default='bwa',
        help="Choose a short-read mapper. Default = bwa")
    parser_mapping.add_argument('--qc_analysis', default=True,
        help="Perform QC analysis. Default = True")

    arguments = parser.parse_args()
    return arguments


#########################################################
if __name__ == '__main__':
    ARGUM = parse_arguments()
    main(ARGUM)

