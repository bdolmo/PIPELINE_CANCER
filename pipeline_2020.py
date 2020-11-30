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

    # Get panel configuration from db
    p.get_panel_configuration(main_dir)

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

    # Create clinical report 
    r.do_report()


def parse_arguments():
    '''parsing input arguments
    '''
    parent_parser = argparse.ArgumentParser(description="Pipeline for NGS analysis v.1.0")
    parent_parser.add_argument("--panel", type=str, required=True,
        help="panel to analyze options today [genes]", dest='panel')
    parent_parser.add_argument("-r", "--reference", required=True, type=str, choices=['hg19', 'hg38'], 
        default='hg19', help="Genome reference to do mapping var calling and annotation",  dest='reference')
    parent_parser.add_argument("-t", "--threads", type=int, default=4,
        help="Number of CPU threads to operate", dest='threads')
    parent_parser.add_argument("--var_class", required=True, type=str, choices=['somatic', 'germline'],
        default='somatic', help="Perform somatic or germline variant analysis", dest='var_class')
    parent_parser.add_argument("-o", "--output_dir", required=True, type=str,
        help="Output directory", dest='output_dir')
    parent_parser.add_argument("-i", "--input_dir", required=True, type=str,
        help="Input directory", dest='input_dir')
    parent_parser.add_argument("--db_dir", type=str, required=True,
        help="Database directory harbouring sqlite files", dest='db_dir')
    parent_parser.add_argument("--ann_dir", type=str, required=True,
        help="Annotation resources directory", dest='ann_dir')
    parent_parser.add_argument("--ref_dir", type=str, required=True,
        help="Genome reference resource directory", dest='ref_dir')
    parent_parser.add_argument("--sample_data", type=str,
        help="Sample data docx", dest='sample_data')
    parent_parser.add_argument("--lab_data", type=str,
        help="Lab data xlsx", dest='lab_data')
    parent_parser.add_argument("--lang", type=str, choices=['cat', 'en', 'esp'], default='cat',
        help="Report language", dest='language')
    parent_parser.add_argument("--min-fusion-size", type=int, default=50000,
        help="Minimum fusion size in bp to be reported", dest='min_fusion_size')
  
    # Now subparsers
    subparsers = parent_parser.add_subparsers(title="sub-commands")
    parser_map = subparsers.add_parser('map', parents=[parent_parser],
        add_help=False, description="Map raw FASTQ files", help='Mapping pipeline')
    parser_map.add_argument("--rm_dup", action='store_true',
        help="Remove duplicates in bam. Default = True")
    parser_map.add_argument('--mapper', choices=['bwa'], default='bwa',
        help="Choose a short-read mapper. Default = bwa")
    parser_map.add_argument('--qc_analysis', default=True,
        help="Perform QC analysis. Default = True")

    parser_call = subparsers.add_parser('call', parents=[parent_parser],
        add_help=False, description="Call variants command", help='Calling pipeline')
    parser_call.add_argument("--rm_dup", action='store_true',
        help="Remove duplicates in bam. Default = True")
    parser_call.add_argument('--mapper', choices=['bwa'], default='bwa',
        help="Choose a short-read mapper. Default = bwa")
    parser_call.add_argument('--qc_analysis', default=True,
        help="Perform QC analysis. Default = True")

    parser_annotate = subparsers.add_parser('annotate', parents=[parent_parser],
        add_help=False, description="Anotate variants command", help='Variant annotatione pipeline')
    parser_annotate.add_argument("--rm_dup", action='store_true',
        help="Remove duplicates in bam. Default = True")
    parser_annotate.add_argument('--mapper', choices=['bwa'], default='bwa',
        help="Choose a short-read mapper. Default = bwa")
    parser_annotate.add_argument('--qc_analysis', default=True,
        help="Perform QC analysis. Default = True")

    parser_report = subparsers.add_parser('report', parents=[parent_parser],
        add_help=False, description="Anotate variants command", help='Reporting pipeline')
    parser_report.add_argument("--rm_dup", action='store_true',
        help="Remove duplicates in bam. Default = True")
    parser_report.add_argument('--mapper', choices=['bwa'], default='bwa',
        help="Choose a short-read mapper. Default = bwa")
    parser_report.add_argument('--qc_analysis', default=True,
        help="Perform QC analysis. Default = True")
    arguments = parent_parser.parse_args()
    return arguments


#########################################################
if __name__ == '__main__':
    ARGUM = parse_arguments()
    main(ARGUM)

