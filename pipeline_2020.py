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

    # Setting analysis variables
    s.init()

    p.set_labdata_env()

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

    s.update_sample_db()
    s.update_summary_db()

def parse_arguments():
    '''parsing input arguments
    '''
    parent_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parent_parser.add_argument("-i", "--input_dir", required=True, type=str,
        help="Input directory with any of the following formats: *.fastq, *.fastq.gz, *.bam, *.vcf", dest='input_dir')
    parent_parser.add_argument("-o", "--output_dir", required=True, type=str,
        help="Output directory", dest='output_dir')
    parent_parser.add_argument("-s", "--sequencing", type=str, required=True,
        help="Sequencing experiment.", dest='sequencing', choices=['targeted', 'wgs'], default="targeted")
    parent_parser.add_argument("-r", "--reference", required=True, type=str, choices=['hg19', 'hg38'], 
        default='hg19', help="Genome reference to do mapping var calling and annotation.",  dest='reference')
    parent_parser.add_argument("-t", "--threads", type=int, default=4,
        help="Number of CPU threads to operate.", dest='threads')
    parent_parser.add_argument("--var_class", required=True, type=str, choices=['somatic', 'germline'],
        default='somatic', help="Perform somatic or germline variant analysis.", dest='var_class')
    parent_parser.add_argument("--panel", type=str, required=False,
        help="Available gene panel", dest='panel')
    parent_parser.add_argument("--lang", type=str, choices=['cat', 'en', 'esp'], default='cat',
        help="Report language.", dest='language')
    parent_parser.add_argument("--min_fusion_size", type=int, default=50000,
        help="Minimum fusion size in bp to be reported.", dest='min_fusion_size')
    parent_parser.add_argument("--sample_data", type=str,
        help="Sample data docx", dest='sample_data')
    parent_parser.add_argument("--lab_data", type=str,
        help="Lab data xlsx", dest='lab_data')
    parent_parser.add_argument("--user_id", type=str, 
        help="User id for database linking", dest='user_id')
    parent_parser.add_argument("--db_dir", type=str, required=True,
        help="Database directory harbouring sqlite files", dest='db_dir')
    parent_parser.add_argument("--ann_dir", type=str, required=True,
        help="Annotation resources directory", dest='ann_dir')
    parent_parser.add_argument("--ref_dir", type=str, required=True,
        help="Genome reference resource directory", dest='ref_dir')

    # https://stackoverflow.com/questions/33645859/how-to-add-common-arguments-to-argparse-subcommands/33646419
    parser = argparse.ArgumentParser(description="Pipeline for NGS analysis v.1.0")
    # Now subparsers
    subparsers = parser.add_subparsers(dest='command', title="sub-commands")

    # All command 
    parser_all = subparsers.add_parser('all', parents=[parent_parser],
        add_help=False, description="Perform all steps (mapping, calling, annotation and report)"
        ,help='Perform all steps (mapping, calling, annotation and report)'
        ,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # mapping 
    parser_all.add_argument("--rm_dup", action='store_true',
        help="Remove duplicates in bam.")
    parser_all.add_argument('--mapper', choices=['bwa'], default='bwa',
        help="Choose a short-read mapper.")
    parser_all.add_argument('--qc_analysis', default=True,
        help="Perform QC analysis.")
    # variant calling
    parser_all.add_argument("--gatk", action='store_true',
        help="Perform variant calling throught GATK.")
    parser_all.add_argument('--freebayes', default=False,
        help="Perform variant calling throught FreeBayes.")
    parser_all.add_argument('--manta', default=True,
        help="Perform SV calling throught Manta.")
    parser_all.add_argument('--grapes', default=True,
        help="Perform SV calling throught GRAPES.")
    # annotation
    parser_all.add_argument('--civic', default=True,
        help="Annotate with CIViC")
    parser_all.add_argument('--cgi', default=True,
        help="Annotate with Cancer Genome Interpreter (CGI)")
    parser_all.add_argument('--gnomad', dest='gnomad', default=True,
        help="Annotate gnomAD frequencies")  
    parser_all.add_argument('--1kg', dest='thousand_genomes', default=True,
        help="Annotate 1000Genomes frequencies")  

    # predictors
    parser_all.add_argument('--missense_predictors', dest='missense_predictors',
        default='sift,polyphen2,mutationtaster2,provean,fathmm,revel,mutpred', type=str,
        help="Choose missense in-silico predictors") 
    parser_all.add_argument('--splicing_predictors', dest='splicing_predictors',
        default='spliceai,maxentscan,dbscsnv', type=str, 
        help="Choose genomewide predictors")  
    parser_all.add_argument('--genomewide_predictors', dest='genomewide_predictors',
        default='cadd,ncer', type=str,
        help="Choose genomewide predictors")  
    parser_all.add_argument('--conservation', dest='conservation',
        default='phastcons,phylop,gerp', type=str,
        help="Choose missense in-silico predictors")  

    # Map command 
    parser_map = subparsers.add_parser('map', parents=[parent_parser],
        add_help=False, description="Map raw FASTQ files", help='Mapping pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_map.add_argument("--rm_dup", action='store_true',
        help="Remove duplicates in bam.")
    parser_map.add_argument('--mapper', choices=['bwa'], default='bwa',
        help="Choose a short-read mapper.")
    parser_map.add_argument('--qc_analysis', default=True,
        help="Perform QC analysis.")

    # Call command 
    parser_call = subparsers.add_parser('call', parents=[parent_parser],
        add_help=False, description="Call variants command", help='Calling pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_call.add_argument("--gatk", action='store_true',
        help="Perform variant calling throught GATK.")
    parser_call.add_argument('--freebayes', default=False,
        help="Perform variant calling throught FreeBayes.")
    parser_call.add_argument('--manta', default=True,
        help="Perform SV calling throught Manta.")
    parser_call.add_argument('--grapes', default=True,
        help="Perform SV calling throught GRAPES.")

    # Annotation command 
    parser_annotate = subparsers.add_parser('annotate', parents=[parent_parser],
        add_help=False, description="Anotate variants command", 
        help='Variant annotatione pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Cancer specific databases  
    parser_annotate.add_argument('--civic', default=True,
        help="Annotate with CIViC")
    parser_annotate.add_argument('--cgi', default=True,
        help="Annotate with Cancer Genome Interpreter (CGI)")
    # Population frequencies 
    parser_annotate.add_argument('--gnomad', default=True,
        help="Annotate gnomAD frequencies")  
    parser_annotate.add_argument('--1kg', default=True,
        help="Annotate 1000Genomes frequencies")  
    # Missense predictors
    parser_annotate.add_argument('--missense_predictors', dest='missense_predictors',
        default='sift,polyphen2,mutationtaster2,provean,fathmm,revel,mutpred', type=str,
        help="Choose missense in-silico predictors")  
    # Splicing predictors
    parser_annotate.add_argument('--splicing_predictors', dest='splicing_predictors',
        default='spliceai,maxentscan,dbscsnv', type=str, 
        help="Choose genomewide predictors")  
    # Genomewide predictors (can apply to non-coding also)
    parser_annotate.add_argument('--genomewide_predictors', dest='genomewide_predictors',
        default='cadd,ncer', type=str,
        help="Choose genomewide predictors")  
    # Conservation scores
    parser_annotate.add_argument('--conservation', dest='conservation',
        default='phastcons,phylop,gerp', type=str,
        help="Choose missense in-silico predictors")  

    parser_report = subparsers.add_parser('report', parents=[parent_parser],
        add_help=False, description="Anotate variants command", help='Reporting pipeline')
    parser_report.add_argument("--rm_dup", action='store_true',
        help="Remove duplicates in bam.")
    parser_report.add_argument('--mapper', choices=['bwa'], default='bwa',
        help="Choose a short-read mapper. Default = bwa")
    parser_report.add_argument('--qc_analysis', default=True,
        help="Perform QC analysis.")
    arguments = parser.parse_args()

    cmd_options = ['all', 'map', 'call', 'annotate', 'report']
    if  len(sys.argv) < 2:
        print(" ERROR: Invalid subcommand option. Please choose between: all, map, call, annotate, report")
        sys.exit()

    else:
        if sys.argv[1] not in cmd_options:
            print(" ERROR: Invalid subcommand option. Please choose between: all, map, call, annotate, report")   
            sys.exit()

    return arguments


#########################################################
if __name__ == '__main__':
    ARGUM = parse_arguments()
    main(ARGUM)

