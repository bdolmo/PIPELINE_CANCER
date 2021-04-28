
import os
import sys
import re
import shutil
import logging
import gzip
from collections import defaultdict
from pathlib import Path
import subprocess
import pybedtools
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
#from matplotlib.patches import Rectangle
import seaborn as sns
import pandas as pd
import numpy as np
from modules import params as p
from modules import trimming as t
from modules import map as m
from modules import var_call as v
from modules import annotate as a
from modules import sqlite as s
from modules import report as r
from modules import lowpass as l
from modules import utils as u

def do_sarscov2():
    '''Main function for SARS-COV-2 analysis
    '''

    fastq_list =  u.get_input_files( p.analysis_env['INPUT_DIR'], "fastq")

    for fq in fastq_list:
        print(fq)
        fq_tmp = os.path.basename(fq).split("_")
        sample  = fq_tmp[0]
        # Create sample directories
        p.sample_env[sample]['SAMPLE_FOLDER'] = \
          p.analysis_env['OUTPUT_DIR'] + "/" + sample

        sample_path = Path(p.sample_env[sample]['SAMPLE_FOLDER'])
        if not sample_path.is_dir():
          os.mkdir(sample_path)

        p.sample_env[sample]['FASTQ_FOLDER'] = \
          p.sample_env[sample]['SAMPLE_FOLDER'] + "/" + "FASTQ_FOLDER"

        if '_R1_' in fq:
            p.sample_env[sample]['RAW_FQ1'] = fq
            p.sample_env[sample]['READY_FQ1'] = fq

        if '_R2_' in fq:
            p.sample_env[sample]['RAW_FQ2'] = fq
            p.sample_env[sample]['READY_FQ2'] = fq

    #trim_primers_fastp()

    # Perform mapping
    m.map_fastq()

    # Remove primer sequences
    trim_primers_ivar()

    # Extract qc, coverage metrics
    qc_metrics()

    # Create a consesnsus sequence
    create_consensus_sequence()

    # Perform variant calling
    variant_calling()

    # Experimental annotation with VEP
    # Not working properly right now
    annotate()

    # Lineage identification using pangolin
    identify_lineage()

    # Create coverage plots
    plot_coverage()

    # Dump a general summary of the results
    summarize_results()

def summarize_results():
    '''Get a results summary in csv format
    '''

    p.analysis_env['SUMMARY_RESULTS'] = p.analysis_env['OUTPUT_DIR'] + \
        "/summary_results.csv"
    o = open(p.analysis_env['SUMMARY_RESULTS'], "w")
    o.write(",,,Call Rate\n")
    o.write("Sample,N Reads,%Gap,1X,10X,30X,100X,1000X,Median Cov,Lineage_QC,Lineage_Probability,Lineage\n")
    for sample in p.sample_env:

        lineage_dict = defaultdict(dict)
        header = []
        with open(p.sample_env[sample]['PANGOLIN_REPORT']) as f:
            for line in f:
                line = line.rstrip("\n")
                tmp = line.split(",")
                if line.startswith("taxon"):
                    for field in tmp:
                        header.append(field)
                        lineage_dict[field] = "."
                else:
                    i = 0
                    for info in tmp:
                        lineage_dict[header[i]]=info
                        i+=1

        qc_dict = p.sample_env[sample]['QC_DICT']
        results = []
        results.append(str(sample))
        results.append(str(qc_dict['total_reads']))
        results.append(str(qc_dict['gaps']))
        results.append(str(qc_dict['1X']))
        results.append(str(qc_dict['10X']))
        results.append(str(qc_dict['30X']))
        results.append(str(qc_dict['100X']))
        results.append(str(qc_dict['1000X']))
        results.append(str(qc_dict['median_coverage']))
        results.append(str(lineage_dict['status']))
        results.append(str(lineage_dict['probability']))
        results.append(str(lineage_dict['lineage']))

        o.write(",".join(results) + "\n")
        f.close()
    o.close()

def identify_lineage():

    for sample in p.sample_env:

        # Linage report files
        p.sample_env[sample]['PANGOLIN_REPORT'] = p.sample_env[sample]['SAMPLE_FOLDER']\
            + "/lineage_report.csv"

        if not os.path.isfile(p.sample_env[sample]['PANGOLIN_REPORT']):
            # Now use pangolin to assign a lineage
            bashCommand = ('{} run -v {}:/sample_folder/ -it {} pangolin '
                ' /sample_folder/{} --outdir /sample_folder/{} '
                ' '.format(
                p.system_env['DOCKER'], p.sample_env[sample]['SAMPLE_FOLDER'],
                p.docker_env['PANGO'],p.sample_env[sample]['CONSENSUS_FASTA_FILENAME'],
                '.'))

            msg = " INFO: Identifying Lineage on sample " + sample
            print (msg)
            logging.info(msg)
            logging.info(bashCommand)

            p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output = p1.stdout.decode('UTF-8')
            error  = p1.stderr.decode('UTF-8')

            if not error:
                msg = " INFO: Lineage identification ended successfully on sample " + sample
                print (msg)
                logging.error(msg)
            else:
                msg = " ERROR: Something went wrong with lineage identification on sample " + sample
                print (msg)
                logging.info(msg)

def plot_coverage():
    '''Plot coverage along the genome, and variants
    '''

    # Setting a custom color palette
    color_palette = ["#96bb7c", "#c64756", "#96bb7c", "#fad586", "#28b5b5",
        "#4b778d", "#8fd9a8", "#d2e69c", "#d8ebe4", "#007580", "#282846",
        "#fed049", "#93329e", "#440a67", "#b4aee8"]

    for sample in p.sample_env:

        p.sample_env[sample]['COVERAGE_PLOT'] = \
            p.analysis_env['OUTPUT_DIR']+ "/" + sample + ".coverage.png"

        sns.set(rc={"figure.dpi":300, 'savefig.dpi':300})
        sns.set_theme()
        sns.set_style('ticks')

        fig, axes = plt.subplots(figsize=(20, 5))
        fig.suptitle(sample, fontsize=20)

        df_cov = pd.read_csv(p.sample_env[sample]['PER_BASE_COV'], sep="\t",
            names=['Chr', 'Position', 'Depth'])

        #axes[0].plt.plot()
        #axes[0].stackplot(df_cov['Position'], (0,10), colors=['grey'],
            #labels=[sample])

        axes.stackplot(df_cov['Position'], df_cov['Depth'], colors=['grey'],
            labels=[sample])

        genes_dict = defaultdict(dict)
        # Open Sars-cov-2 gtf file and create a dict
        i = 0
        with gzip.open(p.aux_env['GTF_GENOME'],'rt') as fin:
            for line in fin:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    continue
                if not 'mRNA' in line:
                    continue
                tmp = line.split("\t")
                info_list = tmp[8].split(";")
                gene_name = ""
                pos = int(tmp[3])
                end = int(tmp[4])
                size = end-pos
                for field in info_list:
                    if field.startswith("Name"):
                        gene_name = field
                        break
                i+=1
                gene_name = gene_name.replace("Name=", "")
                rect = plt.Rectangle((pos, -550), size, 500,
                    edgecolor="black",linewidth=1, color=color_palette[i])
                axes.add_patch(rect)
                #axes.annotate()
        fin.close()

        #axes[1].stackplot(df_cov['Position'], df_cov['Depth'], colors=['grey'],
        #    labels=[sample])
        #cov_plot.figure.savefig(plot)
        #cov_plot.figure.savefig(p.sample_env[sample]['COVERAGE_PLOT'])
        plt.savefig(p.sample_env[sample]['COVERAGE_PLOT'])
        plt.close()
    pass

def annotate():
    '''Annotate with VEP
    '''

    #vep -i myvariants.vcf -gff Sars_cov_2.ASM985889v3.100.primary_assembly.MN908947.3.gff3.gz -fasta Sars_cov_2.ASM985889v3.dna_sm.toplevel.fa.gz -synonyms chr_synonyms.txt
    for sample in p.sample_env:

        p.sample_env[sample]['VEP_VCF'] = \
            p.sample_env[sample]['READY_SNV_VCF'].replace(".vcf", ".vep.vcf")

        p.sample_env[sample]['VEP_VCF_NAME'] = \
            os.path.basename(p.sample_env[sample]['READY_SNV_VCF_NAME']).replace(".vcf", ".vep.vcf")

        vep_output_vcf = p.aux_env['VEP_FOLDER_OUTPUT'] + "/" + p.sample_env[sample]['VEP_VCF_NAME']

        #vep -i myvariants.vcf -gff Sars_cov_2.ASM985889v3.100.primary_assembly.MN908947.3.gff3.gz
        # -fasta Sars_cov_2.ASM985889v3.dna_sm.toplevel.fa.gz -synonyms chr_synonyms.t
        bashCommand = ('{} run -t -i -v {}:/genomedir/ -v {}:/anndir/ -v {}:/opt/vep/.vep {}'
        ' perl vep --input_file /opt/vep/.vep/input/{} --output_file /opt/vep/.vep/output/{}'
        ' --gff /opt/vep/.vep/sars_cov_2/{} --synonyms /opt/vep/.vep/sars_cov_2/{}'
        ' --format vcf --vcf --force_overwrite'
        ' --fasta /genomedir/{} '
        .format(p.system_env['DOCKER'], p.aux_env['GENOME_FOLDER'],\
        p.analysis_env['ANN_DIR'], p.aux_env['VEP_FOLDER'], p.docker_env['VEP'],\
        p.sample_env[sample]['READY_SNV_VCF_NAME'], p.sample_env[sample]['VEP_VCF_NAME'],\
        p.aux_env['GTF_GENOME_NAME'], p.aux_env['CHR_SYNONYMS_NAME'], p.aux_env['GENOME_NAME'] ))

        if not os.path.isfile(p.sample_env[sample]['VEP_VCF']):

            # First, copy vcf to vep input dir
            shutil.copy2(p.sample_env[sample]['READY_SNV_VCF'],
                p.aux_env['VEP_FOLDER_INPUT']+"/"+p.sample_env[sample]['READY_SNV_VCF_NAME'])

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
        p.sample_env[sample]['READY_SNV_JSON'] = p.sample_env[sample]['VEP_VCF'].\
            replace(".vcf", ".json")
        p.sample_env[sample]['READY_SNV_VCF_NAME'] = p.sample_env[sample]['VEP_VCF_NAME']

        # Create JSON file
        u.convert_vcf_2_json(p.sample_env[sample]['READY_SNV_VCF'])

def qc_metrics():
    '''Extract qc, and coverage
    '''
    for sample in p.sample_env:

        # Create QC directory
        p.sample_env[sample]['QC_FOLDER'] = \
            p.sample_env[sample]['SAMPLE_FOLDER'] + "/" + "QC_FOLDER"

        qc_path = Path(p.sample_env[sample]['QC_FOLDER'])
        if not qc_path.is_dir():
            os.mkdir(qc_path)

        #p.sample_env[sample]['SUMMARY_QC'] = \
        #    p.sample_env[sample]['QC_FOLDER'] + "/" + sample + "_summary_qc.csv"

        p.sample_env[sample]['PER_BASE_COV'] = \
            p.sample_env[sample]['QC_FOLDER'] + "/" + sample + ".coverage.bed"

        # Getting cov metrics with Mosdepth
        if not os.path.isfile(p.sample_env[sample]['PER_BASE_COV']):

            bashCommand = ('{} depth -aa {} > {}').format(p.system_env['SAMTOOLS'],\
                 p.sample_env[sample]['READY_BAM'],p.sample_env[sample]['PER_BASE_COV'])
            # bashCommand = ('{} --fast-mode {} {}') \
            #     .format(p.system_env['MOSDEPTH'],  \
            #     p.sample_env[sample]['QC_FOLDER'] + "/" + sample, p.sample_env[sample]['READY_BAM'])
            # print(bashCommand)
            msg = " INFO: Extracting coverage metrics for sample " + sample
            print (msg)
            logging.info(msg)
            logging.info(bashCommand)

            p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output = p1.stdout.decode('UTF-8')
            error  = p1.stderr.decode('UTF-8')

            if not error:
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

    # Merged coverages files
    p.analysis_env['MERGED_COVERAGES'] = p.analysis_env['OUTPUT_DIR'] +\
        "/merged.coverages.bed"

    for sample in p.sample_env:

        qc_dict = {
            'total_bases' : 0,
            'gaps' : 0,
            'median_coverage' : 0,
            'mean_coverage' : 0,
            'std_coverage' : 0,
            '1X' : 0,
            '10X': 0,
            '20X': 0,
            '30X': 0,
            '100X': 0,
            '1000X' : 0,
        }

        coverage_dict = defaultdict(dict)
        sample_header = []
        cov_list = []

        # Getting total number of reads (we could use pysam instead)
        bashCommand = ('{} view -c {}').format(p.system_env['SAMTOOLS'], \
        p.sample_env[sample]['READY_BAM'])
        msg = " INFO: Getting total reads for sample " + sample
        print (msg)
        logging.info(bashCommand)

        p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = p1.stdout.decode('UTF-8')
        error  = p1.stderr.decode('UTF-8')
        if error:
            msg = " ERROR: Could not extract the total number reads for sample "+ sample
            print(msg)
            p.logging.error(error)
            sys.exit()
        else:
            output = output.rstrip("\n")
            qc_dict['total_reads'] = output
            #p.sample_env[sample]['TOTAL_READS'] = output
        total_bases = 0
        gaps = 0
        bases_1X = 0
        bases_10X= 0
        bases_30X= 0
        bases_100X = 0
        bases_1000X = 0
        sample_header.append(sample)
        #with gzip.open(p.sample_env[sample]['PER_BASE_COV'],'rt') as fin:
        with open(p.sample_env[sample]['PER_BASE_COV']) as fin:
            for line in fin:
                line = line.rstrip('\n')
                tmp = line.split("\t")
                coordinate = tmp[0]+"\t"+tmp[1]+"\t"+tmp[1]
                coverage   = tmp[2]
                cov_list.append(coverage)
                qc_dict['total_bases'] = str(int(qc_dict['total_bases'])+ 1)
                total_bases+=1
                if int(coverage) == 0:
                    gaps+=1
                if int(coverage) >= 1:
                    bases_1X+=1
                if int(coverage) >= 10:
                    bases_10X+=1
                #if int(coverage) >= 20:
                #    qc_dict["20X"]+=1
                if int(coverage) >= 30:
                    bases_30X+=1
                if int(coverage) >= 100:
                    bases_100X+=1
                if int(coverage) >= 1000:
                    bases_1000X+=1

                if not coordinate in coverage_dict:
                    coverage_dict[coordinate] = defaultdict(dict)
                    coverage_dict[coordinate][sample] = coverage
                else:
                    coverage_dict[coordinate][sample] = coverage
        fin.close()
        qc_dict["gaps"] = str(round(100*(gaps/total_bases), 2))
        qc_dict["1X"]   = str(round(100*(bases_1X/total_bases), 2))
        qc_dict["10X"]  = str(round(100*(bases_10X/total_bases), 2))
        qc_dict["30X"]  = str(round(100*(bases_30X/total_bases), 2))
        qc_dict["100X"] = str(round(100*(bases_100X/total_bases), 2))
        qc_dict["1000X"]= str(round(100*(bases_1000X/total_bases), 2))

        # create a numpy array from a python list
        cov_arr = np.array(cov_list).astype(np.int)

        # Calculate median coverage
        median_cov = np.median(cov_arr)
        mean_cov   = np.mean(cov_arr)
        std_cov    = np.std(cov_arr)
        qc_dict['median_coverage'] = median_cov
        qc_dict['mean_coverage']   = mean_cov
        qc_dict['std_coverage']    = std_cov

        p.sample_env[sample]['QC_DICT'] = qc_dict

    o = open(p.analysis_env['MERGED_COVERAGES'], "w")
    o.write("Chromosome\tStart\tEnd\t" + "\t".join(sample_header)+"\n" )
    for coordinate in coverage_dict:
        o.write(coordinate)
        for sample in coverage_dict[coordinate]:
            coverage = coverage_dict[coordinate][sample]
            o.write("\t" + coverage)
        o.write("\n")
    o.close()

def variant_calling():

    for sample in p.sample_env:
        # Create BAM_FOLDER if not present
        p.sample_env[sample]['VCF_FOLDER'] = \
         p.sample_env[sample]['SAMPLE_FOLDER'] + "/" + "VCF_FOLDER"

        vcf_folder_path = Path(p.sample_env[sample]['VCF_FOLDER'])
        if not vcf_folder_path.is_dir():
            os.mkdir(vcf_folder_path)

        p.sample_env[sample]['FREEBAYES_VCF'] = p.sample_env[sample]['VCF_FOLDER'] + \
             "/" + sample + ".freebayes.vcf"
        p.sample_env[sample]['FREEBAYES_VCF_NAME'] =  sample +  ".freebayes.vcf"

        p.sample_env[sample]['READY_SNV_VCF'] = p.sample_env[sample]['FREEBAYES_VCF']
        p.sample_env[sample]['READY_SNV_VCF_NAME'] = p.sample_env[sample]['FREEBAYES_VCF_NAME']

        if not os.path.isfile(p.sample_env[sample]['FREEBAYES_VCF']):

            bashCommand = ('{} -p 1 -q 20 -m 60 --min-coverage 10 -V -f {}  {} > {}')\
                .format(p.system_env['FREEBAYES'], \
                p.aux_env['GENOME_FASTA'], p.sample_env[sample]['READY_BAM'],
                p.sample_env[sample]['FREEBAYES_VCF'] )
            logging.info(bashCommand)
            p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            output = p1.stdout.decode('UTF-8')
            error  = p1.stderr.decode('UTF-8')

            if not error:
                msg = " INFO: Freebayes was executed successfully on " + sample
                print (msg)
                logging.info(msg)
            else:
                msg = " ERROR: Something happened when using Freebayes on " + sample
                print (msg)
                logging.error(msg)

def create_consensus_sequence():
    '''Perform variant calling and create a consensus fasta sequence
    '''
    for sample in p.sample_env:
        # Setting consensus fasta name
        p.sample_env[sample]['CONSENSUS_FASTA'] = \
            p.sample_env[sample]['SAMPLE_FOLDER'] + "/" + sample + ".consensus.fa"

        p.sample_env[sample]['CONSENSUS_FASTA_FILENAME'] = sample + ".consensus.fa"

        if not os.path.isfile(p.sample_env[sample]['CONSENSUS_FASTA']):

            bashCommand = ('{} mpileup -A -d 6000000 -B -Q 0 --reference {} {}'
            ' | {} consensus -p {} -n N')\
            .format(p.system_env['SAMTOOLS'], p.aux_env['GENOME_FASTA'],
            p.sample_env[sample]['READY_BAM'], p.system_env['IVAR'],
            p.sample_env[sample]['CONSENSUS_FASTA'])
            print(bashCommand)

            msg = " INFO: Calling variants and creating a consensus sequence on sample " + sample
            print(msg)
            p.logging.info(msg)
            p.logging.info(bashCommand)

            p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output = p1.stdout.decode('UTF-8')
            error  = p1.stderr.decode('UTF-8')

def trim_primers_fastp():
    '''Trim FASTQ files with fastp
    '''
    for sample in p.sample_env:

        # Setting trimmed fq's
        p.sample_env[sample]['READY_FQ1'] =  p.sample_env[sample]['RAW_FQ1'].\
            replace(".fastq.gz", ".trimmed.fastq.gz")

        p.sample_env[sample]['READY_FQ2'] =  p.sample_env[sample]['RAW_FQ2'].\
            replace(".fastq.gz", ".trimmed.fastq.gz")

        # Setting fastq names in the global dictionaries
        p.sample_env[sample]['FASTP_JSON'] = \
          p.sample_env[sample]['FASTQ_FOLDER'] + "/" + "fastp.json"

        # Now trimming with fastp
        bashCommand = ('{} -i  {} -I {} -o {} -O {} -w {} -j {}') \
          .format(p.system_env['FASTP'], p.sample_env[sample]['RAW_FQ1'],
          p.sample_env[sample]['RAW_FQ2'],
          p.sample_env[sample]['READY_FQ1'],
          p.sample_env[sample]['READY_FQ2'],
          p.analysis_env['THREADS'],
          p.sample_env[sample]['FASTP_JSON'] )

        if not os.path.isfile(p.sample_env[sample]['READY_FQ1']) \
          and not os.path.isfile(p.sample_env[sample]['READY_FQ2']):

          msg = " INFO: Trimming sample " +  sample
          print(msg)
          p.logging.info(msg)
          p.logging.info(bashCommand)

          p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
          output = p1.stdout.decode('UTF-8')
          error  = p1.stderr.decode('UTF-8')
          if error:
            if re.search("error", error):
              msg = " ERROR: Something wrong happened with fastp trimming:"
              print(msg)
              p.logging.error(error)
              sys.exit()
            else:
              msg = " INFO: FASTQ Trimming ended successfully for sample " + sample
              p.logging.info(msg)
              print(msg)
          else:
            msg = " INFO: FASTQ Trimming ended successfully for sample " + sample
            p.logging.info(msg)
            print(msg)
        else:
          msg = " INFO: Skipping trimming of sample " +  sample
          print(msg)
          p.logging.info(msg)

def trim_primers_ivar():
    '''Trim primers on BAM files using ivar
    '''

    for sample in p.sample_env:

        p.sample_env[sample]['TRIMMED_BAM_NAME'] = p.sample_env[sample]['RAW_BAM']\
            .replace(".bam", ".trimmed")

        p.sample_env[sample]['TRIMMED_BAM'] = p.sample_env[sample]['RAW_BAM']\
            .replace(".bam", ".trimmed.bam")

        p.sample_env[sample]['READY_BAM'] = p.sample_env[sample]['RAW_BAM']\
            .replace(".bam", ".trimmed.sorted.bam")

        # Trim primer sequences using IVAR
        if not os.path.isfile(p.sample_env[sample]['TRIMMED_BAM']):
            bashCommand = ('{} trim -e -i {} -b {} -p {}')\
            .format(p.system_env['IVAR'],p.sample_env[sample]['RAW_BAM'] ,\
            p.analysis_env['PRIMERS'],p.sample_env[sample]['TRIMMED_BAM_NAME'])

            msg = " INFO: Trimming primers on sample " + sample
            print(msg)
            print(bashCommand)
            p.logging.info(msg)
            p.logging.info(bashCommand)

            p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output = p1.stdout.decode('UTF-8')
            error  = p1.stderr.decode('UTF-8')

        # Now sort BAMs
        if not os.path.isfile(p.sample_env[sample]['READY_BAM']):
            bashCommand = ('{} sort {} -T {} -o {}').format(p.system_env['SAMTOOLS'],
            p.sample_env[sample]['TRIMMED_BAM'],
            sample, p.sample_env[sample]['READY_BAM'])
            msg = " INFO: Sorting trimmed BAM from " + sample
            print(msg)
            p.logging.info(msg)
            p.logging.info(bashCommand)

            p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output = p1.stdout.decode('UTF-8')
            error  = p1.stderr.decode('UTF-8')

        m.index_bam(p.sample_env[sample]['READY_BAM'])
