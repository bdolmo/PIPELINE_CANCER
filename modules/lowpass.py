#!/usr/bin/env python3

import os
import sys
import re
import logging
import gzip
from collections import defaultdict
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.signal import savgol_coeffs, savgol_filter
import subprocess
from natsort import natsorted, index_natsorted, order_by_index
from modules import utils as u
from modules import params as p
import pybedtools
import pyranges as pr
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
pd.options.mode.chained_assignment = None  # default='warn'

class CnvPlot:
    '''Class for plotting cnvs
    ''' 
    def __init__(self, cnr_file, cns_file, calls, sample, output_dir):
        # ratios  
        self.cnr_file   = cnr_file
        # segments
        self.cns_file   = cns_file
        # calls
        self.calls   = calls 

        self.sample     = sample
        self.output_dir = output_dir

    def plot_genome(self, genomewide, by_chr):

        cnr_df = pd.read_csv(self.cnr_file, sep ="\t")

        # Naming the genome plot 
        plot = self.output_dir + "/" + self.sample + ".genomewide.png"

        # Setting chromosome color 
        palette_dict = defaultdict(dict)
        color_list =["#4f6b76", "#b2bbc0"] 
    
        unique_chromosomes = cnr_df['Chromosome'].unique().tolist()
        idx = 0
        for chr in unique_chromosomes:
            if idx==2:
                idx = 0
            palette_dict[chr] = color_list[idx]
            idx+=1
        x_list = []

        # Setting xtick divisions and labels 
        chromosomes = cnr_df['Chromosome'].tolist()
        i = 0
        xtick_list = []
        for chr in chromosomes:
            if not chr in x_list:
                x_list.append(chr)
                xtick_list.append(i)
            else:
                x_list.append("")
            i+=1
        if genomewide == True:

            if not os.path.isfile(plot):
                sns.set(rc={"figure.dpi":300, 'savefig.dpi':300})
                sns.set_theme()
                sns.set_style('ticks')
                fig, axes = plt.subplots(figsize=(20,7))
                fig.suptitle(self.sample, fontsize=20)
                sample_ratio = self.sample + "_ratio"

                ratio_plot = sns.scatterplot(data=cnr_df, x=cnr_df.index, 
                    y=cnr_df[sample_ratio], size=0.015, hue=cnr_df.Chromosome, 
                    palette=palette_dict, alpha=0.4, edgecolor="none")
                
                # Setting y limits 
                ratio_plot.set(ylim=(0,3))
                ratio_plot.set_xticks(xtick_list)
                ratio_plot.set_xticklabels(unique_chromosomes, rotation=30, size=15)
                ratio_plot.set_yticks(ratio_plot.get_yticks())
                ratio_plot.set_yticklabels(ratio_plot.get_yticks(), size=15)
                ratio_plot.get_legend().remove()
                ratio_plot.set_xlabel("", fontsize=20)
                ratio_plot.set_ylabel("Ratio", fontsize=20)

                # Adding vertical lines to separate chromosomes 
                for xc in xtick_list:
                    ratio_plot.axvline(x=xc, color="black")

                # Saving as png  
                ratio_plot.figure.savefig(plot)
                plt.close()

                return plot

    def plot_cnv(self, chr, start, end, gene_list=None, add_genes=False):
        '''Plot a CNV call, add gene labels (optionally)
        '''        
        
        cnr_df = pd.read_csv(self.cnr_file, sep ="\t")
        #cns_df = pd.read_csv(self.cns_file, sep = "\t")
        if gene_list:
            genes_df = pr.read_bed(gene_list)

        tmp_call = self.output_dir + "/" + chr + "." + start + "." + end + ".tmp.bed"
        o = open (tmp_call, 'w')
        o.write(chr + "\t" + start + "\t" + end + "\n")
        o.close()

        segments_no_header = remove_bed_header(self.cns_file, 'Chromosome')

        a = pybedtools.BedTool(tmp_call)
        b = pybedtools.BedTool(segments_no_header)
        c = a.intersect(b, wo=True, stream=True)

        calls_with_segments = tmp_call.replace(".tmp.bed", ".segments.bed")
        o = open(calls_with_segments, 'w')
        for line in iter(c):
            line = str(line)
            line = line.rstrip()
            tmp  = line.split('\t')
            o.write(tmp[0]+"\t"+tmp[1]+"\t"+tmp[2]+"\t"+tmp[7]+"\n")
            o.write(tmp[0]+"\t"+tmp[2]+"\t"+tmp[2]+"\t"+tmp[7]+"\n")
        o.close()
        segment_data = pd.read_csv(calls_with_segments, sep="\t", names=['Chromosome', 'Start', 'End', 'Ratio'])

        plot = self.output_dir + "/" + chr + "." + start + "." + end + ".png"

        #if not os.path.isfile(plot):
        ratio_data   = cnr_df.loc[(cnr_df['Chromosome'] == chr) & 
            (cnr_df['Start'] >= int(start)) & (cnr_df['End'] <= int(end))]
        sns.set(rc={"figure.dpi":300, 'savefig.dpi':300})
        sns.set_theme()
        sns.set_style('ticks')

        plot_title = chr + ":" + start + "-" + end
        sample_ratio = self.sample + "_ratio"
        fig, axes = plt.subplots(figsize=(20,7))
        fig.suptitle(plot_title, fontsize=20)
          
        cnv_plot = sns.scatterplot(data=ratio_data, x=ratio_data.Start, 
            y=ratio_data[sample_ratio], size=0.015, alpha=0.4, 
            edgecolor="none")

        cnv_plot = sns.lineplot(data=segment_data, x=segment_data.Start, 
             y=segment_data.Ratio)
        cnv_plot.figure.savefig(plot)
        plt.close()

        pass

# class Cnv:
#   '''Cnv class. Input format is a pd's dataframe
#   '''
#   def __init__(self, df, sample, output_dir):
#     # Defining attributes
#     self.df         = df
#     self.sample     = sample
#     self.output_dir = output_dir

#   def call_cnv(self, del_threshold, dup_threshold, z_score):
#     '''Cnv class
#     ''' 

def do_lowpass():
    '''Main function for lowpass CNV analysis
    '''

    if p.analysis_env['ANALYSIS_MODE'] == "call": 
      bam_list =  u.get_input_files( p.analysis_env['INPUT_DIR'], "bam")

      for bam in bam_list:
        sample_name = os.path.basename(bam).replace(".bam", "")
        p.sample_env[sample_name]['READY_BAM']      = bam
        p.sample_env[sample_name]['READY_BAM_NAME'] = os.path.basename(bam)
        p.sample_env[sample_name]['BAM_FOLDER']     = p.analysis_env['INPUT_DIR']
        p.sample_env[sample_name]['SAMPLE_FOLDER']  = p.analysis_env['OUTPUT_DIR']\
          + "/" + sample_name

        sample_path = Path(p.sample_env[sample_name]['SAMPLE_FOLDER'])
        if not sample_path.is_dir():
            os.mkdir(sample_path)

    # bin size in base pairs
    bin_size = 50000
    p.analysis_env['LOWPASS_BINS'] = \
        p.aux_env['GENOME_FOLDER'] + "/" + "lowpass." +  str(bin_size) + ".bp.bed"

    gc_dict  = defaultdict(dict)
    map_dict = defaultdict(dict)

    if not os.path.isfile(p.analysis_env['LOWPASS_BINS']):
        
        # Create bin regions with a fi 
        create_bins(bin_size)

        # Annotate GC content
        gc_dict = annotate_gc(p.analysis_env['LOWPASS_BINS'])

        # Annotate mappability
        map_dict = annotate_mappability(p.analysis_env['LOWPASS_BINS'])

        # Add gc and map
        tmp_gc_map = p.analysis_env['LOWPASS_BINS'].replace(".bed", ".tmp.bed")
        o = open(tmp_gc_map, 'w')
        with open (p.analysis_env['LOWPASS_BINS']) as f:
            for line in f:
                line = line.rstrip("\n")
                tmp = line.split("\t")
                coordinate = tmp[0] + "\t" + tmp[1] + "\t" + tmp[2]
                o.write(coordinate + "\t" + str(gc_dict[coordinate]) + "\t" + str(round(map_dict[coordinate], 2)) + "\n" )   
        f.close()
        o.close()
        os.remove(p.analysis_env['LOWPASS_BINS'])
        os.rename(tmp_gc_map, p.analysis_env['LOWPASS_BINS'])
    else:
        with open (p.analysis_env['LOWPASS_BINS']) as f:
            for line in f:
                line = line.rstrip("\n")
                tmp = line.split("\t")
                coordinate = tmp[0] + "\t" + tmp[1] + "\t" + tmp[2]
                gc_dict[coordinate]  = str(round(float(tmp[3] ), 2))  
                map_dict[coordinate] = str(round(float(tmp[4] ), 2)) 

    # Extract coverage and add GC-content 
    extract_coverage(gc_dict, map_dict)

    # Normalization 
    for sample in p.sample_env:
        p.sample_env[sample]['NORMALIZED_COVERAGE'] = \
            p.sample_env[sample]['COV_FOLDER'] + "/" + sample + ".normalized.coverage.bed"

        # Load a dataframe of raw coverage data
        df = pd.read_csv(p.sample_env[sample]['RAW_COVERAGE'], sep="\t")

        # Calculate median coverage
        median_cov = df['coverage'].median()

        # Normalize by GC-content 
        df = normalize(df, sample, median_cov, ['gc', 'map'] )
       
        # Export normalized coverage 
        df.to_csv(p.sample_env[sample]['NORMALIZED_COVERAGE'], sep="\t", mode='w', index=None)

    # Plotting normalized coverage 
    for sample in p.sample_env:

        p.sample_env[sample]['NORMALIZED_COVERAGE'] = \
            p.sample_env[sample]['COV_FOLDER'] + "/" + sample + ".normalized.coverage.bed"

        # Read normalized coverage as a pandas dataframe 
        df = pd.read_csv(p.sample_env[sample]['NORMALIZED_COVERAGE'], sep="\t")

        # Plotting 
        cov_plot = plot_normalization(df, sample, p.sample_env[sample]['COV_FOLDER'])

    # Merge sample coverages into a single file/dataframe
    merged_df = merge_samples_coverage()

    # Reference creation
    #create_heatmap(merged_df)

    calculate_ratios(merged_df)

    segment_coverage(3, 0.0001)

    call_cnvs(0.6, 1.4, 3)

    # Plotting
    for sample in p.sample_env:

        cnp = CnvPlot(
            cnr_file=p.sample_env[sample]['RATIO_FILE'], 
            cns_file=p.sample_env[sample]['SEGMENT_FILE'], 
            calls= p.sample_env[sample]['CALLS_FILE'], 
            sample=sample, 
            output_dir=p.sample_env[sample]['COV_FOLDER']
        )

        # Plotting genomewide CNV profile
        sample_plot = cnp.plot_genome(genomewide=True, by_chr=False)

        # Now plotting isolated CNVs
        with open (p.sample_env[sample]['CALLS_FILE']) as cf:
            for line in cf:
                line = line.rstrip('\n')
                if line.startswith('Chromosome'):
                    continue
                tmp = line.split('\t')
                chr  = tmp[0]
                start= tmp[1]
                end  = tmp[2]   
                cnv_plot = cnp.plot_cnv(chr=chr, start=start, end=end)

def create_heatmap(df):
    '''Plotting coverage normalization 
    '''

    sample_list = []
    for sample in p.sample_env:
        sample_list.append(sample) 
    data = df[sample_list]
    dat_corr = data.corr()
    sns.set_context("talk")

    heatmap = sns.clustermap(dat_corr, metric="correlation", cmap='coolwarm')
    p.analysis_env['HEATMAP'] = \
        p.analysis_env['OUTPUT_DIR'] + "/" + "heatmap.png"    
      
    heatmap.figure.savefig(p.analysis_env['HEATMAP'])
    heatmap.close()

def merge_samples_coverage():
    '''Plotting coverage normalization 
    '''
    p.analysis_env['MERGED_COVERAGES'] = \
        p.analysis_env['OUTPUT_DIR'] + "/" + "merged.normalized.coverage.bed"

    merged_df = pd.read_csv(p.analysis_env['LOWPASS_BINS'], names=['Chromosome', 'Start', 'End', 'gc', 'map'], sep ="\t")

    for sample in p.sample_env:
        sample_df = pd.read_csv(p.sample_env[sample]['NORMALIZED_COVERAGE'], sep ="\t")
        merged_df[sample] = sample_df['normalized_map']

    merged_df.to_csv(p.analysis_env['MERGED_COVERAGES'], sep="\t", mode='w', index=None)
    return merged_df

def calculate_ratios(cnr_df):
    '''Calculate bin ratio
    '''
    ratio_fields = [] 
 
    for sample in p.sample_env:

        baseline_samples = []
        for oth in p.sample_env:
            if oth == sample:
                continue
            baseline_samples.append(oth)

        p.sample_env[sample]['RATIO_FILE'] \
            = p.sample_env[sample]['COV_FOLDER'] + "/" + sample + ".ratios.bed"
        p.sample_env[sample]['RATIO_PLOT'] = \
            p.sample_env[sample]['COV_FOLDER'] + "/" + sample + ".ratios.png" 

        if not os.path.isfile(p.sample_env[sample]['RATIO_FILE']):
            sample_ratio = sample + "_ratio"

            median_autosomes  = cnr_df[(cnr_df['Chromosome'] != "chrX") & (cnr_df['Chromosome'] != "chrY")][sample].median()
            median_chrx       = cnr_df[cnr_df['Chromosome'] != "chrX"][sample].median()
            median_chry       = cnr_df[cnr_df['Chromosome'] != "chrY"][sample].median()

            msg = " INFO: Sample {}  median autosomes={}, median chrX={}, median chrY={}"\
                .format(sample, median_autosomes, median_chrx, median_chry)
            #print (msg)

            cnr_df[sample_ratio] = cnr_df.apply(do_ratio_ref, baseline=baseline_samples, sample=sample, axis=1)

            # Filter low mappability and extreme gc content bins
            cnr_df = cnr_df.loc[(cnr_df['map'] >= 90) & (cnr_df['gc'] >=30) & (cnr_df['gc'] < 75 )]

            # Now substract sample ratios 
            cnr_df = cnr_df[['Chromosome','Start','End','gc','map', sample_ratio]]

            # Write dataframe as bed 
            cnr_df.to_csv(p.sample_env[sample]['RATIO_FILE'], sep="\t", mode='w', index=None)
            p.sample_env[sample]['CNR_DF'] = cnr_df

        else:
            pass

        # if not os.path.isfile(p.sample_env[sample]['RATIO_PLOT']):
        #     cnr_df = pd.read_csv(p.sample_env[sample]['RATIO_FILE'], sep ="\t")
        #     palette_dict = defaultdict(dict)

        #     color_list =["#4f6b76", "#b2bbc0"] 
        #     unique_chromosomes = cnr_df['chromosome'].unique().tolist()
        #     idx = 0
        #     for chr in unique_chromosomes:
        #         if idx==2:
        #             idx = 0
        #         palette_dict[chr] = color_list[idx]
        #         idx+=1
        #     x_list = []
        #     chromosomes = cnr_df['chromosome'].tolist()
        #     i = 0
        #     xtick_list = [] 
        #     for chr in chromosomes:
        #         if not chr in x_list:
        #             x_list.append(chr)
        #             xtick_list.append(i)
        #         else:
        #             x_list.append("")
        #         i+=1

        #     sns.set(rc={"figure.dpi":300, 'savefig.dpi':300})
        #     sns.set_theme()
        #     sns.set_style('ticks')
        #     fig, axes = plt.subplots(figsize=(20,7))
        #     fig.suptitle(sample, fontsize=20)
        #     sample_ratio = sample + "_ratio"

        #     ratio_plot = sns.scatterplot(data=cnr_df, x=cnr_df.index, 
        #         y=cnr_df[sample_ratio], size=0.015, hue=cnr_df.chromosome, 
        #         palette=palette_dict, alpha=0.4, edgecolor="none")
            
        #     # Setting y limits 
        #     ratio_plot.set(ylim=(0,3))
        #     ratio_plot.set_xticks(xtick_list)
        #     ratio_plot.set_xticklabels(unique_chromosomes, rotation=30, size=15)
        #     ratio_plot.set_yticks(ratio_plot.get_yticks())
        #     ratio_plot.set_yticklabels(ratio_plot.get_yticks(), size=15)
        #     ratio_plot.get_legend().remove()
        #     ratio_plot.set_xlabel("", fontsize=20)
        #     ratio_plot.set_ylabel("Ratio", fontsize=20)

        #     # Adding vertical lines to separate chromosomes 
        #     for xc in xtick_list:
        #         ratio_plot.axvline(x=xc, color="black")

        #     # Saving as png  
        #     ratio_plot.figure.savefig(p.sample_env[sample]['RATIO_PLOT'])


def do_ratio_same(row, median_sample, sample):
    ratio = row[sample]/median_sample
    return ratio

def do_ratio_ref (row, baseline, sample):
    a = row[baseline].to_numpy()
    median_baseline = np.median(a)
    if median_baseline == 0:
        ratio = 0
    else:
        ratio = row[sample]/median_baseline
    return ratio

def plot_normalization(df, sample, output_dir):
    '''Plotting coverage normalization 
    '''
    msg = " INFO: Plotting normalized coverage for sample " + sample
    print(msg)
    p.logging.info(msg)

    out_png = output_dir + "/" + sample + ".png"

    if not os.path.isfile(out_png):
        sns.set(font_scale=2)
        fig, axes = plt.subplots(2, 2, figsize=(25,22))
        fig.suptitle('GC-content & Mappability correction', fontsize=50)
        axes[0, 0] = sns.boxplot(ax=axes[0, 0], x="gc_integer", y="coverage", data=df, showfliers=False, palette="Blues")
        axes[0, 0].set_title("Raw coverage vs GC", fontsize=30)
        axes[0, 0].set_xticklabels(axes[0, 0].get_xticklabels(), rotation=30)
        axes[0, 0].set(xlabel='%GC')
        axes[0, 0].set(ylabel='Coverage')
        axes[0, 0].xaxis.set_major_locator(ticker.MultipleLocator(base=5))

        axes[0, 1] = sns.boxplot(ax=axes[0, 1], x="gc_integer", y="normalized_gc", data=df, showfliers=False, palette="Blues")
        axes[0, 1].set_title("GC-content corrected coverage", fontsize=30)
        axes[0, 1].set_xticklabels(axes[0, 1].get_xticklabels(), rotation=30)
        axes[0, 1].set(xlabel='%GC')
        axes[0, 1].set(ylabel='Coverage')
        axes[0, 1].xaxis.set_major_locator(ticker.MultipleLocator(base=5))

        axes[1, 0] = sns.boxplot(ax=axes[1, 0],x="map_integer", y="coverage", data=df, showfliers=False, palette="Blues")
        axes[1, 0].set_title("Raw coverage vs mappability", fontsize=30)
        axes[1, 0].set_xticklabels(axes[1, 0].get_xticklabels(), rotation=30)
        axes[1, 0].set(xlabel='%Mappability')
        axes[1, 0].set(ylabel='Coverage')
        axes[1, 0].set(xlim=(0,100)) 
        axes[1, 0].xaxis.set_major_locator(ticker.MultipleLocator(base=10))

        axes[1, 1] = sns.boxplot(ax=axes[1, 1],x="map_integer", y="normalized_map", data=df, showfliers=False, palette="Blues")
        axes[1, 1].set_title("GC-Mappability corrected coverage", fontsize=30)
        axes[1, 1].set_xticklabels(axes[1, 1].get_xticklabels(), rotation=30)
        axes[1, 1].set(xlabel='%Mappability')
        axes[1, 1].set(ylabel='Coverage')
        axes[1, 1].set(xlim=(0,100)) 
        axes[1, 1].xaxis.set_major_locator(ticker.MultipleLocator(base=10))

        fig.savefig(out_png)
    return out_png
 
def create_bins(bin_size):
    '''Create a BED file with genomewide bins
    '''

    # Set lowpass genomewide bins with a fixed bin size
    bashCommand = ('{} makewindows -g {} -w {}| {} intersect -a stdin -b {} -v > {}')\
    .format(p.system_env['BEDTOOLS'],p.aux_env['CHROM_SIZES'] ,\
    bin_size, p.system_env['BEDTOOLS'], p.aux_env['EXCLUDE_REGIONS'], p.analysis_env['LOWPASS_BINS'])

    msg = " INFO: Creating bins of " + str(bin_size) + " size (bp)"
    print(msg)
    p.logging.info(msg)
    p.logging.info(bashCommand)
    
    p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = p1.stdout.decode('UTF-8')
    error  = p1.stderr.decode('UTF-8')
    if not error:
      msg = " INFO: Bin regions successfully created"
      print(msg)
      logging.info(msg)        
    else:
      msg = " ERROR: Could not create bin regions"
      print (msg)
      logging.error(error)
      logging.error(msg)        

def annotate_mappability(input_bed):

    tmp_input_bed = input_bed.replace(".bed", ".tmp.bed")
    map_input_bed = input_bed.replace(".bed", ".map.bed")

    msg = " INFO: Extracting mappability"
    print(msg)
    p.logging.info(msg)

    a = pybedtools.BedTool(p.aux_env['MAPPABILITY_FILE'])
    b = pybedtools.BedTool(input_bed)
    c = a.intersect(b, wo=True, stream=True)
    map_dict = defaultdict(dict)

    for line in iter(c):
        line = str(line)
        line = line.rstrip()
        tmp  = line.split('\t')
        mappability = float(tmp[3])
        size        = int(tmp[6])-int(tmp[5])
        bases = int(tmp[7])
        marginal = ((mappability*bases)/size)*100
        coordinate = tmp[4]+"\t"+tmp[5]+"\t"+tmp[6]
        if not coordinate in map_dict:
            map_dict[coordinate] = marginal
        else: 
            map_dict[coordinate]+= marginal

    return map_dict

def annotate_gc(input_bed):
    '''Add gc content, return a dict with gc content
    '''    

    gc_input_bed = input_bed.replace(".bed", ".gc.bed")

    bashCommand = ('{} nuc -fi {} -bed {} > {} ')\
    .format(p.system_env['BEDTOOLS'], p.aux_env['GENOME_FASTA'], input_bed, gc_input_bed)

    msg = " INFO: Extracting gc content"
    print(msg)
    p.logging.info(msg)
    p.logging.info(bashCommand)
    
    gc_dict = defaultdict(dict)
    if not os.path.isfile(gc_input_bed):
        p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = p1.stdout.decode('UTF-8')
        error  = p1.stderr.decode('UTF-8')
        print(error)

    with open (gc_input_bed) as f:
        for line in f:
            if line.startswith('#'):
                continue
            tmp = line.split("\t")
            chr   = tmp[0]
            start = tmp[1]
            end   = tmp[2]
            gc = round(float(tmp[4]) * 100, 2)
            coordinate = chr + "\t" + start + "\t" + end
            gc_dict[coordinate] = gc
    return gc_dict

def extract_coverage(gc_dict, map_dict):
    '''Coverage extraction with megadepth 
    '''

    for sample in p.sample_env:

        msg = " INFO: Extracting coverage from sample " + sample
        print(msg)
        logging.info(msg)

        # Create COV_FOLDER if not present
        p.sample_env[sample]['COV_FOLDER'] = \
            p.sample_env[sample]['SAMPLE_FOLDER'] + "/" + "COV_FOLDER"

        cov_folder_path = Path(p.sample_env[sample]['COV_FOLDER'])
        if not cov_folder_path.is_dir():
            os.mkdir(cov_folder_path)

        p.sample_env[sample]['RAW_COVERAGE'] = \
            p.sample_env[sample]['COV_FOLDER'] + "/" + sample + ".raw.coverage.bed"

        bashCommand = ('{} {} --threads {} --annotation {} > {}')\
        .format(p.system_env['MEGADEPTH'],p.sample_env[sample]['READY_BAM'],\
        p.analysis_env['THREADS'],p.analysis_env['LOWPASS_BINS'],\
        p.sample_env[sample]['RAW_COVERAGE'])

        p.logging.info(bashCommand)
        
        p.sample_env[sample]['COV_METRICS'] =  p.sample_env[sample]['COV_FOLDER'] + "/" + \
            sample + ".cov.metrics.txt"  

        if not os.path.isfile(p.sample_env[sample]['RAW_COVERAGE']):
            p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output = p1.stdout.decode('UTF-8')
            error  = p1.stderr.decode('UTF-8')
            total_reads = ""
            if not error:
                msg = " INFO: Coverage extraction was successful"
                print(msg)
                logging.info(msg)        
            else:
                if re.search("error", error):
                    msg = " ERROR: Could not extract coverage"
                    print (msg)
                    logging.error(error)
                    logging.error(msg)
                else:
                    tmp = error.split("\n")
                    for line in tmp:
                        if line.startswith("Read"):
                            m = re.search(r'\d+', line)
                            total_reads = m.group()
            p.sample_env[sample]['TOTAL_READS'] = total_reads
            f = open(p.sample_env[sample]['COV_METRICS'], 'w')
            f.write(sample + "\t" + total_reads + "\n")
            f.close()
        else:
            with open(p.sample_env[sample]['COV_METRICS']) as f:
                for line in f:
                    tmp = line.split("\t")
                    p.sample_env[sample]['TOTAL_READS'] = tmp[1]
            f.close()

        annot_cov = p.sample_env[sample]['RAW_COVERAGE'].replace(".bed", ".annot.bed")
        o = open(annot_cov, 'w')
        o.write("Chromosome"+"\t"+"Start"+"\t"+"End"+"\t"+"bin"+"\t"+"gc"+"\t"+"map"+"\t"+"coverage"+"\n")
        bin_n = 0
        with open (p.sample_env[sample]['RAW_COVERAGE']) as f:
            for line in f:
                if line.startswith("Chromosome"):
                    continue              
                bin_n+=1
                line = line.rstrip("\n")
                tmp  = line.split("\t")
                chr   = tmp[0]
                start = tmp[1]
                end   = tmp[2]
                if len(tmp) == 4:
                    cov = tmp[3]
                else:
                    cov = tmp[6]

                coordinate = chr + "\t" + start + "\t" + end
                gc  = str(gc_dict[coordinate])
                mappability = str(map_dict[coordinate])

                out_bin = "bin_" + str(bin_n)
                o.write(chr+"\t"+start+"\t"+end+"\t"+out_bin+"\t"+gc+"\t"+mappability+"\t"+cov+"\n")
        f.close()
        o.close()
        os.remove(p.sample_env[sample]['RAW_COVERAGE'])
        os.rename(annot_cov, p.sample_env[sample]['RAW_COVERAGE'])

def normalize(df, sample, median_cov, fields):

    msg = " INFO: Normalizing  " + str(fields)
    print(msg)
    logging.info(msg)

    df['cov_bylib'] = (df['coverage']/median_cov)*100
    cov_target = 'cov_bylib'
      
    for field in fields:

        field_int = field + "_integer"

        # Get integer field (gc or map) value 
        df[field_int] = df[field].apply(int)

        median_field_cov = "median_" + field + "_cov"

        # Group coverage by field and calculate the median
        df[median_field_cov]  = df.groupby(field_int)[cov_target].transform("median")
        median_cov = df[(df['Chromosome'] != "chrX") & (df['Chromosome'] != "chrY")][cov_target].median()

        # Apply rolling median 
        normalized_field = "normalized_" + field
        df[normalized_field]  = (df[cov_target]*median_cov)/df[median_field_cov]
        cov_target = normalized_field

    return df

def segment_coverage(n_segments, alpha):
    '''segment with CBS
    '''

    for sample in p.sample_env:

        to_segment = p.sample_env[sample]['RATIO_FILE'].replace(".bed", ".tosegment.bed")
        o = open(to_segment, 'w')
        with open (p.sample_env[sample]['RATIO_FILE']) as f:
            for line in f:
                if line.startswith('Chromosome'):
                    continue
                line = line.rstrip("\n")
                tmp  = line.split("\t")
                o.write(tmp[0]+"\t"+tmp[1]+"\t"+tmp[2]+"\t"+tmp[3]+"\t"+tmp[-1]+"\n")
        o.close()

        rscript = p.sample_env[sample]['RATIO_FILE'].replace(".ratios.bed", ".CBS.R")
        p.sample_env[sample]['SEGMENT_FILE'] = p.sample_env[sample]['RATIO_FILE'].replace(".bed", ".segment.bed")
        r = open(rscript, 'w')
        r.write("library(DNAcopy)"+"\n")
        line = "cn <- read.table(\"{}\", header=F)".format(to_segment)
        r.write(line + "\n")
        line = "CNA.object <-CNA( genomdat = cn[,5], chrom = cn[,1], maploc = cn[,2], data.type = \'logratio\')"
        r.write(line + "\n")
        line = "CNA.smoothed <- smooth.CNA(CNA.object)"
        r.write(line + "\n")
        line = "segs <- segment(CNA.object, verbose=0, min.width={}, alpha = {})".format(n_segments, alpha)
        r.write(line + "\n")
        line = "segs2=segs$output"
        r.write(line + "\n")
        line = "write.table(segs2[,2:6], file=\"{}\",row.names=F, col.names=F, quote=F, sep=\"\t\")".format(p.sample_env[sample]['SEGMENT_FILE'])
        r.write(line + "\n")
        r.close()

        if not os.path.isfile(p.sample_env[sample]['SEGMENT_FILE']):
            bashCommand = ('Rscript {} ').format(rscript)
            p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output = p1.stdout.decode('UTF-8')
            error  = p1.stderr.decode('UTF-8')

        # add header to segment file
        idx = 0
        tmp_segment = p.sample_env[sample]['SEGMENT_FILE'].replace(".bed", ".tmp.bed")
        o = open(tmp_segment, 'w')
        with open (p.sample_env[sample]['SEGMENT_FILE']) as f:
            for line in f:
                line = line.rstrip('\n')
                if idx == 0:
                    if not line.startswith("Chromosome"):
                        o.write("Chromosome\tStart\tEnd\tSegments\tRatio\n")
                o.write(line+"\n")
                idx+=1
        f.close()
        o.close()

        os.remove(p.sample_env[sample]['SEGMENT_FILE'])
        os.rename(tmp_segment, p.sample_env[sample]['SEGMENT_FILE'])


def mean_zscore():
    pass

def remove_bed_header(file, pattern):

    no_header_file = file.replace(".bed", ".noheader.bed")
    nh = open(no_header_file, 'w')
    with open(file) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(pattern):
                continue
            nh.write(line+"\n")
    f.close()        
    nh.close()
    return no_header_file

def call_cnvs(del_threshold, dup_threshold, z_score):
    '''Calling CNVs
    ''' 
    del_threshold = float(del_threshold)
    dup_threshold = float(dup_threshold)

    for sample in p.sample_env:
        p.sample_env[sample]['CALLS_FILE'] \
             = p.sample_env[sample]['COV_FOLDER'] + "/" + sample + ".calls.bed"

        ratio_no_header = remove_bed_header(p.sample_env[sample]['RATIO_FILE'], "Chromosome")
        seg_no_header   = remove_bed_header(p.sample_env[sample]['SEGMENT_FILE'], "Chromosome")

        # We will use plain dicts and bedtools stuff, instead of pandas
        tmp_calls = p.sample_env[sample]['CALLS_FILE'].replace(".bed", ".tmp.bed")
        o = open(tmp_calls, 'w')

        with open (p.sample_env[sample]['SEGMENT_FILE'])  as seg:
            for line in seg:
                line = line.rstrip("\n")
                #chr11	19200000	26550000	133	0.5013	DEL
                if line.startswith('Chromosome'):
                    continue
                tmp = line.split('\t')
                chr   = tmp[0]
                start = tmp[1]
                end   = tmp[2]
                n_bins= tmp[3]
                ratio = tmp[4]
                cnvtype = ""
                if float(ratio) >= dup_threshold:
                    cnvtype = "DUP"
                    o.write(line + "\t" + cnvtype + "\n")
                elif float(ratio) < del_threshold:
                    cnvtype = "DEL"
                    o.write(line + "\t" + cnvtype + "\n")
        seg.close()
        o.close()

        calls_dict = defaultdict(dict)
        a = pybedtools.BedTool(ratio_no_header)
        b = pybedtools.BedTool(tmp_calls)
        c = a.intersect(b, wo=True, stream=True)
        for line in iter(c):
            line = str(line)
            line = line.rstrip()
            tmp  = line.split('\t')

            variant = tmp[6]+"\t"+tmp[7]+"\t"+tmp[8]+"\t"+tmp[9]+"\t"+tmp[10]+"\t"+tmp[11]
            if not variant in calls_dict:
                calls_dict[variant] = list()
            else:
                calls_dict[variant].append(float(tmp[5])) 

        o = open(p.sample_env[sample]['CALLS_FILE'], 'w')
        o.write("Chromosome\tStart\tEnd\tRegions\tRatio\tCnvtype\tStd")

        for variant in calls_dict:
            arr = np.array(calls_dict[variant])
            std = round(np.std(arr), 3)
            o.write(variant + "\t" + str(std) + "\n")
        o.close()

        os.remove(tmp_calls)
        os.remove(ratio_no_header)
        os.remove(seg_no_header)

        # seg_df = pd.read_csv(p.sample_env[sample]['SEGMENT_FILE'], 
        #     sep="\t", names=['chromosome', 'start', 'end', 'n_bins', 'ratio'] )

        # calls_df = seg_df.loc[ (seg_df['ratio'] <= del_threshold) | (seg_df['ratio'] >= dup_threshold)]

        # calls_df['cnvtype'] = np.where(calls_df['ratio'] >= dup_threshold, 'DUP', 'DEL')

        #calls_df.sort_values(by=['chromosome'], inplace=True)
        #calls_df.reindex(index=order_by_index(calls_df.index, index_natsorted(calls_df.chromosome)))
        #calls_df.to_csv(p.sample_env[sample]['CALLS_FILE'], sep="\t", mode='w', index=None)

