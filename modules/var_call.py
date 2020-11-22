#!/usr/bin/env python3

import os
import sys
import re
import logging
import gzip
import shutil
import csv
from collections import defaultdict
from pathlib import Path
from scipy import stats
import numpy as np
import subprocess
from modules import utils as u
from modules import params as p
from modules import trimming as t

def do_var_call():

    if p.analysis_env['VARIANT_CLASS'] == 'somatic':        
       # do_freebayes()
        do_mutect2()
        do_manta()
        do_cnvkit()
        
def do_cnvkit():

  # Create accessible regions
  paneldir = p.analysis_env['PANEL_WORKDIR']
  access_regions_hg19 = paneldir + "/" + "access.hg19.bed"

  # Create target regions 
  target_bed = paneldir + "/" + p.analysis_env['PANEL_NAME'].replace(".bed", ".target.bed")
  target_bed_name = p.analysis_env['PANEL_NAME'].replace(".bed", ".target.bed")

  # Create off-target regions 
  antitarget_bed = paneldir + "/" + p.analysis_env['PANEL_NAME'].replace(".bed", ".antitarget.bed")
  antitarget_bed_name = p.analysis_env['PANEL_NAME'].replace(".bed", ".antitarget.bed")

  if not os.path.isfile(access_regions_hg19):

    bashCommand = ('python3 {} access {} -o {}').format(p.system_env['CNVKIT'], \
      p.aux_env['GENOME_FASTA'], access_regions_hg19)

    msg = " INFO: Creating access regions file for cnvkit"
    print(msg)
    logging.info(msg)
    logging.info(bashCommand)

    process = subprocess.Popen(bashCommand,#.split(),
      shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()

    if not error.decode('UTF-8'):
      msg = " INFO: Access regions successfully created"
      print(msg)
      logging.info(msg)
    else:
      msg = " ERROR: Could not create accessible regions file"
      print (msg)
      logging.error(msg)

  if not os.path.isfile(target_bed):

    # Create on/offtarget bins with autobin command
    bashCommand = ('python3 {} target {} -o {}').format(p.system_env['CNVKIT'], \
     p.analysis_env['PANEL'], target_bed)

    msg = " INFO: Creating target regions file"
    print (msg)
    logging.info(msg)
    logging.info(bashCommand)

    process = subprocess.Popen(bashCommand,#.split(),
      shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
      
    if not error.decode('UTF-8'):
      msg = " INFO: target regions successfully created"
      print(msg)
      logging.info(msg)    
    else:
      msg = " ERROR: Could not create bin regions file"
      print (msg)

  if not os.path.isfile(antitarget_bed):

    # Create on/offtarget bins with autobin command
    bashCommand = ('python3 {} antitarget  {} -o {}').format(p.system_env['CNVKIT'] , \
     access_regions_hg19, antitarget_bed)

    msg = " INFO: Creating off-target regions file"
    print (msg)
    logging.info(msg)     
    logging.info(bashCommand)

    process = subprocess.Popen(bashCommand,#.split(),
      shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
      
    if not error.decode('UTF-8'):
      msg = " INFO: off-target regions successfully created"
      print(msg)
      logging.info(msg)       
    else:
      msg = " ERROR: Could not create off-target regions file"
      print (msg)
      logging.error(msg)

  if not os.path.isfile(p.aux_env['NORMALS_REF_CNA']):

    msg = " ERROR: Could not find reference for normals"
    print (msg)
    logging.error(msg)
    sys.exit()

  coverage_list = [] 
  anticoverage_list = [] 
  
  # Extract coverage
  for sample in p.sample_env:

      p.sample_env[sample]['CNV_FOLDER'] = \
       p.sample_env[sample]['VCF_FOLDER'] + "/" + "CNV_FOLDER"

      if not os.path.isdir(p.sample_env[sample]['CNV_FOLDER']):
        os.mkdir(p.sample_env[sample]['CNV_FOLDER'])

      bam = p.sample_env[sample]['READY_BAM']   

      coverage_file = p.sample_env[sample]['CNV_FOLDER'] + "/" + sample + ".targetcoverage.cnn"
      coverage_list.append(coverage_file)

      bashCommand = ('python3 {} coverage {} {} -o {} ').format(p.system_env['CNVKIT'], \
      bam, target_bed, coverage_file)

      if not os.path.isfile(coverage_file):

        msg = " INFO: Extracting on-target coverage for sample " + sample
        print(msg)
        logging.info(msg)
        logging.info(bashCommand)

        process = subprocess.Popen(bashCommand,#.split(),
          shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()
          
        if not error.decode('UTF-8'):
          msg = " INFO: Coverage was extracted for sample " + sample
        else:
          msg = " ERROR: Could not create bin regions files"
          print (msg)
          logging.error(msg)
      else:
        msg = " INFO: Skipping coverage extraction for sample " + sample
        print (msg)
        logging.info(msg)

      anticoverage_file = p.sample_env[sample]['CNV_FOLDER'] + "/" + sample + ".antitargetcoverage.cnn"
      anticoverage_list.append(anticoverage_file)

      bashCommand = ('python3 {} coverage {} {} -o {}').format(p.system_env['CNVKIT'], \
      bam, antitarget_bed, anticoverage_file)

      if not os.path.isfile(anticoverage_file):
        process = subprocess.Popen(bashCommand,#.split(),
          shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()

        msg = " INFO: Extracting off-target coverage for sample " + sample
        print(msg)
        logging.info(msg)
        logging.info(bashCommand)

        if not error.decode('UTF-8'):
          pass
        else:
          msg = " ERROR: Could not create bin regions files"
          print (msg)
          logging.error(msg)
      else:
        msg = " INFO: Skipping coverage extraction for sample " + sample
        print (msg)
        logging.error(msg)

  # Segment and calling
  for sample in p.sample_env:

      coverage_str = ' '.join(coverage_list)
      coverage_file = p.sample_env[sample]['CNV_FOLDER'] + "/" + sample + ".targetcoverage.cnn"

      anticov_str = ' '.join(anticoverage_list)
      anticov_file = p.sample_env[sample]['CNV_FOLDER'] + "/" + sample + ".antitargetcoverage.cnn"

      p.sample_env[sample]['CNV_VCF'] = p.sample_env[sample]['VCF_FOLDER'] + "/" + sample + ".CNV.vcf"

      # Fix and ratio creation
      ratio_file = p.sample_env[sample]['CNV_FOLDER'] + "/" +  sample + ".cnr"
      bashCommand = ('python3 {} fix {} {} {} -o {}').format(p.system_env['CNVKIT'], \
      coverage_file, coverage_file, p.aux_env['NORMALS_REF_CNA'], ratio_file)

      if not os.path.isfile(ratio_file):
        process = subprocess.Popen(bashCommand,#.split(),
          shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()
        
        msg = " INFO: Fixing (normalizing) sample " + sample
        print(msg)
        logging.info(msg)
        logging.info(bashCommand)

        if not error.decode('UTF-8'):
          msg = " INFO: Fix (normalization) on sample " + sample + " ended successfully"
          print(msg)
          logging.info(msg)
        else:
          msg = " ERROR: Could not fix sample " + sample
          print (msg)
          logging.error(msg)
      else:
        msg = " INFO: Skipping fix of sample " + sample
        print (msg)
        logging.info(msg)

      # Fix cnr weird results
      seen = defaultdict(dict)
      tmp_cnr = ratio_file.replace(".cnr", ".tmp.cnr")
      o = open(tmp_cnr, "w")
      with open (ratio_file) as f:
       for line in f:
         line = line.rstrip("\n")
         if line.startswith("chromosome"):
           o.write(line+ "\n")
         else:
           tmp = line.split("\t")
           coordinates = tmp[0]+tmp[1]
           if not coordinates in seen:
             seen[coordinates] = 0
           else:
             if seen[coordinates] > 0:
               continue
           o.write(line + "\n")
           seen[coordinates]+=1
      f.close()
      o.close()
      os.remove(ratio_file)
      os.rename(tmp_cnr, ratio_file)

      # Segment
      segment_file = p.sample_env[sample]['CNV_FOLDER'] + "/" +  sample + ".cns"
      bashCommand = ('python3 {} segment {} -o {}').format(p.system_env['CNVKIT'], \
      ratio_file, segment_file)
      if not os.path.isfile(segment_file):

        msg = " INFO: Segmenting sample " + sample
        print(msg)
        logging.info(msg)
        logging.info(bashCommand)
        process = subprocess.Popen(bashCommand,#.split(),
          shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()
          
        if not error.decode('UTF-8'):
          msg = " INFO: Segmentation of sample " + sample + " ended successfully"
          print (msg)
          logging.info(msg)
        else:
          msg = " ERROR: Could not segment sample " + sample
          print (msg)
          logging.info(msg)
      else:
        msg = " INFO: Skipping segmentation of sample " + sample
        print (msg)
        logging.info(msg)

      # Export segmented ratios to DNAcopy (.seg) format
      export_segfile = p.sample_env[sample]['CNV_FOLDER'] + "/" +  sample + ".seg"
      o = open(export_segfile, "w")
      o.write("ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean"+"\n")
      with open (segment_file) as f:
        for line in f:
          line = line.rstrip("\n")
          if line.startswith("chromosome"):
            continue
          tmp = line.split("\t")
          #chromosome	start	end	gene	log2	depth	probes	weight
          o.write(sample+"\t"+tmp[0]+"\t"+tmp[1]+"\t"+tmp[2]+"\t"+tmp[6]+"\t"+tmp[4]+"\n")
      o.close()

     # vcf_file = tumors_dir + "/" + sample + ".mutect2.vcf.gz"
      # purity = do_pureCN(vcf_file, export_segfile, ratio_file, sample)
     # purity = 0.9
      # print(purity)

      purity = 1.00
      if sample in p.lab_data:
        purity = p.lab_data[sample]['PURITY']
        purity = purity.replace("%", "")
        if purity == "." or purity == 0:
          purity = 1
        purity = str(float(purity)/100)
        msg = " INFO: " + sample + " purity: " + purity
        print(msg)
        logging.info(msg)

      # # Call CNVs
      call_file = p.sample_env[sample]['CNV_FOLDER'] + "/" +  sample + ".calls.cns"
      bashCommand = ('python3 {} call {} -m clonal --purity {} -o {}').format(p.system_env['CNVKIT'], \
      segment_file, purity, call_file)

      if not os.path.isfile(call_file):

        msg = " INFO: Calling CNAs on sample " + sample
        print(msg)
        logging.info(msg)
        logging.info(bashCommand)
        
        process = subprocess.Popen(bashCommand,#.split(),
          shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()
          
        if not error.decode('UTF-8'):
          msg = " INFO: Calling CNAs on sample ended successfully" + sample
          print(msg)
          logging.info(msg)        
        else:
          msg = " ERROR: Could not call CNA on sample " + sample
          print (msg)
          logging.error(msg)
      else:
        msg = " INFO: Skipping CNA calling of sample " + sample
        print (msg)
        logging.info(msg)

      if not os.path.isfile(p.sample_env[sample]['CNV_VCF']):

        bashCommand = ('python3 {} export vcf {} -y  -o {}').format(p.system_env['CNVKIT'], \
          call_file, p.sample_env[sample]['CNV_VCF'])
        msg = " INFO: Exporting CNA calls to vcf for sample " + sample
        print(msg)
        logging.info(msg)
        logging.info(bashCommand)
        
        process = subprocess.Popen(bashCommand,#.split(),
          shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()
          
        if not error.decode('UTF-8'):
          msg = " INFO: Exporting CNA calls to vcf for sample " +sample + " ended successfully"
          print(msg)
          logging.info(msg)        
        else:
          msg = " ERROR: Could not export CNA on sample " + sample
          print (msg)
          logging.error(msg)
      else:
        msg = " INFO: Skipping CNA exporting to vcf for sample " + sample
        print (msg)
        logging.info(msg)        


      plot_scatter( ratio_file, call_file, sample, p.sample_env[sample]['CNV_FOLDER'])

      for gene in p.cna_plot_env:
        plot = p.sample_env[sample]['CNV_FOLDER'] + "/" + gene + ".png"
        chrom = p.cna_plot_env[gene]['chromosome']
        if not 'chr' in chrom:
          chrom = 'chr' + chrom
        bashCommand = ('python3 {} scatter {} -s {} -c {} -g {} --y-min -2 --segment-color red --by-bin -o {}').format(p.system_env['CNVKIT'], \
          ratio_file, call_file, chrom, gene, plot)
        if not os.path.isfile(plot):
          process = subprocess.Popen(bashCommand,#.split(),
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
          output, error = process.communicate()

def plot_scatter(cnr_file, cns_file, sample, outdir):

  tmp_cnr = outdir + "/" + "tmpcnr.bed"
  r = open(tmp_cnr, "w")
  with open (cnr_file, "r") as f:
    for line in f:
      line = line.rstrip("\n")
      if line.startswith("chromosome"):
        continue
      else:
        r.write(line+"\n")
  r.close()
  f.close()

  tmp_cns = outdir + "/" + "tmpcns.bed"
  s = open(tmp_cns, "w")
  with open (cns_file, "r") as f:
    for line in f:
      line = line.rstrip("\n")
      if line.startswith("chromosome"):
        continue
      else:
        s.write(line+"\n")
  s.close()
  awk = "awk \'{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$12 }\'"
  intersect_file = outdir + "/" + "intersect.bed"
  bashCommand = ('bedtools intersect -a {} -b {} -wa -wb | {}  > {}').format(tmp_cnr, \
    tmp_cns, awk, intersect_file)

  process = subprocess.Popen(bashCommand,#.split(),
    shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  output, error = process.communicate()

  os.remove(tmp_cns)
  os.remove(tmp_cnr)

  genomewide_plot = outdir + "/" + sample + ".genomewide.png"
  rscript = outdir + "/" + "plotScatter.R"
  f = open(rscript, "w")

  f.write("library(ggplot2)"+"\n")
  f.write("library(plyr)"+"\n")
  f.write("library(grid)"+"\n")
  f.write("library(gridExtra)"+"\n")
  f.write("library(gtools)"+"\n")
  f.write("mydata<-read.table(file=\"")
  f.write(intersect_file + "\", sep =\"\t\",check.names = FALSE, header=FALSE)\n")
  f.write("attach(mydata)"+"\n")
  f.write("chrsort <- mixedsort(levels(factor(V1)))"+"\n")
  f.write("mydata$V1<-factor(mydata$V1, levels=chrsort)"+"\n")
  f.write("mydata<-mydata[order(mydata$V1),]"+"\n")
  f.write("xval<-seq(1,length(mydata$V1))"+"\n")
  f.write("mycount<-count(mydata$V1)"+"\n")
  f.write("attach(mycount)"+"\n")
  f.write("mysum<-cumsum(freq)"+"\n")
  f.write("mysum<-append(mysum,0, after=0)"+"\n")
  f.write("mysum<-mysum[-length(mysum)]"+"\n")
  f.write("png(\"" + genomewide_plot + "\", res=137, width=1600, height=600)"+"\n")
  f.write("myplot<-ggplot(mydata,aes(x=xval, y =V5), group=V1) + geom_point(size=0.4, aes(colour=V1))")
  f.write("+ xlab(\"Ratio\") + geom_point(aes(y=V6), size=0.2, colour=\"red\")")
  f.write("+ scale_x_continuous(\"chromosome\",breaks = mysum, labels=chrsort) + theme_bw()")
  f.write("+ ylim(-3,3)")
  f.write("+ theme(panel.grid.major.x = element_line(colour = \"black\"),panel.grid.minor.y=element_blank(), panel.grid.major.y = element_blank(),  panel.background = element_rect(colour = \"black\", size=1))")
  f.write("+ theme(axis.text.x = element_text(size = 12, hjust = 1, angle=45))")
  f.write("+ theme(axis.text.y = element_text(size = 12, hjust = 1))")
  f.write("+ scale_colour_manual(values=c(\"#4f6b76\", \"#b2bbc0\", \"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\",\"#4f6b76\", \"#b2bbc0\", \"#4f6b76\", \"#b2bbc0\", \"#4f6b76\", \"#b2bbc0\", \"#4f6b76\", \"#b2bbc0\", \"#4f6b76\", \"#b2bbc0\")) + theme(legend.position=\"none\")")
  f.write("+ geom_hline(yintercept = -1, colour=\"blue\", linetype=\"dashed\", size=0.3)")
  f.write("+ geom_hline(yintercept = 0.58, colour=\"red\", linetype=\"dashed\", size=0.3) +ylab(\"Ratios\")")
  f.write("+ ggtitle(\"" + sample + "\") + theme(plot.title = element_text(hjust = 0.5))"+"\n")
  f.write("myplot" +"\n")
  f.write("dev.off()" +"\n")
  f.close()

  bashCommand = ('Rscript {} ').format(rscript)
  process = subprocess.Popen(bashCommand,#.split(),
    shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  output, error = process.communicate()

  os.remove(intersect_file)

def do_pureCN(vcf, segfile, ratiofile, sample):

 normaldb = "~/Escriptori/PIPELINE_CANCER/BIN_FOLDER/PureCN/extdata/NormalDB.R"
 purecn   = "~/Escriptori/PIPELINE_CANCER/BIN_FOLDER/PureCN/extdata/PureCN.R"
 bashCommand = ('Rscript {} --out {} --sampleid {} --tumor {} --segfile {} --vcf {}'
 ' --genome hg19 --funsegmentation Hclust --force --postoptimize --seed 123 --sex=\"diploid\"').format(purecn, \
  outdir, sample, ratiofile, segfile, vcf)

 summary_purecn = outdir + "/" + sample + ".csv"

 if not os.path.isfile(summary_purecn):
  process = subprocess.Popen(bashCommand,#.split(),
    shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  output, error = process.communicate()

 purity = ""
 with open(summary_purecn) as f:
   for line in f:
     line= line.rstrip("\n")
     tmp = line.split(",")
     if line.startswith("sampleid"):
       continue
     else:
       purity = tmp[1]
 return purity



def do_freebayes():

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

            bashCommand = ('{} -f {}  {} >  {}').format(p.system_env['FREEBAYES'], \
              p.aux_env['GENOME_FASTA'], p.sample_env[sample]['READY_BAM'],p.sample_env[sample]['FREEBAYES_VCF'] )
            logging.info(bashCommand)
            process = subprocess.Popen(bashCommand,#.split(),
                shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            output, error = process.communicate()

            if not error.decode('UTF-8'):
                msg = " INFO: Freebayes was executed successfully on " + sample
                print (msg)
                logging.info(msg)
            else:
                msg = " ERROR: Something happeend when using Freebayes on " + sample
                print (msg)
                logging.error(msg)

def do_cna():

    for sample in p.sample_env:
        p.sample_env[sample]['CNA_FOLDER'] = \
            p.sample_env[sample]['VCF_FOLDER'] + "/" "CNA_FOLDER"

        cna_folder_path = Path(p.sample_env[sample]['CNA_FOLDER'])
        if not cna_folder_path.is_dir():
            os.mkdir(cna_folder_path)

def do_manta():

    for sample in p.sample_env:

        # Create BAM_FOLDER if not present
        p.sample_env[sample]['VCF_FOLDER'] = \
         p.sample_env[sample]['SAMPLE_FOLDER'] + "/" + "VCF_FOLDER"

        vcf_folder_path = Path(p.sample_env[sample]['VCF_FOLDER'])
        if not vcf_folder_path.is_dir():
            os.mkdir(vcf_folder_path)

        # This is the predefined output from manta
        p.sample_env[sample]['MANTA_TUMOR'] = \
            p.sample_env[sample]['VCF_FOLDER'] + "/results/variants/tumorSV.vcf.gz"
        p.sample_env[sample]['MANTA_TUMOR_NAME'] = "tumorSV.vcf.gz"

        # Renamed manta output for our sample
        p.sample_env[sample]['MANTA_VCF'] = \
            p.sample_env[sample]['VCF_FOLDER'] + "/" + sample + ".manta.vcf"
        p.sample_env[sample]['MANTA_VCF_NAME'] = sample + ".manta.vcf"

        # Ready SV file
        p.sample_env[sample]['READY_SV_VCF'] = \
            p.sample_env[sample]['VCF_FOLDER'] + "/" + sample + ".manta.vcf"
        p.sample_env[sample]['READY_SV_VCF_NAME'] = sample + ".manta.vcf"

        if not os.path.isfile(p.sample_env[sample]['MANTA_VCF']):
            
            bashCommand = ('{} --tumorBam {} --referenceFasta {} --exome --runDir {}').format(
             p.system_env['MANTA_CONFIG'], p.sample_env[sample]['READY_BAM'], \
             p.aux_env['GENOME_FASTA'], p.sample_env[sample]['VCF_FOLDER'] )
            logging.info(bashCommand)
            process = subprocess.Popen(bashCommand,#.split(),
                shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            output, error = process.communicate()

            if not error.decode('UTF-8'):
                msg = " INFO: Manta config was created successfully for sample " + sample
                print (msg)
                logging.info(msg)
            else:
                msg = " ERROR: Manta config could not be created for sample " + sample
                print (msg)
                logging.error(msg)

            run_manta_script = p.sample_env[sample]['VCF_FOLDER'] + "/" + "runWorkflow.py"

            bashCommand = ('{} --j {}').format(run_manta_script, p.analysis_env['THREADS'])
            msg = " INFO: Calling SVs with Manta on sample " + sample
            print(msg)
            logging.info(msg)
            logging.info(bashCommand)

            process = subprocess.Popen(bashCommand,#.split(),
                shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            output, error = process.communicate()

            if not error.decode('UTF-8'):
                msg = " INFO: Manta was executed successfully for sample " + sample
                print (msg)
                logging.info(msg)
            else:
                msg = " ERROR: Manta could not be executed for sample " + sample
                print (msg)
                logging.error(msg)

        if not os.path.isfile(p.sample_env[sample]['MANTA_VCF']):
            input = gzip.GzipFile(p.sample_env[sample]['MANTA_TUMOR'], 'rb')
            s = input.read()
            input.close()

            output = open(p.sample_env[sample]['MANTA_VCF'], 'wb')
            output.write(s)
            output.close()
        else:
            msg = " INFO: Skipping Manta analysis for "+ sample
            print(msg)
            logging.info(msg)

def do_mutect2():

    vcf_list = [] 
    for sample in p.sample_env:
        print(sample)
        # Create BAM_FOLDER if not present
        p.sample_env[sample]['VCF_FOLDER'] = \
         p.sample_env[sample]['SAMPLE_FOLDER'] + "/" + "VCF_FOLDER"

        vcf_folder_path = Path(p.sample_env[sample]['VCF_FOLDER'])
        if not vcf_folder_path.is_dir():
            os.mkdir(vcf_folder_path)

        p.sample_env[sample]['MUTECT2_VCF'] = p.sample_env[sample]['VCF_FOLDER'] + \
            "/" + sample + ".mutect2.vcf"
        p.sample_env[sample]['MUTECT2_VCF_NAME'] =  sample +  ".mutect2.vcf"
        vcf_list.append(p.sample_env[sample]['MUTECT2_VCF'])

        p.sample_env[sample]['READY_SNV_VCF'] = p.sample_env[sample]['MUTECT2_VCF']
        p.sample_env[sample]['READY_SNV_VCF_NAME'] = p.sample_env[sample]['MUTECT2_VCF_NAME']
        
        filtered_vcf = p.sample_env[sample]['MUTECT2_VCF'].replace(".vcf", ".filtered.vcf")
        filtered_vcf_name = p.sample_env[sample]['MUTECT2_VCF_NAME'].replace(".vcf", ".filtered.vcf")

        bashCommand = ('{} run -v {}:/bam_data/ -v {}:/vcf_data/ -v {}:/bundle/ '
            '-v {}:/panel_data/ -v {}:/gnomad_data/ -it {} gatk Mutect2 -I /bam_data/{} '
            '-O /vcf_data/{} -L /panel_data/{} -R /bundle/{} --genotype-germline-sites --germline-resource /gnomad_data/{}'
            ' -max-mnp-distance 0 --native-pair-hmm-threads {}'.format(p.system_env['DOCKER'],
            p.sample_env[sample]['BAM_FOLDER'], p.sample_env[sample]['VCF_FOLDER'], \
            p.defaults['BUNDLE_FOLDER'], p.analysis_env['PANEL_WORKDIR'], p.defaults['GNOMAD_FOLDER'], \
            p.docker_env['GATK'], p.sample_env[sample]['READY_BAM_NAME'], 
            p.sample_env[sample]['MUTECT2_VCF_NAME'],p.analysis_env['PANEL_LIST_NAME'], \
            p.aux_env['GENOME_NAME'], p.analysis_env['GNOMAD_AF_VCF_NAME'], p.analysis_env['THREADS']
            ))

        if not os.path.isfile(p.sample_env[sample]['MUTECT2_VCF']):
            msg = " INFO: Running Mutect2 for sample "+ sample
            print(msg)
            logging.info(msg)
            logging.info(bashCommand)

            process = subprocess.Popen(bashCommand,#.split(),
                shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = process.communicate()

            if not error.decode('UTF-8'):
                if re.search(r'Mutect2 done.', output.decode('UTF-8')):
                    msg = " INFO: Mutect2 varcall ended OK for sample " + sample
                    print (msg)
                    logging.info(msg)
                    ok = True
                else:
                    msg = " ERROR: Something went wrong with Mutect2 varcall for sample " + sample
                    print (msg)
                    logging.error(msg)
            else:
                msg = " ERROR: Something went wrong with Mutect2 varcall for sample " + sample
                print (msg)
                logging.error(msg)

            # Filtering Mutect2 calls
            bashCommand = ('{} run -v {}:/vcf_data/ -v {}:/bundle/ '
                '-it {} gatk FilterMutectCalls -R /bundle/{} '
                '-V /vcf_data/{} -O /vcf_data/{}'.format(p.system_env['DOCKER'],
                p.sample_env[sample]['VCF_FOLDER'], p.defaults['BUNDLE_FOLDER'],
                p.docker_env['GATK'],  p.aux_env['GENOME_NAME'], \
                p.sample_env[sample]['MUTECT2_VCF_NAME'], \
                filtered_vcf_name
                ))           
            msg = " INFO: Running FilterMutectCalls for sample "+ sample
            print(msg)
            logging.info(msg)
            logging.info(bashCommand)

            process = subprocess.Popen(bashCommand,#.split(),
                shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            output, error = process.communicate()
            if not error.decode('UTF-8'):
                if re.search(r'FilterMutectCalls done.', output.decode('UTF-8')):
                    msg = " INFO: FilterMutectCalls varcall ended OK for sample " + sample
                    print (msg)
                    logging.info(msg)
                    ok = True
                else:
                    msg = " ERROR: Something went wrong with FilterMutectCalls varcall for sample " + sample
                    print (msg)
                    logging.error(msg)
            else:
                msg = " ERROR: Something went wrong with FilterMutectCalls varcall for sample " + sample
                print (msg)
                logging.error(msg)
            os.remove(p.sample_env[sample]['MUTECT2_VCF'])
            os.rename(filtered_vcf, p.sample_env[sample]['MUTECT2_VCF'])
        else:
            msg = " INFO: Skipping Mutect2 analysis for "+ sample
            print(msg)
            logging.info(msg)


#     p.analysis_env['COHORT_MUTECT2_VCF'] = \
#         p.analysis_env['OUTDIR'] + "/" + "combined.mutect2.vcf.gz"

#     vcf_str = ' --variant /vcf_data/'.join(vcf_list)
#     bashCommand = ('{} run -v {}:/vcf_data/ -v {}:/bundle/ '
#         '-it {} gatk CombineGVCFs -R {} {}'.format(p.system_env['DOCKER'],
#         p.sample_env[sample]['VCF_FOLDER'], p.aux_env['GENOME_NAME'], ))

#  gatk CombineGVCFs \
#    -R reference.fasta \
#    --variant sample1.g.vcf.gz \
#    --variant sample2.g.vcf.gz \
#    -O cohort.g.vcf.gz