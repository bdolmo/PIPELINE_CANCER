#!/usr/bin/env python3

import sys
import re
import gzip
import binascii
import os.path
import glob
import json
import requests
from os import path
from collections import defaultdict
from pprint import pprint
import json
from pathlib import Path
import subprocess
from modules import params as p
from modules import utils as u
from modules import trimming as t

from flask import Flask
from flask import request, render_template, url_for, redirect, flash
from flask_wtf import FlaskForm
import sqlite3
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import Table, Column, Float, Integer, String, MetaData, ForeignKey

basedir = os.path.abspath(os.path.dirname(__file__))

# global app 
#global db
app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = ""
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)

def init():

    app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///' + p.analysis_env['DB_DIR'] + "/NGS.db" 
    app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
    app.config.update(dict(
        SECRET_KEY="powerful secretkey",
        WTF_CSRF_SECRET_KEY="a csrf secret key"
    ))
    db = SQLAlchemy(app)

class Variants(db.Model):
    __tablename__ = 'VARIANTS'
    id = db.Column(db.Integer, primary_key=True)
    chromosome = db.Column(db.String(20))
    pos = db.Column(db.String(20))
    ref = db.Column(db.String(20))
    alt = db.Column(db.String(20))
    var_type = db.Column(db.String(20))
    genome_version = db.Column(db.String(20))
    gene = db.Column(db.String(20))
    isoform = db.Column(db.String(20))
    hgvsg = db.Column(db.String(20))
    hgvsp = db.Column(db.String(20))
    hgvsc = db.Column(db.String(20))
    count = db.Column(db.Integer)

class SampleVariants(db.Model):
    __tablename__ = 'SAMPLE_VARIANTS'
    id = db.Column(db.Integer, primary_key=True)
    sample_id = db.Column(db.Integer)
    var_id = db.Column(db.Integer)
    ann_id = db.Column(db.Integer)
    lab_confirmation = db.Column(db.String(20))    
    confirmation_technique = db.Column(db.String(20))
    classification = db.Column(db.String(20))
    ann_key = db.Column(db.String(200))
    ann_json = db.Column(db.String(20000))

class VarAnnotation(db.Model):
    __tablename__= 'VAR_ANNOTATION'
    id = db.Column(db.Integer, primary_key=True)
    var_id   = db.Column(db.Integer)
    ann_key  = db.Column(db.String(20))
    ann_json = db.Column(db.String(20))

class SampleTable(db.Model):
    __tablename__ = 'SAMPLES'
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.String(20))
    lab_id  = db.Column(db.String(120))
    ext1_id = db.Column(db.String(80))
    ext2_id = db.Column(db.String(80))
    run_id  = db.Column(db.String(80))
    petition_id  = db.Column(db.String(80))
    extraction_date =  db.Column(db.String(80))
    analysis_date   =  db.Column(db.String(80))
    tumour_purity   =  db.Column(db.String(80))
    sex  = db.Column(db.String(80))
    diagnosis  = db.Column(db.String(80))
    physician_name  = db.Column(db.String(80))
    medical_center  = db.Column(db.String(80))
    medical_address  = db.Column(db.String(80))
    sample_type  = db.Column(db.String(80))
    panel    =  db.Column(db.String(80))
    panel_version = db.Column(db.String(80))
    subpanel =  db.Column(db.String(80))
    roi_bed  =  db.Column(db.String(80))
    software =  db.Column(db.String(80))
    software_version =  db.Column(db.String(80))
    bam = db.Column(db.String(80))
    merged_vcf = db.Column(db.String(80))
    report_pdf = db.Column(db.String(80))

    def __repr__(self):
        return '<Sample %r>' % self.lab_id

class TherapeuticTable(db.Model):
    __tablename__ = 'THERAPEUTIC_VARIANTS'
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.String(20))
    lab_id  = db.Column(db.String(120))
    ext1_id = db.Column(db.String(80))
    ext2_id = db.Column(db.String(80))
    run_id  = db.Column(db.String(80))
    petition_id  = db.Column(db.String(80))
    gene  = db.Column(db.String(120))
    enst_id  = db.Column(db.String(120))
    hgvsp = db.Column(db.String(120))
    hgvsg =  db.Column(db.String(120))
    hgvsc =  db.Column(db.String(120))
    exon  = db.Column(db.String(120))
    variant_type = db.Column(db.String(120))
    consequence =  db.Column(db.String(120))
    depth = db.Column(db.String(120))
    allele_frequency = db.Column(db.String(120))
    read_support = db.Column(db.String(120))
    max_af = db.Column(db.String(120))
    max_af_pop = db.Column(db.String(120))
    therapies = db.Column(db.String(240))
    clinical_trials = db.Column(db.String(240))
    tumor_type = db.Column(db.String(240))
    var_json = db.Column(db.String(5000))

    def __repr__(self):
        return '<TherapeuticVariants %r>' % self.gene

class OtherVariantsTable(db.Model):
    __tablename__ = 'OTHER_VARIANTS'
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.String(20))
    lab_id  = db.Column(db.String(120))
    ext1_id = db.Column(db.String(80))
    ext2_id = db.Column(db.String(80))
    run_id  = db.Column(db.String(80))
    petition_id  = db.Column(db.String(80))
    gene  = db.Column(db.String(120))
    enst_id  = db.Column(db.String(120))
    hgvsp = db.Column(db.String(120))
    hgvsg =  db.Column(db.String(120))
    hgvsc =  db.Column(db.String(120))
    exon  = db.Column(db.String(120))
    variant_type = db.Column(db.String(120))
    consequence =  db.Column(db.String(120))
    depth = db.Column(db.String(120))
    allele_frequency = db.Column(db.String(120))
    read_support = db.Column(db.String(120))
    max_af = db.Column(db.String(120))
    max_af_pop = db.Column(db.String(120))
    therapies = db.Column(db.String(240))
    clinical_trials = db.Column(db.String(240))
    tumor_type = db.Column(db.String(240))
    var_json = db.Column(db.String(5000))

class RareVariantsTable(db.Model):
    __tablename__ = 'RARE_VARIANTS'
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.String(20))
    lab_id  = db.Column(db.String(120))
    ext1_id = db.Column(db.String(80))
    ext2_id = db.Column(db.String(80))
    run_id  = db.Column(db.String(80))
    petition_id  = db.Column(db.String(80))
    gene  = db.Column(db.String(120))
    enst_id  = db.Column(db.String(120))
    hgvsp = db.Column(db.String(120))
    hgvsg =  db.Column(db.String(120))
    hgvsc =  db.Column(db.String(120))
    exon  = db.Column(db.String(120))
    variant_type = db.Column(db.String(120))
    consequence =  db.Column(db.String(120))
    depth = db.Column(db.String(120))
    allele_frequency = db.Column(db.String(120))
    read_support = db.Column(db.String(120))
    max_af = db.Column(db.String(120))
    max_af_pop = db.Column(db.String(120))
    therapies = db.Column(db.String(240))
    clinical_trials = db.Column(db.String(240))
    tumor_type = db.Column(db.String(240))
    var_json = db.Column(db.String(5000))

class SummaryQcTable(db.Model):
    __tablename__ = 'SUMMARY_QC'

    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.String(20))
    lab_id  = db.Column(db.String(120))
    ext1_id = db.Column(db.String(80))
    ext2_id = db.Column(db.String(80))
    run_id  = db.Column(db.String(80))
    petition_id  = db.Column(db.String(80))
    summary_json = db.Column(db.String(12000))
    def __repr__(self):
        return '<SummaryQc %r>' % self.user_id

class Petition(db.Model):
    __tablename__ = 'PETITIONS' 
    Id             = db.Column(db.Integer(), primary_key=True, autoincrement=True)
    Petition_id    = db.Column(db.String(20))
    #User_id        = db.Column(db.String(20))
    Date           = db.Column(db.String(20))
    AP_code        = db.Column(db.String(20))
    HC_code        = db.Column(db.String(20))
    Tumour_pct     = db.Column(db.String(20))
    Volume         = db.Column(db.String(20))
    Conc_nanodrop  = db.Column(db.String(20))
    Ratio_nanodrop = db.Column(db.String(20))
    Medical_doctor = db.Column(db.String(50))
    Billing_unit   = db.Column(db.String(50))

    def __init__(self, Petition_id,  Date, AP_code, HC_code, Tumour_pct, Volume, Conc_nanodrop, 
        Ratio_nanodrop, Medical_doctor, Billing_unit):
        self.Petition_id = Petition_id
        #self.User_id  = User_id
        self.Date   = Date
        self.AP_code= AP_code
        self.HC_code  = HC_code
        self.Tumour_pct= Tumour_pct
        self.Volume = Volume
        self.Conc_nanodrop = Conc_nanodrop
        self.Ratio_nanodrop = Ratio_nanodrop
        self.Medical_doctor = Medical_doctor
        self.Billing_unit = Billing_unit
        
class Biomarker(db.Model):
  __tablename__ = 'BIOMARKERS'
  id = db.Column(Integer, primary_key=True)
  gene = db.Column(String(50))
  variant = db.Column(String(50))
  exon = db.Column(String(50))
  chr  = db.Column(String(50))
  pos  = db.Column(String(50))
  end  = db.Column(String(50))
  panel= db.Column(String(50))
  def __init__(self, ID, gene, variant, exon, chr, pos, end, panel):
      self.id      = ID
      self.gene    = gene
      self.variant = variant
      self.exon    = exon
      self.chr   = chr
      self.pos   = pos
      self.end   = end
      self.panel = panel

class Roi(db.Model):
  __tablename__ = 'ROI_TABLE'
  id = db.Column(Integer, primary_key=True)
  chromosome = db.Column(String(50))
  start = db.Column(String(50))
  end = db.Column(String(50))
  ensg_id = db.Column(String(50))
  enst_id = db.Column(String(50))
  gene_name = db.Column(String(50))
  genome_version = db.Column(String(50))
  panel_name = db.Column(String(50))
  panel_version = db.Column(String(50))

  def __init__(self, ID, CHROMOSOME, START, END, ENSG_ID, ENST_ID, GENE_NAME, GENOME_VERSION, PANEL_NAME, PANEL_VERSION):
      self.id  = ID
      self.chromosome= CHROMOSOME
      self.start     = START
      self.end       = END
      self.ensg_id   = ENSG_ID
      self.enst_id   = ENST_ID
      self.gene_name = GENE_NAME
      self.genome_version = GENOME_VERSION
      self.panel_name = PANEL_NAME
      self.panel_version = PANEL_VERSION

class Disclaimer(db.Model):
  __tablename__ = 'DISCLAIMERS'
  id = db.Column(Integer, primary_key=True)
  panel = db.Column(String(3000))
  genes = db.Column(String(3000))
  methodology = db.Column(String(3000))
  analysis = db.Column(String(3000))
  lab_confirmation = db.Column(String(3000))
  technical_limitations = db.Column(String(3000))
  legal_provisions = db.Column(String(3000))
  language = db.Column(String(3000))

  def __init__(self, id, panel, genes, methodology, analysis, lab_confirmation,technical_limitations, legal_provisions, language):
      self.id  = id
      self.panel     = panel
      self.genes     = genes
      self.methodology     = methodology
      self.analysis  = analysis
      self.lab_confirmation= lab_confirmation
      self.technical_limitations = technical_limitations
      self.legal_provisions= legal_provisions
      self.language = language

class Cna(db.Model):
  __tablename__ = 'CNA'
  id = db.Column(Integer, primary_key=True)
  chromosome = db.Column(String(50))
  start = db.Column(String(50))
  end = db.Column(String(50))
  gene = db.Column(String(50))
  genome_version = db.Column(String(50))
  panel_name = db.Column(String(50))
  panel_version = db.Column(String(50))
  dump_therapeutic = db.Column(String(50))
  min_cn = db.Column(String(50))
  def __init__(self, id, chromosome, start, end, gene, genome_version, panel_name, panel_version, dump_therapeutic, min_cn):
      self.id  = id
      self.chromosome= chromosome
      self.start     = start
      self.end       = end
      self.gene      = gene
      self.genome_version = genome_version
      self.panel_name= panel_name
      self.panel_version= panel_version
      self.dump_therapeutic = dump_therapeutic
      self.min_cn = min_cn

def update_sample_db():

    for sample in p.sample_env:

        result = SampleTable.query.filter_by(user_id=p.analysis_env['USER_ID'])\
            .filter_by(lab_id=sample).filter_by(run_id=p.analysis_env['OUTPUT_NAME']).first()

        if not result:

            ext1_id = '.'
            petition_id = '.'
            if p.lab_data[sample]['AP_CODE']:
                ext1_id = p.lab_data[sample]['AP_CODE']
                petition_id = p.sample_data[ext1_id]['PETITION_ID'] 
            ext2_id = '.'
            if p.lab_data[sample]['HC_CODE']:
                ext2_id = p.lab_data[sample]['HC_CODE']

            Sample = SampleTable(user_id=p.analysis_env['USER_ID'], lab_id=sample, ext1_id=ext1_id, 
                ext2_id=ext2_id, run_id=p.analysis_env['OUTPUT_NAME'],petition_id=petition_id, sex='.', 
                diagnosis='.', physician_name='.', medical_center='.', medical_address='.', sample_type='.',
                extraction_date=p.lab_data[sample]['PETITION_DATE'], analysis_date=p.analysis_env['ANALYSIS_DATE'], 
                tumour_purity=p.lab_data[sample]['PURITY'], panel=p.analysis_env['PANEL_NAME'], panel_version=p.analysis_env['PANEL_VERSION'] , 
                subpanel="ALL", roi_bed=p.analysis_env['PANEL_NAME'], software="varMut", software_version="0.9.0", bam=p.sample_env[sample]['READY_BAM'],
                merged_vcf=p.sample_env[sample]['READY_MERGED_VCF'], report_pdf=p.sample_env[sample]['REPORT_PDF']  )
            
            db.session.add(Sample)
            db.session.commit()

def update_therapeutic_variants_db():
    pass

def update_other_variants_db():
    pass

def update_rare_variants_db():
    pass

def update_summary_db():

    for sample in p.sample_env:
        summary_dict = defaultdict(dict)
        summary_dict['MEAN_COVERAGE']    = p.sample_env[sample]['MEAN_COVERAGE']
        summary_dict['TOTAL_READS']      = p.sample_env[sample]['TOTAL_READS']
        summary_dict['ON_TARGET_READS']  = p.sample_env[sample]['ON_TARGET_READS']
        summary_dict['ROI_PERCENTAGE']   = p.sample_env[sample]['ROI_PERCENTAGE']
        summary_dict['ROI_PERCENTAGE']   = p.sample_env[sample]['PCR_DUPLICATES_PERCENTAGE']
        summary_dict['MEAN_INSERT_SIZE'] = p.sample_env[sample]['MEAN_INSERT_SIZE']
        summary_dict['SD_INSERT_SIZE']   = p.sample_env[sample]['SD_INSERT_SIZE']
      
        for field in p.sample_env[sample]['CALL_RATE']:

            # Setting call rate values 
            summary_dict['CALL_RATE'][field] =  p.sample_env[sample]['CALL_RATE'][field] 

            # Setting lost exon values 
            if field in p.sample_env[sample]['LOST_EXONS']:
                summary_dict['LOST_EXONS'][field] =  p.sample_env[sample]['LOST_EXONS'][field]
            else:
                summary_dict['LOST_EXONS'][field] =  '.'

        ext1_id = '.'
        petition_id = '.'
        if p.lab_data[sample]['AP_CODE']:
            ext1_id = p.lab_data[sample]['AP_CODE']
            petition_id = p.sample_data[ext1_id]['PETITION_ID'] 
        ext2_id = '.'
        if p.lab_data[sample]['HC_CODE']:
            ext2_id = p.lab_data[sample]['HC_CODE']

        summary_json_str = json.dumps(summary_dict)
  
        result = SummaryQcTable.query.filter_by(user_id=p.analysis_env['USER_ID'])\
            .filter_by(lab_id=sample).filter_by(run_id=p.analysis_env['OUTPUT_NAME']).first()

        if not result:

            Summary = SummaryQcTable(user_id=p.analysis_env['USER_ID'], lab_id=sample, ext1_id=ext1_id, 
                ext2_id=ext2_id, run_id=p.analysis_env['OUTPUT_DIR'],petition_id=petition_id, summary_json=summary_json_str )
            db.session.add(Summary)
            db.session.commit()

def load_petitions():

    all_petitions = Petition.query.all()
    return all_petitions

def load_cna(pname):
    pname = pname.replace(".bed", "")
    pname = pname.replace(".v1", "")
    cna_env = defaultdict(dict)
    cna_info = Cna.query.filter_by(panel_name=pname).all()
    if cna_info:
        for g in cna_info:
            #print(g.gene)
            cna_env[g.gene] = defaultdict(dict)
            cna_env[g.gene]['gene'] = g.gene
            cna_env[g.gene]['chromosome'] = g.chromosome
            cna_env[g.gene]['start'] = g.start
            cna_env[g.gene]['end'] = g.end
            cna_env[g.gene]['genome_version'] = g.genome_version
            cna_env[g.gene]['min_cn'] = g.min_cn

    return cna_env


def load_panel_biomarkers(pname):

    pname = pname.replace(".bed", "")
    pname = pname.replace(".v1", "")
    biomarkers_env = defaultdict(dict)
    biomarkers_info = Biomarker.query.filter_by(panel=pname).all()
    if biomarkers_info:
        for biomarker in biomarkers_info:
            biomarkers_env[biomarker.id]['GENE'] = biomarker.gene
            biomarkers_env[biomarker.id]['VARIANT'] = biomarker.variant
            biomarkers_env[biomarker.id]['EXON'] = biomarker.exon
            biomarkers_env[biomarker.id]['CHR'] = biomarker.chr
            biomarkers_env[biomarker.id]['POS'] = biomarker.pos
            biomarkers_env[biomarker.id]['END'] = biomarker.end
    return biomarkers_env

def load_panel_disclaimers(pname, lang):
    pname= pname.replace(".bed", "")
    pname= pname.replace(".v1", "")

    disclaimer_env = defaultdict(dict)
    disclaimer_info = Disclaimer.query.filter_by(panel=pname).filter_by(language=lang).all()
    if disclaimer_info:
        for entry in disclaimer_info:
            #print(entry.language)
            disclaimer_env['GENES']            = entry.genes
            disclaimer_env['METHODOLOGY']      = entry.methodology
            disclaimer_env['ANALYSIS']         = entry.analysis
            disclaimer_env['LAB_CONFIRMATION'] = entry.lab_confirmation
            disclaimer_env['TECHNICAL_LIMITATIONS'] = entry.technical_limitations
            disclaimer_env['LEGAL_PROVISIONS'] = entry.legal_provisions
            disclaimer_env['LANGUAGE'] = entry.language

    return disclaimer_env

def load_panel_transcripts(pname):

    pname= pname.replace(".bed", "")
    pname= pname.replace(".v1", "")

    roi_env = defaultdict(dict)
    roi_info = Roi.query.filter_by(panel_name=pname).all()
    if roi_info:
        for roi in roi_info:
            roi_env[roi.ensg_id]   = roi.enst_id
            roi_env[roi.gene_name] = roi.enst_id
    return roi_env