#!/usr/bin/env python3

import os
import sys
import shutil
import re
import logging
import gzip
import csv
import hashlib
from collections import defaultdict
from pathlib import Path
from scipy import stats
import numpy as np
import subprocess
from civicpy import civic, exports
import requests
import json
from sqlalchemy import create_engine, Column, Integer, String
from sqlalchemy_utils import database_exists, create_database
from sqlalchemy.ext.declarative import declarative_base
from flask import Flask

from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import Table, Column, Float, Integer, String, MetaData, ForeignKey
from googletrans import Translator

from modules import utils as u
from modules import params as p
from modules import trimming as t
from modules import sqlite as s

def do_report():

    clinical_vcf()

    if p.analysis_env['VARIANT_CLASS'] == "somatic":

      annotate_biomarkers()

      create_somatic_report()


def annotate_biomarkers():

  for sample in p.sample_env:
    p.sample_env[sample]['BIOMARKER'] = defaultdict(dict)

    idx = 0
    for biomarker in p.biomarker_env:

      chr = p.biomarker_env[biomarker]['CHR']
      if not 'chr' in chr:
        chr = "chr" + chr
      coordinate = chr + ":" \
        + p.biomarker_env[biomarker]['POS'] + "-" + p.biomarker_env[biomarker]['END'] 

      mean_depth = u.mean_depth_coordinate(p.sample_env[sample]['READY_BAM'], coordinate )

      p.sample_env[sample]['BIOMARKER'][idx] = defaultdict(dict)
      p.sample_env[sample]['BIOMARKER'][idx]['COORDINATE'] = coordinate
      p.sample_env[sample]['BIOMARKER'][idx]['GENE'] =  p.biomarker_env[biomarker]['GENE']
      p.sample_env[sample]['BIOMARKER'][idx]['VARIANT'] =  p.biomarker_env[biomarker]['VARIANT']
      p.sample_env[sample]['BIOMARKER'][idx]['EXON'] =  p.biomarker_env[biomarker]['EXON']
      p.sample_env[sample]['BIOMARKER'][idx]['MEAN_DEPTH'] = mean_depth
      idx+=1

def create_somatic_report():

    for sample in p.sample_env:

        p.sample_env[sample]['REPORT_FOLDER'] = p.sample_env[sample]['SAMPLE_FOLDER'] + "/" + "REPORT_FOLDER"

        report_folder_path = Path(p.sample_env[sample]['REPORT_FOLDER'])
        if not report_folder_path.is_dir():
            os.mkdir(report_folder_path)

        p.sample_env[sample]['REPORT_PDF'] = p.sample_env[sample]['REPORT_FOLDER'] + "/" + sample + ".report"

        p.sample_env[sample]['CLINICAL_REPORT_CSV'] = \
            p.sample_env[sample]['REPORT_FOLDER'] + "/" + sample + ".clinical.csv"
        o =  open (p.sample_env[sample]['CLINICAL_REPORT_CSV'], "w")

        sample_db = p.sample_env[sample]['REPORT_FOLDER'] + "/" + sample + ".db"  
        p.sample_env[sample]['REPORT_DB'] = sample_db

        app = Flask(__name__)
        app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///' + sample_db
        app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

        db = SQLAlchemy(app)
        class Sample(db.Model):
            __tablename__ = 'SAMPLE_INFORMATION'
            id = db.Column(db.Integer, primary_key=True)
            lab_id = db.Column(db.String(120))
            ext1_id = db.Column(db.String(80))
            ext2_id = db.Column(db.String(80))
            sex    =  db.Column(db.String(80))
            diagnoses =  db.Column(db.String(80))
            physician_name =  db.Column(db.String(80))
            medical_center =  db.Column(db.String(80))
            medical_address =  db.Column(db.String(80))
            sample_type =  db.Column(db.String(80))
            extraction_date =  db.Column(db.String(80))
            tumor_purity =  db.Column(db.String(80))
            panel =  db.Column(db.String(80))
            analysis_date =  db.Column(db.String(80))
            def __repr__(self):
                return '<Sample %r>' % self.lab_id

        class TherapeuticVariants(db.Model):
            __tablename__ = 'THERAPEUTIC_VARIANTS'
            id = db.Column(db.Integer, primary_key=True)
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

            def __repr__(self):
                return '<TherapeuticVariants %r>' % self.gene

        class OtherClinicalVariants(db.Model):
            __tablename__ = 'OTHER_VARIANTS'
            id = db.Column(db.Integer, primary_key=True)
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

            def __repr__(self):
                return '<OtherClinicalVariants %r>' % self.gene

        class RareVariants(db.Model):
            __tablename__ = 'RARE_VARIANTS'
            id = db.Column(db.Integer, primary_key=True)
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

            def __repr__(self):
                return '<RareVariants %r>' % self.gene

        class Biomarkers(db.Model):
            __tablename__ = 'BIOMARKERS'
            id = db.Column(db.Integer, primary_key=True)
            gene  = db.Column(db.String(120))
            variant = db.Column(db.String(120))
            exon =  db.Column(db.String(120))
            allele_fraction =  db.Column(db.String(120))
            sequencing_depth = db.Column(db.String(120))
            def __repr__(self):
                return '<Biomarkers %r>' % self.gene

        class SummaryQc(db.Model):
            __tablename__ = 'SUMMARY_QC'

            id = db.Column(db.Integer, primary_key=True)
            total_reads = db.Column(db.String(120))
            mean_coverage = db.Column(db.String(120))
            enrichment =  db.Column(db.String(120))
            call_rate = db.Column(db.String(120))
            lost_exons =  db.Column(db.String(120))
            pct_read_duplicates = db.Column(db.String(120))
            def __repr__(self):
                return '<SummaryQc %r>' % self.total_reads

        class Disclaimers(db.Model):
            __tablename__ = 'DISCLAIMERS'

            id = db.Column(db.Integer, primary_key=True)
            gene_list = db.Column(db.String(3000))
            lab_methodology = db.Column(db.String(3000))
            analysis =  db.Column(db.String(3000))
            lab_confirmation = db.Column(db.String(3000))
            technique_limitations =  db.Column(db.String(3000))
            legal_provisions = db.Column(db.String(3000))
            def __repr__(self):
                return '<Disclaimers %r>' % self.gene_list

        if os.path.isfile(sample_db) and os.path.isfile(p.sample_env[sample]['REPORT_PDF']):
          msg = " INFO: Skipping report creation for sample " + sample
          print (msg)
          logging.info(msg)
          continue

        if os.path.isfile(sample_db) and not os.path.isfile(p.sample_env[sample]['REPORT_PDF']):
           # Create PDF report 
          bashCommand = ('{} pr {} -r {} -f pdf -t generic --db-url jdbc:sqlite:{} --db-driver org.sqlite.JDBC -o {} --jdbc-dir {}') \
            .format(p.system_env['JASPERSTARTER'],p.aux_env['REPORT_JRXML'], 
            p.defaults['JASPERREPORT_FOLDER'], sample_db, p.sample_env[sample]['REPORT_PDF'], p.defaults['JDBC_FOLDER'])
          msg = " INFO: Generating pdf report for sample " + sample
          print (msg)
          print(bashCommand)
          logging.info(bashCommand)
          
          p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
          output = p1.stdout.decode('UTF-8')
          error  = p1.stderr.decode('UTF-8')
          if not error and not output:
            msg = " INFO: Report generation for " + sample + " ended successfully"
            print(msg)
            logging.info(msg)
          else:
            msg = " ERROR: Failed report generation for " + sample
            print(msg)
            logging.error(msg)
          continue     
        db.create_all()       

        ext1_id = '.'
        petition_id = '.'
        if p.lab_data[sample]['AP_CODE']:
          ext1_id = p.lab_data[sample]['AP_CODE']
          print(ext1_id)
          petition_id = p.sample_data[ext1_id]['PETITION_ID'] 
        ext2_id = '.'
        if p.lab_data[sample]['HC_CODE']:
          ext2_id = p.lab_data[sample]['HC_CODE']

        sample_r = Sample(lab_id=sample, ext1_id=p.lab_data[sample]['AP_CODE'], ext2_id=p.lab_data[sample]['HC_CODE'], sex='.', diagnoses='.', 
          physician_name='.', medical_center='.', medical_address='.', sample_type='.', extraction_date=p.lab_data[sample]['PETITION_DATE'],
          tumor_purity=p.lab_data[sample]['PURITY'], panel=p.analysis_env['PANEL_NAME'], analysis_date= p.analysis_env['ANALYSIS_DATE'])
        db.session.add(sample_r)
        db.session.commit()

        sample_found = s.SampleTable.query.filter_by(user_id=p.analysis_env['USER_ID'])\
            .filter_by(lab_id=sample).filter_by(run_id=p.analysis_env['OUTPUT_NAME']).first()

        if not sample_found:

          cnv_dict = defaultdict(dict)
          with open(p.sample_env[sample]['CNV_JSON']) as jf:
            cnv_dict = json.load(jf)
          cnv_json_str = json.dumps(cnv_dict)

          Sample_db = s.SampleTable(user_id=p.analysis_env['USER_ID'], lab_id=sample, ext1_id=ext1_id, 
            ext2_id=ext2_id, run_id=p.analysis_env['OUTPUT_NAME'],petition_id=petition_id, sex='.', diagnosis='.', 
            physician_name='.', medical_center='.', medical_address='.', sample_type='.',
            extraction_date=p.lab_data[sample]['PETITION_DATE'], analysis_date=p.analysis_env['ANALYSIS_DATE'], 
            tumour_purity=p.lab_data[sample]['PURITY'], panel=p.analysis_env['PANEL_NAME'], 
            subpanel="ALL", roi_bed=p.analysis_env['PANEL_NAME'], software="varMut", software_version="0.9.0",
            bam=p.sample_env[sample]['READY_BAM'],merged_vcf=p.sample_env[sample]['READY_MERGED_VCF'], 
            report_pdf=p.sample_env[sample]['REPORT_PDF']+".pdf", cnv_json=cnv_json_str,
            report_db=sample_db, sample_db_dir=p.sample_env[sample]['REPORT_FOLDER'])
          print("COMMITTING")

          s.db.session.add(Sample_db)
          s.db.session.commit()
          result = s.SampleTable.query.filter_by(user_id=p.analysis_env['USER_ID'])\
            .filter_by(lab_id=sample).filter_by(run_id=p.analysis_env['OUTPUT_NAME']).first()
          sample_id = result.id
        else:
          sample_id = sample_found.id

        disclaimer_r = Disclaimers(gene_list=p.disclaimers_env['GENES'], lab_methodology=p.disclaimers_env['METHODOLOGY'], 
          analysis=p.disclaimers_env['ANALYSIS'], lab_confirmation=p.disclaimers_env['LAB_CONFIRMATION'], 
          technique_limitations=p.disclaimers_env['TECHNICAL_LIMITATIONS'],legal_provisions=p.disclaimers_env['LEGAL_PROVISIONS'])
        db.session.add(disclaimer_r)
        db.session.commit()

        # filling biomarker report  
        idx = 0
        for biomarker in p.biomarker_env:
          biomarker_r = Biomarkers(gene=p.sample_env[sample]['BIOMARKER'][idx]['GENE'],
          variant=p.sample_env[sample]['BIOMARKER'][idx]['VARIANT'],
          exon=p.sample_env[sample]['BIOMARKER'][idx]['EXON'],
          allele_fraction='.', sequencing_depth=p.sample_env[sample]['BIOMARKER'][idx]['MEAN_DEPTH'])
          db.session.add(biomarker_r)
          idx+=1
        db.session.commit()

        # filling biomarker database
        idx = 0
        for biomarker in p.biomarker_env:
          biomarker_db = s.BiomarkerTable(user_id=p.analysis_env['USER_ID'], lab_id=sample, ext1_id=ext1_id, 
          ext2_id=ext2_id, run_id=p.analysis_env['OUTPUT_NAME'], gene=p.sample_env[sample]['BIOMARKER'][idx]['GENE'],
          variant=p.sample_env[sample]['BIOMARKER'][idx]['VARIANT'], exon=p.sample_env[sample]['BIOMARKER'][idx]['EXON'],
          chr='.', pos='.', end='.', vaf='.', depth=p.sample_env[sample]['BIOMARKER'][idx]['MEAN_DEPTH'])
          s.db.session.add(biomarker_db)
          idx+=1
        s.db.session.commit()

        summary_r = SummaryQc(total_reads=p.sample_env[sample]['TOTAL_READS'], 
        mean_coverage=str(p.sample_env[sample]['MEAN_COVERAGE']), 
        enrichment=str(p.sample_env[sample]['ROI_PERCENTAGE']), 
        call_rate=str(p.sample_env[sample]['CALL_RATE']['30X']), 
        lost_exons=str(p.sample_env[sample]['LOST_EXONS']['30X']), 
        pct_read_duplicates=str(p.sample_env[sample]['PCR_DUPLICATES_PERCENTAGE']))
        db.session.add(summary_r)
        db.session.commit()

        # report_header = [] 
        # report_header.append("Lab_ID")
        # report_header.append("Ext_ID")
        # report_header.append("Gene")
        # report_header.append("HGVSp")
        # report_header.append("HGVSg")
        # report_header.append("HGVSc")
        # report_header.append("Exon")
        # report_header.append("ENST_ID")
        # report_header.append("VAF")
        # report_header.append("Identifiers")
        # report_header.append("MAX_AF")
        # report_header.append("MAX_AF_POP")
        # report_header.append("Drugs")
        # report_header.append("Clinical_trials")
        # o.write('\t'.join(report_header)+"\n")

        var_dict = defaultdict(dict)
        with open(p.sample_env[sample]['READY_MERGED_JSON']) as f:
            var_dict = json.load(f)
        
        for variant in var_dict['variants']:
          results_list = [] 
          lab_id = sample
          results_list.append(lab_id)
          ext_id = "999999"
          results_list.append(ext_id)

          go_therapeutic = False
          go_other = False
          go_rare  = False
          go_common= False
          chromosome = var_dict['variants'][variant]['CHROM']
          pos     = var_dict['variants'][variant]['POS']
          ref     = var_dict['variants'][variant]['REF'] 
          alt     = var_dict['variants'][variant]['ALT'] 
          p_code  = '.'
          g_code  = '.'
          c_code  = '.'
          gene    = '.'
          depth   = '.'
          exon    = '.'
          enst_id = '.'
          clin_sig= '.'
          rs_id   = '.'
          max_af  = '.'
          max_af_pop   = '.'
          drugs        = '.'
          clin_trials  = '.'
          variant_type = '.'
          isoform      = '.'
          classification = '.'
          svlen = '.'              

          ann_dict = defaultdict(dict)
          for field in var_dict['variants'][variant]:
            ann_dict[field] = var_dict['variants'][variant][field] 
          m = hashlib.md5()
          ann_json_str = json.dumps(ann_dict)
          m.update(ann_json_str.encode('UTF-8'))
          ann_key = m.hexdigest()
          
          var_json_str = json.dumps(var_dict['variants'][variant])

          if 'SVTYPE' in var_dict['variants'][variant]['INFO']:
            if 'CNA_GENES' in var_dict['variants'][variant]['INFO']:
              gene  = var_dict['variants'][variant]['INFO']['CNA_GENES'] 
              if 'CN' in var_dict['variants'][variant]:
                VAF = int(var_dict['variants'][variant]['CN'])
              else:
                GT = var_dict['variants'][variant]['GT']
                VAF = int(GT[0])
              svtype = var_dict['variants'][variant]['INFO']['SVTYPE']

              # Setting CNA type
              if svtype == "DUP":
                variant_type = "Amplification"
              elif svtype == "DEL":
                variant_type = "Loss"

              # Filling evidence info from CIViC
              drugs_list          = []
              clin_trials_list    = []
              ev_direction_list   = []
              diseases_list       = []
              ev_significance_list= []

              if 'CIVIC' in var_dict['variants'][variant]['INFO']:
                for ev_id in var_dict['variants'][variant]['INFO']['CIVIC']:
                  if 'EV_DRUGS' in var_dict['variants'][variant]['INFO']['CIVIC'][ev_id]:
                    drugs_result = var_dict['variants'][variant]['INFO']['CIVIC'][ev_id]['EV_DRUGS'].split('&')
                    for drug in drugs_result:
                      if drug != '.':
                        drugs_list.append(drug)
                  if 'EV_CLINICAL_TRIALS' in var_dict['variants'][variant]['INFO']['CIVIC'][ev_id] :
                    clin_trials_result = var_dict['variants'][variant]['INFO']['CIVIC'][ev_id]['EV_CLINICAL_TRIALS'].split('&')
                    for ctrial in clin_trials_result:
                      if ctrial != '.':
                        clin_trials_list.append(ctrial)
                  if 'EV_DIRECTION' in var_dict['variants'][variant]['INFO']['CIVIC'][ev_id] :
                    direction_result = var_dict['variants'][variant]['INFO']['CIVIC'][ev_id]['EV_DIRECTION'].split('&')
                    for direction in direction_result:
                      if direction != '.':
                        ev_direction_list.append(direction)
                  if 'EV_SIGNIFICANCE' in var_dict['variants'][variant]['INFO']['CIVIC'][ev_id] :
                    significance_result = var_dict['variants'][variant]['INFO']['CIVIC'][ev_id]['EV_SIGNIFICANCE'].split('&')
                    for significance in significance_result:
                      if significance != '.':
                        ev_significance_list.append(significance)                      
                  if 'EV_DISEASE' in var_dict['variants'][variant]['INFO']['CIVIC'][ev_id] :
                    disease_result = var_dict['variants'][variant]['INFO']['CIVIC'][ev_id]['EV_DISEASE'].split('&')
                    for disease in disease_result:
                      if disease != '.':
                        disease = disease.replace("_", " ")                   
                        diseases_list.append(disease)
              if drugs_list :
                drugs_str = ', '.join(list(set(drugs_list)))
              else:
                drugs_str = '.'
              if clin_trials_list:
                clintrials_str = ','.join(list(set(clin_trials_list)))
              else:
                clintrials_str = '.'
              if ev_direction_list:
                direction_str = ','.join(list(set(ev_direction_list)))
              else:
                direction_str = '.'
              if ev_significance_list:
                significance_str = ','.join(list(set(ev_significance_list)))
              else:
                significance_str = '.'
              if diseases_list:
                diseases_str = ','.join(list(set(diseases_list)))
                diseases_str = diseases_str.replace("_", " ")
              else:
                diseases_str = '.'

              # Now classifying each variant into their corresponding clinical class
              if drugs_str != '.' and 'Supports' in direction_str and 'Sensitivity/Response' in significance_str:
                if gene in p.cna_env:
                  if VAF >= int(p.cna_env[gene]['min_cn']):
                    go_therapeutic = True
              elif drugs_str != '.' and not 'Sensitivity/Response' in significance_str:             
                if gene in p.cna_env:
                  if VAF >= int(p.cna_env[gene]['min_cn']): 
                    go_other = True               

              # Filling therapeutic table  
              if go_therapeutic == True:
                classification = "Therapeutic"
                therapeutic_r = TherapeuticVariants(gene=gene, enst_id=enst_id, hgvsp=p_code, hgvsg=g_code,
                hgvsc=c_code, exon=exon, variant_type=variant_type, consequence=consequence, depth=depth,
                allele_frequency=str(VAF),max_af=max_af,max_af_pop=max_af_pop,
                therapies=drugs_str, clinical_trials=clintrials_str, tumor_type=diseases_str)
                db.session.add(therapeutic_r)
                db.session.commit()

                result = s.TherapeuticTable.query.filter_by(user_id=p.analysis_env['USER_ID'])\
                .filter_by(lab_id=sample).filter_by(run_id=p.analysis_env['OUTPUT_NAME']).filter_by(gene=gene)\
                .filter_by(enst_id=enst_id).filter_by(hgvsp=p_code).filter_by(hgvsg=g_code).filter_by(hgvsc=c_code).first()

                if not result:
                  therapeutic_db = s.TherapeuticTable(user_id=p.analysis_env['USER_ID'], lab_id=sample, ext1_id=ext1_id, 
                  ext2_id=ext2_id, run_id=p.analysis_env['OUTPUT_NAME'],petition_id=petition_id,gene=gene, enst_id=enst_id, 
                  hgvsp=p_code, hgvsg=g_code, hgvsc=c_code, exon=exon, variant_type=variant_type, consequence=consequence, 
                  depth=depth, allele_frequency=str(VAF),read_support=".",max_af=max_af,max_af_pop=max_af_pop,therapies=drugs_str,
                  clinical_trials=clintrials_str,tumor_type=diseases_str, var_json=var_json_str, validated_bioinfo="no", 
                  validated_assessor="no")
                  s.db.session.add(therapeutic_db)
                  s.db.session.commit()

              # Filling Other clinical variants  
              if go_other == True:
                classification = "Other clinical"
                other_r = OtherClinicalVariants(gene=gene,  enst_id=enst_id, hgvsp=p_code, hgvsg=g_code,
                hgvsc=c_code, exon=exon, variant_type=variant_type, consequence=consequence, depth=depth,
                allele_frequency=str(VAF),max_af=max_af,max_af_pop=max_af_pop,therapies=drugs_str,
                clinical_trials=clintrials_str,tumor_type=diseases_str)
                db.session.add(other_r)
                db.session.commit()

                result = s.OtherVariantsTable.query.filter_by(user_id=p.analysis_env['USER_ID'])\
                .filter_by(lab_id=sample).filter_by(run_id=p.analysis_env['OUTPUT_NAME']).filter_by(gene=gene)\
                .filter_by(enst_id=enst_id).filter_by(hgvsp=p_code).filter_by(hgvsg=g_code).filter_by(hgvsc=c_code).first()
                if not result:
                  other_db = s.OtherVariantsTable(user_id=p.analysis_env['USER_ID'], lab_id=sample, ext1_id=ext1_id, 
                  ext2_id=ext2_id, run_id=p.analysis_env['OUTPUT_NAME'],petition_id=petition_id,gene=gene, enst_id=enst_id, 
                  hgvsp=p_code, hgvsg=g_code, hgvsc=c_code, exon=exon, variant_type=variant_type, consequence=consequence, 
                  depth=depth, allele_frequency=str(VAF),read_support=".",max_af=max_af,max_af_pop=max_af_pop,therapies=drugs_str,
                  clinical_trials=clintrials_str,tumor_type=diseases_str, var_json=var_json_str,validated_bioinfo="no", 
                  validated_assessor="no", classification="Other")
                  s.db.session.add(other_db)
                  s.db.session.commit()
            else:
              variant_type = "SV"
              if 'SVLEN' in var_dict['variants'][variant]['INFO']:
                svlen = abs(int(var_dict['variants'][variant]['INFO']['SVLEN']))
                if svlen < p.analysis_env['MIN_FUSION_SIZE']:
                  continue 

              if 'PR' in var_dict['variants'][variant]:
                pe_info = var_dict['variants'][variant]['PR'].split(",")
                pe_ref = int(pe_info[0])
                pe_alt = int(pe_info[1])
              else:
                pe_ref = 0
                pe_alt = 0
              if 'SR' in var_dict['variants'][variant]:
                sr_info = var_dict['variants'][variant]['SR'].split(",")
                sr_ref = int(sr_info[0])
                sr_alt = int(sr_info[1])
              else:
                sr_ref = 0
                sr_alt = 0
              fusion_out = False
              if (pe_alt > 0 and sr_alt > 0):
                fusion_out = True
              depth = '.'
              read_support = 0
              if sr_alt+sr_ref > 0:
                read_support += sr_alt
                depth = str(sr_alt+sr_ref)
                VAF = str(round((sr_alt/(sr_alt+sr_ref)),3))
              elif pe_ref+pe_alt >0:
                read_support += pe_alt
                depth = str(pe_alt+pe_ref)
                VAF = str(round((pe_alt/(pe_alt+pe_ref)),3))
              else:
                depth = "0"
                VAF = "0"

              # Filter SVs with leass than 3 reads
              if read_support < 3:
                continue 

              gene = '.'
              if var_dict['variants'][variant]['INFO']['FUSION']['PARTNERS'] != '.' and \
                var_dict['variants'][variant]['INFO']['FUSION']['PARTNERS'] != "" :
                gene = var_dict['variants'][variant]['INFO']['FUSION']['PARTNERS'] + " FUSION"
              if gene != '.' and fusion_out == True:
                if gene == '.':
                  gene = var_dict['variants'][variant]['INFO']['EDGE_GENES']
                results_list.append(gene)
                p_code = '.'
                results_list.append(p_code)
                g_code = '.'
                results_list.append(g_code)
                c_code = '.'
                results_list.append(c_code)
                exon = '.'
                results_list.append(exon)
                enst_id = '.'
                results_list.append(enst_id)
                results_list.append(VAF)
                clin_sig = '.'
                results_list.append(clin_sig)
                rs_id = '.'
                results_list.append(rs_id)
                max_af = '.'
                results_list.append(max_af)
                max_af_pop = '.'
                results_list.append(max_af_pop)
                drugs = '.'
                results_list.append(drugs)
                clin_trials = '.'
                results_list.append(clin_trials)
                o.write('\t'.join(results_list)+ "\n")
                print("SV " + p_code)

                classification = "Therapeutic"
                therapeutic_r = TherapeuticVariants(gene=gene, enst_id=enst_id,hgvsp=p_code, hgvsg=g_code,
                hgvsc=c_code, exon=exon, variant_type=variant_type, consequence='.', depth=depth, allele_frequency=VAF,
                max_af=max_af, max_af_pop=max_af_pop, therapies=drugs_str, clinical_trials=clintrials_str, tumor_type=diseases_str)
                db.session.add(therapeutic_r)
                db.session.commit() 

                result = s.TherapeuticTable.query.filter_by(user_id=p.analysis_env['USER_ID'])\
                .filter_by(lab_id=sample).filter_by(run_id=p.analysis_env['OUTPUT_NAME']).filter_by(gene=gene)\
                .filter_by(enst_id=enst_id).filter_by(hgvsp=p_code).filter_by(hgvsg=g_code).filter_by(hgvsc=c_code).first()
                if not result:
                  therapeutic_db = s.TherapeuticTable(user_id=p.analysis_env['USER_ID'], lab_id=sample, ext1_id=ext1_id, 
                  ext2_id=ext2_id, run_id=p.analysis_env['OUTPUT_NAME'], petition_id=petition_id,gene=gene, enst_id=enst_id, 
                  hgvsp=p_code, hgvsg=g_code, hgvsc=c_code, exon=exon, variant_type=variant_type, consequence='.', 
                  depth=depth, allele_frequency=VAF, read_support=read_support, max_af=max_af, max_af_pop=max_af_pop, 
                  therapies=drugs_str, clinical_trials=clintrials_str, tumor_type=diseases_str, var_json=var_json_str,
                  validated_bioinfo="no", validated_assessor="no", classification="Therapeutic")
                  s.db.session.add(therapeutic_db)
                  s.db.session.commit()
          else:
            variant_type = '.'
            ref = var_dict['variants'][variant]['REF']
            alt = var_dict['variants'][variant]['ALT']

            if len(ref) > len(alt):
              variant_type = "Deletion"
            elif len(ref) < len(alt):
              variant_type = "Insertion"
            else:
              if len(ref) > 1:
                variant_type = "MNV"
              else:
                variant_type = "SNV"

            gene   = var_dict['variants'][variant]['INFO']['CSQ']['SYMBOL']
            results_list.append(gene)
            p_code = var_dict['variants'][variant]['INFO']['CSQ']['HGVSp']
            AD = var_dict['variants'][variant]['AD']
            read_support = "."
            if AD:
              ad_list = AD.split(",")
              read_support = ad_list[1] 

            # replace weird symbols (coming from VEP) for synonymous changes
            p_code = p_code.replace('%3D', '=')
            if p_code:
              tmp_p_code = p_code.split(":")
              if len(tmp_p_code)>1:
                p_code = tmp_p_code[1] 
            results_list.append(p_code)
            g_code = var_dict['variants'][variant]['INFO']['CSQ']['HGVSg'] 
            results_list.append(g_code)
            c_code = var_dict['variants'][variant]['INFO']['CSQ']['HGVSc']
            if c_code:
              tmp_c_code = c_code.split(":")
              if len(tmp_c_code)> 1:
                c_code = tmp_c_code[1]
            results_list.append(c_code)
            exon   = var_dict['variants'][variant]['INFO']['CSQ']['EXON']
            results_list.append(exon)
            enst_id= var_dict['variants'][variant]['INFO']['CSQ']['Feature']
            results_list.append(enst_id)
            isoform     = var_dict['variants'][variant]['INFO']['CSQ']['Feature']
            consequence = var_dict['variants'][variant]['INFO']['CSQ']['Consequence']
            c_list = [] 
            conseq_list = consequence.split("&")
            if len(conseq_list)>1:
              for v in conseq_list:
                v = v.replace("_", " ")
                v = v.capitalize()
                c_list.append(v)
              consequence = ','.join(c_list)
            else:
              consequence = consequence.replace("_", " ")
              consequence = consequence.capitalize()
            depth  = var_dict['variants'][variant]['INFO']['DP']
            VAF    = var_dict['variants'][variant]['AF']
            results_list.append(VAF)
            clin_sig = var_dict['variants'][variant]['INFO']['CSQ']['CLIN_SIG']
            rs_id  = var_dict['variants'][variant]['INFO']['CSQ']['Existing_variation']
            rs_id = rs_id.replace('&', ',')
            results_list.append(rs_id)

            # Gnomad frequencies present at the Broad server do not always match VEP data
            # We will check the presence from both sources, but preserving original GnomAD data when available
            if 'GNOMAD_AF_popmax' in var_dict['variants'][variant]['INFO']['CSQ']:
              max_af = var_dict['variants'][variant]['INFO']['CSQ']['GNOMAD_AF_popmax']
            else:
              max_af = var_dict['variants'][variant]['INFO']['CSQ']['MAX_AF']
            results_list.append(max_af)
            if 'GNOMAD_popmax' in var_dict['variants'][variant]['INFO']['CSQ']:
              max_af_pop = var_dict['variants'][variant]['INFO']['CSQ']['GNOMAD_popmax']
            else: 
              max_af_pop = var_dict['variants'][variant]['INFO']['CSQ']['MAX_AF_POPS']
            max_af_pop = max_af_pop.replace('&', ',')

            results_list.append(max_af_pop)
            drugs_list         = []
            clin_trials_list   = []
            ev_direction_list  = []
            diseases_list      = []
            ev_significance_list=[] 
            if 'CIVIC' in var_dict['variants'][variant]['INFO']:
              for ev_id in var_dict['variants'][variant]['INFO']['CIVIC']:
                if 'EV_DRUGS' in var_dict['variants'][variant]['INFO']['CIVIC'][ev_id]:
                  drugs_result = var_dict['variants'][variant]['INFO']['CIVIC'][ev_id]['EV_DRUGS'].split('&')
                  for drug in drugs_result:
                    if drug != '.':
                      drugs_list.append(drug)
                if 'EV_CLINICAL_TRIALS' in var_dict['variants'][variant]['INFO']['CIVIC'][ev_id] :
                  clin_trials_result = var_dict['variants'][variant]['INFO']['CIVIC'][ev_id]['EV_CLINICAL_TRIALS'].split('&')
                  for ctrial in clin_trials_result:
                    if ctrial != '.':
                      clin_trials_list.append(ctrial)
                if 'EV_DIRECTION' in var_dict['variants'][variant]['INFO']['CIVIC'][ev_id] :
                  direction_result = var_dict['variants'][variant]['INFO']['CIVIC'][ev_id]['EV_DIRECTION'].split('&')
                  for direction in direction_result:
                    print(direction)
                    if direction != '.':
                      ev_direction_list.append(direction)
                if 'EV_SIGNIFICANCE' in var_dict['variants'][variant]['INFO']['CIVIC'][ev_id] :
                  significance_result = var_dict['variants'][variant]['INFO']['CIVIC'][ev_id]['EV_SIGNIFICANCE'].split('&')
                  for significance in significance_result:
                    if significance != '.':
                      ev_significance_list.append(significance)                      
                if 'EV_DISEASE' in var_dict['variants'][variant]['INFO']['CIVIC'][ev_id] :
                  disease_result = var_dict['variants'][variant]['INFO']['CIVIC'][ev_id]['EV_DISEASE'].split('&')
                  for disease in disease_result:
                    if disease != '.':
                      disease = disease.replace("_", " ")                   
                      diseases_list.append(disease)
            if drugs_list :
              drugs_str = ', '.join(list(set(drugs_list)))
            else:
              drugs_str = '.'
            if clin_trials_list:
              clintrials_str = ','.join(list(set(clin_trials_list)))
            else:
              clintrials_str = '.'
            if ev_direction_list:
              direction_str = ','.join(list(set(ev_direction_list)))
            else:
              direction_str = '.'
            if ev_significance_list:
              significance_str = ','.join(list(set(ev_significance_list)))
            else:
              significance_str = '.'
            if diseases_list:
              diseases_str = ','.join(list(set(diseases_list)))
              diseases_str = diseases_str.replace("_", " ")
            else:
              diseases_str = '.'
            results_list.append(drugs_str)
            results_list.append(clintrials_str)

            go_therapeutic = False
            go_other       = False
            go_rare        = False
            go_common      = False
            if significance_str != ".":
              print(significance_str)
            clinsig_list = clin_sig.split("&")

            # Selecting actionable variants from CIViC 
            if drugs_str != '.' and 'Supports' in direction_str \
              and 'Sensitivity/Response' in significance_str \
              or 'Outcome' in significance_str:
              if max_af == '.':
                go_therapeutic = True
                classification = "Therapeutic"
              else:
                if float(max_af) < 0.01:
                  go_therapeutic = True
                  classification = "Therapeutic"

            # Selecting pathogenic variants from ClinVar 
            elif not 'Sensitivity/Response' in significance_str:
              for interpretation in clinsig_list:
                if interpretation.endswith('pathogenic'):
                  if max_af == '.':
                    go_other = True
                    classification = "Other clinical"
                  else:
                    if float(max_af) < 0.01:
                      go_other = True
                      classification = "Other clinical"

              # Selecting rare variants without clinical evidence
              else:
                if max_af == '.' :
                  go_rare = True
                  classification = "Rare"
                else:
                  if float(max_af) < 0.01:
                    go_rare = True
                    classification = "Rare"
                  if float(max_af) >= 0.01:
                    go_common = True
                    classification = "Common"

            # VAF might have multiple values if coming from an INDEL
            VAF_list = VAF.split(",")
            if len(VAF_list) > 1:
              VAF = VAF_list[0]   

            if float(VAF) < 0.02:
              continue

            # Filling therapeutic table  
            if go_therapeutic == True:

              therapeutic_r = TherapeuticVariants(gene=gene, enst_id=enst_id, hgvsp=p_code, hgvsg=g_code,
              hgvsc=c_code, exon=exon, variant_type=variant_type, consequence=consequence, depth=depth,
              allele_frequency=VAF,read_support=read_support,max_af=max_af,max_af_pop=max_af_pop,
              therapies=drugs_str,clinical_trials=clintrials_str,tumor_type=diseases_str)
              db.session.add(therapeutic_r)
              db.session.commit()

              result = s.TherapeuticTable.query.filter_by(user_id=p.analysis_env['USER_ID'])\
                .filter_by(lab_id=sample).filter_by(run_id=p.analysis_env['OUTPUT_NAME']).filter_by(gene=gene)\
                .filter_by(enst_id=enst_id).filter_by(hgvsp=p_code).filter_by(hgvsg=g_code).filter_by(hgvsc=c_code).first()
                
              if not result:
                therapeutic_db = s.TherapeuticTable(user_id=p.analysis_env['USER_ID'], lab_id=sample, ext1_id=ext1_id, 
                ext2_id=ext2_id, run_id=p.analysis_env['OUTPUT_NAME'], petition_id=petition_id,gene=gene, enst_id=enst_id, 
                hgvsp=p_code, hgvsg=g_code, hgvsc=c_code, exon=exon, variant_type=variant_type, consequence='.', 
                depth=depth, allele_frequency=VAF, read_support=read_support, max_af=max_af, max_af_pop=max_af_pop, 
                therapies=drugs_str, clinical_trials=clintrials_str, tumor_type=diseases_str, var_json=var_json_str,
                validated_bioinfo="no", validated_assessor="no", classification="Therapeutic")
                s.db.session.add(therapeutic_db)
                s.db.session.commit()

            # Filling Other clinical variants  
            if go_other == True:

              other_r = OtherClinicalVariants(gene=gene,  enst_id=enst_id, hgvsp=p_code, hgvsg=g_code,
              hgvsc=c_code, exon=exon, variant_type=variant_type, consequence=consequence, depth=depth,
              allele_frequency=VAF,read_support=read_support,max_af=max_af,max_af_pop=max_af_pop,therapies=drugs_str,
              clinical_trials=clintrials_str,tumor_type=diseases_str)
              db.session.add(other_r)
              db.session.commit()

              result = s.OtherVariantsTable.query.filter_by(user_id=p.analysis_env['USER_ID'])\
                .filter_by(lab_id=sample).filter_by(run_id=p.analysis_env['OUTPUT_NAME']).filter_by(gene=gene)\
                .filter_by(enst_id=enst_id).filter_by(hgvsp=p_code).filter_by(hgvsg=g_code).filter_by(hgvsc=c_code).first()

              if not result:
                other_db = s.OtherVariantsTable(user_id=p.analysis_env['USER_ID'], lab_id=sample, ext1_id=ext1_id, 
                ext2_id=ext2_id, run_id=p.analysis_env['OUTPUT_NAME'],petition_id=petition_id,gene=gene, enst_id=enst_id, 
                hgvsp=p_code, hgvsg=g_code, hgvsc=c_code, exon=exon, variant_type=variant_type, consequence=consequence, 
                depth=depth, allele_frequency=str(VAF),read_support=read_support,max_af=max_af,max_af_pop=max_af_pop,therapies=drugs_str,
                clinical_trials=clintrials_str,tumor_type=diseases_str, var_json=var_json_str,validated_bioinfo="no", validated_assessor="no",
                classification="Other")
                s.db.session.add(other_db)
                s.db.session.commit()

            # Filling Rare variants  
            if go_rare == True:

              rare_v = RareVariants(gene=gene, enst_id=enst_id, hgvsp=p_code, hgvsg=g_code,
              hgvsc=c_code, exon=exon, variant_type=variant_type, consequence=consequence, depth=depth,
              allele_frequency=VAF,read_support=read_support,max_af=max_af,max_af_pop=max_af_pop,therapies=drugs_str,
              clinical_trials=clintrials_str,tumor_type=diseases_str)  
              db.session.add(rare_v)
              db.session.commit()

              result = s.RareVariantsTable.query.filter_by(user_id=p.analysis_env['USER_ID'])\
                .filter_by(lab_id=sample).filter_by(run_id=p.analysis_env['OUTPUT_NAME']).filter_by(gene=gene)\
                .filter_by(enst_id=enst_id).filter_by(hgvsp=p_code).filter_by(hgvsg=g_code).filter_by(hgvsc=c_code).first()

              if not result:
                rare_db = s.RareVariantsTable(user_id=p.analysis_env['USER_ID'], lab_id=sample, ext1_id=ext1_id, 
                ext2_id=ext2_id, run_id=p.analysis_env['OUTPUT_NAME'],petition_id=petition_id,gene=gene, enst_id=enst_id, 
                hgvsp=p_code, hgvsg=g_code, hgvsc=c_code, exon=exon, variant_type=variant_type, consequence=consequence, 
                depth=depth, allele_frequency=str(VAF), read_support=read_support, max_af=max_af,max_af_pop=max_af_pop,therapies=drugs_str,
                clinical_trials=clintrials_str,tumor_type=diseases_str, var_json=var_json_str,validated_bioinfo="no", validated_assessor="no",
                classification="Rare")
                s.db.session.add(rare_db)
                s.db.session.commit()

            if 'pathogenic' in clin_sig:
              o.write('\t'.join(results_list)+ "\n")
        
          var_id = ""

          # Checking that the variant has been detected previously
          var = s.Variants.query.filter_by(chromosome=chromosome).filter_by(pos=pos)\
            .filter_by(ref=ref).filter_by(alt=alt).filter_by(var_type=variant_type).first()
          if not var:
            count = 1

            new_variant= s.Variants(chromosome=chromosome, pos=pos, ref=ref, alt=alt,
            var_type=variant_type, genome_version=p.analysis_env['GENOME_VERSION'], gene=gene, 
            isoform=isoform, hgvsg=g_code, hgvsp=p_code, hgvsc=c_code,count=count)
            s.db.session.add(new_variant)
            var_id = new_variant.id
          else:
            if not sample_id:
              var.count = var.count+1
            var_id = var.id

          found_sample_var =  s.SampleVariants.query.filter_by(var_id=var_id).first()
          if not found_sample_var:
          # Recording sample variants
            sample_var = s.SampleVariants(sample_id=sample_id, var_id=var_id, ann_id=1, 
            lab_confirmation="None", confirmation_technique="None", classification=classification,
            ann_json=ann_json_str, ann_key=ann_key)       
            s.db.session.add(sample_var)
            s.db.session.commit()

        o.close()

        # Create PDF report 
        bashCommand = ('{} pr {} -r {} -f pdf -t generic --db-url jdbc:sqlite:{} --db-driver org.sqlite.JDBC -o {} --jdbc-dir {}') \
          .format(p.system_env['JASPERSTARTER'],p.aux_env['REPORT_JRXML'], 
          p.defaults['JASPERREPORT_FOLDER'], sample_db, p.sample_env[sample]['REPORT_PDF'], p.defaults['JDBC_FOLDER'])
        msg = " INFO: Generating pdf report for sample " + sample
        print (msg)
        print(bashCommand)
        logging.info(bashCommand)
        
        p1 = subprocess.run(bashCommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = p1.stdout.decode('UTF-8')
        error  = p1.stderr.decode('UTF-8')

        if not error and not output:
          msg = " INFO: Report generation for" + sample + " ended OK"
          print(msg)
          logging.info(msg)
        else:
          msg = " ERROR: Failed report generation for " + sample
          print(msg)
          logging.error(msg)

def clinical_vcf():

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
                        civic_dict, civic_list = u.create_vep_dict(line)
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

                                        if transcript_list[civic_dict['EV_DIRECTION']] == "Supports" \
                                            or 'pathogenic' in transcript_list[vep_dict['CLIN_SIG']]:
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
        # Create JSON file
        u.convert_vcf_2_json(p.sample_env[sample]['CLINICAL_SNV_VCF'])

        #$jasperst pr $report -f pdf -t generic --db-url jdbc:sqlite:$dbName --db-driver org.sqlite.JDBC -o $pdfReport -r
        

