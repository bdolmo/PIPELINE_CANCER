#!/usr/bin/env python3

import sys
import re
import gzip
import binascii
import os.path
import glob
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

global app 
global db

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///NGS_DB.db'

db = SQLAlchemy(app)
app.config.update(dict(
     SECRET_KEY="powerful secretkey",
    WTF_CSRF_SECRET_KEY="a csrf secret key"
))

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
        self.id        = ID
        self.chromosome= CHROMOSOME
        self.start     = START
        self.end       = END
        self.ensg_id   = ENSG_ID
        self.enst_id   = ENST_ID
        self.gene_name = GENE_NAME
        self.genome_version = GENOME_VERSION
        self.panel_name = PANEL_NAME
        self.panel_version = PANEL_VERSION

def load_panel_transcripts(pname):
    roi_env = defaultdict(dict)
    roi_info = Roi.query.filter_by(panel_name=pname).all()
    if roi_info:
        for roi in roi_info:
            roi_env[roi.ensg_id]   = roi.enst_id
            roi_env[roi.gene_name] = roi.enst_id
    return roi_env