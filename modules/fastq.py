#!/usr/bin/env python3

import sys
import re
import gzip
import binascii
import os.path
from os import path
import argparse
from pathlib import Path

class Fastq:
    '''given a fq file name we check
       if gz
       if fastq
       extract R1 or R2 
       extract SX
       extract code and the format
    '''
    def __init__(self, fq1, fq2, paired):
        # Defining attributes
        self.fq1 = fq1
        self.fq2 = fq2
        self.paired = paired

    #def sampleName(self):
        
    #def sampleNum(self):

    #def readNum(self):

    def readLen(self):
        ''' Get read length
        '''
        length_list_fq = list()
        if is_gz(self.fq1):
            with gzip.open(self.fq1, 'rb') as f:
                for line in f:
                    length_list_fq.append(len(line))
        else:
            with open(self.fq1, 'r') as f:
                for line in f:
                   length_list_f1.append(len(line))         
        return most_frequent(length_list_fq1), 

    def isValid (self): # Returns a bool: True or False
        ''' Validate a fastq file
            Returning: True or False
        '''
        # For paired-end mode
        if self.paired == True:

            # First, Use regex to validate suffix
            if re.search('(fq$|fastq$|fastq.gz$|fa.gz$)', self.fq1):
                pass
            else:
                return False
            if re.search('(fq$|fastq$|fastq.gz$|fa.gz$)', self.fq2):
                pass
            else:
                return False

            # Second validate same names between fq1 and fq2
            fq1_tmp = os.path.basename(self.fq1).split("_")
            fq2_tmp = os.path.basename(self.fq2).split("_")

            fq1_name = fq1_tmp[0]
            fq2_name = fq2_tmp[0]
            
            if fq1_name == fq2_name:
                pass
            else:
                return False

            # Third, checking illumina's nomenclature
            if re.search('_S[0-9]+_L[0-9]+_R[12]_[0-9]+', self.fq1):
                pass
            else:
                return False
            if re.search('_S[0-9]+_L[0-9]+_R[12]_[0-9]+', self.fq2):
                pass
            else:
                return False
            return True
            # Third check if f1 and f2 have equal total lines
            # nlines_fq1 = count_lines(self.fq1)
            # nlines_fq2 = count_lines(self.fq2)
            # if nlines_fq1 == nlines_fq2:
            #     return True
            # else:
            #     return False
        # For single-ende mode
        else:
            if re.search('(fq$|fastq$|fastq.gz$|fa.gz$)', self.fq1) is None:
                return False
            if re.search('_S[0-9]+_L[0-9]+_R[12]_[0-9]+', self.fq1) is None:
                return False
            # Returnning true by default but
            # we should check at least the commont structure
            # of a fastq file (4lines/read, header, etc)
            return True


def most_frequent(List): 
    counter = 0
    num = List[0] 
      
    for i in List: 
        curr_frequency = List.count(i) 
        if(curr_frequency> counter): 
            counter = curr_frequency 
            num = i 
  
    return num 

# Recognize if fastq is gzipped
# https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed
def is_gz(file):
    with open(file, 'rb') as f:
        return binascii.hexlify(f.read(2)) == b'1f8b'

def count_lines(file):
    '''Count total lines of a file.
       detect if gzipped
    '''
    count = 0
    if is_gz(file):
        with gzip.open(file, 'rb') as f:
            for line in f:
                count += 1
    else:
        with open(file, 'r') as f:
            for line in f:
                count += 1
    return count