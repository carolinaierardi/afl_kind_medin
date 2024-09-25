#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 17:04:54 2024

@author: carolinaierardi
"""

import pandas as pd
import numpy as np 
import os
from io import StringIO
from Bio import SeqIO
import requests




#variables to change
current_wd = "/Users/carolinaierardi/Documents/Cambridge/SummerInternship/Data"
url = 'https://rest.uniprot.org/uniprotkb/Q08431.fasta'


os.chdir(current_wd)                              #change wd

input_fasta = requests.get(url).text              #obtain sequences from website
input_fasta = StringIO(input_fasta)               #make into stringIo object

input_fasta = list(SeqIO.parse(input_fasta,'fasta'))

print(f"Length of medin: {len(str(input_fasta[0].seq[267:317]))}")


#CHANGE THIS DEPENDING ON WHICH SEQUENCE IS NEEDED
with open("medin.fasta","w") as f:
        for seq_record in input_fasta:
                f.write(">Medin " + "\n")
                f.write(str(seq_record.seq[267:317]) + "\n")  
