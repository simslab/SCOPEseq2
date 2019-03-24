#! /usr/bin/python
from SCOPEseq2_demultiplexer import get_cbc_umi
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

def parse_user_input():
	parser = argparse.ArgumentParser()
	parser.add_argument('-r','--run-name',required=True,help='Name of sequencing run.')
	parser.add_argument('-d','--dir-name',required=True,help='Name of I/O directory.')
	parser.add_argument('-r1','--read1-fastq',required=True,help='Path to gzipped read 1 fastq file.')
	parser.add_argument('-r2','--read2-fastq',required=True,help='Path to gzipped reqd 2 fastq file.')
	parser.add_argument('-bc1','--barcode1-infile',required=True,help='Path to file containing first cell-identifying barcode segment.')
	parser.add_argument('-bc2','--barcode2-infile',required=True,help='Path to file containing second cell-identifying barcode segment.')
	return parser

parser = parse_user_input()
user_input = parser.parse_args()


