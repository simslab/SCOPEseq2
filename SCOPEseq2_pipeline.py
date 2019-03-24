#! /usr/bin/python
from SCOPEseq2_demultiplexer import get_cbc_umi
from SCOPEseq2_clipper import clipper
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import os

def parse_user_input():
	parser = argparse.ArgumentParser()
	parser.add_argument('-r','--run-name',required=True,help='Name of sequencing run.')
	parser.add_argument('-data','--data-dir',required=True,help='Name of I/O directory.')
	parser.add_argument('-r1','--read1-fastq',required=True,help='Path to gzipped read 1 fastq file.')
	parser.add_argument('-r2','--read2-fastq',required=True,help='Path to gzipped reqd 2 fastq file.')
	parser.add_argument('-bc1','--barcode1-infile',required=True,help='Path to file containing first cell-identifying barcode segment.')
	parser.add_argument('-bc2','--barcode2-infile',required=True,help='Path to file containing second cell-identifying barcode segment.')
	parser.add_argument('-ref','--reference-dir',required=True,help='Name of reference genome directory for STAR.')
	parser.add_argument('-gtf','--gtf',required=True,help='Path to transcriptome annotation gtf file for STAR.')
	parser.add_argument('-dist','--overhang_distance',default=65,help='Sets the sjdbOverhang parameter for STAR.')
	parser.add_argument('-t','--threads',default=1,help='Sets the number of threads for STAR and samtools.')
	return parser

parser = parse_user_input()
user_input = parser.parse_args()

barcode1_list = [line.split()[0] for line in open(user_input.barcode1_infile)]
barcode2_list = [line.split()[0] for line in open(user_input.barcode2_infile)]

cbcs,umis = get_cbc_umi(barcode1_list,barcode2_list,user_input.read1_fastq) # get cell-identifying barcodes (CBCs) and unique molecular identifiers (UMIs)
demux_N = sum([1 for bc in cbcs if bc != '0']) # number of demultiplexed reads, '0' indicates no CBC
reads_N = len(cbcs) # number of reads
demux_P = float(demux_N)/float(reads_N)*100. # percentage of reads demultiplexed

print('Found %(demux_N)d reads with CBC/UMI out of %(reads_N)d reads or %(demux_P)f%% demultiplexed...' % vars())

ref = user_input.reference_dir
gtf = user_input.gtf
t = user_input.threads
fq2file = user_input.read2_fastq
fq2clip = user_input.data_dir+'/'+user_input.run_name+'_R2.clip.fastq.gz'
bamfile = user_input.data_dir+'/'+user_input.run_name
dist = user_input.overhang_distance

clipped_N,total_N = clipper(fq2file,fq2clip)
clipped_P = float(clipped_N)/float(total_N)
print('Found %(clipped_N)d reads with poly(A) out of %(total_N)d reads or %(clipped_P)f%% clipped...' % vars())

cmd = '/home/ubuntu/Software/STAR/bin/Linux_x86_64/STAR --readFilesCommand zcat --genomeDir %(ref)s --sjdbOverhang %(dist)s --sjdbGTFfile %(gtf)s --twopassMode Basic --runThreadN %(t)s --readFilesIn %(fq2clip)s --outFileNamePrefix %(bamfile)s --outSAMtype BAM Unsorted --outSAMunmapped Within' % vars()
print('STAR command...')
print(cmd)
os.system(cmd)
