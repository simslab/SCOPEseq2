#! /usr/bin/python3.4
import numpy as np

# The purpose of this program is to apply an H=1 Hamming filter to the UMIs
# in a DropSeq experiment. Note that this function is partially converted into
# C code using cython and compiled for speed.

def umifilter(addresscts_file,filtercts_file): # input file is the output of DropSeqPipeline6_addresscts.py
	gene_dict = {} # dictionary where keys are gene symbols, values are tuples of cell barcode, UMI sequence, number of reads
	init_cts = 0
	filt_cts = 0
	with open(addresscts_file) as f:
		for line in f:
			llist = line.split()
			gene = llist[2] # gene symbol
			tup = llist[0],llist[1],int(llist[3]) # tuple containing cell barcode, UMI sequence, and number of reads
			if gene in gene_dict.keys():
				gene_dict[gene].append(tup)
			else:
				gene_dict[gene] = [tup]
			init_cts+=1
	with open(filtercts_file,'w') as g:
		for gene in gene_dict.keys(): # only need to compare UMIs with the same gene symbol
			tups = gene_dict[gene]
			cells = {} # pre-sort based on cell barcodes- make a dictionary where the keys are cell barcodes, values are tuple with UMI sequence, number of reads
			for tup in tups:
				if tup[0] in cells.keys():
					cells[tup[0]].append((np.array(list(tup[1])),tup[2]))
				else:
					cells[tup[0]] = [(np.array(list(tup[1])),tup[2])]
			keep = []
			for cell in cells.keys(): # only need to compare UMIs with the same cell barcode
				tups = cells[cell] 
				for tup in tups:
					umi = tup[0]
					cts = tup[1]
					go = 1
					for tup2 in tups: # loop through all UMI-UMI pairs with the same cell barcode and gene symbol
						if tup2[1] > cts: # if you find a second UMI with more reads than the first UMI
							d = np.count_nonzero(umi!=tup2[0]) 
							if d < 2: # and that second UMI is with a Hamming distance of 1 from the first
								go = 0 # the first UMI is discarded
								break
					if go == 1:
						st = cell+'\t'+''.join(umi)+'\t'+gene+'\t'+str(cts)+'\n'
						g.write(st)
						filt_cts+=1
	return filt_cts,init_cts

