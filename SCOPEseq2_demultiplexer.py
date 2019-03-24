#! /usr/bin/python
import gzip
import io

# Enumerate HD=1 sequences for a given barcode sequence
def enumerate_bc(bc):
	bcs = []
	N = len(bc)
	for i in range(N):
		bc1=''
		bc2=''
		bc3=''
		bc4=''
		for j in range(N):
			if i==j:
				bc1+='A'
				bc2+='G'
				bc3+='C'
				bc4+='T'
			else:
				bc1+=bc[j]
				bc2+=bc[j]
				bc3+=bc[j]
				bc4+=bc[j]
		bcs.append(bc1)
		bcs.append(bc2)
		bcs.append(bc3)
		bcs.append(bc4)	
	bcs = list(set(bcs))
	return bcs

# Generate dictionary linking HD=1 sequences (keys) to correct barcode (values)
def get_cbc_dict(bcs):
	cbc_dict = {}
	for bc in bcs:
		allbcs = enumerate_bc(bc)
		for allbc in allbcs:
			cbc_dict[allbc] = bc
	return cbc_dict	

# Extract SCOPEseq2 cell barcodes (CBCs) and unique molecular identifiers (UMIs) for read 1 fastq
def get_cbc_umi(first_bc,second_bc,read1_fastq):
	first_bc_dict = get_cbc_dict(first_bc)
	second_bc_dict = get_cbc_dict(second_bc)
	with io.BufferedReader(gzip.open(read1_fastq,'rb')) as f:
		next(f)
		cbcs = []
		umis = []
		for ct,line in enumerate(f,start=0):
			if ct%4 == 0:
				dline = line.decode()
				bc1 = dline[2:10] # first 8-base barcode segment
				if bc1 in first_bc_dict.keys():
					bc1 = first_bc_dict[bc1]
					bc2 = dline[12:20] # second 8-base barcode segment
					if bc2 in second_bc_dict.keys():
						bc2 = second_bc_dict[bc2] 
						umi = dline[0:2]+dline[10:12]+dline[20:24] # 8-base UMI in 2x2-base and 1x4-base blocks
						if umi.find('N') == -1:
							cbcs.append(bc1+bc2)
							umis.append(umi)
						else:
							cbcs.append('0')
							umis.append('0')
					else:
						cbcs.append('0')
						umis.append('0')
				else:
					cbcs.append('0')
					umis.append('0')
	return cbcs,umis
						
