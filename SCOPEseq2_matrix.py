#! /usr/bin/python
import numpy as np

def total_matrix(bcfile,genefile,filtercts_file,matrix_file,hist_file):
	barcodes = {}
	bclist = []
	i=0
	with open(bcfile) as f:
		for line in f:
			bc = line.split()[0]
			bclist.append(bc)
			barcodes[bc] = i
			i+=1
	gidgene = {}
	matrix = {}
	with open(genefile) as f:
		for line in f:
			llist = line.split()
			gid = llist[0]
			gidgene[gid] = llist[1]
			matrix[gid] = [0 for bc in bclist]	
	with open(filtercts_file) as f:
		for line in f:
			llist = line.split()
			matrix[llist[2]][barcodes[llist[0]]]+=1
	with open(matrix_file,'w') as g:
		for gid in gidgene.keys():
			gene = gidgene[gid]
			st = gid+'\t'+gene+'\t'+'\t'.join([str(pt) for pt in matrix[gid]])+'\n'
			g.write(st)
	matrix = np.array(list(matrix.values()))
	norm = sum(matrix)
	with open(hist_file,'w') as g:
		for bc,n in zip(bclist,norm):
			st = bc+'\t'+str(n)+'\n'
			g.write(st)
	return 0
