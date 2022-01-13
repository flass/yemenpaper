#!/usr/bin/env python2.7
import os, sys
from glob import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

dcompl = {'A':'T','T':'A','G':'C','C':'G'}
addheader = ['VARTYPE', 'LOCTAG', 'PRODUCT', 'NOTE']

refgbk = sys.argv[1]
contrastsrad = sys.argv[2]
lnfsnptab = glob(contrastsrad+'*.fixed.snps.tab')
contrasts = [nfsnptab.replace(contrastsrad+'.', '').replace('.fixed.snps.tab', '') for nfsnptab in lnfsnptab]

refseqrecs = [rec for rec in SeqIO.parse(refgbk, format='genbank')]

dchrsnp = {}
for contrast in contrasts:
	print contrast
	nfsnprad = "%s.%s.fixed.snps"%(contrastsrad, contrast)
	nfsnptab = nfsnprad+".tab"
	nfsnpout = nfsnprad+"annots.tab"
	fout = open(nfsnpout, 'w')
	currchrom = 0
	currfeat = None
	nfeat = None
	with open(nfsnptab, 'r') as fsnptab:
		header = fsnptab.readline().rstrip('\n').split('\t')
		fout.write('\t'.join(header+addheader)+'\n')
		ipos = header.index('POS')
		iref = header.index('REF')
		ialt = header.index('ALT')
		ichr = header.index('X.CHROM')
		for line in fsnptab:
			lsp = line.rstrip('\n').split('\t')
			chrom = int(lsp[ichr].strip())
			if chrom != currchrom:
				currchrom = chrom
				print '', 'chr%d'%chrom
				dchrsnp[chrom] = {}
				refseqrec = refseqrecs[chrom-1]
				nfeat = len(refseqrec.features)
				currfeat = 0
			pos = int(lsp[ipos].strip())
			refbase = lsp[iref].strip()
			altbase = lsp[ialt].strip().split(',')[0] # if several alternative states, arbitrarily take the first
			dchrsnp[chrom][pos] = lsp
			lposfeat = []
			lastfeat = None
			# find annotations
			for f in range(currfeat, nfeat):
				feature = refseqrec.features[f]
				lastfeat = feature
				currfeat = f
				if pos in feature:
					lposfeat.append(feature)
				elif lposfeat:
					if feature.location.end > lposfeat[-1].location.end:
						# already found some feature matches, and we're scanning past their range now; can stop now
						break
				elif feature.location.start > pos:
					# went past the position
					break
			vartype = 'IG'
			note = ''
			for pf, feature in enumerate(lposfeat):
				quals = feature.qualifiers
				if feature.type=='gene':
					nextfeature = lposfeat[pf+1]
					if nextfeature.type=='CDS' and feature.location==nextfeature.location:
						# skip the protein coding gene as the corresponding CDS feture comes next
						continue
					else:
						# non coding gene
						vartype = 'NC'
				elif feature.type=='CDS':
					if feature.strand==1:
						relpos = pos - int(feature.location.start) - 1
						srefbase = refbase
						saltbase = altbase
					else:
						relpos = int(feature.location.end) - pos
						srefbase = dcompl[refbase]
						saltbase = dcompl[altbase]
					cdsseq = feature.extract(refseqrec).seq
					if (relpos < 0) or (relpos >= len(cdsseq)):
						print "Warning, base is out of range"
						print "position on genomic record:", pos, 'position in CDS:', relpos
						print feature
						vartype = 'NA'
					else:
						if not str(cdsseq[relpos])==srefbase:
							print feature
							print cdsseq[relpos-10:relpos+11]
							print '"         ^         "'
							print "position on genomic record:", pos, 'position in CDS:', relpos
							raise IndexError, "extracted base '%s' is different from reference variant '%s'"%(cdsseq[relpos], srefbase)
#							print "Warning, extracted base '%s' is different from reference variant '%s'"%(cdsseq[relpos], srefbase)
						mutseq = cdsseq[:relpos] + Seq(saltbase, IUPACAmbiguousDNA()) + cdsseq[(relpos+1):]
						refprot = cdsseq.translate()
						mutprot = mutseq.translate()
						if mutprot == refprot:
							vartype = 'SY'
						else:
							vartype = 'NS'
							for i in range(len(refprot)):
								refaa = refprot[i]
								mutaa = mutprot[i]
								if mutaa != refaa:
									aapos = i+1
									note = "%s%d%s"%(refaa, aapos, mutaa)
				else:
					vartype = ''
				# write results
				outline = '\t'.join(lsp+[vartype, quals.get('locus_tag', [''])[0], quals.get('product', [''])[0], note])+'\n'
				fout.write(outline)

	fout.close()

