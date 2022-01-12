#!/usr/bin/env Rscript
library(ape)
options(width=150)

yemen = Sys.getenv('yemen')

lenchr1 = 3046494
reflaneid = '33224_2#241'

laneid2col = function(laneids){
	gsub("^([23]{2}[0-9]{3}_[0-9])#", "X\\1.", laneids)
}
col2laneid = function(cols){
	sub('(.+).1_\\1_genomic', '\\1', gsub('([23]{2}[0-9]{3}_[0-9])\\.([0-9])', '\\1#\\2', sub('^X', '', cols)))
}

sangerids = read.table(file.path(yemen, 'metadata/Yemen_and_context_pf-info.sampleid2laneid.csv'), comment.char='', quote='', sep=',', header=T)


sanglane2assid = function(inid, outidset, ref=sangerids, incol='Lane', outcol='Sample'){
	if (inid %in% outidset){
		return(inid)
	}else{
		newid.i = which(inid==ref[,incol])
		if (length(newid.i)==1){
			return(ref[newid.i,outcol])
		}else{
			return(inid)
		}
	}
}

cargs = commandArgs(trailingOnly=T)
dirpanaroo = cargs[1]
nffulltree = cargs[2]
nfmainclade = cargs[3]
nfclade2018 = cargs[4]

# Pangenome abs/pres profiles
panarooannot = read.csv(file.path(dirpanaroo, 'gene_annotations.csv'), comment.char='', header=T)
panaroorefloctag = read.csv(file.path(dirpanaroo, sprintf('gene_presence_absence_roary_%s.csv', reflaneid)), comment.char='', header=T)
refltcol = make.names(paste0('loctag_', reflaneid))
colnames(panaroorefloctag) = c("Gene", refltcol)
panaroopangenomegenes = read.table(file.path(dirpanaroo, 'gene_presence_absence.Rtab'), comment.char='', sep='\t', header=T)
rownames(panaroopangenomegenes) = panaroopangenomegenes[,'Gene']
panaroopangenomestruct = read.table(file.path(dirpanaroo, 'struct_presence_absence.Rtab'), comment.char='', sep='\t', header=T)
rownames(panaroopangenomestruct) = panaroopangenomestruct[,'Gene']
npanacols = ncol(panaroopangenomegenes)
stopifnot(all(colnames(panaroopangenomegenes) == colnames(panaroopangenomestruct)))
stopifnot(all(panarooannot[,'Gene'] == panaroopangenomegenes[,'Gene']))
stopifnot(all(panaroorefloctag[,'Gene'] == panaroopangenomegenes[,'Gene']))
panaroopangenomegenes = merge(panaroopangenomegenes, panarooannot)
panaroopangenomegenes = merge(panaroopangenomegenes, panaroorefloctag)
print(head(colnames(panaroopangenomegenes)))
colpanaroo2isol = col2laneid(colnames(panaroopangenomegenes)[2:npanacols])
isol2colpanaroo = colnames(panaroopangenomegenes)[2:npanacols]
names(isol2colpanaroo) = colpanaroo2isol

fulltree = read.tree(nffulltree)
mainclade = read.tree(nfmainclade) # includes ref
clade2018 = read.tree(nfclade2018)

panacolingenomes = which(colpanaroo2isol %in% fulltree$tip.label)

fulltree.cols = intersect(sapply(fulltree$tip.label, sanglane2assid, outidset=colpanaroo2isol), colpanaroo2isol)
mainclade.cols = intersect(sapply(mainclade$tip.label, sanglane2assid, outidset=colpanaroo2isol), colpanaroo2isol)
clade2018.cols = intersect(sapply(clade2018$tip.label, sanglane2assid, outidset=colpanaroo2isol), colpanaroo2isol)

# previous (clade with mostly 2016-2017 Yemen isolates) vs. new (clade with mostly 2018-2019 Yemen isolates)
bgclade.cols = setdiff(fulltree.cols, mainclade.cols)
prevclade.samples = grepl("^CNRVC", mainclade.cols) | (mainclade.cols %in% c("X33224_2.7", "X33224_2.33"))
prevclade.cols = mainclade.cols[prevclade.samples]
print("prevclade.cols")
print(prevclade.cols)
newclade.cols = mainclade.cols[!prevclade.samples] # includes ref
print("newclade.cols")
print(newclade.cols)
# 2018 (clade with only and almost all 2018 isolates) vs the bg of new
print("clade2018.cols")
print(clade2018.cols)
samples.2018 = mainclade.cols %in% clade2018.cols
bgnewclade.cols = setdiff(newclade.cols, clade2018.cols) # includes ref
print("bgnewclade.cols")
print(bgnewclade.cols)

dataset.lengths = c(length(fulltree.cols), length(mainclade.cols), length(bgclade.cols), length(prevclade.cols), length(newclade.cols), length(bgnewclade.cols), length(clade2018.cols))
names(dataset.lengths) = c("full.tree", "main.Yemen.clone.clade", "background.divergent.clade", "mostly2016-2017.sample.clade", "mostly2018-2019.sample.clade", "mostly2019.sample.clade", "2018.sample.clade")
print("dataset.lengths")
print(dataset.lengths)


panatabs = list(gene=panaroopangenomegenes, struct=panaroopangenomestruct)

for (npanatab in names(panatabs)){
	panatab = panatabs[[npanatab]]
	panamat = data.matrix(panatab[, panacolingenomes+1])
	colnames(panamat) = colpanaroo2isol[panacolingenomes]
	panamatfreq = as.double(apply(panamat, 1, sum, na.rm=T)) / ncol(panamat)
	midfreqs = which((panamatfreq > 0.02) & (panamatfreq < 0.98))
	panamat = panamat[midfreqs,]
	pananames = panatab[midfreqs,1]
	if ("Annotation" %in% colnames(panatab)){
		panannot = panatab[midfreqs,c("Annotation", refltcol)]
	}else{ panannot = NULL }
	print(sprintf('    size of pangenome %s presence/absence matrix: %d, %d', npanatab, nrow(panamat), ncol(panamat)), quote=F)

	# main (main yemen cholera clone) vs bg (Yemen divergent isolates)
	mainclade.accgene.mat = panamat[,mainclade.cols]
	bgclade.accgene.mat = panamat[, bgclade.cols]
	accgene.freq.bg = apply(bgclade.accgene.mat, 1, sum) / length(bgclade.cols)
	accgene.freq.main = apply(mainclade.accgene.mat, 1, sum) / length(mainclade.cols)
	uniquevar.bgvsmainclade = (accgene.freq.bg>0.8 & accgene.freq.main<0.2)
	uniquevar.mainvsbgclade = (accgene.freq.main>0.8 & accgene.freq.bg<0.2)
	contrast.mainvsbgclade = uniquevar.mainvsbgclade | uniquevar.bgvsmainclade
	#var.in.mainclade = apply(mainclade.accgene.mat, 1, sum)>0
	#mainclade.filt.accgene.mat = mainclade.accgene.mat[var.in.mainclade | contrast.mainvsbgclade,]
	#bgclade.filt.accgene.mat = bgclade.accgene.mat[var.in.mainclade | contrast.mainvsbgclade,]

	# previous (clade with mostly 2016-2017 Yemen isolates) vs. new (clade with mostly 2018-2019 Yemen isolates)
	newclade.accgene.mat = panamat[,newclade.cols]
	prevclade.accgene.mat = panamat[, prevclade.cols]
	accgene.freq.prev = apply(prevclade.accgene.mat, 1, sum) / length(prevclade.cols)
	accgene.freq.new = apply(newclade.accgene.mat, 1, sum) / length(newclade.cols)
	uniquevar.prevvsnewclade = (accgene.freq.prev>0.8 & accgene.freq.new<0.2)
	uniquevar.newvsprevclade = (accgene.freq.new>0.8 & accgene.freq.prev<0.2)
	contrast.newvsprevclade = uniquevar.newvsprevclade | uniquevar.prevvsnewclade

	## 2018 vs the bg of main
	#clade2018.accgene.mat = panamat[,clade2018.cols]
	#bgmainclade.accgene.mat = panamat[, setdiff(mainclade.cols, clade2018.cols)]
	#unique.2018vsbgmainclade = (apply(clade2018.accgene.mat, 1, sum)>0.8 & apply(bgmainclade.accgene.mat, 1, sum)<0.2)
	#unique.bgmainvs2018clade = (apply(bgmainclade.accgene.mat, 1, sum)>0.8 & apply(clade2018.accgene.mat, 1, sum)<0.2)
	#contrast.mainvsbgclade = unique.2018vsbgmainclade | unique.bgmainvs2018clade
	#clade2018.filt.accgene.mat = clade2018.accgene.mat[var.in.mainclade | contrast.mainvsbgclade,]
	#bgmainclade.filt.accgene.mat = bgmainclade.accgene.mat[var.in.mainclade | contrast.mainvsbgclade,]

	# 2018 (clade with only and almost all 2018 isolates) vs the bg of new
	clade2018.accgene.mat = panamat[,clade2018.cols]
	bgnewclade.accgene.mat = panamat[, bgnewclade.cols]
	accgene.freq.bgnew = apply(bgnewclade.accgene.mat, 1, sum) / length(bgnewclade.cols)
	accgene.freq.2018 = apply(clade2018.accgene.mat, 1, sum) / length(clade2018.cols)
	uniquevar.2018vsbgnewclade = (accgene.freq.2018>0.8 & accgene.freq.bgnew<0.2)
	uniquevar.bgnewvs2018clade = (accgene.freq.bgnew>0.8 & accgene.freq.2018<0.2)
	contrast.2018vsbgnewclade = uniquevar.2018vsbgnewclade | uniquevar.bgnewvs2018clade
	#clade2018.filt.accgene.mat = clade2018.accgene.mat[var.in.mainclade | contrast.2018vsbgnewclade,]
	#bgnewclade.filt.accgene.mat = bgnewclade.accgene.mat[var.in.mainclade | contrast.2018vsbgnewclade,]

	#fixed.panag = which(accgene.freq.prev>0.8 | accgene.freq.new>0.8 | accgene.freq.2018>0.8 | accgene.freq.bgnew>0.8)
	combined.contrasts = which(contrast.mainvsbgclade | contrast.newvsprevclade | contrast.2018vsbgnewclade)

	constrasts = list(bgVSmain=contrast.mainvsbgclade,
					  prevVSnew=contrast.newvsprevclade,
					  bgnewVS2018=contrast.2018vsbgnewclade,
					  allcontrasts=combined.contrasts)
	for (contrasttag in names(constrasts)){
	  fixed.panag = constrasts[[contrasttag]]
	  fixed.panag.info = cbind(pananames[fixed.panag],
							  accgene.freq.bg[fixed.panag], accgene.freq.main[fixed.panag],
							  accgene.freq.prev[fixed.panag], accgene.freq.new[fixed.panag],
							  accgene.freq.bgnew[fixed.panag], accgene.freq.2018[fixed.panag])
	  outcols = paste0(npanatab, c("", ".freq.bgH.clade", ".freq.H89.clade", ".freq.H9ef.clade", ".freq.H9gh.clade", ".freq.H9g.clade", ".freq.H9h.clade"))
	  colnames(fixed.panag.info) = outcols
      if (!is.null(panannot)){ 
		  fixed.panag.info = cbind(fixed.panag.info, panannot[fixed.panag,]) 
	      colnames(fixed.panag.info) = c(outcols, c("Annotation", refltcol))
		  fixed.panag.info = fixed.panag.info[order(fixed.panag.info[,refltcol]),]
	  }
	  write.table(format(fixed.panag.info, digits=5), file=paste0(dirpanaroo, sprintf("%s.fixed.panaroo_%s.tsv", contrasttag, npanatab)),
				  sep='\t', row.names=F, quote=FALSE)
	}
}