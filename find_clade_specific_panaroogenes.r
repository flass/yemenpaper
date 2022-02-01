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
lnfclades = cargs[3:length(cargs)]
# process arguments to recover dataset tags
refstests = strsplit(lnfclades, split=',')
tagnfclades = lapply(1:length(refstests), function(k){
	s = refstests[[k]]
	if (length(s) != 2){ stop(sprintf("expected a string containing 2 clade files in Newick format separated by a comma; got: '%s'", paste(s,collapse=','))) }
	tagnf = sapply(strsplit(s, split='='), function(tn){
		if (length(tn)==1) return(c(NA, tn))
		else return(tn[1:2]) 
	})
	lnf = tagnf[2,]
	names(lnf) = ifelse(is.na(tagnf[1,]), paste(c("reference.clade", "test.clade"), k, sep='.'), tagnf[1,])
	return(lnf)
})
fulltree = read.tree(nffulltree)
fulltree.cols = laneid2col(fulltree[['tip.label']])

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
panaroopangenomegenes = merge(panaroopangenomegenes, panarooannot, by='Gene')
panaroopangenomegenes = merge(panaroopangenomegenes, panaroorefloctag, by='Gene')
colpanaroo2isol = col2laneid(colnames(panaroopangenomegenes)[2:npanacols])
isol2colpanaroo = colnames(panaroopangenomegenes)[2:npanacols]
names(isol2colpanaroo) = colpanaroo2isol
panacolingenomes = which(colpanaroo2isol %in% fulltree$tip.label)
fulltree.cols = intersect(sapply(fulltree$tip.label, sanglane2assid, outidset=colpanaroo2isol), colpanaroo2isol)

panatabs = list(gene=panaroopangenomegenes, struct=panaroopangenomestruct)

for (npanatab in names(panatabs)){
	panatab = panatabs[[npanatab]]
	panamat = data.matrix(panatab[, panacolingenomes+1])
	colnames(panamat) = colpanaroo2isol[panacolingenomes]
	panamatfreq = as.double(apply(panamat, 1, sum, na.rm=T)) / ncol(panamat)
	midfreqs = which((panamatfreq > 0.02) & (panamatfreq < 0.98))
	if ("Annotation" %in% colnames(panatab)){
		panannot = panatab[midfreqs,c("Annotation", refltcol)]
	}else{ panannot = NULL }
	panamat = panamat[midfreqs,]
	pananames = panatab[midfreqs,1]
	print(sprintf('    size of intermediate frequency (.02-.98) pangenome %s presence/absence matrix: %d, %d', npanatab, nrow(panamat), ncol(panamat)), quote=F)
	
	accgenecontrasts = list()
	accgene.freq.full = panamatfreq[midfreqs]
	accgenefreqvec = list(accgene.freq.fulltree=accgene.freq.full)


	for (tagnf in tagnfclades){
		nfrefclade = tagnf[[1]]
		nftestclade = tagnf[[2]]
		reftag = names(tagnf)[1]
		testtag = names(tagnf)[2]
		bgfulltag = sprintf("bg(full|%s)", reftag)
		bgreftag = sprintf("bg(%s|%s)", reftag, testtag)

		refclade = read.tree(nfrefclade)
		testclade = read.tree(nftestclade)

		refclade.cols = laneid2col(refclade[['tip.label']])
		testclade.cols = laneid2col(testclade[['tip.label']])

		refclade.cols = intersect(sapply(refclade$tip.label, sanglane2assid, outidset=colpanaroo2isol), colpanaroo2isol)
		testclade.cols = intersect(sapply(testclade$tip.label, sanglane2assid, outidset=colpanaroo2isol), colpanaroo2isol)

		# ref clade (may include test clade) vs background (rest of the tree)
		bgfull.cols = setdiff(fulltree.cols, refclade.cols)
		refclade.accgene.mat = panamat[,refclade.cols]
		bgfull.accgene.mat = panamat[, bgfull.cols]
		accgene.freq.bgfull = apply(bgfull.accgene.mat, 1, sum) / length(bgfull.cols)
		accgene.freq.ref = apply(refclade.accgene.mat, 1, sum) / length(refclade.cols)
		uniquevar.bgfullvsref = (accgene.freq.bgfull>0.8 & accgene.freq.ref<0.2)
		uniquevar.refvsbgfull = (accgene.freq.ref>0.8 & accgene.freq.bgfull<0.2)
		accgenecontrast.refvsbgfull = uniquevar.refvsbgfull | uniquevar.bgfullvsref

		# test clade vs background of ref clade (rest of the tree)
		bgref.cols = setdiff(refclade.cols, testclade.cols)
		print('bgref.cols')
		print(bgref.cols)
		testclade.accgene.mat = panamat[,testclade.cols]
		bgref.accgene.mat = panamat[, bgref.cols]
		accgene.freq.bgref = apply(bgref.accgene.mat, 1, sum) / length(bgref.cols)
		accgene.freq.test = apply(testclade.accgene.mat, 1, sum) / length(testclade.cols)
		uniquevar.bgrefvstest = (accgene.freq.bgref>0.8 & accgene.freq.test<0.2)
		uniquevar.testvsbgref = (accgene.freq.test>0.8 & accgene.freq.bgref<0.2)
		accgenecontrast.testvsbgref = uniquevar.bgrefvstest | uniquevar.testvsbgref

		newaccgenecontrasts = list(accgenecontrast.refvsbgfull, accgenecontrast.testvsbgref)
		names(newaccgenecontrasts) = c(sprintf("%s.vs.full", reftag), sprintf("%s.vs.%s", testtag, reftag))
		accgenecontrasts = c(accgenecontrasts, newaccgenecontrasts)

		newaccgenefreqvec = list(accgene.freq.ref, accgene.freq.bgfull, accgene.freq.test, accgene.freq.bgref)
		names(newaccgenefreqvec) = paste("accgene.freq", c(reftag, bgfulltag, testtag, bgreftag), sep='.')
		accgenefreqvec = c(accgenefreqvec, newaccgenefreqvec)

		dataset.sizes = c(length(fulltree.cols), length(refclade.cols), length(bgfull.cols), length(testclade.cols), length(bgref.cols))
		names(dataset.sizes) = c("full.tree", reftag, bgfulltag, testtag, bgreftag)
		print("dataset sizes:")
		print(dataset.sizes)
	}
	
	testgenes = pananames %in% c('group_3132','group_7344','group_2300')
	print(lapply(accgenefreqvec, function(agfv){ agfv[testgenes] }))
	
	for (accgenecontrasttag in names(accgenecontrasts)){
	  fixed.accgenes = accgenecontrasts[[accgenecontrasttag]]
	  fixed.accgenes.info = cbind(as.character(pananames[fixed.accgenes]),
								  format(do.call(cbind, lapply(accgenefreqvec, function(agfv){ agfv[fixed.accgenes] })), digits=5))
	  if (!is.null(panannot)){ 
		  fixed.accgenes.info = cbind(fixed.accgenes.info, panannot[fixed.accgenes,])
		  colnames(fixed.accgenes.info) = c("Gene", names(accgenefreqvec), colnames(panannot)) 
	  }else{
	  	 colnames(fixed.accgenes.info) = c("Gene", names(accgenefreqvec))  
	  }
	  write.table(fixed.accgenes.info, file=paste0(dirpanaroo, sprintf("%s.fixed.panaroo_%s.tsv", accgenecontrasttag, npanatab)),
				  sep='\t', row.names=F, quote=FALSE)
	}
}