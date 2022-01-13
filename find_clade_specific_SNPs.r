#!/usr/bin/env Rscript
library(ape)
options(width=160)

laneid2col = function(laneids){
	gsub("^([23]{2}[0-9]{3}_[0-9])#", "X\\1.", laneids)
}

lenchr1 = 3046494
reflaneid = '33224_2#241'

cargs = commandArgs(trailingOnly=T)
if (length(cargs)<3){
	print("Usage: find_new_isolate_specific_SNPs.r VCFfile full.tree.file [refclade1.tag=]refclade1.file,[testclade1.tag=]testclade1.file [[refclade2.tag=]refclade2.file,[testclade2.tag=]testclade2.file [, ...]]", quote=F)
}
nfvcf = cargs[1]
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

snps = read.table(nfvcf, comment.char='', sep='\t', header=T, skip=3)
if (grepl('chr2', nfvcf)){ 
	snps[,1] = 2
	pos.for.artemis = snps[,2] + lenchr1
	vcf.art = read.table(nfvcf, comment.char='', sep='\t', header=F, skip=4)
	vcf.art[,2] = vcf.art[,2] + lenchr1
	header = readLines(nfvcf, n=4)
	lenchr2 = as.numeric(sub('##contig=<ID=[0-9],length=([0-9]+)>', '\\1', header[2]))
	sumlenchr = lenchr2 + lenchr1
	header[2] = sprintf('##contig=<ID=2,length=%d>', sumlenchr)
	write(header, file=sub('.vcf', '.for.artemis.vcf', nfvcf, fixed=TRUE))
	write.table(vcf.art, file=paste0(nfvcf, '.for.artemis'),
			sep='\t', col.names=F, row.names=F, quote=FALSE, append=TRUE)
}else{
	pos.for.artemis = snps[,2]
	snps[snps[,2]>lenchr1,1] = 2
	snps[,2] = ifelse(snps[,2]>lenchr1, (snps[,2] - lenchr1), snps[,2])
}

snpcontrasts = list()
snp.freq.full = apply(snps[,fulltree.cols], 1, sum) / length(fulltree.cols)
snpfreqvec = list(snp.freq.full.tree=snp.freq.full)

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

# ref clade (may include test clade) vs background (rest of the tree)
bgfull.cols = setdiff(fulltree.cols, refclade.cols)
refclade.snp.mat = snps[,refclade.cols]
bgfull.snp.mat = snps[, bgfull.cols]
snp.freq.bgfull = apply(bgfull.snp.mat, 1, sum) / length(bgfull.cols)
snp.freq.ref = apply(refclade.snp.mat, 1, sum) / length(refclade.cols)
uniquevar.bgfullvsref = (snp.freq.bgfull>0.8 & snp.freq.ref<0.2)
uniquevar.refvsbgfull = (snp.freq.ref>0.8 & snp.freq.bgfull<0.2)
snpcontrast.refvsbgfull = uniquevar.refvsbgfull | uniquevar.bgfullvsref
#var.in.refclade = apply(refclade.snp.mat, 1, sum)>0
#refclade.filt.snp.mat = refclade.snp.mat[var.in.refclade | snpcontrast.refvsbgfull,]
#bgfull.filt.snp.mat = bgfull.snp.mat[var.in.refclade | snpcontrast.refvsbgfull,]

# test clade vs background of ref clade (rest of the tree)
bgref.cols = setdiff(refclade.cols, testclade.cols)
testclade.snp.mat = snps[,testclade.cols]
bgref.snp.mat = snps[, bgref.cols]
snp.freq.bgref = apply(bgref.snp.mat, 1, sum) / length(bgref.cols)
snp.freq.test = apply(testclade.snp.mat, 1, sum) / length(testclade.cols)
uniquevar.bgrefvstest = (snp.freq.bgref>0.8 & snp.freq.test<0.2)
uniquevar.testvsbgref = (snp.freq.test>0.8 & snp.freq.bgref<0.2)
snpcontrast.testvsbgref = uniquevar.testvsbgref | uniquevar.bgrefvstest
#var.in.testclade = apply(testclade.snp.mat, 1, sum)>0
#tesstclade.filt.snp.mat = testclade.snp.mat[var.in.testclade | snpcontrast.testvsbgref,]
#bgref.filt.snp.mat = bgref.snp.mat[var.in.testclade | snpcontrast.testvsbgref,]
newsnpcontrasts = list(snpcontrast.refvsbgfull, snpcontrast.testvsbgref)
names(newsnpcontrasts) = c(sprintf("%s.vs.full", reftag), sprintf("%s.vs.%s", testtag, reftag))
snpcontrasts = c(snpcontrasts, newsnpcontrasts)

newsnpfreqvec = list(snp.freq.ref, snp.freq.bgfull, snp.freq.test, snp.freq.bgref)
names(newsnpfreqvec) = paste("snp.freq", c(reftag, bgfulltag, testtag, bgreftag), sep='.')
snpfreqvec = c(snpfreqvec, newsnpfreqvec)

dataset.sizes = c(length(fulltree.cols), length(refclade.cols), length(bgfull.cols), length(testclade.cols), length(bgref.cols))
names(dataset.sizes) = c("full.tree", reftag, bgfulltag, testtag, bgreftag)
print("dataset sizes:")
print(dataset.sizes)
			   
}

for (snpcontrasttag in names(snpcontrasts)){
  fixed.snps = snpcontrasts[[snpcontrasttag]]
  fixed.snps.i = as.numeric(rownames(snps)[fixed.snps])
  fixed.snps.info = cbind(snps[fixed.snps.i,1:6], 
						  pos.for.artemis[fixed.snps.i],
						  do.call(cbind, lapply(snpfreqvec, function(sfv){ sfv[fixed.snps] })))

  colnames(fixed.snps.info) = c(colnames(fixed.snps.info)[1:6], 'pos.for.artemis', names(snpfreqvec)) 
  write.table(format(fixed.snps.info, digits=5), file=paste0(nfvcf, sprintf(".%s.fixed.snps.tab", snpcontrasttag)),
			  sep='\t', row.names=F, quote=FALSE)
}
