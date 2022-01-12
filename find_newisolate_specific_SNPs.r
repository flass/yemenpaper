#!/usr/bin/env Rscript
library(ape)
options(width=160)

laneid2col = function(laneids){
	gsub("^([23]{2}[0-9]{3}_[0-9])#", "X\\1.", laneids)
}

lenchr1 = 3046494
reflaneid = '33224_2#241'

cargs = commandArgs(trailingOnly=T)
nfvcf = cargs[1]
nffulltree = cargs[2]
nfmainclade = cargs[3]
nfclade2018 = cargs[4]

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
fulltree = read.tree(nffulltree)
mainclade = read.tree(nfmainclade)
clade2018 = read.tree(nfclade2018)

fulltree.cols = laneid2col(fulltree[['tip.label']])
mainclade.cols = laneid2col(mainclade[['tip.label']]) # includes ref
clade2018.cols = laneid2col(clade2018[['tip.label']])

# main (main yemen cholera clone) vs bg (Yemen divergent isolates)
bgclade.cols = setdiff(fulltree.cols, mainclade.cols)
mainclade.snp.mat = snps[,mainclade.cols]
bgclade.snp.mat = snps[, bgclade.cols]
snp.freq.bg = apply(bgclade.snp.mat, 1, sum) / length(bgclade.cols)
snp.freq.main = apply(mainclade.snp.mat, 1, sum) / length(mainclade.cols)
uniquevar.bgvsmainclade = (snp.freq.bg>0.8 & snp.freq.main<0.2)
uniquevar.mainvsbgclade = (snp.freq.main>0.8 & snp.freq.bg<0.2)
contrast.mainvsbgclade = uniquevar.mainvsbgclade | uniquevar.bgvsmainclade
#var.in.mainclade = apply(mainclade.snp.mat, 1, sum)>0
#mainclade.filt.snp.mat = mainclade.snp.mat[var.in.mainclade | contrast.mainvsbgclade,]
#bgclade.filt.snp.mat = bgclade.snp.mat[var.in.mainclade | contrast.mainvsbgclade,]

# previous (clade with mostly 2016-2017 Yemen isolates) vs. new (clade with mostly 2018-2019 Yemen isolates)
prevclade.samples = grepl("^CNRVC", mainclade.cols) | (mainclade.cols %in% c("X33224_2.7", "X33224_2.33"))
prevclade.cols = mainclade.cols[prevclade.samples]
newclade.cols = mainclade.cols[!prevclade.samples] # includes ref
print("prevclade.cols")
print(prevclade.cols)
print("newclade.cols")
print(newclade.cols)

newclade.snp.mat = snps[,newclade.cols]
prevclade.snp.mat = snps[, prevclade.cols]
snp.freq.prev = apply(prevclade.snp.mat, 1, sum) / length(prevclade.cols)
snp.freq.new = apply(newclade.snp.mat, 1, sum) / length(newclade.cols)
uniquevar.prevvsnewclade = (snp.freq.prev>0.8 & snp.freq.new<0.2)
uniquevar.newvsprevclade = (snp.freq.new>0.8 & snp.freq.prev<0.2)
contrast.newvsprevclade = uniquevar.newvsprevclade | uniquevar.prevvsnewclade

## 2018 vs the bg of main
#clade2018.snp.mat = snps[,clade2018.cols]
#bgmainclade.snp.mat = snps[, setdiff(mainclade.cols, clade2018.cols)]
#unique.2018vsbgmainclade = (apply(clade2018.snp.mat, 1, sum)>0.8 & apply(bgmainclade.snp.mat, 1, sum)<0.2)
#unique.bgmainvs2018clade = (apply(bgmainclade.snp.mat, 1, sum)>0.8 & apply(clade2018.snp.mat, 1, sum)<0.2)
#contrast.mainvsbgclade = unique.2018vsbgmainclade | unique.bgmainvs2018clade
#clade2018.filt.snp.mat = clade2018.snp.mat[var.in.mainclade | contrast.mainvsbgclade,]
#bgmainclade.filt.snp.mat = bgmainclade.snp.mat[var.in.mainclade | contrast.mainvsbgclade,]

# 2018 (clade with only and almost all 2018 isolates) vs the bg of new
print("clade2018.cols")
print(clade2018.cols)
samples.2018 = mainclade.cols %in% clade2018.cols
bgnewclade.cols = setdiff(newclade.cols, clade2018.cols) # includes ref
print("bgnewclade.cols")
print(bgnewclade.cols)
clade2018.snp.mat = snps[,clade2018.cols]
bgnewclade.snp.mat = snps[, bgnewclade.cols]
snp.freq.bgnew = apply(bgnewclade.snp.mat, 1, sum) / length(bgnewclade.cols)
snp.freq.2018 = apply(clade2018.snp.mat, 1, sum) / length(clade2018.cols)
uniquevar.2018vsbgnewclade = (snp.freq.2018>0.8 & snp.freq.bgnew<0.2)
uniquevar.bgnewvs2018clade = (snp.freq.bgnew>0.8 & snp.freq.2018<0.2)
contrast.2018vsbgnewclade = uniquevar.2018vsbgnewclade | uniquevar.bgnewvs2018clade
#clade2018.filt.snp.mat = clade2018.snp.mat[var.in.mainclade | contrast.2018vsbgnewclade,]
#bgnewclade.filt.snp.mat = bgnewclade.snp.mat[var.in.mainclade | contrast.2018vsbgnewclade,]

dataset.lengths = c(length(fulltree.cols), length(mainclade.cols), length(bgclade.cols), length(prevclade.cols), length(newclade.cols), length(bgnewclade.cols), length(clade2018.cols))
names(dataset.lengths) = c("full.tree", "main.Yemen.clone.clade", "background.divergent.clade", "mostly2016-2017.sample.clade", "mostly2018-2019.sample.clade", "mostly2019.sample.clade", "2018.sample.clade")
print("dataset.lengths")
print(dataset.lengths)

#fixed.snps = which(snp.freq.prev>0.8 | snp.freq.new>0.8 | snp.freq.2018>0.8 | snp.freq.bgnew>0.8)
combined.contrasts = which(contrast.mainvsbgclade | contrast.newvsprevclade | contrast.2018vsbgnewclade)

constrasts = list(bgVSmain=contrast.mainvsbgclade,
				  prevVSnew=contrast.newvsprevclade,
				  bgnewVS2018=contrast.2018vsbgnewclade,
				  allcontrasts=combined.contrasts)
for (contrasttag in names(constrasts)){
  fixed.snps = constrasts[[contrasttag]]
  fixed.snps.i = as.numeric(rownames(snps)[fixed.snps])
  fixed.snps.info = cbind(snps[fixed.snps.i,c(1:6)], pos.for.artemis[fixed.snps.i],
						  snp.freq.bg[fixed.snps], snp.freq.main[fixed.snps],
						  snp.freq.prev[fixed.snps], snp.freq.new[fixed.snps],
						  snp.freq.bgnew[fixed.snps], snp.freq.2018[fixed.snps])

  colnames(fixed.snps.info) = gsub("\\[fixed\\.snps.*\\]", ".clade", colnames(fixed.snps.info)) 
  write.table(format(fixed.snps.info, digits=5), file=paste0(nfvcf, sprintf(".%s.fixed.snps.tab", contrasttag)),
			  sep='\t', row.names=F, quote=FALSE)
}