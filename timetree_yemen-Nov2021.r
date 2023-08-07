#!/usr/bin/env Rscript

library(ape)
library(BactDating)
library(ade4)

translategenomeids2strainnames = function(genomeids, specimen2lane){
	strainnames = sapply(genomeids, function(gid){ 
	  if (gid %in% specimen2lane$Lane.id){
		return(as.character(specimen2lane$Specimen.id[specimen2lane$Lane.id==gid]))
	  }else{ return(gid) }
	})
	return(strainnames)
}
root_to_tip = function(tre){
	N = Ntip(tre)
	root_node = N + 1
	r2td = dist.nodes(tre)[1:N, root_node]
	names(r2td) = tre$tip.label
	return(r2td)
}

yemen = Sys.getenv('yemen')
#nftreeyemen = file.path(yemen, 'read_mapping', 'consensus_ali', 'hybridref-Nov2021', 'hybridref-Nov2021_all_mapped.snp.aln.rba.raxml.support.subcladeH9.rooted.monoH9gmonoH9c')
nftreeyemen = file.path(yemen, 'read_mapping', 'consensus_ali', 'hybridref-Nov2021', 'hybridref-Nov2021_all_mapped.clonalframeml_out.labelled_tree.newick')
treeyemen = ladderize(read.tree(nftreeyemen))

specimen2lane = read.table(file.path(yemen, 'read_mapping', 'specimeid2laneid_SangerPasteur.txt'), sep='\t', header=T, comment.char='')
treeyemen$tip.label = translategenomeids2strainnames(treeyemen$tip.label, specimen2lane)

nfmetayemen = file.path(yemen, "metadata/hybridref-Nov2021.subcladeH9.unified_metadata_yemen2018-2019_extended.csv")
metayemen = read.csv(nfmetayemen)
metayemen[['Isolate.name']] = translategenomeids2strainnames(metayemen[['Isolate.name']], specimen2lane)

i.samples.yemenonly = metayemen[['Country']]=='Yemen' & !is.na(metayemen[['continuous.date']])
yemenonly.contdates = metayemen[i.samples.yemenonly, c('Isolate.name', 'continuous.date')]
rownames(yemenonly.contdates) = as.character(yemenonly.contdates[['Isolate.name']])
treeyemen.yemenonly = ladderize(drop.tip(treeyemen, treeyemen$tip.label[!(treeyemen$tip.label %in% yemenonly.contdates[['Isolate.name']])]))

i.samples.yemen.withcontdate = !is.na(metayemen[['continuous.date']])
yemen.contdates = metayemen[i.samples.yemen.withcontdate, c('Isolate.name', 'continuous.date')]
rownames(yemen.contdates) = as.character(yemen.contdates[['Isolate.name']])
treeyemen.contdates = ladderize(drop.tip(treeyemen, treeyemen$tip.label[!(treeyemen$tip.label %in% yemen.contdates[['Isolate.name']])]))

write.tree(treeyemen.contdates, file=paste0(nftreeyemen, '.withcontdates.nwk'))
treeyemen.contdates.root2tip = root_to_tip(treeyemen.contdates)
cor.test(treeyemen.contdates.root2tip[rownames(yemen.contdates)], yemen.contdates[,2])
pdf(file=paste0(nftreeyemen, '.withcontdates.root2tip_vs_dates.pdf'))
plot(treeyemen.contdates.root2tip[rownames(yemen.contdates)] ~ yemen.contdates[,2])
abline(lm(treeyemen.contdates.root2tip[rownames(yemen.contdates)] ~ yemen.contdates[,2]), col='red')
dev.off()
write.tree(treeyemen.yemenonly, file=paste0(nftreeyemen, '.yemenonly.nwk'))
treeyemen.yemenonly.root2tip = root_to_tip(treeyemen.yemenonly)
cor.test(treeyemen.yemenonly.root2tip[rownames(yemenonly.contdates)], yemenonly.contdates[,2])
pdf(file=paste0(nftreeyemen, '.yemenonly.root2tip_vs_dates.pdf'))
plot(treeyemen.yemenonly.root2tip[rownames(yemenonly.contdates)] ~ yemenonly.contdates[,2])
abline(lm(treeyemen.yemenonly.root2tip[rownames(yemenonly.contdates)] ~ yemenonly.contdates[,2]), col='red')
dev.off()

samples.withcontdate = as.character(metayemen[i.samples.yemen.withcontdate,'Isolate.name'])
time.dist.withcontdate = dist(metayemen[i.samples.yemen.withcontdate, 'continuous.date'])
patristic.dist.all = cophenetic(treeyemen.contdates)
patristic.dist.withcontdate = as.dist(as.matrix(patristic.dist.all)[samples.withcontdate, samples.withcontdate])
				
print("  Test correlation of phylogenetic vs. time distances (Pearson's correlation; all samples with fine time data)", quote=F)
ctpt2 = cor.test(time.dist.withcontdate, patristic.dist.withcontdate)
print(ctpt2)
#mtpt2 = mantel.randtest(time.dist.withcontdate, patristic.dist.withcontdate, n=100000)
#print(mtpt2)

# for BactDating
nsitesfullaln = 4163920
treeyemen.contdates.scaled = treeyemen.contdates
treeyemen.contdates.scaled$edge.length = treeyemen.contdates.scaled$edge.length * nsitesfullaln
#treeyemen.contdates.roottotip = roottotip(treeyemen.contdates.scaled, yemen.contdates[treeyemen.contdates.scaled$tip.label,2])
#treeyemen.contdates.clusteredTest = clusteredTest(treeyemen.contdates.scaled, yemen.contdates[treeyemen.contdates.scaled$tip.label,2])
seeds = c(1234567, 8901234, 5678901)
for (i in 1:3){
	treeyemen.contdates.dated = bactdate(treeyemen.contdates.scaled, yemen.contdates[treeyemen.contdates.scaled$tip.label,2], nbIts=100000, showProgress=T)
	nfradbactdateout.contdates = paste0(nftreeyemen, sprintf('.withcontdates.bactdating.%d', i))
	save(treeyemen.contdates.dated, file=paste0(nfradbactdateout.contdates, '.RData'))
	pdf(file=paste0(nfradbactdateout.contdates, '.pdf'), width=30, height=30)
	#pdf(file=paste0(nftreeyemen, '.withcontdates.bactdating.pdf'), width=20, height=30)
	a=2005 ; b=2020
	par(xaxp=c(a, b, (b-a)*4))
	plot(treeyemen.contdates.dated, 'treeCI')
	plot(treeyemen.contdates.dated, 'treeRoot')
	plot(treeyemen.contdates.dated, 'trace')
	#plot(treeyemen.contdates.roottotip)
	dev.off()
}

treeyemen.yemenonly.scaled = treeyemen.yemenonly
treeyemen.yemenonly.scaled$edge.length = treeyemen.yemenonly.scaled$edge.length * nsitesfullaln
#treeyemen.yemenonly.roottotip = roottotip(treeyemen.yemenonly.scaled, yemenonly.contdates[treeyemen.yemenonly.scaled$tip.label,2])
#treeyemen.yemenonly.clusteredTest = clusteredTest(treeyemen.yemenonly.scaled, yemenonly.contdates[treeyemen.yemenonly.scaled$tip.label,2])
treeyemen.yemenonly.dated = bactdate(treeyemen.yemenonly.scaled, yemenonly.contdates[treeyemen.yemenonly.scaled$tip.label,2], nbIts=100000, showProgress=T)
nfradbactdateout.yemenonly = paste0(nftreeyemen, '.yemenonly.bactdating')
save(treeyemen.yemenonly.dated, file=paste0(nfradbactdateout.yemenonly, '.RData'))
pdf(file=paste0(nfradbactdateout.yemenonly, '.pdf'), width=30, height=30)
#pdf(file=paste0(nftreeyemen, '.yemenonly.bactdating.pdf'), width=20, height=30)
plot(treeyemen.yemenonly.dated, 'treeCI')
plot(treeyemen.yemenonly.dated, 'treeRoot')
plot(treeyemen.yemenonly.dated, 'trace')
#plot(treeyemen.yemenonly.roottotip)
dev.off()


treecladeD = read.tree('~/yemen2019/read_mapping/consensus_ali/hybridref_mapped.consensus.chr1and2.snp.aln.rba.raxml.support.cladeD.nolowcoverage.assembly_ids')
