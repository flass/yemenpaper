#!/usr/bin/env Rscript
library(phytools)
library(ggplot2)
library(ggtree)
library(ggnewscale)
library(RColorBrewer)
library(ggmap)
library(sf)
library(ade4)

microbiomerepo=file.path(Sys.getenv('HOME'), 'software', 'microbiomes')
source(file.path(microbiomerepo, 'scripts', 'gpscoords.r'))

#source(file.path(Sys.getenv('HOME'), 'scripts', 'vibrio', 'plotting_snp_sites_MGEs.r'))
source(file.path(Sys.getenv('HOME'), 'scripts', 'vibrio', 'plot_MGE_from_blast_results.r'))
source(file.path(Sys.getenv('HOME'), 'scripts', 'vibrio', 'plot_hamburger_Oantigen.r'))
source(file.path(Sys.getenv('HOME'), 'scripts', 'UNsubregions_noOceania.r'))
source(file.path(Sys.getenv('HOME'), 'scripts', 'dataframe2iTol_dataset.r'))

all.clades = c(LETTERS[1:11], paste('H', 1:9, sep='.'), paste('H.9', letters[1:8], sep='.'))
allyearmonth1619 = unlist(lapply(2016:2019, function(y){ paste(y, sprintf("%02d", 1:12), sep='-') }))

options(width=160)
annotlist2mat = function(pab, isol=NULL){
	annots = sort(unique(as.character(pab[,3])))
	if (is.null(isol)){ isolates = sort(unique(as.character(pab[,1])))
	}else{ isolates = sort(unique(isol)) }
	annot.abspres = as.data.frame(sapply(annots, function(annot){
		abr = unique(as.character(pab[pab[,3]==annot,1]))
		abrap = as.character(sapply(isolates, function(isolate){
			ifelse((isolate %in% abr), annot, "")
		}))
		return(factor(abrap, levels=c("", annot)))
	}))
	colnames(annot.abspres) = annots
	rownames(annot.abspres) = isolates
	return(annot.abspres)
}

split.date2DMY = function(sdmy){
	dmy = t(simplify2array(lapply(as.character(sdmy), function(s){
		if (!is.na(s)){ return(strsplit(s, split='/')[[1]])
		}else{ return(rep(NA, 3)) }
	})))
	return(dmy)
}
split.Date2DMY = function(dateymd){
	dmy = t(simplify2array(lapply(as.character(dateymd), function(s){
		if (!is.na(s)){ return(rev(strsplit(s, split='-')[[1]]))
		}else{ return(rep(NA, 3)) }
	})))
	return(dmy)
}

months = c('january', 'february', 'march', 'april', 'may', 'june', 'july', 'august', 'september', 'october', 'november', 'december')
months.lengths = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
months.cumlengths = cumsum(months.lengths)

date.DMY2continous = function(dmy){
	DMY = apply(dmy, 1, function(dmyst){
		dmyt = as.numeric(dmyst)
#		ct = (as.double(dmyt[1] - 1)/365) + (as.double(dmyt[2] - 1)/12) + as.double(dmyt[3])
		ct = (as.double(months.cumlengths[dmyt[2] - 1] + dmyt[1] - 1)/365) + as.double(dmyt[3])
#		print(c(dmyt, "=", c((as.double(dmyt[1] - 1)/31), (as.double(dmyt[2] - 1)/12), as.double(dmyt[3])), "=", ct))
		return(ct)
	})
	return(DMY)
}

yemen = Sys.getenv('yemen')
#yemen = '/Users/fl4/yemen2019'
#ptgdbname = 'yemen2019newcontext'
ptgdbname = 'contextallyemen2019'
pantayemencore = file.path(yemen, 'panta', ptgdbname, '05.core_genome')
newcontext = Sys.getenv('newcontext')
dorman2020 = Sys.getenv('dorman2020')
hybconsrad = Sys.getenv('hybconsrad')
n16961 = file.path(yemen, 'N16961v2_mapped_consensus')
#cladetreetags = c('full.nocontam', 'mainclade.nocontam.nodivergentchineseisol')
cladetreetags = c('full.nocontam', 'mainclade.nocontam')
#treetags = c(cladetreetags, paste(cladetreetags, 'lsd.nwk', sep='.'), paste('assembly_ids', cladetreetags, 'tre', sep='.'))
treetags = paste('assembly_ids', cladetreetags, 'tre', sep='.')

datatags = lapply(treetags, function(treetag){ paste(ptgdbname, treetag, sep='.') })
nftrees = lapply(treetags, function(treetag){ file.path(pantayemencore, sprintf('core-genome-based_reference_tree_%s.%s', ptgdbname, treetag)) })
names(nftrees) = datatags
nftrees[['hybridref_mapped.consensus.full']] = paste0(hybconsrad, '.chr1and2.snp.aln.rba.raxml.support.full.nolowcoverage')
nftrees[['hybridref_mapped.consensus.main']] = paste0(hybconsrad, '.chr1and2.snp.aln.rba.raxml.support.mainclade.nolowcoverage')
nftrees[['hybridref_mapped.consensus.narrow']] = paste0(hybconsrad, '.chr1and2.snp.aln.rba.raxml.support.narrowclade.nolowcoverage')
nftrees[['hybridref_mapped.consensus.main.norecomb']] = paste0(hybconsrad, '.chr1and2.final_tree.mainclade.nolowcoverage.tre')
nftrees[['N16961v2_mapped.consensus.main']] = file.path(n16961, 'N16961v2_mapped_pseudo_genome_concat.snp.aln.rba.raxml.support.main.tre')
nftrees[['hybridref-Nov2021']] = file.path(dirname(dirname(hybconsrad)), 'hybridref-Nov2021', 'hybridref-Nov2021_all_mapped.snp.aln.rba.raxml.support.full.rooted.H9hsisterH9guniteH9c.no225')
nftrees[['hybridref-Nov2021.subcladeH689']] = file.path(dirname(dirname(hybconsrad)), 'hybridref-Nov2021', 'hybridref-Nov2021_all_mapped.snp.aln.rba.raxml.support.subcladeH689.rooted.H9hsisterH9guniteH9c.no225')
nftrees[['hybridref-Nov2021.subcladeH9']] = file.path(dirname(dirname(hybconsrad)), 'hybridref-Nov2021', 'hybridref-Nov2021_all_mapped.snp.aln.rba.raxml.support.subcladeH9.rooted.monoH9gmonoH9c')
nftrees[['hybridref-Nov2021.gubbins']] = file.path(dirname(dirname(hybconsrad)), 'hybridref-Nov2021', 'hybridref-Nov2021_all_mapped.final_tree.tre.no225')
nftrees[['ST555_hybrid']] = file.path(dirname(dirname(hybconsrad)), 'ST555_hybrid', 'ST555_hybrid_mapped.snp.aln.rba.raxml.support.full')

lbnftrees = sapply(nftrees, basename)
# to find the same subclade tree files
lbnftrees[['N16961v2_mapped.consensus.main']] = lbnftrees[['hybridref_mapped.consensus.narrow']] = lbnftrees[['hybridref_mapped.consensus.main.norecomb']] = lbnftrees[['hybridref_mapped.consensus.main']] = lbnftrees[['hybridref_mapped.consensus.full']] = basename(paste0(hybconsrad, '.chr1and2.snp.aln.rba.raxml.support.full'))
lbnftrees[['hybridref-Nov2021.gubbins']] = lbnftrees[['hybridref-Nov2021.subcladeH689']] = lbnftrees[['hybridref-Nov2021']]

# metadata
metaweill = read.table(file.path(yemen, 'metadata/Weill_Nature2018-Yemen_TableS2.txt'), comment.char='', quote='', sep='\t', header=T, strip.white=T)
colnames(metaweill)[colnames(metaweill)=='Source'] = 'Data.source'
colnames(metaweill)[colnames(metaweill)=='Source.1'] = 'Isolate.source'
prevstudies = ifelse(metaweill$Data.source=='This study', 'Weill et al. Nature 2018', as.character(metaweill[['Data.source']]))
metaweill[['Data.source']] = factor(prevstudies, levels=sort(c(unique(prevstudies), 'This study')))
metaweill[['EBI.ENA.accession.number']] = sub("^_", "", metaweill[['EBI.ENA.accession.number']])

#metaghulam = read.table(file.path(yemen, 'metadata/Yemen2018-2019_sequenced_sample_metadata.csv'), comment.char='', quote='', sep=',', header=T)
metaghulam = read.table(file.path(yemen, 'metadata/Yemen2018-2019_sequenced_sample_metadata-cor2.csv'), comment.char='', quote='', sep=',', header=T)
metaghulam[metaghulam==''] = NA
metaghulamPI = metaghulam
metaghulamPI[['specimen.id']] = paste0(as.character(metaghulamPI[['specimen.id']]), '-PI')
metaghulam = rbind(metaghulam, metaghulamPI)
print("head(metaghulam)")
print(head(metaghulam))
print("tail(metaghulam)")
print(tail(metaghulam))

seqid2name = read.table(file.path(yemen, 'metadata/yemen2019_SangerPasteurSampleId2SupplierNames.txt'), comment.char='', quote='', sep='\t', header=T)
#sangerids = read.table(file.path(yemen, 'metadata/Yemen_study_pf-info.csv'), comment.char='', quote='', sep=',', header=T)[,1:2]
sangerids = read.table(file.path(yemen, 'metadata/Yemen_and_context_pf-info.sampleid2laneid.csv'), comment.char='', quote='', sep=',', header=T)
enaaccessions = read.table(file.path(yemen, 'STDY5960_biosample_accessions.txt'), comment.char='', quote='', header=T)[,1:3]
colnames(enaaccessions) = c('Sample', 'EBI.ENA.accession.number', 'Lane')
sangerids = merge(sangerids, enaaccessions, all=TRUE)

# put clean strain names instead of genome assembly/lane ids
specimen2lane = read.table(file.path(yemen, 'read_mapping', 'specimeid2laneid_SangerPasteur.txt'), sep='\t', header=T, comment.char='')
translategenomeids2strainnames = function(genomeids){
	strainnames = sapply(genomeids, function(gid){ 
	  if (gid %in% specimen2lane$Lane.id){
		return(as.character(specimen2lane$Specimen.id[specimen2lane$Lane.id==gid]))
	  }else{ return(gid) }
	})
	return(strainnames)
}
#assid2code = read.table(file.path(dorman2020, sprintf('strain_infos_%s.txt', ptgdbname)), header=T, sep='\t', comment.char='')
assid2code = read.table(file.path(newcontext, 'lastcontext', sprintf('strain_infos_%s.txt', ptgdbname)), header=T, sep='\t', comment.char='')

INSDC.acc.col = paste(c('SRA', 'Assembly', 'BioProject', 'BioSample', 'WGS', 'other'), 'accession', sep='.')

parseAccession = function(indf, accol){
	indf[['SRA.accession']] = ifelse(grepl('^[SED]R[RSX]', indf[[accol]]) , indf[[accol]], NA)
	indf[['Assembly.accession']] = ifelse(grepl('^GC[AF]_', indf[[accol]]) , indf[[accol]], NA)
	indf[['BioProject.accession']] = ifelse(grepl('^PRJ', indf[[accol]]) , indf[[accol]], NA)
	indf[['BioSample.accession']] = ifelse(grepl('^SAM', indf[[accol]]) , indf[[accol]], NA)
	indf[['WGS.accession']] = ifelse(grepl('000000$', indf[[accol]]), indf[[accol]], NA)
	all.na = apply(indf[,INSDC.acc.col[grep('other', INSDC.acc.col, inv=T)]], 1, function(x){ all(is.na(x))})
	indf[['other.accession']] = ifelse(all.na, indf[[accol]], NA)
	return(indf)
}

metadorman = read.table(paste0(dorman2020, '.txt'), header=T, sep='\t', comment.char='')
metadorman = parseAccession(metadorman, 'accession')

metawang = read.table(file.path(newcontext, 'Wang_PLOS-NTD_2019_tableS1.txt'), header=T, sep='\t', comment.char='')
metawang = parseAccession(metawang, 'Accession')
colnames(metawang) = sub('Strain', 'name', colnames(metawang))

metaeastaf = read.csv(file.path(newcontext, 'East_africa_combined_info.csv'), comment.char='')
metaeastaf = merge(metaeastaf, sangerids, by.x='assembly_id', by.y='Sample', all.x=T)
colnames(metaeastaf) = sub('Serotype', 'Serotype_phenotypic', 
					   sub('Country', 'country', 
					   sub('Isolation.year', 'Year', colnames(metaeastaf))))
metaeastaf = parseAccession(metaeastaf, 'BioSample')

metamorita = read.csv(file.path(newcontext, 'Morita_mBio_2020-supp-table-S1.csv'), comment.char='')
colnames(metamorita) = sub('Country', 'country', 
					   sub('Isolation.year', 'Year', 
					   sub('Reference', 'reference', colnames(metamorita))))
metamorita = parseAccession(metamorita, 'Sanger_Lane_id')

metanew = merge(metawang, metadorman, all=T)
colnames(metanew) = sub('Sanger_lane_ID', 'assembly_id', colnames(metanew))
metanew[which(metanew$Year %in% c('no data', 'Unknown')), 'Year'] = NA
metanew[which(metanew$Source.y %in% c('no data', 'Not noted', 'Uncertain')), 'Source.y'] = NA
metanew[which(!(metanew$Serotype_phenotypic %in% c('Inaba', 'Ogawa'))), 'Serotype_phenotypic'] = NA
metanew[which(metanew[['Lineage']]==''),'Lineage'] = NA
metanew[which(metanew[['lineage_after_Domman_et_al']]=='/'),'lineage_after_Domman_et_al'] = NA

metanew = merge(metanew, metaeastaf, all=T)
metanew = merge(metanew, metamorita, all=T)

print("head(metanew)")
print(head(metanew))
print("metanew[which(metanew$assembly_id=='ERR3342505'),]")
print(metanew[which(metanew$assembly_id=='ERR3342505'),])

# where relevant, allow to translate Sanger Lane ids into assembly ids i.e. genome ids used in the rest of the study
sanglaneids = as.character(metanew[,'Sanger_Lane_id'])
sanglane2assid = function(laneorassid){
	newassid.i = which(sanglaneids==laneorassid)
	if (length(newassid.i)==1){
#                print(c("newcontext Lane match:", newassid.i, metanew[newassid.i,'assembly_id']), quote=F)
		return(metanew[newassid.i,'assembly_id'])
	}else{ return(laneorassid) }
}
# ARIBA screen of intact wbeT gene i.e. coding for a Ogawa phenotype ina O1 antigen background
aribasero = read.csv(file.path(yemen, 'serotyping_O1/ariba_wbeT/ariba_summary_wbeT.csv'), row.names=1, stringsAsFactors=TRUE)
#rownames(aribasero) = sapply(basename(dirname(rownames(aribasero))), sanglane2assid)
rownames(aribasero) = basename(dirname(rownames(aribasero)))

# Hamburger prediction of O-antigen biosynthesis gene cluster
print("summary(matOantigenclustmatches)")
print(summary(matOantigenclustmatches))
print("head(matOantigenclustmatches)")
print(head(matOantigenclustmatches))
#print(rownames(matOantigenclustmatches))

# map of yemen governorates
yem_governorates_map_polygons = st_read(file.path(yemen, 'QGIS_maps/from_QuickOSM/Yemen_admin_level_4_justfill.dxf'))
# sf data.frame
#Simple feature collection with 51 features and 6 fields
#geometry type:  GEOMETRY
#dimension:      XY
#bbox:           xmin: 41.60825 ymin: 11.90848 xmax: 54.73894 ymax: 19
#CRS:            NA
# CRS should be WGS 84 / EPSG:4326

#print("head(metanew)")
#print(head(metanew))

# Pangenome abs/pres profiles
panaroopangenomegenes = read.table(file.path(yemen, 'panaroo', ptgdbname, 'gene_presence_absence.Rtab'), comment.char='', sep='\t', header=T)
rownames(panaroopangenomegenes) = panaroopangenomegenes[,'Gene']
panaroopangenomestruct = read.table(file.path(yemen, 'panaroo', ptgdbname, 'struct_presence_absence.Rtab'), comment.char='', sep='\t', header=T)
rownames(panaroopangenomestruct) = panaroopangenomestruct[,'Gene']
stopifnot(all(colnames(panaroopangenomegenes) == colnames(panaroopangenomestruct)))
print("head(colnames(panaroopangenomegenes)[2:ncol(panaroopangenomegenes)])")
print(head(colnames(panaroopangenomegenes)[2:ncol(panaroopangenomegenes)]))
colpanaroo2isol = sub('^([0-9]{5}_[0-8])\\.([0-9])', '\\1#\\2', sub('^(.+)\\.1_\\1_genomic', '\\1', sub('^X', '', colnames(panaroopangenomegenes)[2:ncol(panaroopangenomegenes)])))
print("head(colpanaroo2isol)")
print(head(colpanaroo2isol))

genomelengths = apply(panaroopangenomegenes[2:ncol(panaroopangenomegenes)], 2, sum)
#print(genomelengths[genomelengths>5000])
# X33224_2.12 X33224_2.187 X33224_2.201 X33224_2.261  X33224_2.33  X33224_2.68
#        7456         6771         6856         7423         7445         6489
# interestingly, in accordance to their place in the core genome tree:
# #187, #201 carry no plasmid, but carry the resistance gene cassette other wise associated with IncA/C2 plasmid
# #2, #12, #33, #68, #261 carry the plasmid, with the cassette
longgenomes = colpanaroo2isol[genomelengths>5000]

# AMR/MGE predicted profiles
abricate.all = read.table(file.path(yemen, 'AMR/abricate_all.tab'), comment.char='', quote='', sep='\t')
parse.abricate = lapply(unique(abricate.all$V12), function(abdb){
	sub.abricate = abricate.all[abricate.all$V12==abdb, c(1,2,6)]
	sub.abricate$V1 = as.factor(sapply(strsplit(as.character(sub.abricate$V1), split='\\.fa'), `[`, 1))
	colnames(sub.abricate) = c("Isolate.name", paste(c("contig", "predicted"), abdb, sep='.'))
	return(unique(sub.abricate))
})
names(parse.abricate) = unique(abricate.all$V12)

cond.abricate = lapply(names(parse.abricate), function(abdb){
	pab = parse.abricate[[abdb]]
	preds = t(sapply(levels(pab[["Isolate.name"]]), function(i){
		pred = paste(sort(unique(as.character(pab[pab[["Isolate.name"]]==i, 3]))), collapse=', ')
		return(c(i, pred))
	}))
	colnames(preds) = c("Isolate.name", paste("predicted", abdb, sep='.'))
	return(preds)
})
names(cond.abricate) = names(parse.abricate)

abricate.pheno = setdiff(unlist(lapply(cond.abricate, colnames)), "Isolate.name")
#datatag = '299context_258yemen2019'
#tree.544 = read.tree(file.path(yemen, '299context_258yemen2019_pseudo_genome.snp.aln.raxml.reduced.phy.raxml.rba.raxml.bestTree'))
#outgroup.tip.i = which(tree.544$tip.label=='33224_2#255')
#tree.544.rooted = midpoint.root(tree.544)
#treeyemen = tree.544.rooted

cureatemetadata = function(genomes, det.ncbiamr, det.plasmid, det.vfdb, datatag, altjoincol=NULL){
	print("head(genomes)3")
	print(head(genomes))
	meta.yemen = merge(genomes, merge(merge(sangerids, seqid2name, by='Sample', all=T), metaghulam, by='specimen.id', all=T), by='Lane', all.x=T)
	colnames(meta.yemen)[colnames(meta.yemen)=='Lane'] = "Isolate.name"
	print("dim(meta.yemen) (a)")
	print(dim(meta.yemen))

	meta.yemen = merge(metaweill, meta.yemen, by=c("Isolate.name"), all.y=T)
	colnames(meta.yemen)[colnames(meta.yemen)=='EBI.ENA.accession.number.x'] = "EBI.ENA.accession.number"
	meta.yemen[["EBI.ENA.accession.number"]] = ifelse(is.na(meta.yemen[["EBI.ENA.accession.number"]]), meta.yemen[["EBI.ENA.accession.number.y"]], meta.yemen[["EBI.ENA.accession.number"]])
	meta.yemen = meta.yemen[,-which(colnames(meta.yemen)=='EBI.ENA.accession.number.y')]
#	meta.yemen = merge(metaweill, meta.yemen, by=c("Isolate.name", "EBI.ENA.accession.number"), all.y=T)
	print("dim(meta.yemen) (b)")
	print(dim(meta.yemen))
	
	print("summary(as.data.frame(metanew[,c('assembly_id', 'Sanger_Lane_id')]))")
	print(summary(as.data.frame(metanew[,c('assembly_id', 'Sanger_Lane_id')])))
	print("summary(meta.yemen$Isolate.name)")
	print(summary(meta.yemen[['Isolate.name']]))
	
	i.newcontext.my2mn = sapply(1:nrow(meta.yemen), function(imy){
		imn = which(as.character(metanew[['assembly_id']])==as.character(meta.yemen[['assembly_id']][imy]))
		if (length(imn)!=1){ return(NA) 
		}else{ return(imn) }
	})
    print("length(which(!is.na(i.newcontext.my2mn)))")
    print(length(which(!is.na(i.newcontext.my2mn))))
	new.yemen.i = grep('33224_2#|33816_1#', meta.yemen[['Isolate.name']])
    print("length(new.yemen.i)")
    print(length(new.yemen.i))
	pasteur.yemen.i = grep('CNRVC19', meta.yemen[['Isolate.name']])
    print("length(pasteur.yemen.i)")
    print(length(pasteur.yemen.i))
    
	print("dim(meta.yemen) (INSDC.acc.col)")
	print(dim(meta.yemen))
	meta.yemen = parseAccession(meta.yemen, 'EBI.ENA.accession.number')
	
	for (insdcac in INSDC.acc.col){
		meta.yemen[[insdcac]] = as.character(sapply(1:nrow(meta.yemen), function(imy){
			imn = i.newcontext.my2mn[imy]
			if (!is.na(imn)){
				return(metanew[imn,insdcac])
			}else{
				return(meta.yemen[imy,insdcac])
			}	
		}))
	}
	
	print("dim(meta.yemen) (c)")
	print(dim(meta.yemen))
	meta.yemen[['Country']] = as.factor(sapply(1:nrow(meta.yemen), function(imy){
		imn = i.newcontext.my2mn[imy]
		if (!is.na(imn)){
			return(metanew[imn,'country'])
		}else{
            if (imy %in% c(new.yemen.i, pasteur.yemen.i)){
                return('Yemen')
            }else{
                return(meta.yemen[imy,'Country'])
            }
		}
	}))
	meta.yemen[['UN.Subregion']] = as.factor(sapply(1:nrow(meta.yemen), function(imy){
        if (imy %in% new.yemen.i){
            return('Western Asia')
        }else{
            imn = i.newcontext.my2mn[imy]
            if (!is.na(imn)){
                sr = metanew[imn,'UN.Subregion']
            }else{
                sr = meta.yemen[imy,'UN.Subregion']
            }
            if (!is.na(sr) & sr!=''){
                return(sr)
            }else{
                sr.i = which(sapply(subregions, function(x){ meta.yemen[imy,'Country'] %in% x }))
                if (length(sr.i)==1){
                    return(names(subregions)[sr.i])
                }else{ return(NA) }
            }
		}
	}))
	meta.yemen[['Continent']] = as.factor(sapply(1:nrow(meta.yemen), function(imy){
        if (imy %in% new.yemen.i){
            return('Asia')
        }else{
            imn = i.newcontext.my2mn[imy]
            if (!is.na(imn)){
                co = metanew[imn,'Continent']
            }else{
                co = meta.yemen[imy,'Continent']
            }
            if (!is.na(co) & co!=''){
                return(co)
            }else{
                co.i = which(sapply(continents, function(x){ meta.yemen[imy,'UN.Subregion'] %in% x }))
                if (length(co.i)==1){
                    return(names(continents)[co.i])
                }else{ return(NA) }
            }
		}
	}))
	meta.yemen[['governorate.eng']] = as.factor(sapply(1:nrow(meta.yemen), function(imy){
		if (imy %in% c(new.yemen.i, pasteur.yemen.i)){
#		  if (imy %in% new.yemen.i){
            myge = meta.yemen[imy,'governorate.eng']
#          }else{ #if (imy %in% pasteur.yemen.i)
#			specid = sub('-PI', '', translategenomeids2strainnames(meta.yemen[imy,'Isolate.name']))
#			matchspecid = which(meta.yemen[['specimen.id']]==specid)
#			if (length(matchspecid)>0){
#			myge = meta.yemen[matchspecid[1],'governorate.eng']
#			}else{ return(NA) }
#		  }
		  if (myge == 'Almahwait') return('Al Mahwit')
          else return(myge)  
		}else{
            if (!is.na(meta.yemen[imy,'Country']) & meta.yemen[imy,'Country']=='Yemen'){
                locprov = meta.yemen[imy,'Locality.Province']
                locprovcit = strsplit(locprov, split=' \\(')[[1]]
                if (length(locprovcit)>1){
                    if (locprovcit[1]=='Bani al-hareth'){ return('Amanat Al Asimah')
                    }else{ return(strsplit(locprovcit[2], split=' City)')[[1]][1]) }
                }
                if (startsWith(locprov, "Sana'a")){ return("Sana'a") }
                locprovdis = strsplit(locprov, split=' District')[[1]]
                if (endsWith(locprov, " District") | length(locprovdis)>1){ return(locprovdis[1]) }
            }
        }
        return(NA)
	}))
	meta.yemen[['district']] = as.factor(sapply(1:nrow(meta.yemen), function(imy){
		if (imy %in% c(new.yemen.i, pasteur.yemen.i)){
#		  if (imy %in% new.yemen.i){
            myge = meta.yemen[imy,'district']
#          }else{ #if (imy %in% pasteur.yemen.i)
#			specid = sub('-PI', '', translategenomeids2strainnames(meta.yemen[imy,'Isolate.name']))
#			matchspecid = which(meta.yemen[['specimen.id']]==specid)
#			if (length(matchspecid)>0){
#			myge = meta.yemen[matchspecid[1],'district']
#			}else{ return(NA) }
#		  }
        }else{
            if (!is.na(meta.yemen[imy,'Country']) & meta.yemen[imy,'Country']=='Yemen'){
                locprov = meta.yemen[imy,'Locality.Province']
                locprovcit = strsplit(locprov, split=' \\(')[[1]]
                if (length(locprovcit)>1) return(locprovcit[1])
                if (startsWith(locprov, "Sana'a ")) return(strsplit("Sana'a", split=' ')[[1]][2])
                locprovdis = strsplit(locprov, split=' District')[[1]]
                if (endsWith(locprov, ' District')) return(locprovdis[1])
            }
        }
        return(NA)
	}))
	meta.yemen[['Data.source']] = as.factor(sapply(1:nrow(meta.yemen), function(imy){
		imn = i.newcontext.my2mn[imy]
		if (!is.na(imn)){
			return(as.character(metanew[imn,'reference']))
		}else{
            if (imy %in% new.yemen.i){
                return('This study')
            }else{
                return(as.character(meta.yemen[imy,'Data.source']))
            }
        }
	}))
	meta.yemen[['Introduction.event']] = as.factor(sapply(1:nrow(meta.yemen), function(imy){
		imn = i.newcontext.my2mn[imy]
		if (imy %in% new.yemen.i & !is.na(meta.yemen[imy, 'clade'])){
			if (startsWith(as.character(meta.yemen[imy, 'clade']), 'I')) return('T13')
			else return(NA) 
		}else{
			iev = as.character(meta.yemen[imy,'Introduction.event'])
			if (is.na(iev) | iev=='') return(NA)
			else return(iev)
        }
	}))
	meta.yemen[['Sample']] = as.character(sapply(1:nrow(meta.yemen), function(imy){
		imn = i.newcontext.my2mn[imy]
		if (!is.na(imn)){
			return(as.character(metanew[imn,'Sanger_Lane_id']))
		}else{
            return(as.character(meta.yemen[imy,'Sample']))
        }
	}))
	meta.yemen[['ctxB']] = as.character(sapply(1:nrow(meta.yemen), function(imy){
		imn = i.newcontext.my2mn[imy]
		if (!is.na(imn)){
			return(as.character(metanew[imn,'ctxB.type']))
		}else{
            return(as.character(meta.yemen[imy,'ctxB']))
        }
	}))
	meta.yemen[['SXT']] = as.character(sapply(1:nrow(meta.yemen), function(imy){
		imn = i.newcontext.my2mn[imy]
		if (!is.na(imn)){
			return(as.character(metanew[imn,'ICE']))
		}else{
            return(as.character(meta.yemen[imy,'SXT']))
        }
	}))
	meta.yemen[['VSP.II']] = as.character(sapply(1:nrow(meta.yemen), function(imy){
		imn = i.newcontext.my2mn[imy]
		if (!is.na(imn)){
			return(as.character(metanew[imn,'VSP.II.type']))
		}else{
            return(as.character(meta.yemen[imy,'VSP.II']))
        }
	}))
	meta.yemen[['Isolation.year']] = as.factor(sapply(1:nrow(meta.yemen), function(imy){
		imn = i.newcontext.my2mn[imy]
		if (!is.na(imn)){
			isoy = metanew[imn,'Year']
		}else{
		  if (imy %in% c(new.yemen.i, pasteur.yemen.i)){
#		    if (imy %in% new.yemen.i){
			  supn = meta.yemen[imy, 'Supplier.Name']
			  isoy = ifelse(grepl('V. ch.18-', supn), '2018', ifelse(grepl('V. ch.19-', supn), '2019', NA))
#			}else{ # if (imy %in% pasteur.yemen.i){
#			  specid = sub('-PI', '', translategenomeids2strainnames(meta.yemen[imy,'Isolate.name']))
#			  matchspecid = which(meta.yemen[['specimen.id']]==specid)
#			  if (length(matchspecid)>0){ 
#			    supn = meta.yemen[matchspecid[1], 'Supplier.Name']
#			    isoy = ifelse(grepl('V. ch.18-', supn), '2018', ifelse(grepl('V. ch.19-', supn), '2019', NA))
#			  }else{ isoy = NA }
#			}
		  }else{
			isoy = meta.yemen[imy,'Isolation.year']
		  }
		}
		if (isoy %in% c("", NA)){ isoy = NA }
		return(isoy)
	}))
	meta.yemen[['Isolation.date']] = sapply(1:nrow(meta.yemen), function(imy){
		imn = i.newcontext.my2mn[imy]
		if (!is.na(imn)){
			isod = metanew[imn,'Isolation.date']
		}else{ if (imy %in% c(new.yemen.i, pasteur.yemen.i)){
#		  if (imy %in% new.yemen.i){
			isod = as.character(meta.yemen[imy,'collection.date'])
#          }else{ #if (imy %in% pasteur.yemen.i)
#			specid = sub('-PI', '', translategenomeids2strainnames(meta.yemen[imy,'Isolate.name']))
#			matchspecid = which(metaghulam[['specimen.id']]==specid)
##			matchspecid = which(meta.yemen[['specimen.id']]==specid)
#			if (length(matchspecid)>0){ isod = metaghulam[matchspecid[1],'collection.date']
##			if (length(matchspecid)>0){ isod = meta.yemen[matchspecid[1],'collection.date']
#			}else{ isod = NA }
#		  }
		}else{
			isod = as.character(meta.yemen[imy,'Isolation.date'])
		}}
		if (isod %in% c("", NA)){ isod=NA }
		return(isod)
	})
	isodmy = sapply(1:nrow(meta.yemen), function(imy){
		isod = meta.yemen[imy,'Isolation.date']
		isoy = meta.yemen[imy,'Isolation.year']
		if (!is.na(isod)){
			return(as.character(isod))
		}else{
			isom = metanew[imy,'Isolation.month']
			if (!is.na(isom)){
				nmonth = which(months==tolower(as.character(isom)))
				if (length(nmonth)==1){
					return(paste('15', as.character(isom), as.character(isoy), sep='/'))
				}
			}
			return(paste('30', '06', as.character(isoy), sep='/'))
		}
	})
#	meta.yemen[['Isolation.Date']] = as.Date(isodmy, format=ifelse(grepl("[0-9]{2}/[0-9]{2}/[0-9]{4}", isodmy), "%d/%m/%Y", "%Y"))
	meta.yemen[['Isolation.Date']] = as.Date(isodmy, format="%d/%m/%Y")
	meta.yemen[['Isolate.source']] = as.factor(sapply(1:nrow(meta.yemen), function(imy){
		imn = i.newcontext.my2mn[imy]
		if (!is.na(imn)){
			return(metanew[imn,'Source.y'])
		}else{
			return(meta.yemen[imy,'Isolate.source'])
		}
	}))
	meta.yemen[['serotype']] = as.factor(sapply(1:nrow(meta.yemen), function(imy){
		imn = i.newcontext.my2mn[imy]
		if (!is.na(imn)){
			return(metanew[imn,'Serotype_phenotypic'])
		}else{
			return(meta.yemen[imy,'serotype'])
		}
	}))
	print("dim(meta.yemen) (d)")
	print(dim(meta.yemen))
	
	meta.yemen[['lineage_after_Domman_et_al']] = as.factor(as.character(metanew[i.newcontext.my2mn,'lineage_after_Domman_et_al']))
	meta.yemen[['Lineage']] = as.factor(as.character(metanew[i.newcontext.my2mn,'Lineage']))
#	print(summary(meta.yemen[,c("Lineage", "lineage_after_Domman_et_al")]))
	
#	meta.yemen[['Data.source']][new.yemen.i] = 'This study'
#	meta.yemen[['Continent']][new.yemen.i] = 'Asia'
#	meta.yemen[['UN.Subregion']][new.yemen.i] = 'Western Asia'
#	meta.yemen[['Country']][new.yemen.i] = 'Yemen'
	
    meta.yemen[['UN.Subregion']] = as.factor(ifelse(meta.yemen[['UN.Subregion']]=='South Asia', 'Southern Asia', as.character(meta.yemen[['UN.Subregion']])))
    meta.yemen[['UN.Subregion']] = as.factor(ifelse(meta.yemen[['UN.Subregion']]=='Middle Africa', 'Central Africa', as.character(meta.yemen[['UN.Subregion']])))
    
	isolate.source = as.character(meta.yemen[['Isolate.source']])
	meta.yemen[['Isolate.source']] = factor(isolate.source, levels=unique(c(isolate.source, 'Environmental')))
	meta.yemen[['Isolate.source']][new.yemen.i] = ifelse(grepl('E.V.', meta.yemen$Supplier.Name[new.yemen.i]), 'Environmental', 'Human')
	meta.yemen[['Isolate.source']][meta.yemen[['Isolate.source']]==""] = NA
	meta.yemen[['Isolate.source']][startsWith(as.character(meta.yemen[['Isolate.source']]), 'Human')] = 'Human'
	#meta.yemen$Isolation.year = as.numeric(as.character(meta.yemen$Isolation.year))
	meta.yemen[['serotype']][new.yemen.i] = ifelse(meta.yemen[['Ogawa']][new.yemen.i]=='+', 'Ogawa', ifelse(meta.yemen[['Inaba']][new.yemen.i]=='+', 'Inaba', NA))
    
    meta.yemen[['O']][new.yemen.i] = ifelse(meta.yemen[['serotype']][new.yemen.i] %in% c('Inaba', 'Ogawa'), 'O1', meta.yemen[['O']][new.yemen.i])
	char.years = as.character(meta.yemen[['Isolation.year']])
#	meta.yemen[['Isolation.year']] = factor(char.years, levels=as.character(min(as.numeric(char.years), na.rm=T):2019))
#	meta.yemen[['Isolation.year']][new.yemen.i] = ifelse(grepl('V. ch.18-', meta.yemen[['Supplier.Name']][new.yemen.i]), '2018', ifelse(grepl('V. ch.19-', meta.yemen[['Supplier.Name']][new.yemen.i]), '2019', NA))
#	meta.yemen$Isolation.recent.year = as.factor(ifelse(as.numeric(as.character(meta.yemen[,'Isolation.year']))<2012, '.pre-2012', as.character(meta.yemen[,'Isolation.year'])))
	meta.yemen[['Isolation.recent.year']] = factor(ifelse(as.numeric(char.years)<2012, NA, char.years), levels=2012:2019)
#	meta.yemen[['Isolation.yearmonth']] = factor(sapply(strsplit(as.character(meta.yemen[['Isolation.date']]), split='/'), function(x){ if (!is.na(x[1])) return(paste(x[3], x[2], sep='-')) else return(NA) }), levels=allyearmonth1619)
	print("dim(meta.yemen) (e)")
	print(dim(meta.yemen))

	meta.yemen = merge(meta.yemen, det.ncbiamr, by.x="Isolate.name", by.y='row.names', all.x=T)
	meta.yemen = merge(meta.yemen, det.plasmid, by.x="Isolate.name", by.y='row.names', all.x=T)
	meta.yemen = merge(meta.yemen, det.vfdb, by.x="Isolate.name", by.y='row.names', all.x=T)

	# reconcile AMR profiles
#	noAMR.RS = (is.na(meta.yemen[['ampicillin']]) & (as.numeric(as.character(meta.yemen[['Isolation.year']])) < 2018))
	noAMR.RS = is.na(meta.yemen[['ampicillin']])
	meta.yemen[, 'MIC.NAL..mg.L.'] = ifelse(as.character(meta.yemen[, 'MIC.NAL..mg.L.']) %in% c('Unknown', ''), NA, meta.yemen[, 'MIC.NAL..mg.L.'])
	meta.yemen[, 'MIC.NAL..mg.L.'] = as.numeric(ifelse(meta.yemen[, 'MIC.NAL..mg.L.']=='>256', 256, as.character(meta.yemen[, 'MIC.NAL..mg.L.'])))
	meta.yemen[, 'MIC.CIP..mg.L.'] = ifelse(as.character(meta.yemen[, 'MIC.CIP..mg.L.']) %in% c('Unknown', ''), NA, meta.yemen[, 'MIC.CIP..mg.L.'])
	meta.yemen[noAMR.RS, 'nalidixicacid'] = ifelse('R', 'S', meta.yemen[noAMR.RS, 'MIC.NAL..mg.L.']>=64 )
	meta.yemen[noAMR.RS, 'ciprofloxacin'] = ifelse('R', 'S', meta.yemen[noAMR.RS, 'MIC.CIP..mg.L.']>=64 )
	meta.yemen$colistin = ifelse('R', 'S', meta.yemen[noAMR.RS, 'MIC.COL..mg.L.']>2 )
	meta.yemen$polymyxin.B = ifelse('R', 'S', meta.yemen[noAMR.RS, 'MIC.POL..mg.L.']>2 )
	print("dim(meta.yemen) (f)")
	print(dim(meta.yemen))
	
	for (abdb in names(cond.abricate)){
		meta.yemen = merge(meta.yemen, cond.abricate[[abdb]], by="Isolate.name", all.x=T)
		print(c(dim(meta.yemen), abdb))
	}
	print(table(meta.yemen[['Isolate.name']])[table(meta.yemen[['Isolate.name']])>1])
	rownames(meta.yemen) = meta.yemen[['Isolate.name']]

	meta.yemen = meta.yemen[!grepl('^X\\.', colnames(meta.yemen))]
	print("dim(meta.yemen) (g)")
	print(dim(meta.yemen))
	
#	meta.yemen[, 'continuous.date'] = as.double(date.DMY2continous(split.date2DMY(meta.yemen[,'collection.date'])))
	meta.yemen[, 'continuous.date'] = as.double(date.DMY2continous(split.Date2DMY(meta.yemen[,'Isolation.Date'])))
	print("dim(meta.yemen) (h)")
	print(dim(meta.yemen))
	
	nftxtout = file.path(yemen, 'metadata', paste(datatag, 'unified_metadata_yemen2018-2019_basic.txt', sep='.'))
	write.table(meta.yemen, file=nftxtout, row.names=F, sep='\t', fileEncoding="UTF-8")
	print(sprintf("wrote table to '%s'", nftxtout), quote=F)
	return(meta.yemen)
}


#for (datatag in names(nftrees)){
#for (datatag in 'ST555_hybrid'){
for (datatag in 'hybridref-Nov2021'){
#for (datatag in c('hybridref-Nov2021.subcladeH9', 'hybridref-Nov2021.subcladeH689', 'hybridref-Nov2021')){
#for (datatag in 'hybridref-Nov2021.gubbins'){
	print('# # # # # # #')
	print(datatag)
	nftree = nftrees[[datatag]]
    bnftree = lbnftrees[[datatag]]
	print(nftree)
	treeyemen = read.tree(nftree)
    
	treeyemenass = treeyemen
	treeyemenass$tip.label = sapply(treeyemenass$tip.label, sanglane2assid)
	nftreecorrected = paste0(nftree, '.assembly_ids')
	write.tree(treeyemenass, file=nftreecorrected)
	
    treeyemencode = treeyemenass
    treeyemencode$tip.label = sapply(treeyemencode$tip.label, function(x){ 
        code = assid2code[ assid2code['assembly_id']==x, 'locus_tag_prefix']
        if (length(code)==1) return(code)
        else return(x)
    })
    nftreecode = paste0(nftree, '.codes')
    write.tree(treeyemencode, file=nftreecode)
    
    
	genomes = data.frame(Lane=treeyemen$tip.label, assembly_id=treeyemenass$tip.label, locus_tag=treeyemencode$tip.label)
	print("head(genomes) 1")
	print(head(genomes))
    
	nfcladepat = paste(paste0(bnftree, c('\\.clade[A-Z]$', '\\.clade[A-Z].[0-9]$', '\\.clade[A-Z]\\.[0-9]\\.[a-z]$')), collapse='|')
    nfcladetrees = list.files(dirname(nftree), pattern=nfcladepat, full=TRUE)
    if (length(nfcladetrees)>0){
		print("found clade tree files:", quote=F)
		ltipsinclades = lapply(nfcladetrees, function(nfcladetree){
			print(nfcladetree, quote=F)
            cladetree = read.tree(nfcladetree)
            cla = strsplit(nfcladetree, '\\.clade')[[1]][2]
            return(data.frame(Lane=cladetree$tip.label, clade=cla))
        })
        tipsinclades = do.call(rbind, ltipsinclades)
		print("summary(tipsinclades)")
		print(summary(tipsinclades))
		
		meta.yemen.iTOLcladedefs = as.data.frame(t(sapply(nfcladetrees, function(nfcladetree){
			cladetree = read.tree(nfcladetree)
            cla = strsplit(nfcladetree, '\\.clade')[[1]][2]
			n = length(cladetree$tip.label)
			strnames = translategenomeids2strainnames(cladetree$tip.label)
			return(c(paste(strnames[c(1, n)], collapse='|'), cla))
		})))
		colnames(meta.yemen.iTOLcladedefs) = c('clade_id', 'clade')
		rownames(meta.yemen.iTOLcladedefs) = meta.yemen.iTOLcladedefs[['clade']]
		print("head(meta.yemen.iTOLcladedefs)")
		print(head(meta.yemen.iTOLcladedefs))
		
        genomes = merge(genomes, tipsinclades, all.x=T)
#        if (length(intersect(c('main', 'mainclade'), strsplit(datatag, '\\.')[[1]]))>0){
#            genomes[['clade']][is.na(genomes[['clade']])] = 'I'
#			n = length(treeyemen$tip.label)
#			cladeIrow = data.frame(clade_id=paste(treeyemen$tip.label[c(1, n)], collapse='|'), clade='I')
#			meta.yemen.iTOLcladedefs = rbind(cladeIrow, meta.yemen.iTOLcladedefs)
#        }else{ if (length(intersect(c('narrow', 'narrowclade'), strsplit(datatag, '\\.')[[1]]))>0){
#            genomes[['clade']][is.na(genomes[['clade']])] = 'Ig'
#			n = length(treeyemen$tip.label)
#			cladeIgrow = data.frame(clade_id=paste(treeyemen$tip.label[c(1, n)], collapse='|'), clade='Ig')
#			meta.yemen.iTOLcladedefs = rbind(cladeIgrow, meta.yemen.iTOLcladedefs)
#			
#		}}
        print("head(genomes) 2")
        print(head(genomes))
		print("summary(as.factor(genomes[['clade']]))")
		print(summary(as.factor(genomes[['clade']])))
    }else{
        genomes[['clade']] = NA
		meta.yemen.iTOLcladedefs = data.frame(clade_id=character(0), clade=character(0))
    }
    
#	print("names(parse.abricate)")
#	print(names(parse.abricate))
	det.ncbiamr = annotlist2mat(parse.abricate[['ncbi']], genomes[,1])
	colnames(det.ncbiamr) = paste0('.', colnames(det.ncbiamr))
	det.plasmid = annotlist2mat(parse.abricate[['plasmidfinder']], genomes[,1])
	det.vfdb = annotlist2mat(parse.abricate[['vfdb']], genomes[,1])
	colnames(det.vfdb) = paste0('.', colnames(det.vfdb))
	
	det.pheno = setdiff(unlist(lapply(list(det.ncbiamr, det.plasmid, det.vfdb), colnames)), "Isolate.name")

	if (grepl('mapped', nftree)){
        altjoincol = 'Sample'
	}else{
        altjoincol = NULL
    }
	meta.yemen = cureatemetadata(genomes, det.ncbiamr, det.plasmid, det.vfdb, datatag, altjoincol=altjoincol)
	
#	print("meta.yemen[['Isolation.date']]")
#	print(meta.yemen[['Isolation.date']])
	
	meta.amr.plas.int = data.matrix(as.data.frame(meta.yemen[,c(colnames(det.ncbiamr), colnames(det.plasmid))]))
	meta.amr.plas.int[is.na(meta.amr.plas.int)] = 1
    
    plot.amr.corr = FALSE
    if (plot.amr.corr){
        pdf(file=file.path(yemen, 'metadata', paste(datatag, 'yemen2018-2019_ncbiamr_vs_plasmidfinder.pdf', sep='.')), height=20, width=30)
        h = heatmap(meta.amr.plas.int, scale='none', keep.dendro=T)
        meta.amr.plas.int[,colnames(det.plasmid)] = meta.amr.plas.int[,colnames(det.plasmid)] + 3
        heatmap(meta.amr.plas.int, scale='none', Colv=h$Colv, Rowv=h$Rowv, main='NCBI AMR and PlasmidFinder screens')
        dev.off()

        pdf(file=file.path(yemen, 'metadata', paste(datatag, 'yemen2018-2019_all_amr_vs_plasmidfinder.pdf', sep='.')), height=20, width=30)
        for (abdb in c("argannot", "card", "ncbi", "resfinder", "vfdb")){
            vfamrtag = ifelse(abdb=='vfdb', 'virulence factor', 'AMR')
            det.amr = annotlist2mat(parse.abricate[[abdb]])
            yemen.iso = meta.yemen["Isolate.name"]
            yemen.iso = merge(yemen.iso, det.amr, by.x="Isolate.name", by.y='row.names', all=T)
            yemen.iso = merge(yemen.iso, det.plasmid, by.x="Isolate.name", by.y='row.names', all=T)
            amr.plas.int = data.matrix(as.data.frame(yemen.iso[,c(colnames(det.amr), colnames(det.plasmid))]))
            amr.plas.int[is.na(amr.plas.int)] = 0
            h = heatmap(amr.plas.int, scale='none', keep.dendro=T, main=paste(abdb, 'AMR and plasmidfinder screens'))
            amr.plas.int[,colnames(det.plasmid)] = amr.plas.int[,colnames(det.plasmid)] + 3
            heatmap(amr.plas.int, scale='none', Colv=h$Colv, Rowv=h$Rowv, main=paste(abdb, vfamrtag, 'and plasmidfinder screens'))
        }
        dev.off()
    }
	pheno.labo.col = c("amikacin", "ampicillin", "azithromycin", "cephotaxime", "ciprofloxacin", "co.trimoxazol", "doxycycline", "erythromycin", "gentamicin", "chloramphenicol", "nalidixicacid")
#	meta.isol.col = c("Isolation.recent.year", "governorate", "Continent" , "UN.Subregion", "Isolate.source", "Genomic.wave", "Introduction.event", "Lineage", "lineage_after_Domman_et_al")
#	meta.isol.col2 = c("Isolation.year", "Isolation.recent.year", "governorate", "Continent" , "UN.Subregion", "Country", "Locality.Province", "Isolate.source", "Genomic.wave", "Introduction.event", "governorate", "district")
	meta.isol.col = c("clade", "Lineage", "lineage_after_Domman_et_al", "Introduction.event", "Isolate.source", "Continent", "UN.Subregion", "Locality.Province", "governorate.eng")
	meta.isol.col2 = c("Isolation.year", "Isolate.source", "district", "governorate", "governorate.eng", "Locality.Province", "Country", "lineage_after_Domman_et_al", "Genomic.wave", "Introduction.event")
	meta.pheno.col = c(pheno.labo.col, "serotype", colnames(det.ncbiamr), colnames(det.plasmid))
	meta.basic.col = c("Isolation.recent.year", "governorate.eng")
	meta.mge.col = c('PLE1', 'ICP1_phage', 'SXT_ICE', 'IncAC_plasmid', 'MDRcassette')	

	
	for (pl in pheno.labo.col){
#		print(pl)
		p = factor(as.character(meta.yemen[,pl]), levels=c("S", "R", "M"))
#		print(summary(p))
		meta.yemen[[pl]] = p
	}
	
	blast.mge = MGECovSNPFromBLAST(treeyemenass)
	print("grep('^GCF', rownames(blast.mge[['cov']]), value=T)")
	print(grep('^GCF', rownames(blast.mge[['cov']]), value=T))
	write.csv(blast.mge[['cov']], file=file.path(yemen, 'metadata', paste(datatag, 'MGE_coverage.csv')))
	write.csv(blast.mge[['snp']], file=file.path(yemen, 'metadata', paste(datatag, 'MGE_SNPcount.csv')))
	
	mge.pres =  as.data.frame(apply(blast.mge[['cov']], 2, function(x){ 
		sapply(x, function(k){
			if (k >70){ return('yes') 
			}else{if (k >20){ return('partial') 
			}else{ return('no') }}
	})}))
	colnames(mge.pres) = colnames(blast.mge[['cov']])
	print("summary(blast.mge[['cov']])")
	print(summary(blast.mge[['cov']]))
	print("summary(blast.mge[['snp']])")
	print(summary(blast.mge[['snp']]))
	print("summary(mge.pres)")
	print(class(mge.pres))
	print(dim(mge.pres))
	print(summary(mge.pres))
	
	mge.pres[['assemb_wbeT_intact']] = (blast.mge[['cov']][,'wbeT'] == 100) & (is.na(blast.mge[['snp']][,'wbeT']) | (blast.mge[['snp']][,'wbeT'] == 0))
	mge.pres[['assemb_wbeT_mutated']] = (blast.mge[['cov']][,'wbeT'] == 100) & (blast.mge[['snp']][,'wbeT'] > 0)
	
	meta.yemen = merge(meta.yemen, mge.pres, by.x="assembly_id", by.y='row.names', all.x=TRUE, all.y=FALSE)
	rownames(meta.yemen) = meta.yemen[['Isolate.name']]
#	print(dim(meta.yemen))
#	print(summary(meta.yemen))
	
	aribasero[,'reads_wbeT_intact'] = (aribasero[['CP013319.assembled']] %in% c('yes', 'yes_nonunique')) & (aribasero[['CP013319.known_var']]=='no' & aribasero[['CP013319.novel_var']]=='no')
	aribasero[,'reads_wbeT_truncated'] = (aribasero[['CP013319.assembled']]=='interrupted') & (aribasero[['CP013319.known_var']]=='yes' | aribasero[['CP013319.novel_var']]=='yes')
	aribasero[,'reads_wbeT_pointmut'] = (aribasero[['CP013319.assembled']] %in% c('yes', 'yes_nonunique')) & (aribasero[['CP013319.known_var']]=='yes' | aribasero[['CP013319.novel_var']]=='yes')
	aribasero[,'reads_wbeT_mutated'] = aribasero[,'reads_wbeT_truncated'] | aribasero[,'reads_wbeT_pointmut']
	
	meta.yemen = merge(meta.yemen, aribasero, by.x="Isolate.name", by.y='row.names', all.x=TRUE, all.y=FALSE)
	rownames(meta.yemen) = meta.yemen[['Isolate.name']]
	
	meta.yemen = merge(meta.yemen, matOantigenclustmatches, by.x="Isolate.name", by.y='row.names', all.x=TRUE, all.y=FALSE)
	rownames(meta.yemen) = meta.yemen[['Isolate.name']]
	

	names2ids = data.frame(Isolate.name=translategenomeids2strainnames(meta.yemen[['Isolate.name']]), Lane.id=meta.yemen[['Isolate.name']])
	meta.yemen = merge(names2ids, meta.yemen, by.x='Lane.id', by.y='Isolate.name')
	rownames(meta.yemen) = meta.yemen[['Isolate.name']]
	
	
#    meta.yemen[meta.yemen==''] = NA
	
	nfmetaout = file.path(yemen, 'metadata', paste(datatag, 'unified_metadata_yemen2018-2019_extended.RData', sep='.'))
	save(meta.yemen, file=nfmetaout)
	nfcsvout = sub('\\.RData$', '.csv', nfmetaout)
	rowstoexcludeincsv = which(colnames(meta.yemen) %in% c(
		colnames(metaweill[,c(14:19, 22:ncol(metaweill))]), 
		grep("CP013319\\.", colnames(aribasero), value=T),
		det.pheno[sapply(det.pheno, function(dphn){ all(meta.yemen[,dphn]=="", na.rm=T) })]
		)
	)
	write.csv(meta.yemen[,-rowstoexcludeincsv], file=nfcsvout, row.names=F, fileEncoding="UTF-8")
	print(sprintf("wrote table to '%s'", nfcsvout), quote=F)
	
	all.yemen.i = which(meta.yemen[['Country']]=='Yemen')
	
	nftestdists = file.path(yemen, 'metadata', paste(datatag, 'yemen2018-2019_distances_correlation.txt', sep='.'))
	nfplotdists = sub('\\.txt', '.pdf', nftestdists)
    
	test.cor.dists = TRUE
	plot.cor.dists = FALSE
    if (test.cor.dists){
        # compute patristic distances and the time and space distances between samples
		treeyemennames = treeyemen
		treeyemennames$tip.label = translategenomeids2strainnames(treeyemennames$tip.label)
        patristic.dist.all = cophenetic(treeyemennames)

        i.samples.withGPS = which(!is.na(meta.yemen[,'district.GPS.long']) & (as.character(meta.yemen[,'Lane.id']) %in% genomes[,'Lane']))
		samples.withGPS = rownames(meta.yemen[i.samples.withGPS,])

        geo.dist.mat = sapply(i.samples.withGPS, function(i){
            pos1 = as.numeric(meta.yemen[i,c('district.GPS.long', 'district.GPS.lat')])
            sapply(i.samples.withGPS, function(j){
                pos2 = as.numeric(meta.yemen[j,c('district.GPS.long', 'district.GPS.lat')])
                if(all(pos1==pos2)){ return(0) 
                }else{ return(geodetic.distance(pos1, pos2)) }
            })
        })
        colnames(geo.dist.mat) = rownames(geo.dist.mat) = samples.withGPS
        geo.dist = as.dist(geo.dist.mat)
        time.dist.withGPS = dist(meta.yemen[i.samples.withGPS, 'continuous.date'])
        patristic.dist.withGPS = as.dist(as.matrix(patristic.dist.all)[samples.withGPS, samples.withGPS])
        
        i.samples.withcontdate = which(!is.na(meta.yemen[,'continuous.date']) & (as.character(meta.yemen[,'Lane.id']) %in% genomes[,'Lane']))
		
		samples.withcontdate = rownames(meta.yemen[i.samples.withcontdate,])
		time.dist.withcontdate = dist(meta.yemen[i.samples.withcontdate, 'continuous.date'])
        patristic.dist.withcontdate = as.dist(as.matrix(patristic.dist.all)[samples.withcontdate, samples.withcontdate])
		
		samples.all.yemen = rownames(meta.yemen[all.yemen.i,])
		time.dist.allyemen = dist(meta.yemen[all.yemen.i, 'continuous.date'])
        patristic.dist.allyemen = as.dist(as.matrix(patristic.dist.all)[samples.all.yemen, samples.all.yemen])

		print("Test correlation of phylogenetic vs. time distances vs. geographic distances ")
		
		print("  Test correlation of phylogenetic vs. geographic distances (Pearson's correlation)", quote=F)
		write(c("###", "Test correlation of phylogenetic vs. geographic distances (Pearson's correlation)"), file=nftestdists)
		ctgp = cor.test(geo.dist, patristic.dist.withGPS)
#		print(ctgp)
		capture.output(print(ctgp), file=nftestdists, append=T)
		print("  Test correlation of phylogenetic vs. geographic distances (Mantel's test)", quote=F)
		write(c("", "Test correlation of phylogenetic vs. geographic distances (Mantel's test)"), file=nftestdists, append=T)
		mtgp = mantel.randtest(geo.dist, patristic.dist.withGPS, n=100000)
#		print(mtgp)
		capture.output(print(mtgp), file=nftestdists, append=T)
		
		print("  Test correlation of time vs. geographic distances (Pearson's correlation; only samples with GPS data)", quote=F)
		write(c("", "###", "Test correlation of time vs. geographic distances (Pearson's correlation; only samples with GPS data)"), file=nftestdists, append=T)
		ctgt = cor.test(geo.dist, time.dist.withGPS)
#		print(ctgt)
		capture.output(print(ctgt), file=nftestdists, append=T)
		print("  Test correlation of time vs. geographic distances (Mantel's test; only samples with GPS data)",, quote=F)
		write(c("", "Test correlation of time vs. geographic distances (Mantel's test; only samples with GPS data)"), file=nftestdists, append=T)
		mtgt = mantel.randtest(geo.dist, time.dist.withGPS, n=100000)
#		print(mtgt)
		capture.output(print(mtgt), file=nftestdists, append=T)
		
		print("  Test correlation of phylogenetic vs. time distances (Pearson's correlation; only samples with GPS data)", quote=F)
		write(c("", "###", "Test correlation of phylogenetic vs. time distances (Pearson's correlation; only samples with GPS data)"), file=nftestdists, append=T)
		ctpt1 = cor.test(time.dist.withGPS, patristic.dist.withGPS)
#		print(ctpt1)
		capture.output(print(ctpt1), file=nftestdists, append=T)
		print("  Test correlation of phylogenetic vs. time distances (Mantel's test; only samples with GPS data)", quote=F)
		write(c("", "Test correlation of phylogenetic vs. time distances (Mantel's test; only samples with GPS data)"), file=nftestdists, append=T)
		mtpt1 = mantel.randtest(time.dist.withGPS, patristic.dist.withGPS, n=100000)
#		print(mtpt1)
		capture.output(print(mtpt1), file=nftestdists, append=T)
				
		print("  Test correlation of phylogenetic vs. time distances (Pearson's correlation; all samples with fine time data)", quote=F)
		write(c("", "###", "Test correlation of phylogenetic vs. time distances (Pearson's correlation; all samples with fine time data)"), file=nftestdists, append=T)
		ctpt2 = cor.test(time.dist.withcontdate, patristic.dist.withcontdate)
#		print(ctpt2)
		capture.output(print(ctpt2), file=nftestdists, append=T)
		print("  Test correlation of phylogenetic vs. time distances (Mantel's test; all samples with fine time data)", quote=F)
		write(c("", "Test correlation of phylogenetic vs. time distances (Mantel's test; all samples with fine time data)"), file=nftestdists, append=T)
		mtpt2 = mantel.randtest(time.dist.withcontdate, patristic.dist.withcontdate, n=100000)
#		print(mtpt2)
		capture.output(print(mtpt2), file=nftestdists, append=T)
				
		print("  Test correlation of phylogenetic vs. time distances (Pearson's correlation; all Yemen samples, with fine time data)", quote=F)
		write(c("", "###", "Test correlation of phylogenetic vs. time distances (Pearson's correlation; all Yemen samples, with fine time data)"), file=nftestdists, append=T)
		ctpt3 = cor.test(time.dist.allyemen, patristic.dist.allyemen)
#		print(ctpt3)
		capture.output(print(ctpt3), file=nftestdists, append=T)
		print("  Test correlation of phylogenetic vs. time distances (Mantel's test; all Yemen samples, with fine time data)", quote=F)
		write(c("", "Test correlation of phylogenetic vs. time distances (Mantel's test; all Yemen samples, with fine time data)"), file=nftestdists, append=T)
		mtpt3 = mantel.randtest(time.dist.allyemen, patristic.dist.allyemen, n=100000)
#		print(mtpt3)
		capture.output(print(mtpt3), file=nftestdists, append=T)
		
		if (plot.cor.dists){
			pdf(file=nfplotdists, height=15, width=12)
			plot(geo.dist ~ patristic.dist)
			abline(lm(geo.dist ~ patristic.dist), col='red')

			plot(geo.dist ~ time.dist)
			abline(lm(geo.dist ~ time.dist), col='red')

			plot(time.dist ~ patristic.dist)
			abline(lm(time.dist ~ patristic.dist), col='red')

			# plot evolution of samples location with time
			plot(x=meta.yemen[i.samples.withGPS,'continuous.date'], y=as.numeric(meta.yemen[i.samples.withGPS,'governorate']), xlab='Collection date', ylab='Sampling location in Yemen (governorate)')
			plot(meta.yemen[i.samples.withGPS,'governorate'] ~ meta.yemen[i.samples.withGPS,'continuous.date'], xlab='Collection date', ylab='Sampling location in Yemen (governorate)', col=cols.few)
			plot(x=meta.yemen[i.samples.withGPS,'continuous.date'], y=as.numeric(meta.yemen[i.samples.withGPS,'district']), xlab='Collection date', ylab='Sampling location in Yemen (district)')
			plot(meta.yemen[i.samples.withGPS,'district'] ~ meta.yemen[i.samples.withGPS,'continuous.date'], xlab='Collection date', ylab='Sampling location in Yemen (district)', col=cols.seq)
			plot(geo.dist.allyemen ~ time.dist.allyemen)
			abline(lm(geo.dist.allyemen ~ time.dist.allyemen), col='red')
			dev.off()
		}
	}

    recent.years = levels(meta.yemen[,'Isolation.recent.year'])
    vclineages1 = levels(meta.yemen[,'Lineage'])
    vclineages2 = levels(meta.yemen[,'lineage_after_Domman_et_al'])
    all.years = levels(meta.yemen[,'Isolation.year'])
    mic.levels = unlist(sapply(meta.isol.col, function(mic){ levels(meta.yemen[,mic]) }))
    cols.few = c("#FFFFFFFF", brewer.pal(12, 'Set3'))
    col.some = rainbow(38)
    col.mics = rainbow(length(mic.levels))
    cols.seq = c("#FFFFFFFF", rainbow(100))
    cols.rnd = sample(col.mics, length(mic.levels)) ; names(cols.rnd) = mic.levels
    cols.bin = c('white', 'black', 'red')
    cols.recy = c(paste0('green', 1:4), paste0(c('light',''), 'blue'), 'purple', 'red') ; names(cols.recy) = recent.years
#	cols.recy.pal = function(i){ cols.recy[1:i] }
    cols.mge.bin = c(brewer.pal(8, 'Dark2'), brewer.pal(8, 'Accent'))
    cols.long = rep(c(brewer.pal(8, 'Set1'), brewer.pal(8, 'Dark2'), brewer.pal(12, 'Paired'), brewer.pal(12, 'Set3'), brewer.pal(9, 'Pastel1'), brewer.pal(8, 'Accent'), brewer.pal(8, 'Pastel2'), brewer.pal(8, 'Set2')), 4)
    cols.clades = c(brewer.pal(8, 'Set1'), brewer.pal(3, 'Dark2'), brewer.pal(8, 'Paired'), brewer.pal(8, 'Set2')[c(1,3:8)], brewer.pal(8, 'Set3'), brewer.pal(8, 'Set2')[2])[1:length(all.clades)]
	names(cols.clades) = all.clades
	cols.clades['H.9.h'] = '#FFC2D0'
	cols.clades['H.9.e'] = '#BBA094'
	cols.clades['H.9.f'] = '#B3B343'
	print(cols.clades)
    cols.clades.dataset = cols.clades[levels(as.factor(meta.yemen[, "clade"]))]
    cols.clades.yemen = cols.clades[levels(as.factor(meta.yemen[all.yemen.i, "clade"]))]
	cols.govyem = brewer.pal(12, 'Paired'); names(cols.govyem) = sort(unique(meta.yemen[all.yemen.i, "governorate.eng"]))
    
    create.iTol.datasets = TRUE
    if (create.iTol.datasets){
        itoldir = file.path(yemen, 'metadata', paste(datatag, 'iTol.datasets', sep='.'))
        dir.create(itoldir, recursive=T)
        write.iTol.dataset.table.bin(meta.yemen, "reads_wbeT_intact", cols.mge.bin[3+length(meta.mge.col)], "wbeT_intact_in_reads", field_shapes=2, height_factor=3, outdir=itoldir)
        write.iTol.dataset.table.bin(meta.yemen, meta.mge.col, cols.mge.bin[(1:length(meta.mge.col))+2], "MGEs", field_shapes=2, height_factor=3, outdir=itoldir)
        write.iTol.dataset.table.bin(meta.yemen, c('Ogawa', 'Inaba', 'assemb_wbeT_intact',  'assemb_wbeT_mutated', 'reads_wbeT_intact', 'reads_wbeT_mutated'), cols.mge.bin[rep((1:2)+length(meta.mge.col)+2, 3)], "serotype", field_shapes=c(2,2,3,3,4,4), height_factor=3, outdir=itoldir)
        write.iTol.dataset.table.bin(meta.yemen, c('O1', 'O37', 'nonO1_A', 'nonO1_B', 'nonO1_C', 'nonO1_D', 'nonO1_E'), cols.clades[1:7], "serogroup", field_shapes=2, height_factor=3, outdir=itoldir)
        
        for (mge.col in meta.mge.col){
            write.iTol.dataset.table.bin.colstrip(meta.yemen, mge.col, cols.mge.bin[1:2], c('present', 'absent'), outdir=itoldir)
        }
        
        write.iTol.dataset.table.fac.colstrip(meta.yemen, "clade", cols.clades.dataset, colour_branches=FALSE, outdir=itoldir)
#        write.iTol.dataset.table.fac.colstrip(meta.yemen, "clade", cols.clades.dataset, colour_branches=TRUE, outdir=itoldir)
		write.iTol.dataset.table.fac.clade.treecols(meta.yemen.iTOLcladedefs, "clade", cols.clades.dataset, outdir=itoldir)
        
        for (isol.col in c("Lineage", "governorate.eng")){
            write.iTol.dataset.table.fac.colstrip(meta.yemen, isol.col, brewer.pal(12, 'Paired'), outdir=itoldir)
        }
        for (isol.col in c("lineage_after_Domman_et_al","UN.Subregion")){
            write.iTol.dataset.table.fac.colstrip(meta.yemen, isol.col, brewer.pal(12, 'Set3'), outdir=itoldir)
        }
        for (isol.col in c("Introduction.event", "Isolate.source", "Isolation.recent.year", "Continent", "Introduction.event")){
            write.iTol.dataset.table.fac.colstrip(meta.yemen, isol.col, brewer.pal(9, 'Set1'), outdir=itoldir)
        }
        for (isol.col in c("district", "Country")){
            write.iTol.dataset.table.fac.colstrip(meta.yemen, isol.col, cols.long, outdir=itoldir)
        }
		print(sprintf("wrote iTOL dataset tables into folder: '%s/'", itoldir), quote=F)
    }
    
    plot.clade.vs.meta = TRUE
    if (plot.clade.vs.meta){
        print("  plot subclade occurence per source/year/governorate/district", quote=F)
        pdf(file=file.path(yemen, 'metadata', paste(datatag, 'yemen2018-2019_clade_vs_yearloc.pdf', sep='.')), height=7, width=12)

#        for (isol.col in c("Isolate.source", "Isolation.recent.year", "governorate.eng", "district")){
#            conttab = table(meta.yemen[meta.yemen[['Country']]=='Yemen', c('clade', isol.col)])
#            barplot(conttab, beside=T, col=cols.clades.yemen, las=2)
#            legend("topleft", fill=cols.clades.yemen, legend=rownames(conttab))
#        }
#        dev.off()
		
#		pdf(file=file.path(yemen, 'metadata', paste(datatag, 'yemen2018-2019_clade_vs_yearloc.ggplot.pdf', sep='.')), height=7, width=12)

		ggyp = ggplot(meta.yemen[all.yemen.i,], aes(fill=governorate.eng, x=Isolation.Date, drop=FALSE)) +
  		  geom_histogram(binwidth=30, position = 'stack') + 
  		  xlab("Month, Year") + ylab("Number of Isolates") + 
  		  theme_classic(base_size = 20) +
  		  scale_fill_manual(values = brewer.pal(12, 'Paired')) +
  		  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  		  facet_grid(rows = vars(clade, IncAC_plasmid), scales="free", space="free_y", switch="y") +
#  		  scale_x_discrete(breaks = allyearmonth1619)
		  scale_x_date(breaks=seq(from=as.Date('2002/01/01'), to=as.Date('2019/12/30'), by='2 months'), date_labels="%b-%Y")
		ggsave(file.path(yemen, 'metadata', paste(datatag, 'yemen2018-2019_clade_vs_yearloc-Yemen.ggplot.pdf', sep='.')), plot=ggyp, device='pdf', height=7, width=12)

		ggwp = ggplot(meta.yemen[!is.na(meta.yemen[['Isolation.recent.year']]),], 
			   aes(fill=UN.Subregion, x=Isolation.recent.year, drop=FALSE)) +
  		  geom_bar(stat='count', position = 'stack') + 
  		  xlab("Month, Year") + ylab("Number of Isolates") + 
  		  theme_classic(base_size = 20) +
  		  scale_fill_manual(values = brewer.pal(12, 'Set3')) +
  		  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  		  facet_grid(rows = vars(clade, IncAC_plasmid, MDRcassette),
					 scales = "free", space = "free_y", switch = "y") +
		  scale_x_discrete(breaks = 2012:2019)
		ggsave(file.path(yemen, 'metadata', paste(datatag, 'yemen2018-2019_clade_vs_yearloc-world.ggplot.pdf', sep='.')), plot=ggwp, device='pdf', height=7, width=12)
#        dev.off()
    }
	
	plot.yemen.map = TRUE
	if (plot.yemen.map){
		svg(file.path(yemen, 'metadata', paste(datatag, 'yemen2018-2019_yemengovmap_withcladedots.svg', sep='.')))
		plot(st_geometry(yem_governorates_map_polygons))
		latlong = meta.yemen[all.yemen.i,c('clade','governorate.eng', 'district.GPS.lat','district.GPS.long')]
#		latlongstr = paste(latlong[,2:3])
#		uniqlatlong = unique(latlongstr)
#		lsamelatlong = lapply(uniqlatlong, function(ull){ latlong[latlongstr==ull,] })
#		ljittlatlong = lapply(lsamelatlong, function(sll){ data.frame(clade=sll[,1], governorate.eng=sll[,2], jit.lat=jitter(sll[,3]), jit.long=jitter(sll[,4])) })
#		jittlatlong = do.call(rbind, ljittlatlong)
		latlong$jit.long = jitter(latlong$district.GPS.long)
		latlong$jit.lat = jitter(latlong$district.GPS.lat)
		points(x=latlong$jit.long, y=latlong$jit.lat, pch=1, col=cols.clades[latlong$clade])
		dev.off()
		pdf(file.path(yemen, 'metadata', paste(datatag, 'yemen2018-2019_yemengovmap_withcladedots_separateclades.pdf', sep='.')), width=16, height=16)
		layout(matrix(1:4, 2, 2))
		for (cla in unique(latlong$clade)){
			plot(st_geometry(yem_governorates_map_polygons), main=paste("Clade", cla))
			points(x=latlong$jit.long[latlong$clade==cla], y=latlong$jit.lat[latlong$clade==cla], pch=1, col=cols.clades[cla])
		}
		dev.off()
		svg(file.path(yemen, 'metadata', paste(datatag, 'yemen2018-2019_sanaamap_withcladedots.svg', sep='.')))
		latlong.sanaa = latlong[latlong$governorate.eng %in% c('Amanat Al Asimah', 'Sana\'a'),]
		plot(st_geometry(yem_governorates_map_polygons[c(3,51),])) # zoom on Sana'a governorate
		points(x=latlong.sanaa$jit.long, y=latlong.sanaa$jit.lat, pch=1, cex=2, col=cols.clades[latlong.sanaa$clade])
		dev.off()
		pdf(file.path(yemen, 'metadata', paste(datatag, 'yemen2018-2019_sanaamap_withcladedots_separateclades.pdf', sep='.')))
		for (cla in unique(latlong.sanaa$clade)){
			plot(st_geometry(yem_governorates_map_polygons[c(3,51),]), main=paste("Clade", cla))
			points(x=latlong.sanaa$jit.long[latlong.sanaa$clade==cla], y=latlong.sanaa$jit.lat[latlong.sanaa$clade==cla], pch=1, col=cols.clades[cla])
		}
		dev.off()
	}
    
    ggplot.trees = FALSE
    if (ggplot.trees){
    #	m1 = data.matrix(as.data.frame(meta.yemen[,setdiff(colnames(meta.yemen), c('Isolate.name', meta.pheno.col))]))
    #	m1[is.na(m1)] = 0
        mid = data.frame(id=as.character(meta.yemen[,'Isolate.name']), stringsAsFactors=FALSE)
        my1 = meta.yemen[,setdiff(colnames(meta.yemen), c('Isolate.name', meta.pheno.col, meta.isol.col2))]
    #	m2 = meta.yemen[,c(meta.pheno.col)]))
        my2 = data.matrix(as.data.frame(meta.yemen[,c(meta.pheno.col, meta.isol.col2)]))
        my2 = my2 - 1
    #	my = cbind(mid, my1, my2)
        my = cbind(mid, my1)
        mydf = as.data.frame(my, stringsAsFactors=FALSE)[treeyemen$tip.label,]
    #	mydf[,'id'] = rownames(m)
    #	meta.yemen[,'id'] = rownames(meta.yemen)
    #	tmydf = data.frame(treeyemen$tip.label, val=as.double(mydf[treeyemen$tip.label, 'continuous.date']))


        print("  try pretty plots with ggtree", quote=F)
        nfggplots = file.path(yemen, 'metadata', paste(datatag, 'ggplots.pdf', sep='.'))
    #	ptree = ggtree(treeyemen) %<+% meta.yemen[treeyemen$tip.label,] + geom_tippoint(aes(color=governorate)) 
    #	print(facet_plot(ptree, panel="Collection date", data=tmydf, geom=geom_point, mapping=aes(x=continuous.date)))
    #	ptree = ggtree(treeyemen) %<+% meta.yemen[treeyemen$tip.label,] + geom_tippoint(aes(color=governorate.eng))
    #	ptree = ggtree(treeyemen, layout = "circular") + theme_tree2()  + geom_tiplab(size=3)
        ptree = ggtree(treeyemen) %<+% mydf
        ptree = ptree + scale_colour_manual(values=cols.recy, na.translate=FALSE) + geom_tippoint(aes(color=Isolation.recent.year))

        print(summary(mydf[, meta.mge.col]))

        phm1 = gheatmap(ptree, mydf[, meta.mge.col], width=0.4,
                        colnames_angle=-45, colnames_offset_y=.25, hjust=0) +
                        scale_fill_viridis_d(option="D", name="MGE abs/pres")
        ggsave(plot=phm1, filename=sub('\\.pdf', '0.pdf', nfggplots), height=15, width=20, device='pdf')
        phm1 = phm1 + new_scale_fill()
        print(summary(meta.yemen[, meta.isol.col]))
    #	print(summary(mydf[, meta.isol.col[1:4]]))
        phm2 = gheatmap(phm1, meta.yemen[, meta.isol.col], width=0.5, offset=.08,
                        colnames_angle=-45, colnames_offset_y=.25, hjust=0) +
                        scale_fill_manual(values=col.mics, labels=mic.levels, na.translate=FALSE)
    #                    scale_fill_manual(values=cols.rnd, labels=mic.levels, na.translate=FALSE)
    #	phm2 = phm2 + new_scale_fill() 
    #	phm3 = gheatmap(phm2, as.data.frame(mydf[, meta.isol.col[2]]), width=0.5, colnames_angle=-45, offset=.1) +
    #                    scale_fill_viridis_d(option="B", name="Governorate")
    #	phm3 = phm3 + new_scale_fill() 
    #	phm4 = gheatmap(phm3, as.data.frame(mydf[, meta.isol.col[3]]), width=0.5, colnames_angle=-45, offset=.1) +
    #                    scale_fill_viridis_d(option="C", name="Lineage (Wang et al.)")
    #	phm4 = phm4 + new_scale_fill() 
    #	phm5 = gheatmap(phm4, as.data.frame(mydf[, meta.isol.col[4]])) +
    #                    scale_fill_viridis_d(option="D", name="Lineage (Domman et al.)")
        ggsave(plot=phm2, filename=sub('\\.pdf', '1.pdf', nfggplots), height=15, width=20, device='pdf')


        ctree = ggtree(treeyemen, layout='fan', open.angle=20, branch.length='none') %<+% mydf
        ctree = ctree + scale_colour_manual(values=cols.recy, na.translate=FALSE) + geom_tippoint(aes(color=Isolation.recent.year))

        chm1 = gheatmap(ctree, mydf[, meta.mge.col], width=0.3,
                        colnames_angle=-45, colnames_offset_y=.25, hjust=0) +
                        scale_fill_viridis_d(option="D", name="MGE abs/pres")
        chm1 = chm1 + new_scale_fill() 
        chm2 = gheatmap(chm1, meta.yemen[, meta.isol.col], width=0.3, offset=30,
                        colnames_angle=-45, colnames_offset_y=.25, hjust=0) +
                        scale_fill_manual(values=col.mics, labels=mic.levels, na.translate=FALSE)
    #                    scale_fill_manual(values=cols.rnd, labels=mic.levels, na.translate=FALSE)

        ggsave(plot=chm2, filename=sub('\\.pdf', '2.pdf', nfggplots), height=15, width=20, device='pdf')


    #	pcdpla = ptree + geom_facet(panel="Collection date", data=mydf, geom=geom_point, mapping=aes(x='continuous.date', color='IncAC_plasmid'))
    #	ggsave(plot=pcdpla, filename=sub('\\.pdf', '3.pdf', nfggplots), height=15, width=12, device='pdf')

    #	psnptree = plotTreeSNPsMGE(treeyemen)
    #	plapsnptree = psnptree[["IncAC2_hybrid"]] %<+% meta.yemen[treeyemen$tip.label,] + geom_tippoint(aes(color=governorate.eng))
    #	print(gheatmap(plapsnptree, as.data.frame(mydf[treeyemen$tip.label, 'IncA/C2_1'])))
    #	phapsnptree = psnptree[["ICP1_phage"]] %<+% meta.yemen[treeyemen$tip.label,] + geom_tippoint(aes(color=governorate.eng))
    #	print(gheatmap(phapsnptree, as.data.frame(mydf[treeyemen$tip.label, 'ICP1_phage_pcov'])))
    #	phapsnptree = psnptree[["SXT_ICE"]] %<+% meta.yemen[treeyemen$tip.label,] + geom_tippoint(aes(color=governorate.eng))
    #	print(gheatmap(phapsnptree, as.data.frame(mydf[treeyemen$tip.label, 'SXT_ICE_pcov'])))
    #	psnptree + geom_tile(fill=`IncA/C2_1`)
    #	dev.off()

    #	print("  plot metadata and phenotypes vs. phylogenetic tree", quote=F)
    #	pdf(file=file.path(yemen, 'metadata', paste(datatag, 'yemen2018-2019_phylo_vs_metadata.pdf', sep='.')), height=30, width=45)
    #	print(head(mydf[,meta.isol.col]))
    #	phylo.heatmap(treeyemen, mydf[,meta.isol.col], main='isolation metadata', fsize=c(0.25, 3 , 1), col=cols.few) ; legend('topleft', fill=c(cols.few[2:(length(recent.years)+1)], cols.few[2:(length(vclineages1)+1)], cols.few[2:(length(vclineages2)+1)]), legend=c(recent.years, vclineages1, vclineages2))
    #	print(head(my2[,meta.isol.col2]))
    #	phylo.heatmap(treeyemen, my2[,meta.isol.col2], main='isolation metadata', fsize=c(0.25, 3 , 1), col=cols.seq) ; legend('topleft', fill=cols.seq[1:length(all.years)], legend=all.years)
    #	print(head(mydf[,meta.pheno.col]))
    #	phylo.heatmap(treeyemen, mydf[,meta.pheno.col], main='phenotype metadata and predictions', split=c(.25, .75), fsize=c(0.25, 1 , 1), col=cols.few)
    #	print(head(mydf[,c(meta.basic.col, meta.pheno.col)]))
    #	phylo.heatmap(treeyemen, mydf[,c(meta.basic.col, meta.pheno.col)], main='phenotype metadata and predictions', split=c(.25, .75), fsize=c(0.25, 1 , 1), col=cols.few) ; legend('topleft', fill=cols.few[2:(length(recent.years)+1)], legend=recent.years)
    #	dev.off()
	}
	plot.plasmid.years = FALSE
    if (plot.plasmid.years){
        print("  plot plasmid occurence per year", quote=F)
        pdf(file=file.path(yemen, 'metadata', paste(datatag, 'yemen2018-2019_IncAC2_vs_year.pdf', sep='.')), height=15, width=12)
        layout(matrix(1:2, 2, 1))
        plastag = "IncA/C2_1"
        metaplas = meta.yemen[,c("Isolation.recent.year", plastag)]
        plastitle = paste("Occurence of", plastag, "per year")
        barplot(t(table(metaplas)), col=c('lightgrey', 'red'), main=paste(plastitle, "(all isolates)")) ; legend('top', fill=c('lightgrey', 'red'), legend=paste(plastag, c("absent", "present")))
        barplot(t(table(metaplas[meta.yemen[,"Country"]=="Yemen",])), col=c('lightgrey', 'red'), main=paste(plastitle, "(Yemen isolates)")) ; legend('top', fill=c('lightgrey', 'red'), legend=paste(plastag, c("absent", "present")))
        dev.off()
	}
	isolingenomes = (colpanaroo2isol %in% meta.yemen[,'Lane.id'])
	print("table(isolingenomes)")
	print(table(isolingenomes))
	panaroopggenes = t(data.matrix(panaroopangenomegenes[, which(isolingenomes)+1]))
	rownames(panaroopggenes) = translategenomeids2strainnames(colpanaroo2isol[isolingenomes])
	print("dim(panaroopggenes)")
	print(dim(panaroopggenes))
	print("Test correlation between ICP1 phage prsence and Panaroo pangenome genes", quote=F)
	print("panaroopggenes[1:10,1:10]")
	print(panaroopggenes[1:10,1:10])
	print("table(meta.yemen[rownames(panaroopggenes), 'ICP1_phage'])")
	print(table(meta.yemen[rownames(panaroopggenes), 'ICP1_phage']))
	icp1pres = ifelse(as.character(meta.yemen[rownames(panaroopggenes), 'ICP1_phage'])=='yes', 1, 0)
	print("summary(icp1pres)")
	print(summary(icp1pres))
	corICP1vspgg = as.data.frame(t(apply(panaroopggenes, 2, function(x){
		ct = cor.test(x, icp1pres)
		return(c(ct$estimate, ct$p.value))				 
	})))
	colnames(corICP1vspgg) = c("cor", "p.value")
	print("summary(corICP1vspgg)")
	print(summary(corICP1vspgg))
	print(quantile(corICP1vspgg[!is.na(corICP1vspgg[,'p.value']),'cor'], p=0:100/100))
	write.table(corICP1vspgg[which(corICP1vspgg[,'p.value']<0.00001),], file=file.path(yemen, 'metadata', paste(datatag, 'corr-ICP1_phage-panaroogenes.tab', sep='_')), row.names=T, col.names=T, sep='\t')
	
    plot.panaroo = FALSE
    if (plot.panaroo){
        print("  plot pangenome gene/structural variants presence/absence vs. phylogenetic tree", quote=F)
    #	isolingenomes = (colpanaroo2isol %in% genomes[,1])
        panaroopggenefreq = as.double(apply(panaroopggenes, 2, sum, na.rm=T)) / nrow(panaroopggenes)
        panaroopggenes = panaroopggenes[,which((panaroopggenefreq > 0.02) & (panaroopggenefreq < 0.98))]
        panaroopggenes = t(data.matrix(sapply(treeyemen$tip.label, function(isol){ 
            if (isol %in% rownames(panaroopggenes)){ return(panaroopggenes[isol,]) }
            else{ return(rep(NA, ncol(panaroopggenes))) }
        }))) ; rownames(panaroopggenes) = treeyemen$tip.label
        print(sprintf('    size of pangenome gene presence/absence matrix: %d, %d', nrow(panaroopggenes), ncol(panaroopggenes)), quote=F)
        panaroopgstruct = t(data.matrix(panaroopangenomestruct[, which(isolingenomes)+1])) ; rownames(panaroopgstruct) = colpanaroo2isol[isolingenomes]
        panaroopgstrucfreq = as.double(apply(panaroopgstruct, 2, sum, na.rm=T)) / nrow(panaroopgstruct)
        panaroopgstruct = panaroopgstruct[,which((panaroopgstrucfreq > 0.02) & (panaroopgstrucfreq < 0.98))]
        panaroopgstruct = t(data.matrix(sapply(treeyemen$tip.label, function(isol){ 
            if (isol %in% rownames(panaroopgstruct)){ return(panaroopgstruct[isol,]) }
            else{ return(rep(NA, ncol(panaroopgstruct))) }
        }))) ; rownames(panaroopgstruct) = treeyemen$tip.label
        print(sprintf('    size of pangenome structural variants presence/absence matrix: %d, %d', nrow(panaroopgstruct), ncol(panaroopgstruct)), quote=F)

        years1819 = as.numeric(as.character(meta.yemen[treeyemen$tip.label, 'Isolation.recent.year']) %in% c('2018', '2019'))
        print(table(years1819))
        pggenes.cor.years1819 = apply(panaroopggenes, 2, function(x){ cor(x, years1819, use="na.or.complete") })
        pgstruct.cor.years1819 = apply(panaroopgstruct, 2, function(x){ cor(x, years1819, use="na.or.complete") })
        print("  find pangenome genes that correlate with years 2018-2019")
        print(summary(pggenes.cor.years1819))
        print(pggenes.cor.years1819[which(pggenes.cor.years1819 >  0.8)])
        print(pggenes.cor.years1819[which(pggenes.cor.years1819 < -0.8)])
        print("  find pangenome structural variants that correlate with years 2018-2019")
        print(summary(pgstruct.cor.years1819))
        print(pgstruct.cor.years1819[which(pgstruct.cor.years1819 >  0.8)])
        print(pgstruct.cor.years1819[which(pgstruct.cor.years1819 < -0.8)])

    #	ykfG~~~ykfG_1-group_3154-group_6121    group_2600-group_6121-group_3154
    #                          0.8724461                           0.8390600
    ## differential insertion of ISVch4 transposase (group_2600-group_6121) in phage (group_3154 = phage tape tail protein)

    #   speG_1-group_4248-group_4681   ybjQ-group_4681-group_4248
    #                  -0.9285687                   -0.9357462
    ## potential synteny break; at contig end in 33224_2#66

    #> a = panaroopangenomestruct[unique(unlist(sapply(c('ykfG~~~ykfG_1', 'group_3154', 'group_6121', 'group_2600', 'speG_1', 'group_4248', 'group_4681', 'ybjQ', 'group_4681', 'group_4248'), function(x){ grep(x, panaroopangenomestruct[,'Gene'], value=T) }))), c('CNRVC140150', 'X33224_2.66')]
    #> a[apply(a, 1, sum)>0,]
    #                                                CNRVC140150 X33224_2.66
    #ykfG~~~ykfG_1-group_3154-group_6121                       0           1
    #group_3154-ykfG~~~ykfG_1-group_4984                       1           1
    #ykfG~~~ykfG_1-group_4984-nanH_2~~~nanH~~~nanH_1           1           1
    #group_2600-group_6121-group_3154                          0           1
    #group_2600-group_6121-group_1878                          1           0
    #speG_1-group_4248-group_4681                              1           0
    #ybjQ-group_4681-group_4248                                1           0
    #ybjQ-csgD_3~~~csgD_2~~~csgD_1-ppiC                        1           1

        pdf(file=file.path(yemen, 'metadata', paste(datatag, 'yemen2018-2019_phylo_vs_panaroo-pangenome.pdf', sep='.')), height=30, width=45)
        print("  plot pangenome gene variants presence/absence vs. phylogenetic tree", quote=F)
        hist(pggenes.cor.years1819)
    #	mpgg = merge(mydf[,c(meta.basic.col)], panaroopggenes, by='row.names')
    #	rownames(mpgg) = mpgg[,1] ; mpgg = as.data.frame(mpgg[2:ncol(mpgg)])
    #	phylo.heatmap(treeyemen, mpgg, main='Pangenome gene presence/absence', fsize=c(0.25, 1 , 1), col=cols.few) ; legend('topleft', fill=cols.few[1:length(recent.years)], legend=recent.years)
        phylo.heatmap(treeyemen, panaroopggenes, main='Pangenome gene presence/absence', fsize=c(0.25, 1 , 1), col=cols.bin)
        print("  plot pangenome structural variants presence/absence vs. phylogenetic tree", quote=F)
        hist(pgstruct.cor.years1819)
    #	mpgs = merge(mydf[,c(meta.basic.col)], panaroopgstruct, by='row.names')
    #	rownames(mpgs) = mpgs[,1] ; mpgs = as.data.frame(mpgg[2:ncol(mpgs)])
    #	phylo.heatmap(treeyemen, mpgs, main='Pangenome structural variants presence/absence', fsize=c(0.25, 1 , 1), col=cols.few) ; legend('topleft', fill=cols.few[1:length(recent.years)], legend=recent.years)
        phylo.heatmap(treeyemen, panaroopgstruct, main='Pangenome structural variants presence/absence', fsize=c(0.25, .1 , 1), col=cols.bin)
        dev.off()

        pgghm = gheatmap(ptree, panaroopggenes, width=0.65, colnames=FALSE, low='white', high='red')
        ggsave(plot=pgghm, filename=sub('\\.pdf', '_panaroo-pangenome-genes.pdf', nfggplots), height=15, width=20, device='pdf')
        pgshm = gheatmap(ptree, panaroopgstruct, width=0.65, colnames=FALSE, low='white', high='red')
        ggsave(plot=pgshm, filename=sub('\\.pdf', '_panaroo-pangenome-struct.pdf', nfggplots), height=15, width=20, device='pdf')
		
	}
}

lanesamples = merge(sangerids, seqid2name, by='Sample')
mapped450 = read.csv('~/yemen2019/metadata/hybridref_mapped.consensus.full.unified_metadata_yemen2018-2019_extended.csv')
assemb882 = read.csv('~/yemen2019/metadata/contextallyemen2019.assembly_ids.full.nocontam.tre.unified_metadata_yemen2018-2019_extended.csv')
lanesamples[,'in.882.assembled.genomes'] = ifelse(lanesamples$Lane %in% assemb882$Isolate.name, 'yes', 'no')
lanesamples[,'in.450.mapped.genomes'] = ifelse(lanesamples$Lane %in% mapped450$Isolate.name, 'yes', 'no')
write.csv(lanesamples, file=file.path(yemen, 'metadata', 'Yemen20218-2019_sample-lane_ids.csv'))
