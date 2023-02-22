#!/usr/bin/env Rscript
blastoutfileds = c('qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
uniqcovlen = function(startendtable){
    a = apply(startendtable, 1, function(x){ return(x[1]:x[2]) })
    if (class(a)[1]=='list'){ return(length(unique(do.call('c', a))))
    }else{ return(length(a)) }
}

ctgvsrefdir = Sys.getenv('ctgvsrefdir')
lnfblastout = list.files(ctgvsrefdir)
genomes = sub('(^.*)_contigs_vs_OantigenJVI.ref', '\\1', lnfblastout)
genomes = sub('(^.*)\\.1_\\1', '\\1', genomes)
genomes = sub('(^.*\\.[1-9])_ASM[0-9]{5,7}v[1-9]', '\\1', genomes)
names(lnfblastout) = genomes
nfcovscoloci = paste(ctgvsrefdir, 'allcoveragescore.RData', sep='_')
if (file.exists(nfcovscoloci)){
  print(paste("loading:", nfcovscoloci), quote=F)
  load(nfcovscoloci)
}else{
  covscoloci = lapply(lnfblastout, function(nfblastout){
    print(nfblastout)
    # first check file has any content
    scanblastout = scan(file=file.path(ctgvsrefdir, nfblastout), nline=1, what='character', quiet=T)
    if (length(scanblastout)==0){ return(data.frame(uniqcov=0, totcov=0, totsco=0, normsco=0, maxsco=0, numcontig.maxsco=NA, O.type=NA)) }
    # blast output for each genome (query: genome contigs; subject: O-antigen LPS loci db)
    blastout = read.table(file.path(ctgvsrefdir, nfblastout), comment.char='', sep='\t', header=F)
    colnames(blastout) = blastoutfileds
#    # determine syntax for numbering of contigs, between type "ERS1234567|SC|contig000001" or "[NZ_]WXYZ01000001.1" or "Anything.1"
#	# case "ERS1234567|SC|contig000001"
#    sccontig = grepl('\\|SC\\|contig', blastout[,1])
#    blastout[['contignum']] = as.numeric(ifelse(sccontig, sub('^.*\\|SC\\|contig([0-9]+)', '\\1', blastout[,1]), sub('^.*\\.|([0-9]+)', '\\1', blastout[,1])))
#    blastout[['biosample.or.wgs']] = ifelse(sccontig, sub('(^.*)\\|SC\\|contig[0-9]+', '\\1', blastout[,1]), sub('(^.*)\\.[0-9]+', '\\1', blastout[,1]))
#    isol = unique(blastout[['biosample.or.wgs']])
#    if (length(isol)>1){ 
#      # case "[NZ_]WXYZ01000001.1"
#	  ncbicontig = grepl('^(NZ_)*[A-Z]{4}01', blastout[['biosample.or.wgs']])
#      blastout[['contignum']] = ifelse(ncbicontig, as.numeric(sub('^(NZ_)*[A-Z]{4}01', '', blastout[['biosample.or.wgs']])), blastout[['contignum']])
#      blastout[['biosample.or.wgs']] = ifelse(ncbicontig, substr(blastout[['biosample.or.wgs']], 1, 6), blastout[['biosample.or.wgs']])
#      isol = unique(blastout[['biosample.or.wgs']])
#      if (length(isol)>1){
#		  blastout[['biosample.or.wgs']] = sub('(^[^\\.]+).*', '\\1', blastout[,1])
#	  }
#    }
    # define coverage per subject, first broken down per contig/query sequence
    uscblout = unique(blastout[,c('saccver','qaccver')])
    covscopersubjcont = as.data.frame(t(sapply(1:nrow(uscblout), function(k){
      s = uscblout[k,'saccver']
      q = uscblout[k,'qaccver']
      qstartendscotable = blastout[blastout[,'saccver']==s & blastout[,'qaccver']==q, c('qstart', 'qend', 'pident', 'bitscore')]
      # account for multiple HSPs (discontiguous match and/or repeats)
#      if (nfblastout=="GCF_003205555.1_ASM320555v1_contigs_vs_OantigenJVI.ref"){
#          print(c('saccver','qaccver'))
#          print(s)
#          print(c)
#          print("qstartendscotable")
#          print(qstartendscotable)
#      }
      uniqcov = uniqcovlen(qstartendscotable)
	  alnlens = abs(qstartendscotable[,2] - qstartendscotable[,1] + 1)
      totcov = sum(alnlens)
      totsco = sum(qstartendscotable[,4])
      normsco = totsco * uniqcov / totcov
	  normpid = sum(qstartendscotable[,3] * alnlens) / totcov
      return(c(uniqcov, totcov, totsco, normsco, normpid))
    })))
    colnames(covscopersubjcont) = c('uniqcov', 'totcov', 'totsco', 'normsco', 'normpid')
    covscopersubjcont = cbind(uscblout, covscopersubjcont)
	
    covscopersubj = as.data.frame(t(sapply(unique(covscopersubjcont[,'saccver']), function(s){
      contcovsco = covscopersubjcont[covscopersubjcont[,'saccver']==s, c('qaccver','uniqcov', 'totcov', 'totsco', 'normsco', 'normpid')]
      subjcovsco = apply(contcovsco[,2:5], 2, sum)
      maxsco = max(contcovsco[,'normsco'])
	  imaxsco = which(contcovsco[,'normsco']==maxsco)[1]
      contmaxsco = contcovsco[imaxsco,'qaccver']
	  normpid = sum(contcovsco[,6] * contcovsco[,3]) / sum(contcovsco[,3])
      return(c(subjcovsco, normpid, maxsco, contmaxsco))
    })))
    colnames(covscopersubj) = c('uniqcov', 'totcov', 'totsco', 'normsco', 'normpid', 'maxsco', 'contig.maxsco')
    covscopersubj[,'O.type'] = unique(covscopersubjcont[,'saccver'])
    return(covscopersubj)
  })
  save(covscoloci, file=nfcovscoloci)
  print(paste("saved:", nfcovscoloci), quote=F)
}
#print(lapply(covscoloci, head))
covscoloci.ncol = sapply(covscoloci, ncol)
print("ignore genomes with problematic results:")
print(covscoloci.ncol[covscoloci.ncol != 8 ])
covscoloci = covscoloci[-which(covscoloci.ncol != 8 )]

topnormsco = do.call(rbind, lapply(covscoloci, function(covscopersubj){
    covscopersubj[covscopersubj[, 'normsco']==max(covscopersubj[, 'normsco']),]
}))
write.table(topnormsco, paste(ctgvsrefdir, 'topnormscore.txt', sep='_'), sep='\t', col.names=T, row.names=T)

topmaxsco = do.call(rbind, lapply(covscoloci, function(covscopersubj){
    covscopersubj[covscopersubj[, 'maxsco']==max(covscopersubj[, 'maxsco']),]
}))
write.table(topmaxsco, paste(ctgvsrefdir, 'topmaxscore.txt', sep='_'), sep='\t', col.names=T, row.names=T)

topcov = do.call(rbind, lapply(covscoloci, function(covscopersubj){
    covscopersubj[covscopersubj[, 'uniqcov']==max(covscopersubj[, 'uniqcov']),]
}))
write.table(topcov, paste(ctgvsrefdir, 'topcoverage.txt', sep='_'), sep='\t', col.names=T, row.names=T)

topcov = do.call(rbind, lapply(covscoloci, function(covscopersubj){
    covscopersubj[covscopersubj[, 'normpid']==max(covscopersubj[, 'normpid']),]
}))
write.table(topcov, paste(ctgvsrefdir, 'toppercentid.txt', sep='_'), sep='\t', col.names=T, row.names=T)
