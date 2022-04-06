[![DOI](https://zenodo.org/badge/447364930.svg)](https://zenodo.org/badge/latestdoi/447364930)

# yemenpaper
Supporting data and scripts for paper on the Yemen cholera outbreak in 2018-2019

### find lineage-specific SNPs 
```sh
fullvcf='hybridref-Nov2021_all_mapped.snp.vcf'
fulltree='hybridref-Nov2021_all_mapped.snp.aln.rba.raxml.support.full.rooted.H9hsisterH9guniteH9c'
mainclade='hybridref-Nov2021_all_mapped.snp.aln.rba.raxml.support.subcladeH9.rooted.monoH9gmonoH9c'

~/scripts/vibrio/find_clade_specific_SNPs.r ${fullvcf} ${fulltree} "H.9=${mainclade},H.9.Yemen=${mainclade}.cladeH.9.efgh" "H.9.Yemen=${mainclade}.cladeH.9.efgh,H.9.gh=${mainclade}.cladeH.9.gh" "H.9.gh=${mainclade}.cladeH.9.gh,H.9.g=${mainclade}.cladeH.9.g" "H.9.gh=${mainclade}.cladeH.9.gh,H.9.h=${mainclade}.cladeH.9.h" 

t132019refannot='Vibrio_cholerae_CNRVC190243.gbk'
python2.7 ~/scripts/vibrio/extract_annot_clade_specific_SNPs.py ${t132019refannot} ${fullvcf}
```

### find lineage-specific genes
Relies on a prior run of [Panaroo](https://github.com/gtonkinhill/panaroo).
```sh
panaroodir='panaroo/contextallyemen2019'

~/scripts/vibrio/find_clade_specific_panaroogenes.r ${panaroodir} ${fulltree} "H.9=${mainclade},H.9.Yemen=${mainclade}.cladeH.9.efgh" "H.9.Yemen=${mainclade}.cladeH.9.efgh,H.9.gh=${mainclade}.cladeH.9.gh" "H.9.gh=${mainclade}.cladeH.9.gh,H.9.g=${mainclade}.cladeH.9.g" "H.9.gh=${mainclade}.cladeH.9.gh,H.9.h=${mainclade}.cladeH.9.h" 
```
