#!/usr/bin/env bash
### map reads to the MGEs

export yemen=${PWD}/yemen2019

param="${@}"

resume="$(echo ${param} | grep -o 'resume')"
trim="$(echo ${param} | grep -o 'trim')"
ncpus="$(echo ${param} | grep -o 'ncpus=[0-9]\+' | cut -d'=' -f2)"
[ -z "${ncpus}" ] && ncpus=$(nproc)
MGEs="$(echo ${param} | grep -o "MGEs=.\+"	 | cut -d'=' -f2 | tr ',' ' ')"
[ -z "${MGEs}" ] && MGEs="SXT_ICE ICP1_phage ICP1_VMJ710 ICP1_2012_A IncAC2_plasmid IncAC2_hybrid ST555_hybrid"
#MGEs="SXT_ICE ICP1_phage IncAC2_plasmid IncAC2_Pasteur IncAC2_hybrid"
#MGEs="IncAC2_hybrid"

echo "readsetpath2datalist: ${readsetpath2datalist}"
echo "MGEs: ${MGEs}"
echo "trim: ${trim}"
echo "ncpus: ${ncpus}"
echo "resume: ${resume}"

mkdir -p ${yemen}/read_mapping/
cd ${yemen}/read_mapping/
mkdir -p trimmed_reads/ mapped_reads/ references/ compare/
for d in mapped_reads references ; do
  for mge in ${MGEs} ; do
    mkdir -p $d/$mge
  done
done

module load trimmomatic/0.39--1
module load bwa/0.7.17-r1188
module load samtools/1.9

#ln -s "${incac}/plasmidfinder_IncAC_plasmid_33224_2#261_contig23.fna" ${yemen}/read_mapping/references/IncAC2_plasmid/IncAC2_plasmid_reference.fa
#ln -s "${yemen}/Pasteur_seq/RND_CNRVC190243_consensus.fasta" ${yemen}/read_mapping/references/IncAC2_Pasteur/IncAC2_Pasteur_reference.fa
#ln -s "${yemen}/SXT_screen/ICEVchind5_GQ463142.1.fasta" ${yemen}/read_mapping/references/SXT_ICE/SXT_ICE_reference.fa
#ln -s $(grep '33224_2#251' ${yemen}/${vdiv}_contigs_fasta_list) ${yemen}/read_mapping/references/ICP1_phage/ICP1_phage_reference.fa
#ln -s ${yemen}/Pasteur_seq/Vibrio_cholerae_CNRVC190243-7PET-T13/unicyler_hybrid_assembly/assembly.fasta ${yemen}/read_mapping/references/IncAC2_hybrid/IncAC2_hybrid_reference.fa
#ln -s ${yemen}/Pasteur_seq/Vibrio_cholerae_CNRVC190247-ST555/unicyclerHybrid_CNRVC190247_OK/assembly_chr2rotated.fasta ${yemen}/read_mapping/references/ST555_hybrid/ST555_hybrid_reference.fa

for mge in ${MGEs} ; do
  cd ${yemen}/read_mapping/references/${mge}/
  bwa index -p ${mge} ${mge}_reference.fa
  samtools faidx ${mge}_reference.fa
done

cd ${yemen}/read_mapping/

for path2data in $(cat ${readsetpath2datalist}) ; do
  lane=$(basename ${path2data})
  echo -e "\n# # # # # # # # # #\n${lane}\n"
  rawreadsrad=${path2data}/${lane}
  rawreadsin="${rawreadsrad}_1.fastq.gz ${rawreadsrad}_2.fastq.gz"
  echo "input reads: ${rawreadsin}"
  if [ -z "${trim}" ] ; then
    readstomap="${rawreadsin}"
  else
    trimmedreadsdir=${yemen}/read_mapping/trimmed_reads/${lane}
    mkdir -p ${trimmedreadsdir}/
    echo "trim reads; output to: ${trimmedreadsdir}/"
    trimmedreadsrad=${trimmedreadsdir}/${lane}
    trimmedreadout="${trimmedreadsrad}_1.fastq.paired.gz ${trimmedreadsrad}_1.fastq.unpaired.gz ${trimmedreadsrad}_2.fastq.paired.gz ${trimmedreadsrad}_2.fastq.unpaired.gz"
    testtrimmedreadout="-s $(echo ${trimmedreadout} | sed -e 's/ / \&\& -s /g')"
    trimmocmd="trimmomatic PE -threads ${ncpus} ${rawreadsin} ${trimmedreadout} LEADING:3 TRAILING:3 MINLEN:36"
    [[ ${testtrimmedreadout} ]] || (echo "# ${trimmocmd}" && eval "${trimmocmd}")
    readstomap="${trimmedreadsrad}_1.fastq.paired.gz ${trimmedreadsrad}_2.fastq.paired.gz"
  fi
  for mge in ${MGEs} ; do
    echo -e "\n# # #\n${mge}\n"
    cd ${yemen}/read_mapping/references/${mge}/
	mgetag=mapped_${mge}.paired
    outd=${yemen}/read_mapping/mapped_reads/${mge}/${lane}
    outfrad=${outd}/${lane}_${mgetag}
    outf1=${outfrad}.sam
    outf2=${outfrad}.bam
    outdp=${outfrad}.depth.txt
    outf3=${outfrad}.bcf
    outf4=${outfrad}.vcf.gz
    outf41=${outfrad}.noindel.vcf.gz
    outf5=${outfrad}.filteredDP5.vcf
    outf6=${outfrad}.filteredDP10.vcf
    outf61=${outfrad}.noindel.filteredDP10.vcf
    outf7=${outfrad}.consensus.fa
    outf71=${outfrad}.noindel.consensus.fa
    bwamemcmd="bwa mem -t ${ncpus} ${mge} ${readstomap} > ${outf1}"
    samtoolscmd="samtools view -h -F260 -f 2 -q 25 -@ ${ncpus} ${outf1} | awk 'substr(\$0,1,1)==\"@\" || (\$9>= 75 && \$9<=1000) || (\$9<=-75 && \$9>=-1000)' | samtools view -@ ${ncpus} -hb - | samtools sort -@ ${ncpus} - > ${outf2} && rm ${outf1}"
	samtoolsdpcmd="samtools depth ${outf2} > ${outdp}"
    bcftoolscmd="bcftools mpileup -C 50 -L 1000 -d 1000 -m 5 -f ${mge}_reference.fa -O b ${outf2} > ${outf3}"
    bcfvcfcmd="bcftools call -m --ploidy 1 -O z --threads ${ncpus} ${outf3} > ${outf4} && bcftools index --threads ${ncpus} ${outf4}"
    bcfvcfcmdnoindel="bcftools call -m --ploidy 1 -O z -V indels --threads ${ncpus} ${outf3} > ${outf41} && bcftools index --threads ${ncpus} ${outf41}"
    vcffilt5cmd="bgzip -cd -@ ${ncpus} ${outf4} | vcfutils.pl varFilter -d 5 - > ${outf5} && bgzip -f ${outf5} && bcftools index --threads ${ncpus} ${outf5}.gz"
    vcffilt10cmd="bgzip -cd -@ ${ncpus} ${outf4} | vcfutils.pl varFilter -d 10 - > ${outf6} && bgzip -f ${outf6} && bcftools index --threads ${ncpus} ${outf6}.gz"
    vcffilt10cmdnoindel="bgzip -cd -@ ${ncpus} ${outf41} | vcfutils.pl varFilter -d 10 - > ${outf61} && bgzip -f ${outf61} && bcftools index --threads ${ncpus} ${outf61}.gz"
    consensudcmd="bcftools consensus -f ${mge}_reference.fa ${outf6}.gz | sed -e 's/^>\(.\+\)$/>${lane} \[mapped to ${mge} \1\]/g' | gzip > ${outf7}.gz"
    consensudcmdnoindel="bcftools consensus -f ${mge}_reference.fa ${outf61}.gz | sed -e 's/^>\(.\+\)$/>${lane} \[mapped to ${mge} \1\]/g' | gzip > ${outf71}.gz"
	mkdir -p ${outd}/
    [[ ${?} -eq 0 && $(cat ${outf1} 2> /dev/null | wc -c) -lt 400 ]] && (echo "# ${bwamemcmd}" && eval "${bwamemcmd}")
    [[ ${?} -eq 0 && $(cat ${outf2} 2> /dev/null | wc -c) -lt 100 ]] && (echo "# ${samtoolscmd}" && eval "${samtoolscmd}")
    [[ ${?} -eq 0 && $(cat ${outdp} 2> /dev/null | wc -c) -lt 10 ]] && (echo "# ${samtoolsdpcmd}" && eval "${samtoolsdpcmd}")
    [[ ${?} -eq 0 && $(cat ${outf3} 2> /dev/null | wc -c) -lt 1000 ]] && (echo "# ${bcftoolscmd}" && eval "${bcftoolscmd}")
#    [[ ${?} -eq 0 && $(cat ${outf4} 2> /dev/null | wc -c) -lt 1000 ]] && (echo "# ${bcfvcfcmd}" && eval "${bcfvcfcmd}")
    [[ ${?} -eq 0 && $(cat ${outf41} 2> /dev/null | wc -c) -lt 1000 ]] && (echo "# ${bcfvcfcmdnoindel}" && eval "${bcfvcfcmdnoindel}")
#    [[ ${?} -eq 0 && $(cat ${outf5}.gz 2> /dev/null | wc -c) -lt 1000 ]] && (echo "# ${vcffilt5cmd}" && eval "${vcffilt5cmd}")
#    [[ ${?} -eq 0 && $(cat ${outf6}.gz 2> /dev/null | wc -c) -lt 1000 ]] && (echo "# ${vcffilt10cmd}" && eval "${vcffilt10cmd}")
    [[ ${?} -eq 0 && $(cat ${outf61}.gz 2> /dev/null | wc -c) -lt 1000 ]] && (echo "# ${vcffilt10cmdnoindel}" && eval "${vcffilt10cmdnoindel}")
#    [[ ${?} -eq 0 && $(cat ${outf7} 2> /dev/null | wc -c) -lt 10000 ]] && (echo "# ${consensudcmd}" && eval "${consensudcmd}")
    [[ ${?} -eq 0 && $(cat ${outf71}.gz 2> /dev/null | wc -c) -lt 10000 ]] && (echo "# ${consensudcmdnoindel}" && eval "${consensudcmdnoindel}")
	
	if [ "${mge}" == 'IncAC2_hybrid' ] ; then
	  if [[ -s ${outf71}.gz ]] ; then
	    nlineplasmid=$(zcat ${outf71}.gz | grep -n '^>.\+ \[mapped to IncAC2_hybrid 3 ' | cut -d':' -f1)
		zcat ${outf71}.gz | tail -n +${nlineplasmid} | gzip > ${outf71}.plasmidonly.gz
      fi
	fi
  done
done

#for mge in ${MGEs} ; do
#  cd ${yemen}/read_mapping/mapped_reads/${mge}/
#  mgetag=mapped_${mge}.paired
#  bcftools merge --threads ${ncpus} ./*/*_${mgetag}.vcf.gz -Oz -o ./${mge}_merged.vcf.gz
#  bcftools merge --threads ${ncpus} ./*/*_${mgetag}.filteredDP5.vcf.gz -Oz -o ./${mge}_merged.filteredDP5.vcf.gz
#  bcftools merge --threads ${ncpus} ./*/*_${mgetag}.filteredDP10.vcf.gz -Oz -o ./${mge}_merged.filteredDP10.vcf.gz
#  bcftools merge --threads ${ncpus} ./*/*_${mgetag}.noindel.filteredDP10.vcf.gz -Oz -o ./${mge}_merged.noindel.filteredDP10.vcf.gz
#done

cd ${yemen}/read_mapping/

mkdir -p ${yemen}/read_mapping/consensus_ali/
for mge in ${MGEs} ; do
  echo $mge
  if [ "${mge}" == 'IncAC2_hybrid' ] ; then
    for lane in $(ls mapped_reads/IncAC2_hybrid/) ; do
      zcat mapped_reads/${mge}/${lane}/${lane}_mapped_${mge}.paired.noindel.consensus.fa.plasmidonly.gz
    done | gzip > consensus_ali/all_mapped_${mge}.paired.noindel.consensus.plasmidonly.aln.gz
  else
    for lane in $(ls mapped_reads/$mge/) ; do
      zcat mapped_reads/${mge}/${lane}/${lane}_mapped_${mge}.paired.noindel.consensus.fa.gz
    done | gzip > consensus_ali/all_mapped_${mge}.paired.noindel.consensus.aln.gz
  fi
done

for mge in ${MGEs} ; do
  cd ${yemen}/read_mapping/mapped_reads/${mge}/
  for lanedir in $(ls -A) ; do
    if [ -d ${yemen}/read_mapping/mapped_reads/${mge}/${lanedir} ] ; then
	  tar -czf ${lanedir}.tar.gz ${lanedir}/ && rm -r ${yemen}/read_mapping/mapped_reads/${mge}/${lanedir}/
	fi
  done
done


# regenerate the merged VCF files
for mge in ${MGEs} ; do
  mgetag=mapped_${mge}.paired
  cd ${yemen}/read_mapping/mapped_reads/${mge}/
  for lanedirtgz in $(ls -d *tar.gz) ; do
    echo ${lanedirtgz}
    lanedir=${lanedirtgz%.tar.gz}
	lanevcfgz=${lanedir}/${lanedir}_${mgetag}.noindel.filteredDP10.vcf.gz
	tar -xzf ${lanedirtgz} ${lanevcfgz} ${lanevcfgz}.csi
	ls -l ${lanevcfgz}
  done
  for lanevcfgz in $(ls ./*/*_${mgetag}.noindel.filteredDP10.vcf.gz) ; do
   if [ ! -e ${lanevcfgz}.csi ] ; then
     bcftools index --threads ${ncpus} ${lanevcfgz}
   fi
  done
  bcftools merge --threads ${ncpus} ./*/*_${mgetag}.noindel.filteredDP10.vcf.gz -Oz -o ./${mge}_merged.noindel.filteredDP10.vcf.gz
done

