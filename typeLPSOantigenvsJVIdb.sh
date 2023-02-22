#!/usr/bin/env bash
module load blast/2.7.1=h96bfa4b_5-c1
module load r/3.6.0

Oantigenrefcontigs='/nfs/pathogen/cholera_club/O-antiseq_JVIserogroupRefStrainOantigenLPSloci/allOserogroups_blastdbs/whole_locus_contigs/allOserogroups.fna'

thisscript=${0}
querycontigs=${1}
outdir=${2}
nbthreads=${3}

if [ -z "${outdir}" ] ; then
  echo "Usage: typeLPSOantigenvsJVIdb.sh query_list path_output [num_threads]"
  echo "  query_list: path to file containing the list of all query genomic contig fasta files to screen, one full file path per line."
  echo "  path_output: path where to create output folderes and files"
  echo "  num_threads: threads for Blast search; default is 4 or the available number of cores ($(nproc)), whichever the lower"
  exit 1
fi
if [ -z "${nbthreads}" ] ; then
  if [ $(nproc) -lt 4 ] ; then
    nbthreads=$(nproc)
  else
    nbthreads=4
  fi
fi

scriptdir=$(dirname ${thisscript})
export ctgvsrefdir=${outdir}/$(basename ${querycontigs})_vs_OantigenJVIrefdb.blastout

mkdir -p ${ctgvsrefdir}/

for contig in $(cat ${querycontigs}) ; do
  if [ ! -z "$(echo ${contig} | grep '.gz')" ] ; then
    bncontig=$(basename ${contig%.gz})
	genomeid=${bncontig%.*}
    blastn -query <(zcat ${contig}) -db ${Oantigenrefcontigs} -num_threads ${nbthreads} \ 
	-outfmt 6 -out ${ctgvsrefdir}/${genomeid}_contigs_vs_OantigenJVIrefdb.blastout
  else
    bncontig=$(basename ${contig})
	genomeid=${bncontig%.*}
    blastn -query ${contig} -db ${Oantigenrefcontigs} -num_threads ${nbthreads} \
	-outfmt 6 -out ${ctgvsrefdir}/${genomeid}_contigs_vs_OantigenJVIrefdb.blastout
  fi
done
# parse BLAST results to find top hits
${scriptdir}/topblasthit_Oantigen.r && echo -e "best hit results listed in\n: $(ls ${ctgvsrefdir}_topnormscore.txt)"
