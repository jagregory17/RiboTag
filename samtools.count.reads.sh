#!/bin/bash

#checking for 3 arguments
if [ $# -ne 3 ];
then
    echo "Usage: $0 <SAMPLE> <alignment> <genome>(hg38/mm20/hybrid)"
    exit 1
fi


sample=$1  # sample name
alignment=$2  # directory with appropriate STAR alignment
genome=$3 #what genome; this retrives the list of chromosomes from a STAR index for the appropriate genome

## Print bash current base directory and script directory:
script_dir=$(dirname $0)
BASEDIR=$(pwd)
echo "Script directory = $script_dir"
echo "Script was exectured from = $BASEDIR"

if [ ! -d ${BASEDIR}/${sample}/${alignment}/ ];
then
    echo "+++ ERROR +++ ${sample}: $(date): alignment '${alignment}' does not exist in '${BASEDIR}'. "
    exit 1
fi


BAM=$BASEDIR/${sample}/${alignment}/${sample}_Aligned.out.sorted.bam
if [ ! -f $BAM ];then
    echo "+++ERROR+++ could not find BAM file in $BASEDIR/${sample}/${alignment}/STAR/STAR_alignment/"
    exit 1
else
    echo "found BAM file!"
fi

source /etc/profile.d/modules.sh 
module load samtools/1.7 

# set path to STAR correct STAR chrName.txt file
if [ $genome == "mm20" ]; then
    chrName=/gpfs/commons/home/jgregory/Annotations/Genomes/mm20/sequences/Star.oh124.Idx/chrName.txt
elif [ $genome == "hg38" ]; then
    chrName=/gpfs/commons/home/jgregory/Annotations/Genomes/hg38/Gencode29/sequences/Star.oh149.Idx/chrName.txt
elif [ $genome == "hybrid" ]; then
    chrName=/gpfs/commons/home/jgregory/Annotations/Genomes/mm20_gencode29_hybrid/Star.trick.oh124.Idx/chrName.txt
elif [ $genome == "Exo" ]; then
    chrName=/gpfs/commons/home/jgregory/Annotations/Genomes/ExoGenes/Star.exo.Hs_Ensembl94.Mm_Ensembl94.oh.Idx/chrName.txt
else
    echo "ERROR, genome not set"
    exit 1
fi

# count uniquely mapping reads
# used info from Alex Dobin https://groups.google.com/forum/#!topic/rna-star/ipHc9FrgaPU
while IFS='' read -r line || [[ -n "$line" ]]; do
        echo "searching for reads in chromosome $line"
        echo "$line" >> $BASEDIR/${sample}/${alignment}/${sample}_${alignment}_samtools.count.dobin.txt
        Nsingle=($(samtools view -c -q 255 -f 9 -F 4 $BAM $line))
        echo "$Nsingle single reads"
        Npaired=($(samtools view -c -q 255 -f 3 $BAM $line))
        echo "$Npaired paired reads"
        unique=$(( $Nsingle + (( $Npaired/2 )) ))
        echo "$unique uniquely mapping reads"
        echo "$unique" >> $BASEDIR/${sample}/${alignment}/${sample}_${alignment}_samtools.count.dobin.txt
done < $chrName
echo "done"
