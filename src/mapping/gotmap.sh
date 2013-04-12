#!/bin/bash
# author: Artem Tarasov
# based on: https://github.com/ngscomparison/NGS-Benchtop-Comparison/tree/master/scripts

min_args=3
execDir=$(dirname $0)

if [ $# -eq $min_args ]; then        
	tag=$1
	reads=$2
	ref=$3
	threads=16
	
	if [ ! -r ${ref}.tmap.anno ]; then
	    echo "Start indexing of reference"
	    tmap index -f ${ref}
	    echo "finish indexing of reference"	   
	else
	    echo "index file exists"
	fi

  # -o 1: output compressed BAM
  # -a 0: "returns the mapping with the best score only if all other mappings 
  #        had worse score, otherwise the read is unmapped" (from TMAP book)
  tmap mapall -f $ref -r $reads -n $threads -o 1 -s $tag.bam -a 0 -v stage1 map1 map2 map3 map4

	samtools sort $tag.bam $tag.uni.sorted
	rm $tag.bam
  samtools index $tag.uni.sorted.bam
else
  echo "Run tmap plus samtools to create a sorted mapped read file."
  echo "1: Output prefix"
  echo "2: Raw read file"
  echo "3: Reference file"
fi
