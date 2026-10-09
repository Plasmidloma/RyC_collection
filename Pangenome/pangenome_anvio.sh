#!/bin/bash

for bamfile in ./Anvio/Ecoli/Mapping/*.sorted.bam
	do
		samplename=$( basename $bamfile)
		samplename=${samplename::-11 } 	
		anvi-init-bam $bamfile
		anvi-profile -T 15 -i $bamfile -c ./Anvio/Ecoli/complete_contigs.db -o ./Anvio/Ecoli/profiles/Ecoli-profiledb/$samplename -S $samplename --min-contig-length 0



	done

for bamfile in ./Anvio/Klebsiella/Mapping/*.sorted.bam
	do
		samplename=$( basename $bamfile)
		samplename=${samplename::-11 } 	
		anvi-init-bam $bamfile
		anvi-profile -T 15 -i $bamfile -c ./Anvio/Klebsiella/complete_contigs.db -o ./Anvio/Klebsiella/profiles/Klebsiella-profiledb/$samplename -S $samplename ---min-contig-length 0
	done
	 
