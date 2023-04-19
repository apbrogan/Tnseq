#!/bin/bash
FolderLocation='/Users/jsher/Desktop/Generic_TN/B_scratch_v3/'
RscriptLocation='/Library/Frameworks/R.framework/Versions/3.3/Resources/bin/Rscript'

#MUST HAVE BOWTIE INDEXES, FASTA FILE AND GENELIST FOR YOUR ORGANISM
BOWTIE_NAME='BC_CP034551'
FASTA_NAME='BC_CP034551.fa'
GENELIST='BC_CP034551_genelist.txt'

#CAN CHANGE THE NUMBER OF READS TO NORMALIZE ALL SAMPLES TO
Normalized_reads='20000000'

#DON'T CHANGE
TA_SITES_LIST='TAsites'
# python2 Codes/TAfinder.py Codes/$FASTA_NAME > Codes/$TA_SITES_LIST

for SAMPLE in T0_A T10_A
do
	echo
	echo $SAMPLE
	bowtie-1.0.0/bowtie -v 3 -a --best --strata -m 1 -q --threads 7 bowtie-1.0.0/indexes/$BOWTIE_NAME Data/$SAMPLE-trimmed.fastq Data/$SAMPLE-trimmed-bowtieMap;
	python2 Codes/ReadCounter.py Codes/$TA_SITES_LIST Data/$SAMPLE-trimmed-bowtieMap Codes/$GENELIST Codes/$FASTA_NAME > Data/$SAMPLE-trimmed-bowtieMap-TA.txt
	$RscriptLocation Codes/R_normal_artemis_DNAplotter.R $FolderLocation $SAMPLE $Normalized_reads
done

#IF YOU WANT TO COMPARE SAMPLES, EDIT THE CSV "SAMPLES_TO_COMPARE.CSV" TO ADD TO THE LIST OF COMPARISIONS THE FOLLOW CODE WILL MAKE
$RscriptLocation Codes/R_ComparingSamples.R $FolderLocation Codes/$GENELIST

#IF YOU WANT TO COMPARE SAMPLES USING THE SLIDING TA WINDOW, EDIT THE CSV "SAMPLES_TO_COMPARE.CSV" TO ADD TO THE LIST OF COMPARISIONS THE FOLLOW CODE WILL MAKE
$RscriptLocation Codes/R_SlidingTAwindow.R $FolderLocation Codes/$GENELIST