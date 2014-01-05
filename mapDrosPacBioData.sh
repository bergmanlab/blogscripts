#!/bin/bash

#add blasr to your working environment
source install/smrtanalysis-2.1.1.128549/etc/setup.sh

#map .bax.h5 files from individual PacBio cells to dm3 using blasr
echo "mapping .bax.h5 files to dm3"
samples=()
for input in `ls Dro*/*/*/*bax.h5` 
do
 sample=`basename $input .bax.h5`
 output1="${sample}.sam"
 output2="${sample}.unaligned"
 output3="${sample}"
 qsub -l core=12 -V -b y -cwd -N blasr_${sample} "install/smrtanalysis-2.1.1.128549/analysis/bin/blasr $input dm3.fa -bestn 1 -nproc 8 -minPctIdentity 80 -sam -out $output1 -unaligned $output2"
 qsub -V -b y -cwd -N sort_${sample} -hold_jid blasr_${sample} "samtools-0.1.19/samtools view -bS $output1 | samtools-0.1.19/samtools sort - $output3"
 samples+=(sort_${sample})
done

#merge bam files from individual cells and index merged bam file
echo "merging and indexing bam files"
holdlist=`printf -- '%s,' "${samples[@]}"`
input="m131*_p0*.bam"
output="dm3PacBio.bam"
qsub -V -b y -cwd -N merge_pacbio -hold_jid $holdlist "samtools-0.1.19/samtools merge $output $input"
qsub -V -b y -cwd -N index_pacbio -hold_jid merge_pacbio "samtools-0.1.19/samtools index $output"

#merge unmapped reads and remove .sam, .bam and .unaligned files from individual cells
echo "merging unmapped reads and cleaning up"
qsub -V -b y -cwd -N merge_unmapped -hold_jid index_pacbio "cat m131*unaligned > dm3PacBio.unaligned"
qsub -V -b y -cwd -N rm_cell_files -hold_jid merge_unmapped "rm m131*unaligned m131*bam m131*sam"


