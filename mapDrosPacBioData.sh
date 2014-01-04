#!/bin/bash

#map .bax.h5 files from individual PacBio cells to dm3 using blasr
echo "mapping .bax.h5 files to dm3"
samples=()
for input in `ls Dro*/*/*/*bax.h5` 
do
 sample=`basename $input .bax.h5`
 output1="${sample}.sam"
 output2="${sample}.unaligned"
 output3="${sample}"
 qsub -V -b y -cwd -N blasr_${sample} -hold_jid tar_Dro1_24NOV2013_398.tgz,tar_Dro2_25NOV2013_399.tgz,tar_Dro3_26NOV2013_400.tgz,tar_Dro4_28NOV2013_401.tgz,tar_Dro5_29NOV2013_402.tgz,tar_Dro6_1DEC2013_403.tgz "install/smrtanalysis-2.1.1.128549/analysis/bin/blasr $input dm3.fa -bestn 1 -nproc 8 -minPctIdentity 80 -sam -out $output1 -unaligned $output2"
 qsub -V -b y -cwd -N sortpb_${sample} -hold_jid blasr_${sample} "samtools-0.1.19/samtools view -bS $output1 | samtools-0.1.19/samtools sort - $output3"
 samples+=(sortpb_${sample})
done

#merge bam files from individual cells and index merged bam file
echo "merging and indexing bam files"
holdlist=`printf -- '%s,' "${samples[@]}"`
input="$outputdir/m131*_p0*.bam"
output="$outputdir/dm3PacBio.bam"
qsub -V -b y -cwd -N merge_pacbio -hold_jid $holdlist "samtools-0.1.19/samtools merge $output $input"
qsub -V -b y -cwd -N index_pacbio -hold_jid merge_pacbio "samtools-0.1.19/samtools index $output"