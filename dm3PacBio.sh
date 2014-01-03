#!/bin/bash

#get compressed .tar archives of D. melanogaster PacBio runs
echo "get tar archives of PacBio runs"
qsub -b y -cwd -N wget_Dro1_24NOV2013_398.tgz "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/Dro1_24NOV2013_398.tgz"
qsub -b y -cwd -N wget_Dro2_25NOV2013_399.tgz "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/Dro2_25NOV2013_399.tgz"
qsub -b y -cwd -N wget_Dro3_26NOV2013_400.tgz "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/Dro3_26NOV2013_400.tgz"
qsub -b y -cwd -N wget_Dro4_28NOV2013_401.tgz "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/Dro4_28NOV2013_401.tgz"
qsub -b y -cwd -N wget_Dro5_29NOV2013_402.tgz "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/Dro5_29NOV2013_402.tgz"
qsub -b y -cwd -N wget_Dro6_1DEC2013_403.tgz  "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/Dro6_1DEC2013_403.tgz"

#get md5sums for PacBio .tar archives
echo "getting md5sums for PacBio tar archives"
qsub -b y -cwd -N wget_Dro1_24NOV2013_398.md5 "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/Dro1_24NOV2013_398.md5"
qsub -b y -cwd -N wget_Dro2_25NOV2013_399.md5 "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/Dro2_25NOV2013_399.md5"
qsub -b y -cwd -N wget_Dro3_26NOV2013_400.md5 "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/Dro3_26NOV2013_400.md5"
qsub -b y -cwd -N wget_Dro4_28NOV2013_401.md5 "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/Dro4_28NOV2013_401.md5"
qsub -b y -cwd -N wget_Dro5_29NOV2013_402.md5 "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/Dro5_29NOV2013_402.md5"
qsub -b y -cwd -N wget_Dro6_1DEC2013_403.md5  "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/Dro6_1DEC2013_403.md5"

#verify .tar archives have been downloaded properly
echo "verifying archives have been downloaded properly"
qsub -b y -cwd -N md5_Dro1_24NOV2013_398.md5 -hold_jid wget_Dro1_24NOV2013_398.tgz,wget_Dro1_24NOV2013_398.md5 "md5sum -c Dro1_24NOV2013_398.md5 > Dro1_24NOV2013_398.md5.checkresult"
qsub -b y -cwd -N md5_Dro2_25NOV2013_399.md5 -hold_jid wget_Dro2_25NOV2013_399.tgz,wget_Dro2_25NOV2013_399.md5 "md5sum -c Dro2_25NOV2013_399.md5 > Dro2_25NOV2013_399.md5.checkresult"
qsub -b y -cwd -N md5_Dro3_26NOV2013_400.md5 -hold_jid wget_Dro3_26NOV2013_400.tgz,wget_Dro3_26NOV2013_400.md5 "md5sum -c Dro3_26NOV2013_400.md5 > Dro3_26NOV2013_400.md5.checkresult"
qsub -b y -cwd -N md5_Dro4_28NOV2013_401.md5 -hold_jid wget_Dro4_28NOV2013_401.tgz,wget_Dro4_28NOV2013_401.md5 "md5sum -c Dro4_28NOV2013_401.md5 > Dro4_28NOV2013_401.md5.checkresult"
qsub -b y -cwd -N md5_Dro5_29NOV2013_402.md5 -hold_jid wget_Dro5_29NOV2013_402.tgz,wget_Dro5_29NOV2013_402.md5 "md5sum -c Dro5_29NOV2013_402.md5 > Dro5_29NOV2013_402.md5.checkresult"
qsub -b y -cwd -N md5_Dro6_1DEC2013_403.md5  -hold_jid wget_Dro6_1DEC2013_403.tgz,wget_Dro6_1DEC2013_403.md5  "md5sum -c Dro6_1DEC2013_403.md5  > Dro6_1DEC2013_403.md5.checkresult"

#extract PacBio runs
echo "extracting PacBio runs"
qsub -b y -cwd -N tar_Dro1_24NOV2013_398.tgz -hold_jid wget_Dro1_24NOV2013_398.tgz "tar -xvzf Dro1_24NOV2013_398.tgz"
qsub -b y -cwd -N tar_Dro2_25NOV2013_399.tgz -hold_jid wget_Dro2_25NOV2013_399.tgz "tar -xvzf Dro2_25NOV2013_399.tgz"
qsub -b y -cwd -N tar_Dro3_26NOV2013_400.tgz -hold_jid wget_Dro3_26NOV2013_400.tgz "tar -xvzf Dro3_26NOV2013_400.tgz"
qsub -b y -cwd -N tar_Dro4_28NOV2013_401.tgz -hold_jid wget_Dro4_28NOV2013_401.tgz "tar -xvzf Dro4_28NOV2013_401.tgz"
qsub -b y -cwd -N tar_Dro5_29NOV2013_402.tgz -hold_jid wget_Dro5_29NOV2013_402.tgz "tar -xvzf Dro5_29NOV2013_402.tgz"
qsub -b y -cwd -N tar_Dro6_1DEC2013_403.tgz  -hold_jid wget_Dro6_1DEC2013_403.tgz  "tar -xvzf Dro6_1DEC2013_403.tgz"

# get dm3 and reformat into single fasta file
echo "getting dm3 and reformatting into single fasta file"
wget -q http://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/chromFa.tar.gz
tar -xvzf chromFa.tar.gz
cat chr*.fa > dm3.fa
rm chromFa.tar.gz chr*.fa

#install PacBio smrtanalysis suite (including blasr PacBio long-read mapping engine)
echo "installing PacBio smrtanalysis suite"
wget http://programs.pacificbiosciences.com/l/1652/2013-11-05/2tqk4f
bash smrtanalysis-2.1.1-centos-6.3.run --extract-only --rootdir ./
source install/smrtanalysis-2.1.1.128549/etc/setup.sh

#install samtools
echo "installing samtools"
wget http://downloads.sourceforge.net/project/samtools/samtools/0.1.19/samtools-0.1.19.tar.bz2
bunzip2 samtools-0.1.19.tar.bz2
tar -vxf samtools-0.1.19.tar
cd samtools-0.1.19
make
cd ..

#map .bax.h5 files from individual PacBio cells to dm3 using blasr
echo "mapping .bax.h5 files to dm3"
samples=()
for input in `ls Dro*/*/*/*bax.h5` 
do
 input=`basename $input .bax.h5`
 output1="${sample}.sam"
 output2="${sample}.unaligned"
 output3="${sample}"
 qsub -V -b y -cwd -N blasr_${input} -hold_jid tar_Dro1_24NOV2013_398.tgz,tar_Dro2_25NOV2013_399.tgz,tar_Dro3_26NOV2013_400.tgz,tar_Dro4_28NOV2013_401.tgz,tar_Dro5_29NOV2013_402.tgz,tar_Dro6_1DEC2013_403.tgz "install/smrtanalysis-2.1.1.128549/analysis/bin/blasr $input dm3.fa -bestn 1 -nproc 8 -minPctIdentity 80 -sam -out $output1 -unaligned $output2"
 qsub -V -b y -cwd -N sortpb_${input} -hold_jid blasr_${input} "samtools-0.1.19/samtools view -bS $output1 | samtools-0.1.19/samtools sort - $output3"
 samples+=(sortpb_${input})
done

#merge bam files from individual cells and index merged bam file
echo "merging and indexing bam files"
holdlist=`printf -- '%s,' "${samples[@]}"`
input="$outputdir/m131*_p0*.bam"
output="$outputdir/dm3PacBio.bam"
qsub -V -b y -cwd -N merge_pacbio -hold_jid $holdlist "samtools-0.1.19/samtools merge $output $input"
qsub -V -b y -cwd -N index_pacbio -hold_jid merge_pacbio "samtools-0.1.19/samtools index $output"