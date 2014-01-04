#!/bin/bash

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
rm smrtanalysis-2.1.1-centos-6.3.run

#install samtools
echo "installing samtools"
wget http://downloads.sourceforge.net/project/samtools/samtools/0.1.19/samtools-0.1.19.tar.bz2
bunzip2 samtools-0.1.19.tar.bz2
tar -vxf samtools-0.1.19.tar
cd samtools-0.1.19
make
cd ..
rm samtools-0.1.19.tar