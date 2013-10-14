#!/bin/bash 

################################################
# Script to construct reference "hologenome"   #
# for D. melanogaster and associated microbes  #
# Casey M. Bergman (University of Manchester)  #
################################################

# download D. melanogaster (dm3) reference genome from UCSC genome database
# unpack dm3, remove composite mtDNA, extract mtDNA from chrU & combine into
# one .fa file excluding chrU & chrUextra

wget -q http://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/chromFa.tar.gz
tar -xvzf chromFa.tar.gz
rm chrM.fa
echo ">chrM_iso1" > chrM.fa
faFrag chrU.fa 5288527 5305749 stdout | grep -v chrU >> chrM.fa
cat chr2L.fa chr2LHet.fa chr2R.fa chr2RHet.fa chr3L.fa chr3LHet.fa chr3R.fa chr3RHet.fa chr4.fa chrX.fa chrXHet.fa chrYHet.fa chrM.fa > dm3.fa 
rm chromFa.tar.gz chr2L.fa chr2LHet.fa chr2R.fa chr2RHet.fa chr3L.fa chr3LHet.fa chr3R.fa chr3RHet.fa chr4.fa chrX.fa chrXHet.fa chrYHet.fa chrU.fa chrUextra.fa chrM.fa

# download yeast (sacCer3) reference genome from UCSC genome database
# unpack and rename sacCer3 chroms & combine into one .fa file

wget -q http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz
tar -xvzf chromFa.tar.gz
for i in chr*fa; do cat $i | sed 's/>chr/>sacCer3_chr/g' >> sacCer3.fa; done
rm chromFa.tar.gz chrI.fa chrII.fa chrIII.fa chrIV.fa chrIX.fa chrM.fa chrV.fa chrVI.fa chrVII.fa chrVIII.fa chrX.fa chrXI.fa chrXII.fa chrXIII.fa chrXIV.fa chrXV.fa chrXVI.fa

# download genomes for microbes known to be associated with D. melanogaster
# from NCBI, combine multi-contig draft assemblies into single fasta file 
# per species & rename fasta headers to human readable form

# Wolbachia pipientis
wget -q ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/Wolbachia_endosymbiont_of_Drosophila_melanogaster_uid272/AE017196.fna
echo ">w_pipientis" >> w_pipientis.fa
grep -v \> AE017196.fna >> w_pipientis.fa
rm AE017196.fna

# Pseudomonas entomophila
wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Pseudomonas_entomophila_L48_uid58639/NC_008027.fna
echo ">p_entomophila" >> p_entomophila.fa
grep -v \> NC_008027.fna >> p_entomophila.fa
rm NC_008027.fna

# Commensalibacter intestini
wget -q ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria_DRAFT/Commensalibacter_intestini_A911_uid73359/AGFR00000000.contig.fna.tgz
tar -xvzf AGFR00000000.contig.fna.tgz
cat AGFR01*fna | sed 's/>.*|/>/g' | sed 's/> />/' | awk -F, '{ print $1; }' | sed 's/ /_/g' | sed 's/>Commensalibacter_intestini_A911_74_/>c_intestini_/g' > c_intestini.fa
rm AGFR01*fna AGFR00000000.contig.fna.tgz

# Acetobacter pomorum
wget -q ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria_DRAFT/Acetobacter_pomorum_DM001_uid60787/AEUP00000000.contig.fna.tgz
tar -xvzf AEUP00000000.contig.fna.tgz
cat AEUP01*fna | sed 's/>.*|/>/g' | sed 's/> />/' | awk -F, '{ print $1; }' | sed 's/ /_/g' | sed 's/>Acetobacter_pomorum_DM001_Contig00/>a_pomorum_/g' > a_pomorum.fa
rm AEUP01*fna AEUP00000000.contig.fna.tgz

# Gluconobacter morbifer
wget -q ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria_DRAFT/Gluconobacter_morbifer_G707_uid73361/AGQV00000000.contig.fna.tgz
tar -xvzf AGQV00000000.contig.fna.tgz
cat AGQV01*fna | sed 's/>.*|/>/g' | sed 's/> />/' | awk -F, '{ print $1; }' | sed 's/ /_/g' | sed 's/>Gluconobacter_morbifer_G707_75_/>g_morbifer_/g'  > g_morbifer.fa
rm AGQV01*fna AGQV00000000.contig.fna.tgz

# Providencia burhodogranariea
wget -q ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria_DRAFT/Providencia_burhodogranariea_DSM_19968_uid82573/AKKL00000000.contig.fna.tgz
tar -xvzf AKKL00000000.contig.fna.tgz
cat AKKL01*fna | sed 's/>.*|/>/g' | sed 's/> />/' | awk -F, '{ print $1; }' | sed 's/ /_/g' | sed 's/Providencia_burhodogranariea_DSM_19968_contig000/p_burhodogranariea_/g' > p_burhodogranariea.fa
rm AKKL01*fna AKKL00000000.contig.fna.tgz

# download genomes for microbes known to be associated with D. melanogaster
# from Ensembl & rename fasta headers to human readable form

# Providencia alcalifaciens
wget -q ftp://ftp.ensemblgenomes.org/pub/release-19/bacteria/fasta/bacteria_18_collection/providencia_alcalifaciens_dmel2/dna/Providencia_alcalifaciens_dmel2.GCA_000314875.2.19.dna.toplevel.fa.gz
gunzip Providencia_alcalifaciens_dmel2.GCA_000314875.2.19.dna.toplevel.fa.gz
cat Providencia_alcalifaciens_dmel2.GCA_000314875.2.19.dna.toplevel.fa | sed 's/dna.*//g' | sed 's/>contig000/>p_alcalifaciens_/g' > p_alcalifaciens.fa
rm Providencia_alcalifaciens_dmel2.GCA_000314875.2.19.dna.toplevel.fa*

# Providencia rettgeri
wget -q ftp://ftp.ensemblgenomes.org/pub/release-19/bacteria/fasta/bacteria_23_collection/providencia_rettgeri_dmel1/dna/Providencia_rettgeri_dmel1.GCA_000314835.1.19.dna.toplevel.fa.gz
gunzip Providencia_rettgeri_dmel1.GCA_000314835.1.19.dna.toplevel.fa.gz
cat Providencia_rettgeri_dmel1.GCA_000314835.1.19.dna.toplevel.fa | sed 's/dna.*//g' | sed 's/>Contig/>p_rettgeri_/g' > p_rettgeri.fa
rm Providencia_rettgeri_dmel1.GCA_000314835.1.19.dna.toplevel.fa* 

# Enterococcus faecalis
wget -q ftp://ftp.ensemblgenomes.org/pub/release-19/bacteria/fasta/bacteria_15_collection/enterococcus_faecalis_fly1/dna/Enterococcus_faecalis_fly1.GCA_000157415.1.19.dna.toplevel.fa.gz
gunzip Enterococcus_faecalis_fly1.GCA_000157415.1.19.dna.toplevel.fa.gz
cat Enterococcus_faecalis_fly1.GCA_000157415.1.19.dna.toplevel.fa | sed 's/dna.*//g' | sed 's/>GG/>e_faecalis_GG/g' > e_faecalis.fa
rm Enterococcus_faecalis_fly1.GCA_000157415.1.19.dna.toplevel.fa*

# create D. melanogaster "hologenome" 
# from individual species fasta files
cat dm3.fa sacCer3.fa w_pipientis.fa p_entomophila.fa c_intestini.fa a_pomorum.fa g_morbifer.fa p_burhodogranariea.fa p_alcalifaciens.fa p_rettgeri.fa e_faecalis.fa > dm3_hologenome.fa
