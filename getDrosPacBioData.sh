#!/bin/bash

#get compressed .tar archives of D. melanogaster PacBio runs
echo "get tar archives of PacBio runs"
qsub -l core=12 -b y -cwd -N wget_Dro1_24NOV2013_398.tgz "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro1_24NOV2013_398.tgz"
qsub -l core=12 -b y -cwd -N wget_Dro2_25NOV2013_399.tgz "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro2_25NOV2013_399.tgz"
qsub -l core=12 -b y -cwd -N wget_Dro3_26NOV2013_400.tgz "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro3_26NOV2013_400.tgz"
qsub -l core=12 -b y -cwd -N wget_Dro4_28NOV2013_401.tgz "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro4_28NOV2013_401.tgz"
qsub -l core=12 -b y -cwd -N wget_Dro5_29NOV2013_402.tgz "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro5_29NOV2013_402.tgz"
qsub -l core=12 -b y -cwd -N wget_Dro6_1DEC2013_403.tgz  "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro6_1DEC2013_403.tgz"

#get md5sums for PacBio .tar archives
echo "getting md5sums for PacBio tar archives"
qsub -b y -cwd -N wget_Dro1_24NOV2013_398.md5 "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro1_24NOV2013_398.md5"
qsub -b y -cwd -N wget_Dro2_25NOV2013_399.md5 "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro2_25NOV2013_399.md5"
qsub -b y -cwd -N wget_Dro3_26NOV2013_400.md5 "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro3_26NOV2013_400.md5"
qsub -b y -cwd -N wget_Dro4_28NOV2013_401.md5 "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro4_28NOV2013_401.md5"
qsub -b y -cwd -N wget_Dro5_29NOV2013_402.md5 "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro5_29NOV2013_402.md5"
qsub -b y -cwd -N wget_Dro6_1DEC2013_403.md5  "wget https://s3.amazonaws.com/datasets.pacb.com/2014/Drosophila/raw/Dro6_1DEC2013_403.md5"

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

#delete .tar archives and md5 files
qsub -V -b y -cwd -N rm_tar_md -hold_jid tar_Dro1_24NOV2013_398.tgz,tar_Dro2_25NOV2013_399.tgz,tar_Dro3_26NOV2013_400.tgz,tar_Dro4_28NOV2013_401.tgz,tar_Dro5_29NOV2013_402.tgz,tar_Dro6_1DEC2013_403.tgz "rm Dro*.tgz Dro*.md5"
