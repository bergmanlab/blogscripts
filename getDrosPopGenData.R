##################################################
# Script to generate SRA metadata for Drosophila #
# melaogaster population genomics projects       #
# Casey M. Bergman (University of Manchester)    #
##################################################

#install the SRAdb bioconductor package
source("http://www.bioconductor.org/biocLite.R")
biocLite("SRAdb")
library(SRAdb)
 
#download and connect to the SRA SQLlite database
sqlfile <- getSRAdbFile()
sra_connection <- dbConnect(SQLite(), "SRAmetadb.sqlite")
 
# make linking table between SRP, SRA, SRS, SRX & SRR for DPGP & DGRP projects
conversion <- sraConvert(c("SRP000694","SRP000224", "SRP005599"), sra_con = sra_connection)
 
# open outfile & print header
outFile <- file("DrosophilaPopulationGenomicData.tsv", "w")
cat("Project\tAccession\tSample\tExperiment\tRun\tSampleAlias\tSampleDescription\n", file=outFile)
 
# loop through conversion table and look-up meta-data for each sample
for (i in 1:nrow(conversion)) {
sample_metadata_sql <- paste("select sample_alias, description from sample where sample_accession like '", conversion[i,3], "'", sep="")
sample_metadata_results <- dbGetQuery(sra_connection, sample_metadata_sql)
cat(conversion[i,1], conversion[i,2], conversion[i,3], conversion[i,4], conversion[i,5], sample_metadata_results[1,1], sample_metadata_results[1,2], "\n", file=outFile, sep = "\t")
}
 
#close filehandle
close(outFile)

#clean up
system("rm SRAmetadb.sqlite")
