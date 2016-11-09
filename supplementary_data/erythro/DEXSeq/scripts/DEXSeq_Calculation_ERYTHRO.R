#DEXSeq_Calculation_ERYTHRO.R
source("https://bioconductor.org/biocLite.R")
biocLite("DEXSeq")
library("DEXSeq")

#Load the count files created with /genomics/users/jcarlos/dexseq/scripts/dexseq_count.py
countFiles = list.files(path="~/erythro/DEXSeq", full.names = TRUE)
sampleTable = data.frame (row.names = c("SAMPLE_1_REP_1","SAMPLE_1_REP_2","SAMPLE_1_REP_3","SAMPLE_2_REP_1","SAMPLE_2_REP_2","SAMPLE_2_REP_3"),
                          condition = c("proE","proE","proE","orthoE","orthoE","orthoE"),
                          libType = c("single-end","single-end","single-end","single-end","single-end","single-end"))    
dxd = DEXSeqDataSetFromHTSeq(countFiles,sampleData = sampleTable,design = ~ sample + exon + condition:exon)

#Normalize the samples, because of the different depths
dxd = estimateSizeFactors( dxd )
#estimate the dispersion estimates
dxd = estimateDispersions( dxd )
dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
dxr1 = DEXSeqResults( dxd )
jpeg("~/erythro/DEXSeq/MA_plot.jpeg")
plotMA( dxr1, cex=0.8 )
dev.off()

#Save the file
write.table(dxr1,file="~/erythro/DEXSeq/Results.txt",sep="\t",quote=FALSE)
