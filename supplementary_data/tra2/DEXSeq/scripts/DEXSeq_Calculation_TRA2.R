#DEXSeq_Calculation_TRA2.R

library("DEXSeq")

#Load the count files created with dexseq_count.py
countFiles = list.files(path="/~/tra2/DEXSeq", full.names = TRUE)
sampleTable = data.frame (row.names = c("SAMPLE_1_REP_1","SAMPLE_1_REP_2","SAMPLE_1_REP_3","SAMPLE_2_REP_1","SAMPLE_2_REP_2","SAMPLE_2_REP_3"),
                          condition = c("KD","KD","KD","CTRL","CTRL","CTRL"),
                          libType = c("paired-end","paired-end","paired-end","paired-end","paired-end","paired-end"))    
dxd = DEXSeqDataSetFromHTSeq(countFiles,sampleData = sampleTable,design = ~ sample + exon + condition:exon)

#Normalize the samples, because of the different dpeths
dxd = estimateSizeFactors( dxd )
#estimate the dispersion estimates
dxd = estimateDispersions( dxd )
dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
dxr1 = DEXSeqResults( dxd )
jpeg("~/tra2/DEXSeq/MA_plot.jpeg")
plotMA( dxr1, cex=0.8 )
dev.off()

#Save the file
write.table(dxr1,file="~/tra2/DEXSeq/Results.txt",sep="\t",quote=FALSE)
