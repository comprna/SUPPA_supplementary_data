#get_positive_negatives_events.R

#We need to get those genes that has just 2 trancripts. 

#Load the refseq annotation
# the refseq gtf available in UCSC changes from time to time. Be sure to choose the GTF provided
# This is the one that we uploaded to the GitHub
gtf_original <- read.table(file="~/SUPPA_supplementary_data/annotation/refseq_hg19.formatted.gtf")
refseq_ids <- unique(gtf_original[c("V10","V13")])
refseq_table <- table(refseq_ids$V10)
just_2_transcripts <- refseq_table[which(refseq_table==2)]
#Get the transcripts associated to these genes
refseq_ids2 <- refseq_ids[which(as.character(refseq_ids$V10)%in%rownames(just_2_transcripts)),]
refseq_ids2_sorted <- refseq_ids2[order(refseq_ids2[,1]),]
colnames(refseq_ids2_sorted) <- c("gene_id","transcript_id")

#2. Get the events of SUPPA annotation that are associated to these genes
#Load the event annotation 
events_annotation <- read.table(file="~/SUPPA_supplementary_data/annotation/refseq_hg19.events_formatted.pooled.ioe",sep="\t",header = TRUE)
unique_genes <- as.data.frame(unique(as.character(refseq_ids2_sorted$gene_id)))
colnames(unique_genes) <- c("gene_id")
events_genes_annotation <- merge(unique_genes,events_annotation,by="gene_id")
events_annotation_filtered <- events_genes_annotation[order(events_genes_annotation[,1]),]
#3.188 events in the genes with just 2 transcripts

#############
#   RSEM    #
#############

#SAMPLE1
TPM1 <- read.table(file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513329.isoforms.results",
                   sep="\t",header = TRUE)
#format the transcript id
id <- substr(as.character(TPM1$transcript_id),14,nchar(as.character(TPM1$transcript_id)))
TPM1_f <- cbind(id,TPM1$TPM)
rownames(TPM1_f) <- TPM1_f[,1]
TPM1_f <- as.data.frame(TPM1_f[,-1])
colnames(TPM1_f) <- "SRR1513329"

PSI1 <- read.table(file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513329.psi",
                   sep="\t",header = TRUE)


#Get all the TPMs associated to the events_annotation_filtered and generate another table
# with the transcripts, the TPM inverted and the associated event
results <- matrix(data=0,nrow=nrow(events_annotation_filtered),ncol=5)
#DeltaPSI sera la diferencia entre PSI_switched - PSI_original
colnames(results) <- c("tpm1","tpm2","PSI_original","PSI_switched","DeltaPSI")
tpms_inverted1 <- matrix(data=0,nrow=250000,ncol=4)
colnames(tpms_inverted1) <- c("transcript_id","orignal_tpm","inverted_tpm","event")
cont <- 1
i <- 1
for(i in 1:nrow(events_annotation_filtered)){
  # print(i)
  transcripts <- unlist(strsplit(as.character(events_annotation_filtered[i,5]),","))
  results[i,1] <- as.numeric(as.character(TPM1_f[transcripts[1],]))
  results[i,2] <- as.numeric(as.character(TPM1_f[transcripts[2],]))
  if(as.character(events_annotation_filtered[i,4])==transcripts[1]){
    results[i,3] <- as.numeric(as.character(TPM1_f[transcripts[1],])) / 
      (as.numeric(as.character(TPM1_f[transcripts[1],])) + as.numeric(as.character(TPM1_f[transcripts[2],])))
    results[i,4] <- as.numeric(as.character(TPM1_f[transcripts[2],])) / 
      (as.numeric(as.character(TPM1_f[transcripts[1],])) + as.numeric(as.character(TPM1_f[transcripts[2],])))
  }
  else{
    results[i,3] <- as.numeric(as.character(TPM1_f[transcripts[2],])) / 
      (as.numeric(as.character(TPM1_f[transcripts[1],])) + as.numeric(as.character(TPM1_f[transcripts[2],])))
    results[i,4] <- as.numeric(as.character(TPM1_f[transcripts[1],])) / 
      (as.numeric(as.character(TPM1_f[transcripts[1],])) + as.numeric(as.character(TPM1_f[transcripts[2],])))
  }
  results[i,5] <- results[i,4] - results[i,3]
  
  #Invert the TPMs
  tpms_inverted1[cont,1] <- transcripts[1]
  tpms_inverted1[cont,2] <- as.numeric(as.character(TPM1_f[transcripts[1],]))
  tpms_inverted1[cont,3] <- as.numeric(as.character(TPM1_f[transcripts[2],]))
  tpms_inverted1[cont,4] <- as.character(events_annotation_filtered[i,3])
  cont <- cont + 1 
  tpms_inverted1[cont,1] <- transcripts[2]
  tpms_inverted1[cont,2] <- as.numeric(as.character(TPM1_f[transcripts[2],]))
  tpms_inverted1[cont,3] <- as.numeric(as.character(TPM1_f[transcripts[1],]))
  tpms_inverted1[cont,4] <- as.character(events_annotation_filtered[i,3])
  cont <- cont + 1 
}
rownames(results) <- as.character(events_annotation_filtered$event_id)
PSI1_filtered <- PSI1[which(PSI1$event_id%in%as.character(events_annotation_filtered$event_id)),]
PSI1_filtered2 <- merge(PSI1_filtered,results,by.x="event_id",by.y="row.names")  
#Remove the empty rows on tpms_inverted1
tpms_inverted1_f <- tpms_inverted1[which(tpms_inverted1[,1]!="0"),]


#SAMPLE2
TPM2 <- read.table(file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513330.isoforms.results",
                   sep="\t",header = TRUE)
#format the transcript id
id <- substr(as.character(TPM2$transcript_id),14,nchar(as.character(TPM2$transcript_id)))
TPM2_f <- cbind(id,TPM2$TPM)
rownames(TPM2_f) <- TPM2_f[,1]
TPM2_f <- as.data.frame(TPM2_f[,-1])
colnames(TPM2_f) <- "SRR1513330"

PSI2 <- read.table(file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513330.psi",
                   sep="\t",header = TRUE)


#Get all the TPMs associated to the events_annotation_filtered and generate another table
# with the transcripts, the TPM inverted and the associated event
results <- matrix(data=0,nrow=nrow(events_annotation_filtered),ncol=5)
#DeltaPSI sera la diferencia entre PSI_switched - PSI_original
colnames(results) <- c("tpm1","tpm2","PSI_original","PSI_switched","DeltaPSI")
tpms_inverted2 <- matrix(data=0,nrow=250000,ncol=4)
colnames(tpms_inverted2) <- c("transcript_id","orignal_tpm","inverted_tpm","event")
cont <- 1
i <- 2
for(i in 1:nrow(events_annotation_filtered)){
  # print(i)
  transcripts <- unlist(strsplit(as.character(events_annotation_filtered[i,5]),","))
  results[i,1] <- as.numeric(as.character(TPM2_f[transcripts[1],]))
  results[i,2] <- as.numeric(as.character(TPM2_f[transcripts[2],]))
  if(as.character(events_annotation_filtered[i,4])==transcripts[1]){
    results[i,3] <- as.numeric(as.character(TPM2_f[transcripts[1],])) / 
      (as.numeric(as.character(TPM2_f[transcripts[1],])) + as.numeric(as.character(TPM2_f[transcripts[2],])))
    results[i,4] <- as.numeric(as.character(TPM2_f[transcripts[2],])) / 
      (as.numeric(as.character(TPM2_f[transcripts[1],])) + as.numeric(as.character(TPM2_f[transcripts[2],])))
  }
  else{
    results[i,3] <- as.numeric(as.character(TPM2_f[transcripts[2],])) / 
      (as.numeric(as.character(TPM2_f[transcripts[1],])) + as.numeric(as.character(TPM2_f[transcripts[2],])))
    results[i,4] <- as.numeric(as.character(TPM2_f[transcripts[1],])) / 
      (as.numeric(as.character(TPM2_f[transcripts[1],])) + as.numeric(as.character(TPM2_f[transcripts[2],])))
  }
  results[i,5] <- results[i,4] - results[i,3]
  
  #Invert the TPMs
  tpms_inverted2[cont,1] <- transcripts[1]
  tpms_inverted2[cont,2] <- as.numeric(as.character(TPM2_f[transcripts[1],]))
  tpms_inverted2[cont,3] <- as.numeric(as.character(TPM2_f[transcripts[2],]))
  tpms_inverted2[cont,4] <- as.character(events_annotation_filtered[i,3])
  cont <- cont + 1 
  tpms_inverted2[cont,1] <- transcripts[2]
  tpms_inverted2[cont,2] <- as.numeric(as.character(TPM2_f[transcripts[2],]))
  tpms_inverted2[cont,3] <- as.numeric(as.character(TPM2_f[transcripts[1],]))
  tpms_inverted2[cont,4] <- as.character(events_annotation_filtered[i,3])
  cont <- cont + 1 
}
rownames(results) <- as.character(events_annotation_filtered$event_id)
PSI2_filtered <- PSI2[which(PSI2$event_id%in%as.character(events_annotation_filtered$event_id)),]
PSI2_filtered2 <- merge(PSI2_filtered,results,by.x="event_id",by.y="row.names")  
#Remove the empty rows on tpms_inverted2
tpms_inverted2_f <- tpms_inverted2[which(tpms_inverted2[,1]!="0"),]


#SAMPLE3
TPM3 <- read.table(file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513331.isoforms.results",
                   sep="\t",header = TRUE)
#format the transcript id
id <- substr(as.character(TPM3$transcript_id),14,nchar(as.character(TPM3$transcript_id)))
TPM3_f <- cbind(id,TPM3$TPM)
rownames(TPM3_f) <- TPM3_f[,1]
TPM3_f <- as.data.frame(TPM3_f[,-1])
colnames(TPM3_f) <- "SRR1513331"

PSI3 <- read.table(file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513331.psi",
                   sep="\t",header = TRUE)

#Get all the TPMs associated to the events_annotation_filtered and generate another table
# with the transcripts, the TPM inverted and the associated event
results <- matrix(data=0,nrow=nrow(events_annotation_filtered),ncol=5)
#DeltaPSI sera la diferencia entre PSI_switched - PSI_original
colnames(results) <- c("tpm1","tpm2","PSI_original","PSI_switched","DeltaPSI")
tpms_inverted3 <- matrix(data=0,nrow=250000,ncol=4)
colnames(tpms_inverted3) <- c("transcript_id","orignal_tpm","inverted_tpm","event")
cont <- 1
i <- 2
for(i in 1:nrow(events_annotation_filtered)){
  # print(i)
  transcripts <- unlist(strsplit(as.character(events_annotation_filtered[i,5]),","))
  results[i,1] <- as.numeric(as.character(TPM3_f[transcripts[1],]))
  results[i,2] <- as.numeric(as.character(TPM3_f[transcripts[2],]))
  if(as.character(events_annotation_filtered[i,4])==transcripts[1]){
    results[i,3] <- as.numeric(as.character(TPM3_f[transcripts[1],])) / 
      (as.numeric(as.character(TPM3_f[transcripts[1],])) + as.numeric(as.character(TPM3_f[transcripts[2],])))
    results[i,4] <- as.numeric(as.character(TPM3_f[transcripts[2],])) / 
      (as.numeric(as.character(TPM3_f[transcripts[1],])) + as.numeric(as.character(TPM3_f[transcripts[2],])))
  }
  else{
    results[i,3] <- as.numeric(as.character(TPM3_f[transcripts[2],])) / 
      (as.numeric(as.character(TPM3_f[transcripts[1],])) + as.numeric(as.character(TPM3_f[transcripts[2],])))
    results[i,4] <- as.numeric(as.character(TPM3_f[transcripts[1],])) / 
      (as.numeric(as.character(TPM3_f[transcripts[1],])) + as.numeric(as.character(TPM3_f[transcripts[2],])))
  }
  results[i,5] <- results[i,4] - results[i,3]
  
  #Invert the TPMs
  tpms_inverted3[cont,1] <- transcripts[1]
  tpms_inverted3[cont,2] <- as.numeric(as.character(TPM3_f[transcripts[1],]))
  tpms_inverted3[cont,3] <- as.numeric(as.character(TPM3_f[transcripts[2],]))
  tpms_inverted3[cont,4] <- as.character(events_annotation_filtered[i,3])
  cont <- cont + 1 
  tpms_inverted3[cont,1] <- transcripts[2]
  tpms_inverted3[cont,2] <- as.numeric(as.character(TPM3_f[transcripts[2],]))
  tpms_inverted3[cont,3] <- as.numeric(as.character(TPM3_f[transcripts[1],]))
  tpms_inverted3[cont,4] <- as.character(events_annotation_filtered[i,3])
  cont <- cont + 1 
}
rownames(results) <- as.character(events_annotation_filtered$event_id)
PSI3_filtered <- PSI3[which(PSI3$event_id%in%as.character(events_annotation_filtered$event_id)),]
PSI3_filtered2 <- merge(PSI3_filtered,results,by.x="event_id",by.y="row.names")  
#Remove the empty rows on tpms_inverted3
tpms_inverted3_f <- tpms_inverted3[which(tpms_inverted3[,1]!="0"),]

#Merge all the previous files
PSI_all1 <- merge(PSI1_filtered2,PSI2_filtered2,by="event_id",suffixes = c("_sample1","_sample2"))
PSI_all2 <- merge(PSI_all1,PSI3_filtered2,by="event_id",suffixes = c("","_sample3"))

#Filter out all the rows with -1 in one of the samples
PSI_all3 <- PSI_all2[which(PSI_all2$SRR1513329_Sample1!=-1 & PSI_all2$SRR1513330_Sample1!=-1 & PSI_all2$SRR1513331_Sample1!=-1),]
#2124 events

#Mark those events that have at least a mean DeltaPSI of 0.2
meanDeltaPSI <- apply(PSI_all3,1,function(x)mean(as.numeric(x[c(7,13,19)])))
significant <- unlist(lapply(meanDeltaPSI,function(x)ifelse(abs(x)>=0.2,TRUE,FALSE)))
table(significant)
# FALSE  TRUE 
# 225  1898 
gene <- unlist(lapply(PSI_all3$event_id,function(x)strsplit(as.character(x),";")[[1]][1]))
PSI_all4 <- cbind(PSI_all3,gene,meanDeltaPSI,significant)

#########
# A5/A3 #
#########

#How many of this are A5/A3?
event_type <- unlist(lapply(unlist(lapply(as.character(PSI_all4$event_id),function(x)strsplit(x,":")[[1]][1])),function(x)strsplit(x,";")[[1]][2]))
PSI_all4_f <- cbind(PSI_all4,event_type)
PSI_all4_A5_A3 <- PSI_all4_f[grep("A5|A3",PSI_all4_f$event_type,perl = TRUE),]
#663 events
#How many are significant?
table(PSI_all4_A5_A3$significant)
# FALSE  TRUE 
# 83   580 

event_type <- unlist(lapply(unlist(lapply(as.character(PSI_all4_A5_A3$event_id),function(x)strsplit(x,":")[[1]][1])),function(x)strsplit(x,";")[[1]][2]))
table(event_type)
#We are gonna take randomly a number of events (318) from the set of significant as the positive set. Only for this we will switch
# their tpms. All the rest (bw significant and not significat) will be the negative ones, we wont switch their tpms
PSI_all4_A5_A3_TRUE <- PSI_all4_A5_A3[which(PSI_all4_A5_A3$significant=="TRUE"),]

#318 for switching
sampling <- sample(1:nrow(PSI_all4_A5_A3_TRUE), 318, replace=FALSE)
positives <- PSI_all4_A5_A3_TRUE[sampling,]

#318 for non switching. We will generate a number at a time and check if we haven't already taken
# that position for the positives
cont <- 0
sampling2 <- integer(0)
negatives <- NULL
while(cont<318){
  # print(cont)
  random <- sample(1:nrow(PSI_all4_A5_A3), 1, replace=FALSE)
  random_event <- as.character(PSI_all4_A5_A3[random,1])
  random_gene <- strsplit(random_event,";")[[1]][1]
  # If the selected event is already in the positive set or the selected event has been already chosen or
  # the new event is in an event included in the positive set or there is already another event previously 
  # chosen in the same gene, generate a new one
  if(random_event %in% as.character(positives$event_id) | random %in% sampling2 | 
     random_gene %in% as.character(positives$gene) | random_gene %in% as.character(negatives$gene)){
    next   
  }
  #If dont, add it to the negative
  sampling2 <- c(sampling2,random)
  cont <- cont + 1
  negatives <- PSI_all4_A5_A3[sampling2,]
}
negatives <- PSI_all4_A5_A3[sampling2,]

#Save this positive and negative events
write.table(positives,file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/positives_A5_A3.txt",sep="\t",quote=FALSE,row.names = FALSE)
write.table(negatives,file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/negatives_A5_A3.txt",sep="\t",quote=FALSE,row.names = FALSE)

# Switch events
#Take the TPM files from RSEM and switch the TPMs for the 318 positive events (we'll use this switched TPM files for 
# running the simulation step) 
tpms_inverted1_f2 <- as.data.frame(tpms_inverted1_f[which(as.character(tpms_inverted1_f[,4])%in%as.character(positives$event_id)),])
tpms_inverted2_f2 <- as.data.frame(tpms_inverted2_f[which(as.character(tpms_inverted2_f[,4])%in%as.character(positives$event_id)),])
tpms_inverted3_f2 <- as.data.frame(tpms_inverted3_f[which(as.character(tpms_inverted3_f[,4])%in%as.character(positives$event_id)),])
#Format this files and save them
TPM1_RSEM <- read.table(file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513329.isoforms.results",
                        sep="\t",header = TRUE)
#format the transcript id
id <- substr(as.character(TPM1_RSEM$transcript_id),14,nchar(as.character(TPM1_RSEM$transcript_id)))
rownames(TPM1_RSEM) <- id
TPM2_RSEM <- read.table(file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513330.isoforms.results",
                        sep="\t",header = TRUE)
#format the transcript id
id <- substr(as.character(TPM2_RSEM$transcript_id),14,nchar(as.character(TPM2_RSEM$transcript_id)))
rownames(TPM2_RSEM) <- id
TPM3_RSEM <- read.table(file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513331.isoforms.results",
                        sep="\t",header = TRUE)
#format the transcript id
id <- substr(as.character(TPM3_RSEM$transcript_id),14,nchar(as.character(TPM3_RSEM$transcript_id)))
rownames(TPM3_RSEM) <- id

TPM1_inverted <- as.matrix(TPM1_RSEM)
TPM2_inverted <- as.matrix(TPM2_RSEM)
TPM3_inverted <- as.matrix(TPM3_RSEM)

for(i in 1:nrow(tpms_inverted1_f2)){
  # print(i)
  TPM1_inverted[as.character(tpms_inverted1_f2[i,1]),6] <- as.numeric(as.character(tpms_inverted1_f2[i,3]))
  TPM2_inverted[as.character(tpms_inverted2_f2[i,1]),6] <- as.numeric(as.character(tpms_inverted2_f2[i,3]))
  TPM3_inverted[as.character(tpms_inverted3_f2[i,1]),6] <- as.numeric(as.character(tpms_inverted3_f2[i,3]))
}
#Transform the numeric columns
TPM1_inverted_df <- as.data.frame(TPM1_inverted)
TPM1_inverted_df$length <- as.numeric(as.character(TPM1_inverted_df$length))
TPM1_inverted_df$effective_length <- as.numeric(as.character(TPM1_inverted_df$effective_length))
TPM1_inverted_df$expected_count <- as.numeric(as.character(TPM1_inverted_df$expected_count))
TPM1_inverted_df$TPM <- as.numeric(as.character(TPM1_inverted_df$TPM))
TPM1_inverted_df$FPKM <- as.numeric(as.character(TPM1_inverted_df$FPKM))
TPM1_inverted_df$IsoPct <- as.numeric(as.character(TPM1_inverted_df$IsoPct))

TPM2_inverted_df <- as.data.frame(TPM2_inverted)
TPM2_inverted_df$length <- as.numeric(as.character(TPM2_inverted_df$length))
TPM2_inverted_df$effective_length <- as.numeric(as.character(TPM2_inverted_df$effective_length))
TPM2_inverted_df$expected_count <- as.numeric(as.character(TPM2_inverted_df$expected_count))
TPM2_inverted_df$TPM <- as.numeric(as.character(TPM2_inverted_df$TPM))
TPM2_inverted_df$FPKM <- as.numeric(as.character(TPM2_inverted_df$FPKM))
TPM2_inverted_df$IsoPct <- as.numeric(as.character(TPM2_inverted_df$IsoPct))

TPM3_inverted_df <- as.data.frame(TPM3_inverted)
TPM3_inverted_df$length <- as.numeric(as.character(TPM3_inverted_df$length))
TPM3_inverted_df$effective_length <- as.numeric(as.character(TPM3_inverted_df$effective_length))
TPM3_inverted_df$expected_count <- as.numeric(as.character(TPM3_inverted_df$expected_count))
TPM3_inverted_df$TPM <- as.numeric(as.character(TPM3_inverted_df$TPM))
TPM3_inverted_df$FPKM <- as.numeric(as.character(TPM3_inverted_df$FPKM))
TPM3_inverted_df$IsoPct <- as.numeric(as.character(TPM3_inverted_df$IsoPct))

write.table(TPM1_inverted_df,file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513329_switched_A5_A3.tpm",sep="\t",quote=FALSE,row.names = FALSE)
write.table(TPM2_inverted_df,file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513330_switched_A5_A3.tpm",sep="\t",quote=FALSE,row.names = FALSE)
write.table(TPM3_inverted_df,file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513331_switched_A5_A3.tpm",sep="\t",quote=FALSE,row.names = FALSE)

######
# SE #
######

#How many of this are SE?
event_type <- unlist(lapply(unlist(lapply(as.character(PSI_all4$event_id),function(x)strsplit(x,":")[[1]][1])),function(x)strsplit(x,";")[[1]][2]))
PSI_all4_f <- cbind(PSI_all4,event_type)
PSI_all4_SE <- PSI_all4_f[grep("SE",PSI_all4_f$event_type,perl = TRUE),]
#935 events
#How many are significant?
table(PSI_all4_SE$significant)
# FALSE  TRUE 
# 84   850 

event_type <- unlist(lapply(unlist(lapply(as.character(PSI_all4_SE$event_id),function(x)strsplit(x,":")[[1]][1])),function(x)strsplit(x,";")[[1]][2]))
#We are gonna take randomly a number of events (277) from the set of significant as the positive set. Only for this we will switch
# their tpms. All the rest (bw significant and not significat) will be the negative ones, we wont switch their tpms
PSI_all4_SE_TRUE <- PSI_all4_SE[which(PSI_all4_SE$significant=="TRUE"),]

#277 for switching
sampling <- sample(1:nrow(PSI_all4_SE_TRUE), 277, replace=FALSE)
positives <- PSI_all4_SE_TRUE[sampling,]
#277 for non switching. We will generate a number at a time and check if we haven't already taken
# that position for the positives
cont <- 0
sampling2 <- integer(0)
while(cont<277){
  # print(cont)
  random <- sample(1:nrow(PSI_all4_A5_A3), 1, replace=FALSE)
  random_event <- as.character(PSI_all4_A5_A3[random,1])
  random_gene <- strsplit(random_event,";")[[1]][1]
  # If the selected event is already in the positive set or the selected event has been already chosen or
  # the new event is in an event included in the positive set or there is already another event previously 
  # chosen in the same gene, generate a new one
  if(random_event %in% as.character(positives$event_id) | random %in% sampling2 | 
     random_gene %in% as.character(positives$gene) | random_gene %in% as.character(negatives$gene)){
    next   
  }
  #If dont, add it to the negative
  sampling2 <- c(sampling2,random)
  cont <- cont + 1
  negatives <- PSI_all4_A5_A3[sampling2,]
}
negatives <- PSI_all4_A5_A3[sampling2,]

#Save this positive and negative events
write.table(positives,file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/positives_SE.txt",sep="\t",quote=FALSE,row.names = FALSE)
write.table(negatives,file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/negatives_SE.txt",sep="\t",quote=FALSE,row.names = FALSE)

# switched events 
#Take the TPM files from RSEM and switch the TPMs for the 318 positive events (we'll use this switched TPM files for 
# running the simulation step) 
tpms_inverted1_f2 <- as.data.frame(tpms_inverted1_f[which(as.character(tpms_inverted1_f[,4])%in%as.character(positives$event_id)),])
tpms_inverted2_f2 <- as.data.frame(tpms_inverted2_f[which(as.character(tpms_inverted2_f[,4])%in%as.character(positives$event_id)),])
tpms_inverted3_f2 <- as.data.frame(tpms_inverted3_f[which(as.character(tpms_inverted3_f[,4])%in%as.character(positives$event_id)),])
#Format this files and save them
TPM1_RSEM <- read.table(file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513329.isoforms.results",
                        sep="\t",header = TRUE)
#format the transcript id
id <- substr(as.character(TPM1_RSEM$transcript_id),14,nchar(as.character(TPM1_RSEM$transcript_id)))
rownames(TPM1_RSEM) <- id
TPM2_RSEM <- read.table(file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513330.isoforms.results",
                        sep="\t",header = TRUE)
#format the transcript id
id <- substr(as.character(TPM2_RSEM$transcript_id),14,nchar(as.character(TPM2_RSEM$transcript_id)))
rownames(TPM2_RSEM) <- id
TPM3_RSEM <- read.table(file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513331.isoforms.results",
                        sep="\t",header = TRUE)
#format the transcript id
id <- substr(as.character(TPM3_RSEM$transcript_id),14,nchar(as.character(TPM3_RSEM$transcript_id)))
rownames(TPM3_RSEM) <- id

TPM1_inverted <- as.matrix(TPM1_RSEM)
TPM2_inverted <- as.matrix(TPM2_RSEM)
TPM3_inverted <- as.matrix(TPM3_RSEM)

for(i in 1:nrow(tpms_inverted1_f2)){
  # print(i)
  TPM1_inverted[as.character(tpms_inverted1_f2[i,1]),6] <- as.numeric(as.character(tpms_inverted1_f2[i,3]))
  TPM2_inverted[as.character(tpms_inverted2_f2[i,1]),6] <- as.numeric(as.character(tpms_inverted2_f2[i,3]))
  TPM3_inverted[as.character(tpms_inverted3_f2[i,1]),6] <- as.numeric(as.character(tpms_inverted3_f2[i,3]))
}
#Transform the numeric columns
TPM1_inverted_df <- as.data.frame(TPM1_inverted)
TPM1_inverted_df$length <- as.numeric(as.character(TPM1_inverted_df$length))
TPM1_inverted_df$effective_length <- as.numeric(as.character(TPM1_inverted_df$effective_length))
TPM1_inverted_df$expected_count <- as.numeric(as.character(TPM1_inverted_df$expected_count))
TPM1_inverted_df$TPM <- as.numeric(as.character(TPM1_inverted_df$TPM))
TPM1_inverted_df$FPKM <- as.numeric(as.character(TPM1_inverted_df$FPKM))
TPM1_inverted_df$IsoPct <- as.numeric(as.character(TPM1_inverted_df$IsoPct))

TPM2_inverted_df <- as.data.frame(TPM2_inverted)
TPM2_inverted_df$length <- as.numeric(as.character(TPM2_inverted_df$length))
TPM2_inverted_df$effective_length <- as.numeric(as.character(TPM2_inverted_df$effective_length))
TPM2_inverted_df$expected_count <- as.numeric(as.character(TPM2_inverted_df$expected_count))
TPM2_inverted_df$TPM <- as.numeric(as.character(TPM2_inverted_df$TPM))
TPM2_inverted_df$FPKM <- as.numeric(as.character(TPM2_inverted_df$FPKM))
TPM2_inverted_df$IsoPct <- as.numeric(as.character(TPM2_inverted_df$IsoPct))

TPM3_inverted_df <- as.data.frame(TPM3_inverted)
TPM3_inverted_df$length <- as.numeric(as.character(TPM3_inverted_df$length))
TPM3_inverted_df$effective_length <- as.numeric(as.character(TPM3_inverted_df$effective_length))
TPM3_inverted_df$expected_count <- as.numeric(as.character(TPM3_inverted_df$expected_count))
TPM3_inverted_df$TPM <- as.numeric(as.character(TPM3_inverted_df$TPM))
TPM3_inverted_df$FPKM <- as.numeric(as.character(TPM3_inverted_df$FPKM))
TPM3_inverted_df$IsoPct <- as.numeric(as.character(TPM3_inverted_df$IsoPct))

write.table(TPM1_inverted_df,file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513329_switched_SE.tpm",sep="\t",quote=FALSE,row.names = FALSE)
write.table(TPM2_inverted_df,file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513330_switched_SE.tpm",sep="\t",quote=FALSE,row.names = FALSE)
write.table(TPM3_inverted_df,file="~/SUPPA_supplementary_data/supplementary_data/simulated_reads/SRR1513331_switched_SE.tpm",sep="\t",quote=FALSE,row.names = FALSE)
