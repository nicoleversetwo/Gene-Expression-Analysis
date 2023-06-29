# -------------   Main objective -------------------
# create a heatmap of gene expression
# adding mutation, copy number variation information
#---------------------------------------------------
setwd("~/Cedars")
library(readxl)
library(TCGAbiolinks)
library(RTCGAToolbox)
library(rowr)
library(ComplexHeatmap)
library(SummarizedExperiment)
library(circlize)
library(stringr)
library(gplots)
#---------------------------------------------------------
# STEP 1: GET CNV results from GDAC firehose
#---------------------------------------------------------
gistic = getFirehoseData(dataset="LUSC", gistic2_Date= getFirehoseAnalyzeDates(last=1))
gistic.Thresholed <- gistic@GISTIC@ThresholedByGene
gistic.all <- gistic@GISTIC@AllByGene
cnv.nfe2l2 <- t(gistic.Thresholed[gistic.Thresholed$Gene.Symbol == "NFE2L2",-c(1:3)])
cnv.keap1  <- t(gistic.Thresholed[gistic.Thresholed$Gene.Symbol == "KEAP1",-c(1:3)])
cnv.cul3  <- t(gistic.Thresholed[gistic.Thresholed$Gene.Symbol == "CUL3",-c(1:3)])
cnv.annotation <- cbind(cnv.nfe2l2,cnv.keap1,cnv.cul3)
colnames(cnv.annotation) <- c("CNV NFE2L2","CNV KEAP1", "CNV CUL3")
rownames(cnv.annotation) <- substr(gsub("\\.","-",rownames(cnv.annotation)),1,15)
#---------------------------------------------------------
# STEP 2: GET mutations from GDC MAF
#---------------------------------------------------------
# Get mutation annotation file
maf.lusc <- GDCquery_Maf("LUSC")
# Lets consider any type of mutation as the same for the moment
mut.nf2l2 <- data.frame(patient = substr(unique(maf.lusc[maf.lusc$Hugo_Symbol == "NFE2L2",]$Tumor_Sample_Barcode),1,15),
                        mut.nf2l2 = TRUE, row.names = 1)
mut.keap1 <- data.frame(patient = substr(unique(maf.lusc[maf.lusc$Hugo_Symbol == "KEAP1",]$Tumor_Sample_Barcode),1,15),
                        mut.keap1 = TRUE, row.names = 1)
mut.annotation <- merge(mut.nf2l2, mut.keap1, by = 0, all = TRUE)
rownames(mut.annotation) <- mut.annotation$Row.names; mut.annotation$Row.names <- NULL
mut.cul3 <- data.frame(patient = substr(unique(maf.lusc[maf.lusc$Hugo_Symbol == "CUL3",]$Tumor_Sample_Barcode),1,15),
                       mut.cul3 = TRUE, row.names = 1)
mut.annotation <- merge(mut.annotation,mut.cul3, by = 0, all = TRUE)
rownames(mut.annotation) <- mut.annotation$Row.names; mut.annotation$Row.names <- NULL

#---------------------------------------------------------
# STEP 3: Gene expression data
#      3.1 GET gene expression data aligned to hg19 from GDC legacy archive
#      3.2 Do the same normalization made before
#---------------------------------------------------------
query <- GDCquery(project = "TCGA-LUSC",
                  data.category = "Gene expression",
                  legacy = TRUE ,
                  experimental.strategy = "RNA-Seq",
                  data.type = "Gene expression quantification",
                  file.type = "normalized_results",
                  platform = "Illumina HiSeq",
                  sample.type = "Primary solid Tumor")
GDCdownload(query)
exp.file <- "LUSC_gene_expression_hg19_pst.rna"
if(!file.exists(exp.file)){
  lusc.exp <- GDCprepare(query, save = TRUE, save.filename = exp.file)
} else {
  lusc.exp <- get(load(exp.file))
  lusc.exp.tumor <- lusc.exp[,lusc.exp$shortLetterCode == "TP"]
}
# Select same genes for motif NFE2
file <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-015-0668-3/MediaObjects/13059_2015_668_MOESM7_ESM.xlsx"
if(!file.exists(basename(file))) downloader::download(file,basename(file))
genes <- read_excel(basename(file),sheet = 2)
genes.motif.nfe <- genes[genes$motif == "NFE2" & genes$CT == "LUSC",]$Gene
genes.motif.nfe <- genes.motif.nfe[genes.motif.nfe %in% values(lusc.exp)$gene_id]
lusc.exp.tumor <- lusc.exp.tumor[genes.motif.nfe,]


# There is one replicated samples: TCGA-21-1076-01A-01R-0692-07 and  TCGA-21-1076-01A-02R-0692-07
# We will remove the second one
lusc.exp.tumor <- lusc.exp.tumor[,!grepl("TCGA-21-1076-01A-02R-0692-07", colnames(lusc.exp.tumor))]

# Normalizing the gene expression data as before
NormalizeMedian <- function (x, col=FALSE, row=FALSE, na.rm=FALSE){
  if(col){
    Median <- apply(x,2,median,na.rm=na.rm)
    x <- t((t(x)-Median))
  }
  if(row){
    Median <- apply(x,1,median,na.rm=na.rm)
    x <- x-Median
  }
  return(x)
}


mat <- assay(lusc.exp.tumor)
norm <- NormalizeMedian(log2(mat+1), row = T)
norm[norm > 3] <- 3
norm[norm < -3] <- -3

#---------------------------------------------------------
# STEP 4: Create heatmap
#      4.1: Prepare annotation file w/ the following information
#           - mutation for KEAP1 and NFE2L2 genes
#           - CNV for KEAP1 and NFE2L2
#           - Was in old pancan dataset?
#      4.2: order inner clusters
#      4.3: plot
#---------------------------------------------------------

annotation <- colData(lusc.exp)[,c("sample","definition","gender","race","tumor_grade")]
rownames(annotation) <- substr(annotation$sample,1,15)
annotation <- merge(annotation, mut.annotation, by = 0 , sort = FALSE, all=TRUE)
rownames(annotation) <- annotation$Row.names; annotation$Row.names <- NULL
annotation <- merge(annotation, cnv.annotation, by = 0 , sort = FALSE,all=TRUE)
rownames(annotation) <- annotation$Row.names; annotation$Row.names <- NULL
annotation <- annotation[na.omit(match(substr(colnames(norm),1,16),annotation$sample)),]

all(substr(colnames(norm),1,15) == rownames(annotation)) # Should be TRUE
file <- "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-015-0668-3/MediaObjects/13059_2015_668_MOESM1_ESM.xlsx"
if(!file.exists(basename(file))) downloader::download(file, basename(file))
paper.samples <- read_excel(basename(file),sheet = 3) # 3 is RNA-seq samples ID
paper.samples <- paper.samples[paper.samples$TN.cat == "LUSC",]

annotation$oldSamples <- FALSE
annotation[annotation$sample %in% substr(paper.samples$TCGA.ID,1,16),]$oldSamples <- TRUE
# Do intern cluster to be easier to see patterns
# TODO: split the gorups in 4, annotation needs to be ordered.
# Sort 1 oldSamples, mutation NFE2L2
cluster.list <- colnames(annotation)[6:8]
for(i in cluster.list){
  print(i)
  samples <- which(annotation$oldSamples == T & annotation[,i] == TRUE)
  aux1 <- norm[,samples]
  dist1 <-  t(aux1) %>% dist %>% hclust(method = "ward.D2")
  annot1 <- annotation[samples,]
  
  
  samples <- which(annotation$oldSamples == T & is.na(annotation[,i]))
  aux2 <- norm[,samples]
  dist2 <-  t(aux2) %>% dist %>% hclust(method = "ward.D2")
  annot2 <- annotation[samples,]
  
  samples <- which(annotation$oldSamples == F & annotation[,i] == TRUE)
  aux3 <- norm[,samples]
  dist3 <-  t(aux3) %>% dist %>% hclust(method = "ward.D2")
  annot3 <- annotation[samples,]
  
  
  samples <- which(annotation$oldSamples == F & is.na(annotation[,i]))
  aux4 <- norm[,samples]
  dist4 <-  t(aux4) %>% dist %>% hclust(method = "ward.D2")
  annot4 <- annotation[samples,]
  
  # Put the samples in expression matrix in the right order
  norm.order <- cbind(aux1[,dist1$order],
                      aux2[,dist2$order],
                      aux3[,dist3$order],
                      aux4[,dist4$order])
  # Put the annotation in the right order
  annot.order <- rbind(annot1[dist1$order,],
                       annot2[dist2$order,],
                       annot3[dist3$order,],
                       annot4[dist4$order,])
  
  annot.order <- annot.order[,c(12:6,2)]
  
  ha = HeatmapAnnotation(annot.order,
                         col = list("CNV KEAP1" = c("-2" = "red", "-1" = "orange","0"= "grey","1"="green", "2"="black"),
                                    "CNV CUL3" = c("-2" = "red", "-1" = "orange","0"= "grey","1"="green", "2"="black"),
                                    "CNV NFE2L2" = c("-2" = "red", "-1" = "orange","0"= "grey","1"="green", "2"="black")
                         ))
  
  heatmap <- Heatmap(norm.order, col = greenred(255),
                     show_column_names = F, top_annotation = ha,show_row_names = F, cluster_columns = F, cluster_rows = T,
                     heatmap_legend_param = list(color_bar = "continuous"),row_names_gp =  gpar(fontsize = 6))
  
  pdf(paste0("heatmap_order_by_",i,".pdf"), width = 10, height = 9)
  draw(heatmap)
  dev.off()
}

cluster.list <- colnames(annotation)[9:11]
for(i in cluster.list){
  print(i)
  
  samples <- which(annotation$oldSamples == T & annotation[,i] %in% c("-1","-2"))
  aux1 <- norm[,samples]
  dist1 <-  t(aux1) %>% dist %>% hclust(method = "ward.D2")
  annot1 <- annotation[samples,]
  
  
  samples <- which(annotation$oldSamples == T & annotation[,i]  %in% c("1","2"))
  aux2 <- norm[,samples]
  dist2 <-  t(aux2) %>% dist %>% hclust(method = "ward.D2")
  annot2 <- annotation[samples,]
  
  
  samples <- which(annotation$oldSamples == T & annotation[,i]  %in% c("0"))
  aux3 <- norm[,samples]
  dist3 <-  t(aux3) %>% dist %>% hclust(method = "ward.D2")
  annot3 <- annotation[samples,]
  
  samples <- which(annotation$oldSamples == F & annotation[,i]  %in% c("-1","-2"))
  aux4 <- norm[,samples]
  dist4 <-  t(aux4) %>% dist %>% hclust(method = "ward.D2")
  annot4 <- annotation[samples,]
  
  
  samples <- which(annotation$oldSamples == F & annotation[,i]  %in% c("1","2"))
  aux5 <- norm[,samples]
  dist5 <-  t(aux5) %>% dist %>% hclust(method = "ward.D2")
  annot5 <- annotation[samples,]
  
  samples <- which(annotation$oldSamples == F & annotation[,i]  %in% c("0"))
  aux6 <- norm[,samples]
  dist6 <-  t(aux6) %>% dist %>% hclust(method = "ward.D2")
  annot6 <- annotation[samples,]
  
  # Put the samples in expression matrix in the right order
  norm.order <- cbind(aux1[,dist1$order],
                      aux2[,dist2$order],
                      aux3[,dist3$order],
                      aux4[,dist4$order],
                      aux5[,dist5$order],
                      aux6[,dist6$order]
  )
  # Put the annotation in the right order
  annot.order <- rbind(annot1[dist1$order,],
                       annot2[dist2$order,],
                       annot3[dist3$order,],
                       annot4[dist4$order,],
                       annot5[dist5$order,],
                       annot6[dist6$order,])
  
  annot.order <- annot.order[,c(12:6,2)]
  
  ha = HeatmapAnnotation(annot.order,
                         col = list("CNV KEAP1" = c("-2" = "red", "-1" = "orange","0"= "grey","1"="green", "2"="black"),
                                    "CNV CUL3" = c("-2" = "red", "-1" = "orange","0"= "grey","1"="green", "2"="black"),
                                    "CNV NFE2L2" = c("-2" = "red", "-1" = "orange","0"= "grey","1"="green", "2"="black")
                         ))
  
  heatmap <- Heatmap(norm.order, col = greenred(255),
                     show_column_names = F, top_annotation = ha,show_row_names = F, cluster_columns = F, cluster_rows = T,
                     heatmap_legend_param = list(color_bar = "continuous"),row_names_gp =  gpar(fontsize = 6))
  
  pdf(paste0("heatmap_order_by_",i,".pdf"), width = 10, height = 9)
  draw(heatmap)
  dev.off()
}


samples <- which(annotation$oldSamples == T)
aux1 <- norm[,samples]
dist1 <-  t(aux1) %>% dist %>% hclust(method = "ward.D2")
annot1 <- annotation[samples,]


samples <- which(annotation$oldSamples == F)
aux2 <- norm[,samples]
dist2 <-  t(aux2) %>% dist %>% hclust(method = "ward.D2")
annot2 <- annotation[samples,]

# Put the samples in expression matrix in the right order
norm.order <- cbind(aux1[,dist1$order],
                    aux2[,dist2$order])
# Put the annotation in the right order
annot.order <- rbind(annot1[dist1$order,],
                     annot2[dist2$order,])
annot.order <- annot.order[,c(12:6,2)]

ha = HeatmapAnnotation(annot.order,
                       col = list("CNV KEAP1" = c("-2" = "red", "-1" = "orange","0"= "grey","1"="green", "2"="black"),
                                  "CNV CUL3" = c("-2" = "red", "-1" = "orange","0"= "grey","1"="green", "2"="black"),
                                  "CNV NFE2L2" = c("-2" = "red", "-1" = "orange","0"= "grey","1"="green", "2"="black")
                       ))

heatmap <- Heatmap(norm.order, col = greenred(255),
                   show_column_names = F, top_annotation = ha,show_row_names = F, cluster_columns = F, cluster_rows = T,
                   heatmap_legend_param = list(color_bar = "continuous"),row_names_gp =  gpar(fontsize = 6))

pdf(paste0("heatmap_order_oldvsnew.pdf"), width = 10, height = 9)
draw(heatmap)
dev.off()




save.image(file="LUSC_heatmap.Rda")


