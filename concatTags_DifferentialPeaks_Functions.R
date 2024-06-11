## ---------------------------
##
## Copyright (c) Rhodes Ford, 2024
## Author: Dr. B. Rhodes Ford
## Email: rhodes.ford@ucsf.edu
## Date Created: May 2020
## Date Edited: 2024-06-11
## 
## 
## Script name: Concatenate Tag Counts and Differential Peak Analysis 
## (including Differential Accessible Regions)
##
## 
##
## ---------------------------
## Notes
## 
## Purpose of script: To take the tag count values from bedtools coverage of 
## all identified peaks (IDR or otherwise), concatenate tag count values for 
## each sample together as columns in a data.frame. Perform differential 
## analysis with DESeq2, and perform normalization using tag counts per million.
##
## Other functions include: getting heatmaps for the changes in gene expression 
## for differential peaks ()
## 
## README detailing any pre-R processing steps for ATAC, CUT&RUN, and CUT&Tag 
## data can be found within parent directory.
##
## Input: Directory with tag count values in bed format (chr, start, end, tag 
##  counts columns) for each sample
##        info file for differential analysis: format 
##
## Outputs: Seurat objects for each sample included in 10x dataset
##   
##    
##
## ---------------------------
packages <- c(
  "DESeq2", "GenomicRanges", "rtracklayer", "ChIPpeakAnno", "ChIPseeker", 
  "biomaRt", "org.Mm.eg.db", "TxDb.Mmusculus.UCSC.mm10.knownGene", "pheatmap", 
  "tidyverse", "dplyr", "factoextra"
)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    if(!require(pkg, character.only = TRUE)) {
      if(pkg %in% c("DESeq2", "GenomicRanges", "rtracklayer", "ChIPpeakAnno", "ChIPseeker", "biomaRt", "org.Mm.eg.db", "TxDb.Mmusculus.UCSC.mm10.knownGene")) {
        BiocManager::install(pkg)
      }
      if(!require(pkg, character.only = TRUE)) {
        stop(paste("Package", pkg, "could not be installed"))
      }
    }
  }
}

library(DESeq2)
library(GenomicRanges)
library(rtracklayer)
library(ChIPpeakAnno)
library(ChIPseeker)
library(biomaRt)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(factoextra)

##################### FUNCTIONS #####################

PerMillion <- function(x){
  x/(sum(x)/1000000)
}

GetBedFromPeakID <- function(peak.id){
  chr <- strsplit(peak.id, ':')[[1]][1]
  start.end <- strsplit(peak.id, ':')[[1]][2]
  print(start.end)
  start <- strsplit(start.end,  '-')[[1]][1]
  end <- strsplit(start.end,  '-')[[1]][2]
  df.return <- data.frame(chr,start,end)
  return(df.return)
}

# functions to concatenate individual tagCount files into data frame for a 
# provided mark of interest (tagCounts files for multiple marks can be contained 
# in one directory and will be filtered for provided mark)
MergeTags <- function(tag.dir, info.df, mark=NULL){ # Where directory
                                                      # infoFile is the information file for DESeq2
                                                      # and mark is the provided histone modification/data modality included in the file name (default is NULL)
  # Check if output directory exists and make directory if not
  print('Merging Tag Count Files')
  starting.dir <- getwd()
  print('Provided Info File has the following samples included:')
  print(info.df)
  
  # Set directory to the tag counts directory and get list of file names
  setwd(tag.dir)
  file.names <- list.files('./',pattern = '.bed')
  file.names <- file.names[grepl(file.names,pattern = 'allTag')==FALSE]
  if(is.null(mark) == TRUE){
    print("No mark provided. All files will be considered for analysis.")
    file.names.filtered <- file.names
  }else{
    file.names.filtered <- file.names[grepl(file.names,pattern = mark)==TRUE]
  }
  
  # Get sample names and order embedded in info file for each sample
  info.df$order <- 1:nrow(info.df)
  
  file.names.filtered <- as.data.frame(file.names.filtered)
  file.names.filtered$order <- NA
  
  # Re-order tag count file list ot match order in info file
  print("Re-ordering tag count file list to match info file.")
  for(i in 1:nrow(info.df)){
    for(j in 1:nrow(file.names.filtered)){
      if(grepl(file.names.filtered[j,1], pattern=info.df[i, 1])){
        file.names.filtered[j, 'order'] <- info.df[i, 3]
      }
    }
  }
  file.names.filtered <- file.names.filtered[order(file.names.filtered$order),]
  file.names.sample.names <- cbind(file.names.filtered[1], info.df[1])
  
  # read in all files as data frames in a list
  file.list <- list()
  
  for(i in 1:nrow(file.names.sample.names)){
    temp.df <- read.table(file.names.sample.names[i,1], 
                            header = F, 
                            sep = "\t",
                            col.names = c('chr', 'start', 'end', file.names.sample.names[i,2]))
    file.list[[file.names.sample.names[i,2]]] <- temp.df
  }
  
  # merge all the data frames by 'chr', 'start', 'end' columns
  tag.counts.merged <- file.list[[1]]
  for(i in 2:length(file.list)){
    tag.counts.merged <- merge(tag.counts.merged,file.list[[i]],by = c('chr','start','end'))
  }
  
  tag.counts.merged$PeakID <- paste(tag.counts.merged$chr,':',tag.counts.merged$start,'-',tag.counts.merged$end,sep='')
  tag.counts.merged <- tag.counts.merged[c(1:3,ncol(tag.counts.merged),4:(ncol(tag.counts.merged)-1))]
  
  # Extract the last character and add '/' if directory doesn't end with '/'
  if(substr(output.dir, nchar(output.dir), nchar(output.dir)) != '/'){
    output.dir <- paste(output.dir, '/', sep = '')
  }
  
  # determine whether or not mark information should be included in file name
  if(is.null(mark) == TRUE){
    out.file.path <- './allTagCounts.bed'
  }else{
    out.file.path <- paste('./', mark, '_allTagCounts.bed', sep = '')
  }
  
  # write bed file to designated path
  write.table(tag.counts.merged, 
              out.file.path, 
              sep='\t', 
              quote=FALSE, 
              col.names = TRUE, 
              row.names = FALSE)
  cat('The merged tag count file can be found at: ', out.file.path, '\n')
  # return starting directory
  setwd(starting.dir)
  
  # Assign "mark" attribute to the returned data frame
  attr(tag.counts.merged,'mark') <- mark
  return(tag.counts.merged)
}

# function to get differential peaks from file path to merged tag counts file 
# and info file
GetDifPeaks <- function(merged.tags.df, info.df, color.list = c('indianred3','lightgoldenrod2','purple3','steelblue3', 'lightpink1','lightgreen', 'grey85', 'black')){
  # Make differential peaks directory in the working directory 
  if(!dir.exists('./differentialPeaks/')){
    dir.create('./differentialPeaks/', recursive = TRUE)
    cat("Directory created:", "./differentialPeaks/", "\n")
  }else{
    cat("Directory already exists:", "./differentialPeaks/", "\n")
  }
  
  starting.dir <- getwd()
  setwd('./differentialPeaks/')
  
  # Get mark information from MergeTags
  mark <- attr(merged.tags.df, 'mark')
  
  number.of.groups <- length(unique(info.df$Group))
  info.df$Group <- as.factor(info.df$Group)

  # Example check for invalid ranges in merged.tags.df
  invalid.ranges <- merged.tags.df[merged.tags.df$end < merged.tags.df$start - 1, ]
  if (nrow(invalid.ranges) > 0) {
    stop("There are invalid ranges in merged.tags.df where end is less than start minus one.")
    return(invalid.ranges)
  }
  
  row.names(merged.tags.df) <- merged.tags.df$PeakID
  dds <- DESeqDataSetFromMatrix(countData = merged.tags.df[,5:ncol(merged.tags.df)], colData = info.df, design = ~Group)
  FeatureData <- row.names(merged.tags.df)
  dds <- estimateSizeFactors(dds)
  dds_counts <- counts(dds)
  genes <- dds@rowRanges@elementMetadata$gene_ids
  run <- DESeq(dds, fitType='local', useT = TRUE)
  
  lev <- levels(dds@colData$Group)
  
  
  all.peaks <- data.frame()
  
  if(!dir.exists('./PairwiseComparisons')){
    dir.create(file.path('./PairwiseComparisons'), showWarnings = FALSE)
    cat('PairwiseComparisons directory created: ', file.path('./PairwiseComparisons'))
  }else{
    cat('PairwiseComparisons directory already exists: ', file.path('./PairwiseComparisons'))
  }
  
  for(x in 1:length(lev)){
    for(y in 1:length(lev)){
      if(x != y & y > x){
        name <- paste(lev[x], lev[y], sep = 'vs')
        result <- results(run, contrast = c('Group', lev[x], lev[y]), format='DataFrame', independentFiltering = TRUE)
        write.csv(result, paste("PairwiseComparisons/", name, "_unfiltered.csv", sep = ''), row.names = TRUE)
        result <- result[abs(result$log2FoldChange) > 2 & result$padj < 0.05,]
        write.csv(result, paste("PairwiseComparisons/", name, ".csv", sep=''), row.names = TRUE)
        # get bed file
        result.temp <- as.data.frame(row.names(result))
        colnames(result.temp)[1] <- 'PeakID'
        result.temp.split <- separate(result.temp, PeakID, into = c("chr", "start_end"), sep = ":", remove = FALSE)
        result.temp.final <- separate(result.temp.split, start_end, into = c("start", "end"), sep = "-")
        write.table(result.temp.final[,c('chr','start','end')], paste("PairwiseComparisons/", name, ".bed", sep=''), sep='\t', row.names = FALSE, col.names=FALSE, quote = FALSE)
        result.peaks <- readPeakFile(paste("PairwiseComparisons/", name, ".bed", sep=''))
        result.anno <- annotatePeak(result.peaks, annoDb='org.Mm.eg.db', TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene)
        result.anno <- as.data.frame(result.anno)
        result.anno.out <- cbind(result.anno[c('seqnames', 'start', 'end', 'annotation', 'SYMBOL')],result)
        print(head(result.anno.out)[1])
        write.csv(result.anno.out, paste("PairwiseComparisons/", name, "_annotated.csv", sep=''), row.names = FALSE)
        all.peaks <- rbind(all.peaks, as.data.frame(result))
      }
    }
  }
  
  all.peaks$PeakID <- rownames(all.peaks)
  all.peaks.tag.counts <- merged.tags.df[merged.tags.df$PeakID %in% all.peaks$PeakID,]
  print(head(all.peaks.tag.counts)[1])
  ## normalize per million
  all.peaks.tag.counts.permillion <- lapply(all.peaks.tag.counts[,5:ncol(all.peaks.tag.counts)], PerMillion)
  all.peaks.tag.counts.permillion <- as.data.frame(all.peaks.tag.counts.permillion)
  all.peaks.tag.counts.permillion <- cbind(all.peaks.tag.counts[,1:4], all.peaks.tag.counts.permillion)
  rownames(all.peaks.tag.counts.permillion) <- all.peaks.tag.counts.permillion$PeakID
  
  # Make Pearson correlation plot
  corr <- cor(all.peaks.tag.counts.permillion[,5:ncol(all.peaks.tag.counts.permillion)], method = 'pearson')
  
  # Set output file path by whether or not a marker is in the metadata of provided df
  if(is.null(mark)){
    pearson.output.path <- 'pearsonCorrelation.pdf'
    pca.output.path <- 'PCA_plot.pdf'
    pretty.pca.output.path <- 'prettyPCA_plot.pdf'
    permillion.output <- 'differentialPeaks_allTIL_perMillion'
    anno.output <- 'differentialPeaks_annotated.csv'
  }else{
    pearson.output.path <- paste('pearsonCorrelation_', mark, '.pdf', sep='')
    pca.output.path <- paste('PCA_plot_', mark, '.pdf', sep='')
    pretty.pca.output.path <- paste('prettyPCA_plot_', mark, '.pdf', sep='')
    permillion.output <- paste(mark, '_differentialPeaks_allTIL_perMillion', sep='')
    anno.output <- paste(mark, '_differentialPeaks_annotated.csv', sep='')
  }
  
  # Write out Pearson correlation map
  pdf(file = pearson.output.path, onefile = TRUE, paper = 'USr')
  print(pheatmap(corr, fontsize = 6))
  dev.off()
  
  # Make PCA plot
  pca <- all.peaks.tag.counts.permillion[,5:ncol(all.peaks.tag.counts.permillion)]
  log <- log(pca + 1)
  t.log <- t(log)
  ind <- rowVars((t.log)) < .Machine$double.eps
  t.log.v <- t.log[!ind,]
  pca <- prcomp(t.log.v)
  
  pdf(file = pca.output.path, onefile = TRUE, paper = 'USr')
  print(fviz_pca_ind(pca, repel = TRUE))
  dev.off()
  
  pdf(file = pretty.pca.output.path)
  color.list.used <- color.list[1:number.of.groups]
  print(fviz_pca_ind(pca, geom.ind = 'point', pointsize = 7, habillage = info.df$Group, pointshape = 19, palette = color.list.used, mean.point = FALSE))
  dev.off()
  
  # Average the replicates
  index <- ncol(all.peaks.tag.counts.permillion) + 1
  for(i in 5:(index-2)){
    group <- as.character(info.df$Group)
    group <- group[!duplicated(group)]
    for(j in 1:length(group)){
      if(grepl(colnames(all.peaks.tag.counts.permillion)[i], pattern = group[j]) & grepl(colnames(all.peaks.tag.counts.permillion)[i+1], pattern = group[j])){
        all.peaks.tag.counts.permillion[index] <- (all.peaks.tag.counts.permillion[i] + all.peaks.tag.counts.permillion[i+1])/2
        colnames(all.peaks.tag.counts.permillion)[index] <- group[j]
        index <- index + 1
      }
    }
  }
  write.csv(all.peaks.tag.counts.permillion, paste(permillion.output, '.csv', sep=''), row.names = FALSE)
  write.table(all.peaks.tag.counts[,1:4], paste(permillion.output, '.bed', sep=''), sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  peaks <- readPeakFile(paste(permillion.output, '.bed', sep=''))
  all.peaks.anno <- annotatePeak(peaks, annoDb='org.Mm.eg.db', TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene)
  all.peaks.annotated <- as.data.frame(all.peaks.anno)
  colnames(all.peaks.annotated)[6] <- 'PeakID'
  
  # Get ChIPSeeker Outputs
  promoter <- getPromoters(TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, upstream=3000, downstream=3000)
  tag.matrix <- getTagMatrix(peaks, windows=promoter)
  
  pdf(file = paste('ChIPSeeker_plots_', mark, '.pdf', sep=''), onefile = TRUE, paper = 'USr')
  print(tagHeatmap(tag.matrix))
  print(plotAvgProf(tag.matrix, xlim = c(-3000, 3000), xlab = "Genomic Region (5'->3')", ylab = "Read Count Frequency"))
  print(plotAnnoPie(all.peaks.anno))
  dev.off()
  
  write.csv(all.peaks.annotated, anno.output, row.names = FALSE)
  returned.value <- all.peaks.tag.counts.permillion[,(ncol(all.peaks.tag.counts.permillion)-length(group)+1):ncol(all.peaks.tag.counts.permillion)]
  attr(returned.value, 'mark') <- mark
  setwd(starting.dir)
  return(returned.value)
}

getHeatmap <- function(perMillion,n){ # where perMillion is a data frame with the averaged perMillion values
                                      #n is the number of clusters
                                      #anno is the location of the annotated peak file read out in the getDifPeaks function
  if(n>7){
    stop('Number of clusters cannot be greater than 7. Please try running manually')
  }
  startingDir <- getwd()
  mark <- attr(perMillion,'mark')
  print(mark)
  setwd(paste('./',mark,'/differentialPeaks',sep=''))
  allPeaks_annotated <- read.csv(paste(mark,'_differentialPeaks_annotated.csv',sep=''))
  p<- pheatmap(as.matrix(perMillion), # matrix location of average values
               cluster_rows = T, # cluster peaks
               cluster_cols = F, # don't cluster samples
               main = paste(mark,'Differential Peaks \nTag Counts per Million ',sep=' '), # title of geneExp
               cutree_rows = n, # number of clusters to separate geneExp into
               scale='row') # scale each row to the row mean
  
  cluster <- cutree(p$tree_row, k=n)
  clusterNums <- 1:n
  clusterNames <- as.character(clusterNums)
  cluster <- as.data.frame(cluster)
  cluster$PeakID <- rownames(cluster)
  annotation <- as.factor(cluster$cluster)
  annotation <- as.data.frame(annotation)
  colnames(annotation)<-'cluster'
  rownames(annotation)<-rownames(cluster)
  cluster_col <- c('purple3','palegreen3','tomato','khaki1','dodgerblue','chocolate1','gray0')
  cluster_col <- cluster_col[1:n]
  names(cluster_col)<- clusterNames
  anno_colors <- list(cluster=cluster_col)
  
  p <- pheatmap(as.matrix(perMillion), 
                show_rownames=F, # don't show rownames
                show_colnames = T, # show samples names
                cluster_rows = T, # cluster peaks
                cluster_cols = F, # don't cluster samples
                annotation_row = annotation, # annotate rows with clusters or anything else you specifiy
                annotation_colors = anno_colors, # colors that are going with each cluster
                border_color = F, # no border for the boxes
                main = paste(mark,'Differential Peaks \nTag Counts per Million ',sep=' '), # title of geneExp
                treeheight_row=3, # height of dendogram
                width=6, # width of each box in geneExp
                cutree_rows = n, # number of clusters to separate geneExp into
                scale='row', # scale each row to the row mean
                filename = paste(mark,'_differentialPeaks_Heatmap_anno.pdf',sep='')
                )
  
  perMillion$PeakID <- rownames(perMillion)
  reordered <- rownames(perMillion[p$tree_row[["order"]],])
  reordered <- as.data.frame(reordered)
  colnames(reordered)<- c('PeakID')
  cluster <- cutree(p$tree_row, k=n)
  cluster <- as.data.frame(cluster)
  cluster$PeakID <- rownames(cluster)
  reordered_cut <- merge(reordered, cluster, by = 'PeakID', sort = F)
  annot_cut_reordered <- merge(reordered_cut, perMillion, by = "PeakID", sort = F)
  rownames(annot_cut_reordered)  <- annot_cut_reordered$PeakID
  differentialPeaks_annotated_clustered <- merge(allPeaks_annotated,annot_cut_reordered,by='PeakID',sort=F)
  write.table(differentialPeaks_annotated_clustered, paste(mark,'_annotatedPeaks_inHeatmapOrder.txt',sep=''), quote = F, sep = '\t',col.names = T, row.names = F)
  # make bed files for Homer motif analysis
  for(i in clusterNums){
    anno_sorted <- differentialPeaks_annotated_clustered[differentialPeaks_annotated_clustered$cluster == i,]
    bed <- anno_sorted[,c('seqnames','start','end','PeakID')]
    bed <- cbind(bed,'.')
    bed <-cbind(bed,0)
    write.table(anno_sorted, paste(mark,'_cluster',i,'all_fromHeatmap.txt',sep=''), quote = F, sep = '\t',col.names = T, row.names = F)
    write.table(bed, paste(mark,'_cluster',i,'_forHomer.bed',sep=''), quote = F, sep = '\t',col.names = T, row.names = F)
  }
  setwd(startingDir)
  attr(differentialPeaks_annotated_clustered,'mark') <- mark
  return(differentialPeaks_annotated_clustered)
}


getGeneExp <- function(heatmapOutput,tpm,dge,difPeakClusterNames,clusterNum){ # where heatmapOutput is the output of the getHeatmap function
                                              #tpm is the file path to the TPM csv file
                                              #dge is thefile path to the dge txt file
                                              #difPeakClusterNames is the names of clusters in order eg c(cluster1name, cluster2name, cluster3name...)
                                              #clusterNum is the number of clusters for output gene expression heatmaps
  returnedValues <- list()
  mark <- attr(heatmapOutput,'mark')
  # read in tpm and dge
  tpm <- read.csv(tpm)
  dge <- read.table(dge, sep ='\t',col.names='GeneID')  # Load in differential peaks with tags per million
  if(is.null(mark)){
    stop('Sorry the heatmapOutput input has no mark information')
  }
  
  startingDir <- getwd()
  setwd(paste('./',mark,'/differentialPeaks',sep = ''))
  if(length(difPeakClusterNames) !=2){
    stop('Sorry this script can only create gene expression heatmaps using two differential peak clusters right now')
  }
  
  for(n in 1:length(difPeakClusterNames)){
    print(paste("Getting gene expression information for the ",difPeakClusterNames[n],' cluster',sep=''))
    up_bed <- heatmapOutput[heatmapOutput$cluster == as.character(n),c('seqnames','start','end','PeakID')]
    up_bed$Blank <- '.'
    up_bed$Strand <- 0
    write.table(up_bed, paste('./',mark,'_differentialPeaks_',difPeakClusterNames[n],'EnrichedPeaks.bed',sep=''),sep='\t',col.names=F,row.names=F,quote=F)
    
    up_annotated <- heatmapOutput[heatmapOutput$cluster == as.character(n),'SYMBOL']
    up_annotated <- up_annotated[!duplicated(up_annotated)==TRUE]
    up_tpms <- tpm[tpm$GeneID %in% up_annotated,]
    rownames(up_tpms)<-up_tpms$GeneID
    up_tpms_var <- up_tpms[apply(up_tpms[,2:ncol(up_tpms)], MARGIN = 1, FUN= var) !=0,2:ncol(up_tpms)]
    
    annotation <- data.frame('DGE' = 1:nrow(up_tpms_var))
    rownames(annotation) <- rownames(up_tpms_var)
    annotation$DGE <- as.factor(ifelse(rownames(annotation) %in% dge$GeneID,'yes','no'))
    
    if(nrow(up_tpms_var) > 100){
      fontsize = 2
    }else{
      fontsize = 5
    }
    
    p_dp <- pheatmap(as.matrix(up_tpms_var), show_rownames=T,show_colnames = T, 
                     cluster_rows = T, cluster_cols = F, border_color = F, fontsize_row=fontsize, 
                     annotation_row = annotation, cutree_rows = clusterNum,
                     main = paste('Genes Associated with \n',difPeakClusterNames[n],'-Specific ',mark,sep=''), treeheight_row=4, width=4, scale = 'row',
                     filename = paste('./',difPeakClusterNames[n],'Specific',mark,'_associatedGenes.pdf',sep='')
                     )
    up_tpms_var$GeneID <- rownames(up_tpms_var)
    reordered <- rownames(up_tpms_var[p_dp$tree_row[["order"]],])
    reordered <- as.data.frame(reordered)
    colnames(reordered)<- c('GeneID')
    cluster <- cutree(p_dp$tree_row, k=clusterNum)
    cluster <- as.data.frame(cluster)
    cluster$GeneID <- rownames(cluster)
    reordered_cut <- merge(reordered, cluster, by = 'GeneID', sort = F)
    annot_cut_reordered <- merge(reordered_cut, up_tpms_var, by = "GeneID", sort = F)
    rownames(annot_cut_reordered)  <- annot_cut_reordered$GeneID
    write.table(annot_cut_reordered, paste('./',difPeakClusterNames[n],'Specific',mark,'_associatedGenes_inHeatmapOrder_test.txt',sep=''), quote = F, sep = '\t',col.names = T, row.names = F)
    returnedValues[[as.character(difPeakClusterNames[n])]] <- annot_cut_reordered
  }
  setwd(startingDir)
  attr(returnedValues,'mark') <- mark
  return(returnedValues)
} 


getGeneExpBed <- function(geneExp,expInfo,heatmapOutput){
  startingDir <- getwd()
  expInfo <- read.csv(expInfo,row.names=1)
  counter <- numeric()
  for(a in 1:ncol(expInfo)){
    if(any(is.na(expInfo[,a]))){
      print(paste(colnames(expInfo)[a],' was not processed because it contained NA values. If you want ',colnames(expInfo)[a],' to be considered, add cluster numbers to expInfo file',sep=''))
      expInfo_processed <- subset(expInfo, select = -c(a))
    }else{
      expInfo_processed <- expInfo
      counter <- c(counter,a)
    }
  }
  colnames(expInfo_processed)[1:ncol(expInfo_processed)] <- colnames(expInfo)[counter]
  mark <- attr(geneExp,'mark')
  setwd(paste('./',mark,'/differentialPeaks/',sep = ''))
  peakDataColNames <- as.character(info[!duplicated(info$Group),'Group'])
  if(is.null(mark)){
    stop('Sorry the geneExp input has no mark information')
  }
  
  for(i in 1:length(geneExp)){ # Loops for each of the data frames from the differential peaks clusters
    if(length(geneExp)==ncol(expInfo_processed)){
      clusters <- as.data.frame(expInfo_processed[,(attr(geneExp[i],'name'))])
      row.names(clusters) <- row.names(expInfo_processed)
      print(paste('Making files for ',mark,' ',as.character(attr(geneExp[i],'name')),sep=''))
      geneExp_working <- as.data.frame(geneExp[i])
      colnames(geneExp_working) <- sub(colnames(geneExp_working),pattern = paste(attr(geneExp[i],'name'),'.',sep=''),replacement = '')
      for(n in 1:nrow(clusters)){
        print(paste('Making bed files for ',mark,' ',as.character(attr(geneExp[i],'name')),' ',row.names(clusters)[n],sep=''))
        geneExp_subset <- geneExp_working[geneExp_working$cluster == as.character(clusters[[n,1]]),]
        geneExp_total <- heatmapOutput[heatmapOutput$cluster == as.character(i),c('seqnames','start','end','PeakID','SYMBOL','cluster',peakDataColNames)]
        geneExp_subset_annotated <- geneExp_total[geneExp_total$SYMBOL %in% geneExp_subset$GeneID,]
        View(geneExp_subset_annotated)
        #write.table(geneExp_subset_annotated, paste('./',mark,'_',attr(geneExp[i],'name'),'SpecificPeaks_',row.names(clusters)[n],'GeneExpression_annotated_test.txt',sep=''), 
        #            sep='\t',col.names=T, row.names=F,quote=F)
        # make files that can be used in homer
        geneExp_subset_annotated$Blank <- '.'
        geneExp_subset_annotated$Strand <- 0
        print(mark)
        #write.table(geneExp_subset_annotated[,c('seqnames','start','end','PeakID','Blank','Strand')], paste('./',mark,'_',attr(geneExp[i],'name'),'SpecificPeaks_',row.names(clusters)[n],'GeneExpression.bed',sep=''), 
        #            sep='\t',col.names=F, row.names=F,quote=F)
      }
    }else if(grepl(colnames(expInfo_processed),pattern = attr(geneExp[i],'name'))){
      clusters <- as.data.frame(expInfo_processed[,attr(geneExp[i],'name')])
      row.names(clusters) <- row.names(expInfo_processed)
      print(paste('Making files for ',mark,' ',as.character(attr(geneExp[i],'name')),sep=''))
      geneExp_working <- as.data.frame(geneExp[i])
      colnames(geneExp_working) <- sub(colnames(geneExp_working),pattern = paste(attr(geneExp[i],'name'),'.',sep=''),replacement = '')
      for(n in 1:nrow(clusters)){
        print(paste('Making bed files for ',mark,' ',as.character(attr(geneExp[i],'name')),' ',row.names(clusters)[n],sep=''))
        geneExp_subset <- geneExp_working[geneExp_working$cluster == as.character(clusters[[n,1]]),]
        geneExp_total <- heatmapOutput[heatmapOutput$cluster == as.character(i),c('seqnames','start','end','PeakID','SYMBOL','cluster',peakDataColNames)]
        geneExp_subset_annotated <- geneExp_total[geneExp_total$SYMBOL %in% geneExp_subset$GeneID,]
        #write.table(geneExp_subset_annotated, paste('./',attr(geneExp[i],'name'),'SpecificPeaks_',row.names(clusters)[n],'GeneExpression_annotated.txt',sep=''), 
        #            sep='\t',col.names=T, row.names=F,quote=F)
        # make files that can be used in homer
        geneExp_subset_annotated$Blank <- '.'
        geneExp_subset_annotated$Strand <- 0
        #write.table(geneExp_subset_annotated[,c('seqnames','start','end','PeakID','Blank','Strand')], paste('./',attr(geneExp[i],'name'),'SpecificPeaks_',row.names(clusters)[n],'GeneExpression.bed',sep=''), 
        #            sep='\t',col.names=F, row.names=F,quote=F)
      }
    }
  }
  setwd(startingDir)
}

  
# make differential gene list
dgeGo <-function(directory){
  origDir <- getwd()
  setwd(directory)
  files <- list.files('./',pattern = '.csv')
  y <- lapply(files, # y is a list of tag count data frames
              FUN = function(files) {
                read.csv(files,col.names = c('GeneID','baseMean','log2FC','lfcSE','stat','pvalue','padj'))
              })
  data <- as.data.frame(y[1])
  for(i in 2:length(y)){
    data <- rbind(data,y[[i]])
  }
  data <- data[duplicated(data$GeneID),]
  write.table(data$GeneID,'../../CSH_dge_list.txt',sep='\t',row.names=F,col.names = T)
  setwd(origDir)
  return(data)
}

