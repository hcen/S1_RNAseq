#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
library("DESeq2")

# install.packages("tidyverse")
library("tidyverse")
library(readxl)
#install.packages("xlsx")
#library("xlsx")

#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
#BiocManager::install("clusterProfiler")
#BiocManager::install("ReactomePA")
library("clusterProfiler")
library("ReactomePA")

#install.packages("pheatmap")
library(pheatmap)
#install.packages("UpSetR")
library(UpSetR)

#BiocManager::install("biomaRt")
library("biomaRt")

#install.packages("devtools") # most recent version of complexheatmap
#library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(circlize) 
#install.packages("gridtext")
library(gridtext)
library(scales)

#install.packages("msigdbr")
library(msigdbr) # load MSigDB gene sets, v7.5.1 (released January 2022) 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) ##Set working directory to where this file is.
getwd()

h.cor=read_excel("output/correlated_genes_H1_new_all5.xlsx")[-1]%>%
  select(-c(entrez:description))
View(h.cor)
m.cor=read_excel("output/correlated_gene_heatmap.xlsx")[-1] %>%
  select(symbol:description)
View(m.cor)
mh.cor.consistent=m.cor%>%left_join(h.cor,by="symbol",suffix=c("_M","_H"))%>%
  mutate(both.correlated=ifelse(adj.p_H<0.05, "yes", "no"))
View(mh.cor.consistent)
write.csv(mh.cor.consistent,"output/MH_both_correlated.csv",row.names = F)

# prepare raw counts and meta data ===================

M1.G.68 = read.table("input/For Howard - All H1 data package/Diff68-69 H1 cells_raw counts/H1-M1G-S1D3-Diff68.counts.genes.coverage.txt") %>% 
  dplyr::rename(M1.G.68=V2)
M1.G.69 = read.table("input/For Howard - All H1 data package/Diff68-69 H1 cells_raw counts/H1-M1G-S1D3-Diff69.counts.genes.coverage.txt") %>% 
  dplyr::rename(M1.G.69=V2)
C3.G.68 = read.table("input/For Howard - All H1 data package/Diff68-69 H1 cells_raw counts/H1-C3G-S1D3-Diff68.counts.genes.coverage.txt") %>% 
  dplyr::rename(C3.G.68=V2)
C3.G.69 = read.table("input/For Howard - All H1 data package/Diff68-69 H1 cells_raw counts/H1-C3G-S1D3-Diff69.counts.genes.coverage.txt") %>% 
  dplyr::rename(C3.G.69=V2)
W50.G.68 = read.table("input/For Howard - All H1 data package/Diff68-69 H1 cells_raw counts/H1-W50G-S1D3-Diff68.counts.genes.coverage.txt") %>% 
  dplyr::rename(W50.G.68=V2)
W50.G.69 = read.table("input/For Howard - All H1 data package/Diff68-69 H1 cells_raw counts/H1-W50G-S1D3-Diff69.counts.genes.coverage.txt") %>% 
  dplyr::rename(W50.G.69=V2)
C3.A.74 = read.table("input/For Howard - All H1 data package/Diff74 H1 cells_raw counts/74-C3A.counts.genes.coverage.txt") %>% 
  dplyr::rename(C3.A.74=V2)
C3.G.74 = read.table("input/For Howard - All H1 data package/Diff74 H1 cells_raw counts/74-C3G.counts.genes.coverage.txt") %>% 
  dplyr::rename(C3.G.74=V2)
M1.A.74 = read.table("input/For Howard - All H1 data package/Diff74 H1 cells_raw counts/74-M1A.counts.genes.coverage.txt") %>% 
  dplyr::rename(M1.A.74=V2)
M1.G.74 = read.table("input/For Howard - All H1 data package/Diff74 H1 cells_raw counts/74-M1G.counts.genes.coverage.txt") %>% 
  dplyr::rename(M1.G.74=V2)
W50.A.74 = read.table("input/For Howard - All H1 data package/Diff74 H1 cells_raw counts/74-W50A.counts.genes.coverage.txt") %>% 
  dplyr::rename(W50.A.74=V2)
W50.G.74 = read.table("input/For Howard - All H1 data package/Diff74 H1 cells_raw counts/74-W50G.counts.genes.coverage.txt") %>% 
  dplyr::rename(W50.G.74=V2)
C3.A.77 = read.table("input/For Howard - All H1 data package/Diff77 H1 cells_raw counts/77-C3A.counts.genes.coverage.txt") %>% 
  dplyr::rename(C3.A.77=V2)
C3.G.77 = read.table("input/For Howard - All H1 data package/Diff77 H1 cells_raw counts/77-C3G.counts.genes.coverage.txt") %>% 
  dplyr::rename(C3.G.77=V2)
M1.A.77 = read.table("input/For Howard - All H1 data package/Diff77 H1 cells_raw counts/77-M1A.counts.genes.coverage.txt") %>% 
  dplyr::rename(M1.A.77=V2)
M1.G.77 = read.table("input/For Howard - All H1 data package/Diff77 H1 cells_raw counts/77-M1G.counts.genes.coverage.txt") %>% 
  dplyr::rename(M1.G.77=V2)
W50.A.77 = read.table("input/For Howard - All H1 data package/Diff77 H1 cells_raw counts/77-W50A.counts.genes.coverage.txt") %>% 
  dplyr::rename(W50.A.77=V2)
W50.G.77 = read.table("input/For Howard - All H1 data package/Diff77 H1 cells_raw counts/77-W50G.counts.genes.coverage.txt") %>% 
  dplyr::rename(W50.G.77=V2)
C3.A.80 = read.table("input/For Howard - All H1 data package/Diff80 H1 cells_raw counts/80-C3A.counts.genes.coverage.txt") %>% 
  dplyr::rename(C3.A.80=V2)
C3.G.80 = read.table("input/For Howard - All H1 data package/Diff80 H1 cells_raw counts/80-C3G.counts.genes.coverage.txt") %>% 
  dplyr::rename(C3.G.80=V2)
M1.A.80 = read.table("input/For Howard - All H1 data package/Diff80 H1 cells_raw counts/80-M1A.counts.genes.coverage.txt") %>% 
  dplyr::rename(M1.A.80=V2)
M1.G.80 = read.table("input/For Howard - All H1 data package/Diff80 H1 cells_raw counts/80-M1G.counts.genes.coverage.txt") %>% 
  dplyr::rename(M1.G.80=V2)
W50.A.80 = read.table("input/For Howard - All H1 data package/Diff80 H1 cells_raw counts/80-W50A.counts.genes.coverage.txt") %>% 
  dplyr::rename(W50.A.80=V2)
W50.G.80 = read.table("input/For Howard - All H1 data package/Diff80 H1 cells_raw counts/80-W50G.counts.genes.coverage.txt") %>% 
  dplyr::rename(W50.G.80=V2)



raw.counts <- W50.G.68 %>% full_join(W50.G.69,by="V1") %>%
  full_join(W50.G.74,by="V1") %>%full_join(W50.G.77,by="V1") %>%
  full_join(W50.G.80,by="V1") %>% full_join(M1.A.74,by="V1") %>% 
  full_join(M1.A.77,by="V1") %>% full_join(M1.A.80,by="V1") %>% 
  full_join(W50.A.74,by="V1") %>% full_join(W50.A.77,by="V1") %>% 
  full_join(W50.A.80,by="V1") %>% full_join(C3.A.74,by="V1") %>% 
  full_join(C3.A.77,by="V1") %>% full_join(C3.A.80,by="V1") %>%
  full_join(C3.G.68,by="V1") %>% full_join(C3.G.69,by="V1") %>% 
  full_join(C3.G.74,by="V1") %>% full_join(C3.G.77,by="V1") %>%
  full_join(C3.G.80,by="V1") %>% full_join(M1.G.68,by="V1") %>%
  full_join(M1.G.69,by="V1") %>% full_join(M1.G.74,by="V1") %>% 
  full_join(M1.G.77,by="V1") %>% full_join(M1.G.80,by="V1")

raw.counts <- column_to_rownames(raw.counts,"V1")
view(raw.counts)
dim(raw.counts) #26364
raw.counts <- raw.counts[rowSums(raw.counts)!=0,]
dim(raw.counts) #22220
write.table(raw.counts, sep="\t",file="input/rawcount_H1_new.txt", row.names=TRUE,col.names=NA,quote=FALSE)

#raw.counts = raw.counts[!duplicated(raw.counts$Gene),]
#dup=raw.counts %>% group_by(Gene) %>% dplyr::filter(n() > 1)
#view(dup)
#raw.counts = read.table(file="input/rawcount.txt", row.names=1) # load counts of MEL1 alone without H1
#raw.counts = read.table(file="input/rawcount_MEL1_H1.txt", row.names=1) # or load counts of MEL1 + H1
#raw.counts = read.table(file="input/rawcount_H1_new.txt", row.names=1) # load new H1 n=3 or 5

raw.counts <- round(raw.counts,0) #Some values are not integers. Rounded up.


#
sample <- colnames(raw.counts)
group <- c(rep("W50.G",5),rep("M1.A",3),rep("W50.A",3),rep("C3.A",3),rep("C3.G",5),rep("M1.G",5))

meta.data <- data.frame(sample, group)
row.names(meta.data) <- meta.data$sample
all(colnames(raw.counts)==rownames(meta.data))
meta.data$group <- factor(meta.data$group, levels=c("W50.G",
                                                    "M1.A",
                                                    "W50.A",
                                                    "C3.A",
                                                    "C3.G",
                                                    "M1.G"))

# end of raw counts and meta data ===========


# DESeq2 and PCA ===============
count.data.set <- DESeqDataSetFromMatrix(countData=raw.counts, 
                                         colData=meta.data, design= ~ group) 
# count.data.set <- DESeqDataSetFromMatrix(countData=raw.counts, colData=meta.data, design= ~ treatment_cell) # also checked DEGs between cell lines

# Filter low count
nrow(count.data.set) #22220

keep <- rowSums(counts(count.data.set)>5) >= nrow(meta.data)*0.25  # 12729 genes
keep <- rowSums(counts(count.data.set)>5) >= nrow(meta.data)*0.9 # 10926 genes
keep <- rowSums(counts(count.data.set)>5) == nrow(meta.data) # 10484 genes

count.filter <- count.data.set[keep,]
nrow(count.filter) 

# create DESeq object
count.data.set.object <- DESeq(count.filter)

vsd <- vst(count.data.set.object)
norm.data = assay(vsd)
head(norm.data)
write.table(norm.data, sep="\t",file="output/Norm_data_H1_new_all5.txt", row.names=TRUE,col.names=NA,quote=FALSE)

# cluster
sampleDists <- dist(t(norm.data),  method = "euclidean") # "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters=hclust(sampleDists)
plot(clusters)

# PCA
# modify the PCA function to change format --------
getMethod("plotPCA","DESeqTransform")
plotPCA.format <- function (object, ...) 
{
  .local <- function (object, intgroup = "condition", 
                      ntop = 500, returnData = FALSE) 
  {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
      colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                    intgroup.df, name = colnames(object))
    if (returnData) {
      attr(d, "percentVar") <- percentVar[1:2]
      return(d)
    }
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", 
                                fill = "group")) + geom_point(size = 5, shape=21,alpha=0.6) + 
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
      coord_fixed() #+  geom_label_repel((aes(label=sample)))
  }
  .local(object, ...)
}
myplot = plotPCA.format(vsd,intgroup=c("Protocols"))
names = colData(vsd)$sample
myplot+  geom_label_repel((aes(label=names)))

myplot = update(myplot, panel = function(x, y, ...) {
  panel.xyplot(x, y, ...);
  ltext(x=x, y=y, labels=names, pos=1, offset=1, cex=0.8)
})
#------------------------------------------------------------------------

# use color blind friendly palette

cbPalette <- c("#E69F00", #lightorange
               "#56B4E9", #blue
               "#D55E00", #darkorange
               "#009E73", #green
               "#CC79A7", #magenta
               "#0072B2", #darkblue
               "#F0E442", #yellow
               "#999999" #grey
               )

#
library(ggrepel) # https://ggrepel.slowkow.com/articles/examples.html
p.pca <- plotPCA.format(vsd, intgroup=c("group"))+ 
  geom_text_repel(aes(label=colData(vsd)$sample),size=3,
                  color="grey50",
                  box.padding   = 0.4,
                  point.padding = 0,
                  #force=1,
                  #force_pull=10,
                  max.overlaps = Inf, # always show all label, regardless of overlap
                  #min.segment.length = 0, # always draw line
                  segment.color = 'grey50')+
  scale_fill_manual(values=cbPalette) +
  theme_bw()+
  theme(panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x  = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_blank()
        )+
  theme(aspect.ratio=1/1)
p.pca

ggsave(filename="figures/PCA_H1_new_all5.png",width=18,height=15,units="cm",dpi=400)

#
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pDist<-pheatmap(sampleDistMatrix,
                clustering_distance_rows = sampleDists,
                clustering_distance_cols = sampleDists,
                col = colors)
save_pheatmap_png <- function(x, filename, width=800, height=700, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_png(pDist, "figures/pDist_H1_new_all5.png")
dev.off()

# end of DESeq2 and PCA ===========



# DEG ============

#resultsNames(count.data.set.object)

results_pairwise <- function(exp,ctrl){
  res <- results(count.data.set.object, contrast=c("group",exp,ctrl),alpha=0.05)
  # save result summary
  summary(res) 
  out <- capture.output(summary(res))
  cat(paste0(exp," vs ",ctrl), out, file="output/results_summary_H1_new.txt", sep="\n", append=TRUE)
  
  # save results - all genes
  res = na.omit(res) # omit NA
  res = res[order(res$padj),] 
  write.table(res, sep="\t",file=paste0("output/Results_",exp,"_",ctrl,"_H1_new.txt"), row.names=TRUE,col.names=NA,quote=FALSE)
  
  # save results - significant genes (DEG)
  res.sig <<- res[res$padj <= 0.05,]
  write.table(res.sig, sep="\t",file=paste0("output/Results_DE_",exp,"_",ctrl,"_H1_new.txt"), row.names=TRUE,col.names=NA,quote=FALSE)
  
  colnames(res.sig) <- paste0(colnames(res.sig),"_",exp,".",ctrl)
  res.sig.sub <<- res.sig[,c(2,6)]
}
results_pairwise("M1.G","C3.G")
res.de.65 <- res.sig.sub %>% as.data.frame()
head(res.de.65)
results_pairwise("M1.G","W50.G")
res.de.64 <- res.sig.sub %>% as.data.frame()
results_pairwise("M1.G","M1.A")
res.de.63 <- res.sig.sub %>% as.data.frame()
results_pairwise("M1.G","C3.A")
res.de.62 <- res.sig.sub %>% as.data.frame()
results_pairwise("M1.G","W50.A")
res.de.61 <- res.sig.sub %>% as.data.frame()

##
results_pairwise("C3.G","W50.G")
res.de.54 <- res.sig.sub %>% as.data.frame()
results_pairwise("C3.G","M1.A")
res.de.53 <- res.sig.sub %>% as.data.frame()
results_pairwise("C3.G","C3.A")
res.de.52 <- res.sig.sub %>% as.data.frame()
results_pairwise("C3.G","W50.A")
res.de.51 <- res.sig.sub %>% as.data.frame()
results_pairwise("W50.G","M1.A")
res.de.43 <- res.sig.sub %>% as.data.frame()
results_pairwise("W50.G","C3.A")
res.de.42 <- res.sig.sub %>% as.data.frame()
results_pairwise("W50.G","W50.A")
res.de.41 <- res.sig.sub %>% as.data.frame()
results_pairwise("M1.A","C3.A")
res.de.32 <- res.sig.sub %>% as.data.frame()
results_pairwise("M1.A","W50.A")
res.de.31 <- res.sig.sub %>% as.data.frame()
results_pairwise("C3.A","W50.A")
res.de.21 <- res.sig.sub %>% as.data.frame()

##
res.merge<- res.de.65 %>% merge(res.de.64, by.x='row.names', by.y="row.names",all=TRUE) %>% 
  merge(res.de.63, by.x= "Row.names", by.y="row.names", all=TRUE) %>% 
  merge(res.de.62, by.x= "Row.names", by.y="row.names", all=TRUE) %>% 
  merge(res.de.61, by.x= "Row.names", by.y="row.names", all=TRUE) %>%
  merge(res.de.54, by.x= "Row.names", by.y="row.names", all=TRUE) %>%
  merge(res.de.53, by.x= "Row.names", by.y="row.names", all=TRUE) %>%
  merge(res.de.52, by.x= "Row.names", by.y="row.names", all=TRUE) %>%
  merge(res.de.51, by.x= "Row.names", by.y="row.names", all=TRUE) %>%
  merge(res.de.43, by.x= "Row.names", by.y="row.names", all=TRUE) %>%
  merge(res.de.42, by.x= "Row.names", by.y="row.names", all=TRUE) %>%
  merge(res.de.41, by.x= "Row.names", by.y="row.names", all=TRUE) %>%
  merge(res.de.32, by.x= "Row.names", by.y="row.names", all=TRUE) %>%
  merge(res.de.31, by.x= "Row.names", by.y="row.names", all=TRUE) %>%
  merge(res.de.21, by.x= "Row.names", by.y="row.names", all=TRUE) %>%
  rename(symbol=Row.names)

View(res.merge)
dim(res.merge) #5759 genes

# end of DEG ===========



# convert gene names and IDs for all DE genes ========

# using biomart -------

#http://www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#introduction
# ###  strip components from the probe names
# # IDlist <- sub('X', '',  featureNames(IDlist)) 
# # ensembl_genes <- gsub("_at","",fdr_df$gene)
# head(listAttributes(useDataset(dataset = "hsapiens_gene_ensembl", mart= useMart("ENSEMBL_MART_ENSEMBL", host    = "www.ensembl.org"))), 120)
# head(listFilters(useDataset(dataset = "hsapiens_gene_ensembl", mart= useMart("ENSEMBL_MART_ENSEMBL", host    = "www.ensembl.org"))), 50)

library(biomaRt)

# read in file composed of exon, trans or noncode ID's
#IDlist <- as.matrix(read.table("gestalt_biomart.csv", sep="\t", header=FALSE))
#filename <- 'gestalt_biomart_annotation_check'

# connect to ensembl
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart= useMart("ENSEMBL_MART_ENSEMBL", host= "https://www.ensembl.org"))

annotated_IDs <- getBM (attributes = c("ensembl_gene_id", "external_gene_name", 
                                       "description", "gene_biotype","entrezgene_id"
                                       #"name_1006","namespace_1003"
                                       ),  #"ensembl_transcript_id",
                        filters = "external_gene_name", 
                        values = res.merge$symbol, 
                        mart= ensembl, 
                        useCache = FALSE) 
View(annotated_IDs)

#View(listAttributes(ensembl)) # check attribute names

dup=annotated_IDs %>% group_by(external_gene_name) %>% dplyr::filter(n() > 1) 
View(dup) # some gene symbols have multiple Ensembl. Related to this question https://www.biostars.org/p/119540/

annotated_IDs.uniuqe <- annotated_IDs[!duplicated(annotated_IDs$external_gene_name),]
#write.table (annotated_IDs, file=paste(filename, ".csv"), sep=",", quote=T, row.names = FALSE, col.names =TRUE)

# using org.Hs.eg.db-----------

res.merge$entrez.s <- mapIds(org.Hs.eg.db, keys=res.merge$symbol, column="ENTREZID", keytype="SYMBOL", multiVals="first")
res.merge$entrez.a <- mapIds(org.Hs.eg.db, keys=res.merge$symbol, column="ENTREZID", keytype="ALIAS", multiVals="first")
res.merge<- res.merge %>% mutate(entrez=coalesce(entrez.s,entrez.a))

res.merge$ensembl.s <- mapIds(org.Hs.eg.db, keys=res.merge$symbol, column="ENSEMBL", keytype="SYMBOL", multiVals="first")
res.merge$ensembl.a <- mapIds(org.Hs.eg.db, keys=res.merge$symbol, column="ENSEMBL", keytype="ALIAS", multiVals="first")
res.merge<- res.merge %>% mutate(ensembl=coalesce(ensembl.s,ensembl.a)) 

dim(res.merge) #5759
length(unique(res.merge$symbol)) #5759
length(unique(annotated_IDs$external_gene_name)) #5462
length(unique(annotated_IDs$entrezgene_id)) #5438

length(unique(annotated_IDs$ensembl_gene_id)) #5883
length(unique(annotated_IDs.uniuqe$ensembl_gene_id)) #5462

length(unique(res.merge$ensembl.s)) #5442
length(unique(res.merge$ensembl.a)) #5674
length(unique(res.merge$ensembl)) # 5715

length(unique(res.merge$entrez.s)) #5476
length(unique(res.merge$entrez.a)) #5711
length(unique(res.merge$entrez)) #5752
# org.Hs.eg.db can map more entrez and ensembl than biomart. Used org.Hs.eg.db 

res.merge$description <- mapIds(org.Hs.eg.db, keys=res.merge$symbol, column="GENENAME", keytype="SYMBOL", multiVals="first")
res.merge$description <- gsub(x = res.merge$description, pattern = "\\ ",replacement = "_") 

res.merge <- res.merge %>% dplyr::select(-entrez.a,-entrez.s,-ensembl.a,-ensembl.s)

View(res.merge)
write.table(res.merge, sep="\t",file="output/Results_DE_merge_H1_new_all5.txt", row.names=TRUE,col.names=NA,quote=FALSE)

# end of converting IDs ========



# overall upset plot of all DE genes ===============

#res.65 <- res.merge %>% select(log2FoldChange_M1.G.C3.G, symbol) %>% filter(!is.na(log2FoldChange_M1.G.C3.G))
#res.64 <- res.merge %>% select(log2FoldChange_M1.G.W50.G, symbol) %>% filter(!is.na(log2FoldChange_M1.G.W50.G))
#res.63 <- res.merge %>% select(log2FoldChange_M1.G.M1.A, symbol) %>% filter(!is.na(log2FoldChange_M1.G.M1.A))
#res.62 <- res.merge %>% select(log2FoldChange_M1.G.C3.A, symbol) %>% filter(!is.na(log2FoldChange_M1.G.C3.A))
#res.61 <- res.merge %>% select(log2FoldChange_M1.G.W50.A, symbol) %>% filter(!is.na(log2FoldChange_M1.G.W50.A))
#--------------
listInput0 <- list("M1.G vs C3.G"=row.names(res.de.65),
                   "M1.G vs W50.G"=row.names(res.de.64),
                   "M1.G vs M1.A"=row.names(res.de.63),
                   "M1.G vs C3.A"=row.names(res.de.62),
                   "M1.G vs W50.A"=row.names(res.de.61)
                   )
View(listInput0)
upset(fromList(listInput0),nsets=10, 
      nintersects = NA,
      mb.ratio = c(0.6, 0.4), 
      text.scale = c(2,2, 2, 1.5, 2.5, 2), point.size = 6, line.size = 2,  
      order.by = c("degree"),keep.order = T, sets = c("M1.G vs W50.G", "M1.G vs M1.A",
                                                      "M1.G vs W50.A","M1.G vs C3.A",
                                                      "M1.G vs C3.G"))

#ggsave(filename="figures/upset_DE_H1_new.png",width=6.4,height=4,units="in",dpi=400)

png(file = "figures/upset_DE_H1_new_all5.png",
    width = 25, 
    height = 18, 
    units = "cm", res = 400)
# plot upset above
dev.off()

# end ==========



# overall heatmap of all DE genes ===============

# old code using pheatmap -------------
#norm.de <- norm[,c(5,6,7,8,11,12,9,10,3,4,1,2)] %>% filter(rownames(norm) %in% c(res.merge$symbol)) 
norm.M.de <- norm.M %>% filter(rownames(norm.M) %in% c(res.merge.M$symbol)) 
norm.MH.de <- norm.MH %>% filter(rownames(norm.MH) %in% c(res.merge.MH$symbol)) 
norm.no.cell <- norm.MH %>% filter(!rownames(norm.MH) %in% c(res.merge.cell$symbol)) 
norm.cell <- norm %>% filter(rownames(norm.MH) %in% c(res.merge.cell$symbol)) 
dim(norm.no.cell)
dim(norm.MH.de)
dim(norm.cell)
dim(norm.M.de)
head(norm.MH.de)
df_num = as.matrix(norm.M.de)
df_num = as.matrix(norm.MH.de)
df_num = as.matrix(norm.no.cell)
df_num = as.matrix(norm.cell[,c(1:12)])
pheatmap(df_num, #main="",
         #labels_row = norm.cor$symbol,
         scale = "row",
         cluster_cols = T,
         cluster_rows = T,
         show_rownames = F,
         show_colnames = T,
         legend = T,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         #cellheight = 0.12, 
         cellwidth = 14,
         fontsize = 8,
         border_color = NA ,filename = "heatmap_all_genes_cluster_MEL1.pdf",width = 4,height = 7 # with = 4 or 5, height = 6 or 7
         #gaps_col = c(2,4,6,8,10)
         #clustering_method = "ward.D2", #"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
         #clustering_distance_rows= #"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
)
dev.off()
# end -------------------

# old code for PCA?-----------
#install.packages("ggfortify")
#library(ggfortify)
#df <- iris[1:4]
#view(df)
#pca_res <- prcomp(df, scale. = TRUE)
#autoplot(pca_res)
#pca_res <- prcomp(norm.MH.de, scale. = TRUE)
#autoplot(pca_res)
# end -----------------

# complex heatmap for detected and DE genes -------

View(norm.data)
m<- norm.data

m<- as.data.frame(norm.data) %>% filter(row.names(norm.data) %in% res.merge$symbol)
dim(m)
#m <- as.matrix(as.data.frame(lapply(m.anno, as.numeric),check.names=F))
View(m)
m.z <- t(scale(t(m))) #%>% as.data.frame()
View(m.z)
colnames(m.z)

m.z[m.z>4] <- NA
m.z <- t(scale(t(m.z)))
m.z[is.na(m.z)] <- 4
m.z[m.z>4] <- 4
max(abs(m.z))

colnames(m.z)
#summary(m.z)
number_of_W50G <- length(grep("W50.G", colnames(m.z)))
number_of_M1A <- length(grep("M1.A", colnames(m.z)))
number_of_W50A <- length(grep("W50.A", colnames(m.z)))
number_of_C3A <- length(grep("C3.A", colnames(m.z)))
number_of_C3G <- length(grep("C3.G", colnames(m.z)))
number_of_M1G <- length(grep("M1.G", colnames(m.z)))

end_index_W50G <- grep("W50.G", colnames(m.z))[number_of_W50G]
end_index_M1A <- grep("M1.A", colnames(m.z))[number_of_M1A]
end_index_W50A <- grep("W50.A", colnames(m.z))[number_of_W50A]
end_index_C3A <- grep("C3.A", colnames(m.z))[number_of_C3A]
end_index_C3G <- grep("C3.G", colnames(m.z))[number_of_C3G]
end_index_M1G <- grep("M1.G", colnames(m.z))[number_of_M1G]

number_of_M1G
end_index_M1G

###

# check colors of ggplot2 http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
library(scales)
#extract hex color codes for a plot with 6 elements in ggplot2 
hex <- hue_pal()(6)
#overlay hex color codes on actual colors
show_col(hex)

# use color blind friendly palette

cbPalette <- c("#E69F00", #lightorange
               "#56B4E9", #blue
               "#D55E00", #darkorange
               "#009E73", #green
               "#CC79A7", #magenta
               "#0072B2", #darkblue
               "#F0E442", #yellow
               "#999999" #grey
)
show_col(cbPalette)
show_col(c("thistle1", "orchid4","purple3"))
#
col_fun = colorRamp2(c(1, 6), c("thistle1", "orchid4"))
brewer.pal(n = 3, name = "Purple")

ha = HeatmapAnnotation(rank = c(rep(1,5),rep(2,3),rep(3,3),rep(4,5),rep(5,3),rep(6,5)),
                       col = list(rank = col_fun)
                       )


heatmap.H1 <- Heatmap(m.z, #matrix_Z_score_total,
                          name = "Z score",
                          show_row_names = FALSE,
                          show_column_names = FALSE,
                          
                          # row_labels = gt_render(m.anno$protein.symbol),
                          row_names_gp = gpar(fontsize = 8),
                          column_names_side = "top",
                          column_dend_side = "bottom",
                      
                      cluster_rows = TRUE,
                      cluster_row_slices = TRUE,
                      clustering_distance_rows = "euclidean",
                      clustering_method_rows = "ward.D2",
                      row_dend_side = "left",
                      row_dend_width = unit(20, "mm"),
                      show_row_dend = TRUE,
          
                      cluster_columns = FALSE,
                      cluster_column_slices = FALSE,
                      # column_split = col.split,
                      column_title = NULL,
                      
                      
                      #row_km = 7,
                      
                      #split=3,
                      
                      bottom_annotation = ha,
                      
                          layer_fun = function(j, i, x, y, width, height, fill) {
                            mat = restore_matrix(j, i, x, y)
                            ind = unique(c(mat[, c(end_index_W50G, 
                                                   end_index_M1A, 
                                                   end_index_W50A, 
                                                   end_index_C3A,
                                                   end_index_C3G
                                                   )
                                               ]
                                           )
                                         )
                            grid.rect(x = x[ind] + unit(0.5/ncol(m.z), "npc"), 
                                      y = y[ind], 
                                      width = unit(0.02, "inches"), 
                                      height = unit(1/nrow(m.z), "npc"),
                                      gp = gpar(col = "white")
                            )
                          },
                          col = colorRamp2(c(-3,0,3), 
                                           c("blue", "white", "red")),
                          top_annotation = columnAnnotation(empty = anno_empty(border = FALSE
                                                                               , height = unit(12, "mm")
                          )),
                          #       left_annotation = 
                          #         rowAnnotation(`Significant CW vs CS` = m.anno$CS_CW,
                          #                       
                          #                       `Significant KW vs KS` = m.anno$KS_KW,
                          
                          #                       `Significant CW vs KW` = m.anno$KW_CW,
                          #                       
                          #                       `Significant CS vs KS` = m.anno$KS_CS,
                          
                          #                       col = list(`Significant CW vs CS` = c("TRUE" = hue_pal()(4)[1]),
                          #                                  `Significant KW vs KS` = c("TRUE" = hue_pal()(4)[2]),
                          #                                  `Significant CW vs KW` = c("TRUE" = hue_pal()(4)[3]),
                          #                                  `Significant CS vs KS` = c("TRUE" = hue_pal()(4)[4])
                          #                       ),
                          #                       annotation_legend_param = list(
                          #                         
                          #                         `Significant CW vs CS` = list(
                          #                           title = c(""),
                          #                           labels = expression(italic("Ins1")^italic("+/+")*" Sucrose vs Water") #expression() issue - https://github.com/jokergoo/ComplexHeatmap/issues/678
                          #                         ),
                          #                         
                          #                         `Significant KW vs KS` = list(
                          #                           title = c(""),
                          #                           labels = expression(italic("Ins1")^italic("+/-")*" Sucrose vs Water")
                          #                         ),
                          
                          #                         `Significant CW vs KW` = list(
                          #                           title = c(""),
                          #                           labels = expression("Water "*italic("Ins1")^italic("+/-")*" vs "*italic("Ins1")^italic("+/+")*"")
                          #                         ),
                          
                          #                         `Significant CS vs KS` = list(
                          #                           title = c(""),
                          #                           labels = expression("Sucrose "*italic("Ins1")^italic("+/-")*" vs "*italic("Ins1")^italic("+/+")*"")
                          #                         )
                          #                       ),
                          #                       show_annotation_name = FALSE,
                          #                       show_legend=FALSE,
                          #                       na_col = "white"
                          #                       
                          #         ),
                          column_order = 1:ncol(m.z),
                          height = 
                            
                            
                            unit(150, "mm"), 
                          
                          width = ncol(m.z)*unit(5, "mm"),
                          border_gp = gpar(col = "grey"),
                          show_heatmap_legend = TRUE,
                          heatmap_legend_param = list(
                            title = "Z-score",
                            title_position = "topleft",
                            legend_height = unit(4, "cm")))

draw(heatmap.H1)
#

png(file = "figures/heatmap_H1_correlated.png",
     width = 7, 
     height = 7, 
     units = "in", res = 600)
pdf(file = "figures/heatmap_H1_correlated.pdf",
    width = 7, 
    height = 7)

png(file = "figures/heatmap_H1_DE.png",
    width = 7, 
    height = 7, 
    units = "in", res = 600)

draw(heatmap.H1)

dev.off()
#

seekViewport("annotation_empty_1")
loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))


####Condition labels

#Condition label 1

grid.rect(x = (loc2$x - loc1$x)*(end_index_W50G)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_W50G)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[1],0.5),
                    col = "white"
          )
)
grid.text(expression("W50.G"), 
          x = (loc2$x - loc1$x)*(end_index_W50G)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11,
                    col = "black"))

#Condition label 2
grid.rect(x = (loc2$x - loc1$x)*(end_index_W50A + 
                                   end_index_W50G)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_W50A)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[3],0.5),
                    col ="white"
          )
)
grid.text(expression("W50.A"), 
          x = (loc2$x - loc1$x)*(end_index_W50A + 
                                   end_index_W50G)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))


#Condition label 3
grid.rect(x = (loc2$x - loc1$x)*(end_index_M1A + 
                                   end_index_W50A)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_M1A)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[2],0.5),
                    col = "white"
          )
)
grid.text(expression("M1.A"), 
          x = (loc2$x - loc1$x)*(end_index_M1A + 
                                   end_index_W50A)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))


#Condition label 4
grid.rect(x = (loc2$x - loc1$x)*(end_index_C3G +
                                   end_index_M1A)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_C3G)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[5],0.5),
                    col = "white"
          )
)
grid.text(expression("C3.G"), 
          x = (loc2$x - loc1$x)*(end_index_C3G +
                                   end_index_M1A)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))

#Condition label 5


grid.rect(x = (loc2$x - loc1$x)*(end_index_C3A +
                                   end_index_C3G)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_C3A)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[4],0.5),
                    col = "white"
          )
)
grid.text(expression("C3.A"), 
          x = (loc2$x - loc1$x)*(end_index_C3A +
                                   end_index_C3G)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))

#Condition label 6

grid.rect(x = (loc2$x - loc1$x)*(end_index_M1G + 
                                   end_index_C3A)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_M1G)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[6],0.5),
                    col = "white"
          )
)
grid.text(expression("M1.G"), 
          x = (loc2$x - loc1$x)*(end_index_M1G + 
                                   end_index_C3A)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))

###Top label gaps


#Vertical lines
grid.rect(x = end_index_W50G/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y)/2,
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_M1A/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y)/2,
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_W50A/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y)/2,
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_C3A/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y)/2,
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_C3G/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y)/2,
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

# add white vertical lines to bottom annotation (rank)
#list_components()
seekViewport("annotation_rank_1")
loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))

#Vertical lines
grid.rect(x = end_index_W50G/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y),
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_M1A/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y),
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_W50A/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y),
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_C3A/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y),
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_C3G/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y),
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))
dev.off()

# end of heatmap =======



#=============================================================================================

# Correlation analysis ==========

#norm.M <- read.table(file="output/Norm_data_MEL1.txt",row.names=1)

#res.merge.M <- read.table(file="output/Results_DE_merge_MEL1.txt",row.names=1)
#norm.M.de <- norm.M %>% filter(rownames(norm.M) %in% c(res.merge.M$symbol)) 

norm.data <- read.table("output/Norm_data_H1_new_all5.txt")
View(norm.data)
cn <- colnames(norm.data)
norm.data <- norm.data %>% select(c(grep("W50.G", cn), grep("W50.A", cn), grep("M1.A",cn), 
                       grep("C3.G",cn), grep("C3.A",cn), grep("M1.G",cn)))

correlation_function<- function(df){
  r <- (1:ncol(df))
  df[ , r] <- apply(df[ , r], 2,    # Specify own function within apply, character to numeric
                    function(x) as.numeric(as.character(x)))
  df.t = as.data.frame(t(df)) #convert row to column
  #rank <- as.numeric(as.character(c(6,6,5,5,1,1,2,2,4,4,3,3)))
  #rank <- as.numeric(as.character(c(6,6,6,6,5,5,5,5,1,1,1,1,2,2,4,4,3,3)))
  rank <- c(rep(1,5),rep(2,3),rep(3,3),rep(4,5),rep(5,3),rep(6,5))
  n <- ncol(df.t)
  correlations <-lapply(1:n, function(x) (cor.test(rank,df.t[,x],method='spearman')))
  #r value
  Rho <- sapply(correlations,function(x) x[["estimate"]])
  #p value
  p<- sapply(correlations,function(x) x[["p.value"]])
  corr<-data.frame(Rho=Rho,p=p,Rho2=Rho^2)
  rownames(corr) <- colnames(df.t)
  corr$adj.p <- p.adjust(corr$p, method = 'BH', n = length(corr$p)) # p.adjust.methods c("holm", "hochberg", "hommel", "bonferroni", "BH" (same as "fdr"), "BY", none")
  #corr$symbol <- row.names(corr)
  #corr$entrez = mapIds(org.Hs.eg.db, keys=rownames(corr), column="ENTREZID", keytype="ALIAS", multiVals="first")
  #corr$description = mapIds(org.Hs.eg.db, keys=rownames(corr), column="GENENAME", keytype="ALIAS", multiVals="first")
  corr <<- corr[order(corr$Rho2,decreasing=T),]
  head(corr)
  #write.xlsx2(corr, file="output/correlated_gene_test.xlsx", sheetName = paste(deparse(substitute(df))),col.names = TRUE, row.names = TRUE, append = TRUE)
  #write.xlsx2(corr, file=paste0("output/correlated_gene.xlsx"),col.names = TRUE, row.names = TRUE, append = TRUE)
  #write.xlsx2(corr, file=paste0("output/correlated_gene_MH_noCell.xlsx"),col.names = TRUE, row.names = TRUE, append = TRUE)
  
  corr.sig <<- corr %>% filter(adj.p<0.05)
  corr.sig.ua <<- corr %>% filter(p<0.05)
  #df.order <<- df[,c(5,6,7,8,11,12,9,10,3,4,1,2)]
  #df.cor <- corr.sig %>% merge(df.order,by.x='symbol',by.y='row.names',all.x=TRUE)
  #df.cor <<- df.cor[order(df.cor$Rho,decreasing = T),]
  
  #df.cor.ua <- corr.sig.ua %>% merge(df.order,by.x='symbol',by.y='row.names',all.x=TRUE)
  #df.cor.ua <<- df.cor.ua[order(df.cor.ua$Rho,decreasing = T),]
  
  #write.xlsx2(df.cor, file="output/correlated_gene_heatmap.xlsx",col.names = TRUE, row.names = TRUE, append = TRUE)
}

correlation_function(norm.data)
rank.cor <- corr
rank.sig <- corr.sig
rank.sig.ua <- corr.sig.ua
View(rank.sig)

#dup=rank.sig.ua.M %>% group_by(entrez) %>% dplyr::filter(n() > 1) #has duplicated Entrez ID!!
#head(dup)
#write.txt(as.data.frame(dup), file="output/dup_cor_ua_M.xlsx",col.names = TRUE, row.names = TRUE, append = F)

dim(rank.sig) # 1767 genes adj. p<0.05
dim(rank.sig.ua) # 3320 genes unadj. p<0.05
dim(res.merge) # 5759 DE genes
dim(norm.data) # 10484 reliably detected genes (>5 counts all samples)

#write.csv(rank.sig, "output/H1_correlated_genes.csv", row.names = F)

# end of correlation analysis =========


# heatmap of correlated genes ===========

m <- rownames_to_column(rank.sig[order(rank.sig$Rho,decreasing=T),]) %>% 
  left_join(rownames_to_column(as.data.frame(norm.data))) %>%
  column_to_rownames() %>% select(!c(1:4))
View(m)

# heatmap for selected markers
corr.marker <-  rank.cor %>% filter(row.names(rank.cor) %in% c( "FOXA2", "SOX17", "CXCR4", "CD177", 
                                               "KIT", "CER1", "GSC", "HHEX", "MIXL1", 
                                               "EOMES", "NODAL", "GATA4", "GATA6", 
                                               "LGR5", "SMAD2", "HNF1B", "HNF4A", "GATA3", 
                                               "KRT18", "KRT8", "FN1", "LINC00458")) 
View(corr.marker)

#write.csv(corr.marker, "output/H1_markers_correlation.csv", row.names = F)

corr.marker <- corr.marker[order(corr.marker$Rho,decreasing = T),]

m <- rownames_to_column(corr.marker) %>% 
  left_join(rownames_to_column(as.data.frame(norm.data))) %>%
  column_to_rownames() %>% select(!c(1:4))

#m <- as.matrix(as.data.frame(lapply(m.anno, as.numeric),check.names=F))
View(m)
m.z <- t(scale(t(m))) #%>% as.data.frame()
View(m.z)
colnames(m.z)

#m.z[m.z>4] <- NA
#m.z <- t(scale(t(m.z)))
#m.z[is.na(m.z)] <- 4
#m.z[m.z>4] <- 4
max(abs(m.z))

colnames(m.z)
#summary(m.z)
number_of_W50G <- length(grep("W50.G", colnames(m.z)))
number_of_M1A <- length(grep("M1.A", colnames(m.z)))
number_of_W50A <- length(grep("W50.A", colnames(m.z)))
number_of_C3A <- length(grep("C3.A", colnames(m.z)))
number_of_C3G <- length(grep("C3.G", colnames(m.z)))
number_of_M1G <- length(grep("M1.G", colnames(m.z)))

end_index_W50G <- grep("W50.G", colnames(m.z))[number_of_W50G]
end_index_M1A <- grep("M1.A", colnames(m.z))[number_of_M1A]
end_index_W50A <- grep("W50.A", colnames(m.z))[number_of_W50A]
end_index_C3A <- grep("C3.A", colnames(m.z))[number_of_C3A]
end_index_C3G <- grep("C3.G", colnames(m.z))[number_of_C3G]
end_index_M1G <- grep("M1.G", colnames(m.z))[number_of_M1G]

number_of_M1G
end_index_M1G

###

# use color blind friendly palette

cbPalette <- c("#E69F00", #lightorange
               "#56B4E9", #blue
               "#D55E00", #darkorange
               "#009E73", #green
               "#CC79A7", #magenta
               "#0072B2", #darkblue
               "#F0E442", #yellow
               "#999999" #grey
)
show_col(cbPalette)
show_col(c("thistle1", "orchid4","purple3"))
#
col_fun = colorRamp2(c(1, 6), c("thistle1", "orchid4"))
brewer.pal(n = 3, name = "Purple")

ha = HeatmapAnnotation(rank = c(rep(1,5),rep(2,3),rep(3,3),rep(4,5),rep(5,3),rep(6,5)),
                       col = list(rank = col_fun)
)


heatmap.H1 <- Heatmap(m.z, #matrix_Z_score_total,
                      name = "Z score",
                      
                      #show_row_names = FALSE,
                      
      # change here !!!
                      show_row_names = #FALSE, 
                                      TRUE,
                      
                      show_column_names = FALSE,
                      
                      #row_labels = gt_render(rownames(m.z)),
                      row_names_gp = gpar(fontsize = 13),
                      column_names_side = "top",
                      column_dend_side = "bottom",
                      
                      cluster_rows = FALSE,
                      cluster_row_slices = TRUE,
                      clustering_distance_rows = "euclidean",
                      clustering_method_rows = "ward.D2",
                      row_dend_side = "left",
                      row_dend_width = unit(20, "mm"),
                      show_row_dend = TRUE,
                      
                      cluster_columns = FALSE,
                      cluster_column_slices = FALSE,
                      # column_split = col.split,
                      column_title = NULL,
                      
                      
                      #row_km = 7,
                      
                      #split=3,
                      
                      bottom_annotation = ha,
                      
                      layer_fun = function(j, i, x, y, width, height, fill) {
                        mat = restore_matrix(j, i, x, y)
                        ind = unique(c(mat[, c(end_index_W50G, 
                                               end_index_M1A, 
                                               end_index_W50A, 
                                               end_index_C3A,
                                               end_index_C3G
                        )
                        ]
                        )
                        )
                        grid.rect(x = x[ind] + unit(0.5/ncol(m.z), "npc"), 
                                  y = y[ind], 
                                  width = unit(0.02, "inches"), 
                                  height = unit(1/nrow(m.z), "npc"),
                                  gp = gpar(col = "white")
                        )
                      },
                      col = colorRamp2(c(-3,0,3), 
                                       c("blue", "white", "red")),
                      top_annotation = columnAnnotation(empty = anno_empty(border = FALSE
                                                                           , height = unit(12, "mm")
                      )),
                      
                      column_order = 1:ncol(m.z),
                      height = 
                        
     # change here !!!                   
                        unit(
                      #    150,
                          100, 
                             "mm"), 
                      
                      width = ncol(m.z)*unit(5, "mm"),
                      border_gp = gpar(col = "grey"),
                      show_heatmap_legend = TRUE,
                      heatmap_legend_param = list(
                        title = "Z-score",
                        title_position = "topleft",
                        legend_height = unit(4, "cm")))

draw(heatmap.H1)
#

png(file = "figures/heatmap_H1_correlated.png",
    width = 7, 
    height = 7, 
    units = "in", res = 600)

pdf(file = "figures/heatmap_H1_correlated.pdf",
    width = 7, 
    height = 7)


png(file = "figures/heatmap_H1_markers.png",
    width = 7, 
    height = 5, 
    units = "in", res = 600)
pdf(file = "figures/heatmap_H1_markers.pdf",
    width = 7, 
    height = 5)

draw(heatmap.H1)

dev.off()
#

seekViewport("annotation_empty_1")
loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))


####Condition labels


####Condition labels

#Condition label 1

grid.rect(x = (loc2$x - loc1$x)*(end_index_W50G)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_W50G)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[1],0.5),
                    col = "white"
          )
)
grid.text(expression("W50.G"), 
          x = (loc2$x - loc1$x)*(end_index_W50G)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11,
                    col = "black"))

#Condition label 2
grid.rect(x = (loc2$x - loc1$x)*(end_index_W50A + 
                                   end_index_W50G)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_W50A)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[3],0.5),
                    col ="white"
          )
)
grid.text(expression("W50.A"), 
          x = (loc2$x - loc1$x)*(end_index_W50A + 
                                   end_index_W50G)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))


#Condition label 3
grid.rect(x = (loc2$x - loc1$x)*(end_index_M1A + 
                                   end_index_W50A)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_M1A)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[2],0.5),
                    col = "white"
          )
)
grid.text(expression("M1.A"), 
          x = (loc2$x - loc1$x)*(end_index_M1A + 
                                   end_index_W50A)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))


#Condition label 4
grid.rect(x = (loc2$x - loc1$x)*(end_index_C3G +
                                   end_index_M1A)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_C3G)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[5],0.5),
                    col = "white"
          )
)
grid.text(expression("C3.G"), 
          x = (loc2$x - loc1$x)*(end_index_C3G +
                                   end_index_M1A)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))

#Condition label 5


grid.rect(x = (loc2$x - loc1$x)*(end_index_C3A +
                                   end_index_C3G)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_C3A)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[4],0.5),
                    col = "white"
          )
)
grid.text(expression("C3.A"), 
          x = (loc2$x - loc1$x)*(end_index_C3A +
                                   end_index_C3G)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))

#Condition label 6

grid.rect(x = (loc2$x - loc1$x)*(end_index_M1G + 
                                   end_index_C3A)/2/ncol(m.z),
          y = 0,
          width = (loc2$x - loc1$x)*(number_of_M1G)/ncol(m.z), 
          height = (loc2$y - loc1$y)/2,
          just = c("center", "bottom"),
          gp = gpar(fill = alpha(cbPalette[6],0.5),
                    col = "white"
          )
)
grid.text(expression("M1.G"), 
          x = (loc2$x - loc1$x)*(end_index_M1G + 
                                   end_index_C3A)/2/ncol(m.z),
          y = 0.25,
          just = c("center", "center"),
          gp = gpar(fontsize = 11))
###Top label gaps


#Vertical lines
grid.rect(x = end_index_W50G/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y)/2,
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_M1A/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y)/2,
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_W50A/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y)/2,
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_C3A/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y)/2,
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_C3G/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y)/2,
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

# add white vertical lines to bottom annotation (rank)
#list_components()
seekViewport("annotation_rank_1")
loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))

#Vertical lines
grid.rect(x = end_index_W50G/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y),
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_M1A/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y),
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_W50A/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y),
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_C3A/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y),
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))

grid.rect(x = end_index_C3G/ncol(m.z),
          y = 0,
          height = (loc2$y - loc1$y),
          width = unit(0.02, "inches"),
          just = c("center", "bottom"),
          gp = gpar(fill = "white", 
                    col = "white", 
                    lwd = 1
          ))
dev.off()

# end of heatmap for correlated genes ==========

# convert gene IDs for correlated genes ========

rank.cor$entrez.s <- mapIds(org.Hs.eg.db, keys=row.names(rank.cor), column="ENTREZID", keytype="SYMBOL", multiVals="first")
rank.cor$entrez.a <- mapIds(org.Hs.eg.db, keys=row.names(rank.cor), column="ENTREZID", keytype="ALIAS", multiVals="first")
rank.cor<- rank.cor %>% mutate(entrez=coalesce(entrez.s,entrez.a))
View(rank.cor)
rank.cor$ensembl.s <- mapIds(org.Hs.eg.db, keys=row.names(rank.cor), column="ENSEMBL", keytype="SYMBOL", multiVals="first")
rank.cor$ensembl.a <- mapIds(org.Hs.eg.db, keys=row.names(rank.cor), column="ENSEMBL", keytype="ALIAS", multiVals="first")
rank.cor<- rank.cor %>% mutate(ensembl=coalesce(ensembl.s,ensembl.a)) 

length(unique(row.names(rank.cor))) #10484
length(unique(rank.cor$ensembl)) # 10345
length(unique(rank.cor$entrez)) #10441


rank.cor$description <- mapIds(org.Hs.eg.db, keys=row.names(rank.cor), column="GENENAME", keytype="SYMBOL", multiVals="first")

rank.cor <- rank.cor %>% dplyr::select(-entrez.a,-entrez.s,-ensembl.a,-ensembl.s)

rank.cor$symbol=row.names(rank.cor)

View(rank.cor)

rank.sig <- rank.cor %>% filter(adj.p<0.05)
rank.pos <- rank.sig %>% filter(Rho>0)
rank.neg <- rank.sig %>% filter(Rho<0)

write.csv(rank.cor, file="output/correlated_genes_H1_new.csv", row.names=F)

corr.marker <-  rank.cor %>% filter(row.names(rank.cor) %in% c( "FOXA2", "SOX17", "CXCR4", "CD177", 
                                                                "KIT", "CER1", "GSC", "HHEX", "MIXL1", 
                                                                "EOMES", "NODAL", "GATA4", "GATA6", 
                                                                "LGR5", "SMAD2", "HNF1B", "HNF4A", "GATA3", 
                                                                "KRT18", "KRT8", "FN1", "LINC00458")) 

rank.cor.marker <- rank.cor %>%
  filter(symbol %in% c( "FOXA2", "SOX17", "CXCR4", "CD177", 
                        "KIT", "CER1", "GSC", "HHEX", "MIXL1", 
                        "EOMES", "NODAL", "GATA4", "GATA6", 
                        "LGR5", "SMAD2", "HNF1B", "HNF4A", "GATA3", 
                        "KRT18", "KRT8", "FN1", "LINC00458"))
View(rank.cor.marker)
write.csv(rank.cor.marker, file="output/correlated_genes_H1_new_marker.csv", row.names=F)

# end of converting IDs for correlated genes ==========


# annotate TF and surface markers for correlated genes ============

# read the TF database, downloaded from < http://humantfs.ccbr.utoronto.ca/ > Lambert SA, Jolma A, Campitelli LF, Das PK, Yin Y, Albu M, Chen X, Taipale J, Hughes TR, Weirauch MT.(2018) The Human Transcription Factors. Cell. 172(4):650-665. doi: 10.1016/j.cell.2018.01.029. Review.

#tf.db <- read.table(file = "input/TF_DatabaseExtract_v_1.01.txt", sep = "\t",header = TRUE, row.names = 1)
tf.db <- read_excel("input/TF_DatabaseExtract_v_1.01.xlsx", )[,-1] 
View(tf.db)
tf.yes <- tf.db %>% filter(`Is TF?`=="Yes")

# read the cell compartment database, downloaded from <  https://compartments.jensenlab.org/  >  Binder,J.X., Pletscher-Frankild,S., Tsafou,K. et al. COMPARTMENTS: unification and visualization of protein subcellular localization evidence. Database (2014) Vol. 2014: article ID bau012; doi:10.1093/database/bau012.

#cc.db <- read.table(file = "input/human_compartment_integrated_full.tsv", sep = "\t")
#colnames(cc.db) <- c("id","symbol", "GO", "compartment","confidence")
cc.db <- read.table(file = "input/human_compartment_knowledge_full.tsv", sep = "\t")
View(cc.db)
colnames(cc.db) <- c("id","symbol","GO","compartment","source","evidence","confidence")
cc <- unique(cc.db$compartment)
cc
membrane <- grep("membrane", cc, value=TRUE)
cat("Membrane-related compartments:", membrane, file="output/membrane_related_compartment.txt", sep="\n")

#cc.pm <- cc.db %>% filter(compartment %in% c("Plasma membrane",
#                                             "Integral component of plasma membrane",
#                                             "Intrinsic component of plasma membrane",
#                                             "Plasma membrane region")) # these other plasma membrane compartments are subsets of "plasma membrane"
cc.pm <- cc.db %>% filter(compartment %in% c("Plasma membrane"))
length(unique(cc.pm$symbol)) # 2076 plasma membrane proteins

rank.cor <- rank.cor %>% 
  mutate(TF = ifelse(symbol %in% tf.yes$`HGNC symbol`, "TF","No")) %>%
  mutate(surface = ifelse(symbol %in% cc.pm$symbol, "surface","No"))
View(rank.cor)
rank.cor.TF.surface <- rank.cor %>% 
  filter(TF != "No" | surface != "No") %>%
  filter(adj.p<0.05)
View(rank.cor.TF.surface)
dim(rank.cor.TF.surface) #263
write.csv(rank.cor, file="output/correlated_genes_H1_new.txt", row.names=F)
write.csv(rank.cor.TF.surface, file="output/correlated_genes_H1_new_TF_surface.txt", row.names=F)

# end of annotte TF and surface markers ==========


# GSEA of gene correlation =============

# Option1 - use Rho as rank score
genelist <- rank.cor$Rho

# Option2 - use p value * Rho sign as rank score
genelist <- sign(rank.cor$Rho) * (-log10(rank.cor$adj.p))

# entrez id as names of the gene list
names(genelist) <- rank.cor$entrez
genelist = sort(genelist, decreasing = TRUE)
head(genelist)


# Run GSEA

x <- genelist

gsea <- function(x){
  
  gse.kegg <- gseKEGG(
    geneList=x,
    organism = "hsa",
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    keyType = "kegg"
  )
  gse.kegg <<- setReadable(gse.kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gse.kegg.df <<- as.data.frame(gse.kegg)
  View(gse.kegg.df)
  #write.xlsx2(gse.kegg.df, file="output/GSEA_.xlsx", sheetName = "KEGG",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  gseReactome <- gsePathway(
    geneList=x,
    organism = "human",
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE)
  gseReactome <<- setReadable(gseReactome, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseReactome.df <<- as.data.frame(gseReactome)
  View(gseReactome.df)
  #write.xlsx2(gseReactome.df, file="output/GSEA_.xlsx", sheetName = "Reactome",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  
  
  gseGO.bp <- gseGO(
    geneList=x,
    ont = "BP",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.bp <<- setReadable(gseGO.bp, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.bp.df <<- as.data.frame(gseGO.bp)
  head(gseGO.bp.df)
  #write.xlsx2(gseGO.bp.df, file="output/GSEA_.xlsx", sheetName = "GO_BP",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  
  gseGO.mf <- gseGO(
    geneList=x,
    ont = "MF",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.mf <<- setReadable(gseGO.mf, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.mf.df <<- as.data.frame(gseGO.mf)
  head(gseGO.mf.df)
  #write.xlsx2(gseGO.mf.df, file="output/GSEA_.xlsx", sheetName = "GO_MF",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  
  
  gseGO.cc <- gseGO(
    geneList=x,
    ont = "CC",
    OrgDb = org.Hs.eg.db,
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    verbose = TRUE,
    keyType = "ENTREZID"
  )
  gseGO.cc <<- setReadable(gseGO.cc, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gseGO.cc.df <<- as.data.frame(gseGO.cc)
  head(gseGO.cc.df)
  #write.xlsx2(gseGO.cc.df, file="output/GSEA_.xlsx", sheetName = "GO_CC",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  
  gse.mkegg <- gseMKEGG(
    geneList=x,
    organism = "hsa",
    minGSSize = 15, 
    maxGSSize = 500,
    pvalueCutoff = 1,
    keyType = "kegg"
  )
  gse.mkegg <<- setReadable(gse.mkegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gse.mkegg.df <<- as.data.frame(gse.mkegg)
  head(gse.mkegg)
  #write.xlsx2(gse.mkegg.df, file="output/GSEA_.xlsx", sheetName = "mKEGG",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  
  gse.wp <- gseWP(x, organism = "Homo sapiens",pvalueCutoff = 1)
  gse.wp <<- setReadable(gse.wp, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  gse.wp.df <<- as.data.frame(gse.wp)
  head(gse.wp)
  #write.xlsx2(gse.wp.df, file="output/GSEA_.xlsx", sheetName = "mKEGG",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  
  # mSigDB Hallmark gene sets  
  #msigdbr_show_species()
  m_df <- msigdbr(species = "Homo sapiens")
  head(m_df, 2) %>% as.data.frame
  
  # H: hallmark gene sets
  # C1: positional gene sets
  # C2: curated gene sets
  # C3: motif gene sets
  # C4: computational gene sets
  # C5: GO gene sets
  # C6: oncogenic signatures
  # C7: immunologic signatures
  
  # MSigDb GSEA
  h_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
    dplyr::select(gs_name, entrez_gene)
  head(h_t2g)
  
  gse.hallmark <- GSEA(x, #geneList 
                       TERM2GENE = h_t2g)
  
  gse.hallmark <<- setReadable(gse.hallmark, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  
  gse.hallmark.df <- as.data.frame(gse.hallmark)
  
  View(gse.hallmark.df)
  
  
}

# option 1 rank by Rho
write.csv(gse.kegg.df, file = "output/GSEA_Rho_kegg_H1_new.csv")
write.csv(gseReactome.df, file = "output/GSEA_Rho_reactome_H1_new.csv")
write.csv(gseGO.bp.df, file = "output/GSEA_Rho_GObp_H1_new.csv")
write.csv(gseGO.mf.df, file = "output/GSEA_Rho_GOmf_H1_new.csv")
write.csv(gseGO.cc.df, file = "output/GSEA_Rho_GOcc_H1_new.csv")
write.csv(gse.mkegg.df, file = "output/GSEA_Rho_mkegg_H1_new.csv")
write.csv(gse.wp.df, file = "output/GSEA_Rho_wiki_H1_new.csv")

# option 2 rank by p value sign
write.csv(gse.kegg.df, file = "output/GSEA_p_kegg_H1_new.csv")
write.csv(gseReactome.df, file = "output/GSEA_p_reactome_H1_new.csv")
write.csv(gseGO.bp.df, file = "output/GSEA_p_GObp_H1_new.csv")
write.csv(gseGO.mf.df, file = "output/GSEA_p_GOmf_H1_new.csv")
write.csv(gseGO.cc.df, file = "output/GSEA_p_GOcc_H1_new.csv")
write.csv(gse.mkegg.df, file = "output/GSEA_p_mkegg_H1_new.csv")
write.csv(gse.wp.df, file = "output/GSEA_p_wiki_H1_new.csv")
write.csv(gse.hallmark.df, file = "output/GSEA_p_hallmark_H1_new.csv")

# end of GSEA for correlated genes ==================


# ORA for correlated genes ==============

geneList <- na.omit(rank.cor)$Rho
names(geneList) <- as.character(na.omit(rank.cor)$entrez)
geneList = sort(geneList, decreasing = TRUE)
head(geneList)

all_genes <- as.character(na.omit(rank.cor)$entrez)

x <- na.omit(rank.sig)
x <- na.omit(rank.pos)
x <- na.omit(rank.neg)

ora<- function(x){
  
  #KEGG ORA
  KEGG <- enrichKEGG(gene         = x$entrez,
                     organism     = 'hsa',
                     universe      = all_genes,
                     pvalueCutoff=1, pAdjustMethod="BH", 
                     qvalueCutoff=1)
  KEGG <<- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  KEGG.df <<- as.data.frame(KEGG)
  View(KEGG.df)
  #KEGG.df2 <- KEGG.df
  #KEGG.df2$Description <- gsub(x = KEGG.df2$Description, pattern = "\\ ",replacement = "_") 
  #KEGG.df2$Description <- gsub(x = KEGG.df2$Description, pattern = "\\,",replacement = ".")
  #write.xlsx2(KEGG.df, file="output/correlated_genes_ORA.xlsx", sheetName = "KEGG",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  #write.csv(KEGG.df2, sep="\t",file="output/correlated_genes_KEGG.csv", row.names=TRUE,col.names=NA,quote=FALSE)
  
  #Reactome ORA
  react <- enrichPathway(gene         = x$entrez,
                         organism     = 'human',
                         universe      = all_genes,
                         minGSSize = 10,
                         maxGSSize = 500,
                         pvalueCutoff=1, pAdjustMethod="BH", 
                         qvalueCutoff=1)
  react <<- setReadable(react, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  react.df <<- as.data.frame(react)
  View(react.df)
  #react.df2 <- react.df
  #react.df2$Description <- gsub(x = react.df2$Description, pattern = "\\ ",replacement = "_") 
  #react.df2$Description <- gsub(x = react.df2$Description, pattern = "\\,",replacement = ".")
  #write.xlsx2(react.df, file="output/correlated_genes_ORA.xlsx", sheetName = "Reactome",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  
  gobp <<- enrichGO(gene        = x$entrez,
                    universe      = all_genes,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    readable      = TRUE)
  gobp.df <<- as.data.frame(gobp)
  View(gobp.df)
  #write.xlsx2(gobp.df, file="output/correlated_genes_ORA.xlsx", sheetName = "GO_BP",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  
  gomf <<- enrichGO(gene       = x$entrez,
                    universe      = all_genes,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    readable      = TRUE)
  gomf.df <<- as.data.frame(gomf)
  head(gomf.df)
  #write.xlsx2(gomf.df, file="output/correlated_genes_ORA.xlsx", sheetName = "GO_MF",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  
  gocc <<- enrichGO(gene      = x$entrez,
                    universe      = all_genes,
                    OrgDb         = org.Hs.eg.db,
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    readable      = TRUE)
  gocc.df <<- as.data.frame(gocc)
  head(gocc.df)
  #write.xlsx2(gocc.df, file="output/correlated_genes_ORA.xlsx", sheetName = "GO_CC",
  #            col.names = TRUE, row.names = TRUE, append = TRUE)
  
  # MSigDb ORA
  hallmark <- enricher(x$entrez, 
                       universe = all_genes,
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       TERM2GENE=h_t2g)
  hallmark <- setReadable(hallmark, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  hallmark.df <- as.data.frame(hallmark)
  View(hallmark.df)
  
}

# save ORA all correlated genes
write.csv(KEGG.df, file = "output/ORA_cor_kegg_H1_new.csv")
write.csv(react.df, file = "output/ORA_cor_reactome_H1_new.csv")
write.csv(gobp.df, file = "output/ORA_cor_gobp_H1_new.csv")
write.csv(gomf.df, file = "output/ORA_cor_gomf_H1_new.csv")
write.csv(gocc.df, file = "output/ORA_cor_gocc_H1_new.csv")

#save ORA positively correlated genes
write.csv(KEGG.df, file = "output/ORA_cor_positive_kegg_H1_new.csv")
write.csv(react.df, file = "output/ORA_cor_positive_reactome_H1_new.csv")
write.csv(gobp.df, file = "output/ORA_cor_positive_gobp_H1_new.csv")
write.csv(gomf.df, file = "output/ORA_cor_positive_gomf_H1_new.csv")
write.csv(gocc.df, file = "output/ORA_cor_positive_gocc_H1_new.csv")
write.csv(hallmark.df, file = "output/ORA_cor_positive_hallmark_H1_new.csv")

#save ORA negatively correlated genes
write.csv(KEGG.df, file = "output/ORA_cor_negative_kegg_H1_new.csv")
write.csv(react.df, file = "output/ORA_cor_negative_reactome_H1_new.csv")
write.csv(gobp.df, file = "output/ORA_cor_negative_gobp_H1_new.csv")
write.csv(gomf.df, file = "output/ORA_cor_negative_gomf_H1_new.csv")
write.csv(gocc.df, file = "output/ORA_cor_negative_gocc_H1_new.csv")

write.csv(hallmark.df, file = "output/ORA_cor_negative_hallmark_H1_new.csv")

# end ORA =========



# old code, compared correlation of all genes vs only de genes===============
listInput1 <- list("Correlated genes (adj. p<0.05)"=rank.sig.M$symbol,
                   "Correlated DE genes (adj. p<0.05)"=rank.sig.M.de$symbol)
upset(fromList(listInput1),nsets=10, nintersects = NA,mb.ratio = c(0.6, 0.4), 
      text.scale = c(2,2, 2, 1.5, 2.5, 2), point.size = 6, line.size = 2,  order.by = c("degree"))
      #queries = list(list(query = intersects, params = list("Insulin_SMP", "Insulin_FUSION", "BS_DEGs"), color = "blue", active = T), 
      #               list(query = intersects,  params = list("Insulin_SMP", "Insulin_FUSION"), color = "purple", active = T), 
      #               list(query = intersects, params = list("Insulin_SMP", "BS_DEGs"),color = "purple", active = T),
      #               list(query = intersects, params = list("Insulin_FUSION", "BS_DEGs"), color = "purple", active = T)))
listInput2 <- list(#"Correlated genes (adj. p<0.05)"=rank.sig.M$symbol,
                   "Correlated genes (p<0.05)"=rank.sig.ua.M$symbol,
                   "DE genes "=res.merge.M$symbol)
upset(fromList(listInput2),nsets=10, nintersects = NA,mb.ratio = c(0.6, 0.4), 
      text.scale = c(2,2, 2, 1.5, 2.5, 2), point.size = 6, line.size = 2,  order.by = c("degree"))

listInput3 <- list("Correlated genes (adj. p<0.05)"=rank.sig.M$symbol,
  "Correlated genes (p<0.05)"=rank.sig.ua.M$symbol,
  #"DE genes "=res.merge.M$symbol,
  "DE MEL1 vs H1 "=res.merge.cell$symbol)
upset(fromList(listInput3),nsets=10, nintersects = NA,mb.ratio = c(0.6, 0.4), 
      text.scale = c(2,2, 2, 1.5, 2.5, 2), point.size = 6, line.size = 2,  order.by = c("degree"))

listInput4 <- list("Correlated genes MEL1 (adj. p<0.05)"=rank.sig.M$symbol,
                   "Correlated genes MEL1 (p<0.05)"=rank.sig.ua.M$symbol,
                   "Correlated genes  MEL1 & H1 (adj. p<0.05)"=rank.sig.no.cell$symbol,
                   "Correlated genes MEL1 & H1 (p<0.05)"=rank.sig.ua.no.cell$symbol
                   )
upset(fromList(listInput4),nsets=10, nintersects = NA,mb.ratio = c(0.6, 0.4), 
      text.scale = c(2,2, 2, 1.5, 2.5, 2), point.size = 6, line.size = 2,  order.by = c("degree"))

listInput5 <- list("Correlated genes MEL1 (adj. p<0.05)"=rank.sig.M$symbol,
                   "Correlated genes MEL1 (p<0.05)"=rank.sig.ua.M$symbol,
                   "Correlated genes H1 (adj. p<0.05)"=rank.sig.H$symbol,
                   "Correlated genes H1 (p<0.05)"=rank.sig.ua.H$symbol
)
upset(fromList(listInput5),nsets=10, nintersects = NA,mb.ratio = c(0.6, 0.4), 
      text.scale = c(2,2, 2, 1.5, 2.5, 2), point.size = 6, line.size = 2,  order.by = c("degree"))

#---------------------------------------------------------
library(VennDiagram)
## get Venn partition
venn1 <- get.venn.partitions(list(Co_MEL1_adj=rank.sig.M$symbol,
                   Co_MEL1_unadj=rank.sig.ua.M$symbol,
                   Co_H1_adj=rank.sig.H$symbol,
                   Co_H1_unadj=rank.sig.ua.H$symbol)) 

#### separate nested lists, tranfer wide to long format
library(splitstackshape)

df_venn1 <- listCol_w(venn1, "..values..", fill = NA)
View(df_venn1)
venn_long1 <- df_venn1 %>% 
  tidyr::pivot_longer(cols=..values.._fl_0001:..values.._fl_2404,
                      values_to = "symbol")%>% 
  na.omit()%>% dplyr::select(!name)
venn_long1$entrez <- mapIds(org.Hs.eg.db, keys=venn_long1$symbol, column="ENTREZID", keytype="ALIAS", multiVals="first")
#venn_long1$degree <- rowSums(venn_long1 == "TRUE")
#venn_long1$description <- mapIds(org.Hs.eg.db, keys=venn_long1$entrez, column="GENENAME", keytype="ENTREZID", multiVals="first")
#venn_long1$description <- gsub(x = venn_long1$description, pattern = "\\ ",replacement = "_")
#venn_long1$description <- gsub(x = venn_long1$description, pattern = "\\,",replacement = ".")
#venn_long1 <- venn_long1[order(-venn_long1$degree),]
View(venn_long1)
colnames(venn_long1)

co.gene1 <- venn_long1 %>% filter(Co_MEL1_adj=="TRUE" & Co_MEL1_unadj=="TRUE" &
                                   Co_H1_adj=="TRUE" & Co_H1_unadj=="TRUE")
co.gene2 <- venn_long1 %>% filter(Co_MEL1_adj=="TRUE" & Co_MEL1_unadj=="TRUE" &
                                    Co_H1_adj=="FALSE" & Co_H1_unadj=="TRUE")
co.gene3 <- venn_long1 %>% filter(Co_MEL1_adj=="FALSE" & Co_MEL1_unadj=="TRUE" &
                                    Co_H1_adj=="TRUE" & Co_H1_unadj=="TRUE")
co.gene4 <- venn_long1 %>% filter(Co_MEL1_adj=="FALSE" & Co_MEL1_unadj=="TRUE" &
                                    Co_H1_adj=="FALSE" & Co_H1_unadj=="TRUE")
#venn_long1 <- read.table("data/venn_long1.txt", row.names = 1)
co.gene <- rbind(co.gene1,co.gene2,co.gene3,co.gene4)
dim(co.gene)
dim(co.gene1)
dim(co.gene2)
dim(co.gene3)
dim(co.gene4)
View(co.gene)
colnames(rank.sig.ua.H)
head(rank.sig.ua.M[,c(1:4)])
co.gene <- co.gene %>% left_join(rank.sig.ua.M[,c(1:5)], by=c("symbol"="symbol")) %>%
  left_join(rank.sig.ua.H[,c(1:5)], by=c("symbol"="symbol")) %>% 
  rename(Rho.M=Rho.x, Rho.H=Rho.y)
co.gene$description = mapIds(org.Hs.eg.db, keys=co.gene$symbol, column="GENENAME", keytype="ALIAS", multiVals="first")
dup= co.gene %>% group_by(entrez) %>% dplyr::filter(n() > 1) # no duplicated genes
View(dup)

co.gene$direction <- co.gene$Rho.M*co.gene$Rho.H
co.gene.up <- co.gene %>% filter(!direction<0) %>% filter(Rho.M>0|Rho.H>0)
co.gene.down <- co.gene %>% filter(!direction<0) %>% filter(Rho.M<0|Rho.H<0)
dim(co.gene.up)
dim(co.gene.down)
write.xlsx2(co.gene.down, file="output/correlated_genes_list_MH_common.xlsx",col.names = TRUE, row.names = TRUE, sheetName = "common neg correlated", append = T)
write.xlsx2(co.gene.up, file="output/correlated_genes_list_MH_common.xlsx",col.names = TRUE, row.names = TRUE, sheetName = "common pos correlated",append = T)
write.xlsx2(co.gene, file="output/correlated_genes_list_MH_common.xlsx",col.names = TRUE, row.names = TRUE, sheetName = "all common correlated",append = T)
write.xlsx2(as.data.frame(dup), file="output/correlated_genes_list_MH_common.xlsx",col.names = TRUE, row.names = TRUE, sheetName = "duplicated or missing EntrezID",append = T)

#end ==============

#clustering_method = "ward.D2", #"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#clustering_distance_rows= #"euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"


# convert gene IDs for all norm genes (not needed) ======
norm.data <- as.data.frame(norm.data)
norm.data$entrez.s <- mapIds(org.Hs.eg.db, keys=row.names(norm.data), column="ENTREZID", keytype="SYMBOL", multiVals="first")
norm.data$entrez.a <- mapIds(org.Hs.eg.db, keys=row.names(norm.data), column="ENTREZID", keytype="ALIAS", multiVals="first")
norm.data<- norm.data %>% mutate(entrez=coalesce(entrez.s,entrez.a))

norm.data$ensembl.s <- mapIds(org.Hs.eg.db, keys=row.names(norm.data), column="ENSEMBL", keytype="SYMBOL", multiVals="first")
norm.data$ensembl.a <- mapIds(org.Hs.eg.db, keys=row.names(norm.data), column="ENSEMBL", keytype="ALIAS", multiVals="first")
norm.data<- norm.data %>% mutate(ensembl=coalesce(ensembl.s,ensembl.a)) 



length(unique(row.names(norm.data))) #10484
length(unique(norm.data$ensembl)) # 10345
length(unique(norm.data$entrez)) #10441

norm.data <- norm.data %>% dplyr::select(-entrez.a,-entrez.s,-ensembl.a,-ensembl.s)
View(norm.data)

#write.table(rank.cor, sep="\t",file="output/correlated_genes_H1_new.txt", row.names=TRUE,col.names=NA,quote=FALSE)

# end ==========

#====================================================
# Plot pathways 
react.df1 <-  react.df %>% mutate(Log10adj.P=log10(p.adjust)*(-1)) #rank.down
react.df1$GeneRatio_num <- sapply(react.df1$GeneRatio, function(x) eval(parse(text=x)))

react.df.select <- react.df1 %>% filter(qvalue<0.005) %>% 
  filter(ID %in% c("R-HSA-195253", "R-HSA-71291","R-HSA-9010553","R-HSA-9675108","R-HSA-983705",
                                                  "R-HSA-8957322","R-HSA-195721","R-HSA-422475","R-HSA-5339716","R-HSA-2454202",
                                                  "R-HSA-9604323","R-HSA-196299","R-HSA-5610780","R-HSA-5610783","R-HSA-5621481","R-HSA-5607764")) 
View(react.df)
View(react.df.select)

p <- ggplot(react.df.select, 
            aes(y = fct_reorder(Description, -pvalue))) + 
  geom_point(aes(size = Count,
                 x=GeneRatio_num,colour = Log10adj.P)) +
  theme_bw(base_size = 16) +
  scale_color_gradient2(#limits=c(3,60),
    midpoint = 0, mid="pink",
    low = "pink", high = "red4" )+
  ylab(NULL)+
  xlab("Gene Ratio")+
  theme(#axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black",size=14),
    axis.text.x = element_text(colour = "black",size=12),
    legend.title = element_text(color = "black", size = 13))+
  labs(color = "-log10 adj. p",size="Gene counts") +
  theme(panel.grid.minor.y = element_blank())+
  theme(panel.grid.minor.x  = element_blank())
#facet_grid(~group, scales='free',space = "free")
#theme(panel.grid.major = element_blank()) #10x7.5inch
# legend.text = element_blank())+
p

p <- ggplot(react.df.select, 
            aes(y = fct_reorder(Description, -pvalue))) + 
  geom_point(aes(size = Count,
                 x=GeneRatio_num,colour = Log10adj.P)) +
  theme_bw(base_size = 16) +
  scale_color_gradient2(limits=c(2,6.5), 
                        #midpoint = 2.1, 
                        low = "lightblue", high = "blue4", space = "Lab" )+
  ylab(NULL)+
  xlab("Gene Ratio")+
  theme(#axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black",size=14),
    axis.text.x = element_text(colour = "black",size=12),
    legend.title = element_text(color = "black", size = 13))+
  labs(color = "-log10 adj. p",size="Gene counts") +
  theme(panel.grid.minor.y = element_blank())+
  theme(panel.grid.minor.x  = element_blank())
p


# concept network plot==========================================
head(rank.down)
geneList.down = rank.down[,1]
## feature 2: named vector
names(geneList.down) = as.character(rank.down$symbol)
## feature 3: decreasing orde
geneList.down = sort(geneList.down, decreasing = TRUE)
head(geneList.down)
#
list(react.df.select$Description)
p <- cnetplot(react, showCategory=c("Signaling by the B Cell Receptor (BCR)", 
                                    "Degradation of GLI2 by the proteasome",                 
                                    "Fc epsilon receptor (FCERI) signaling",                 
                                    "Degradation of beta-catenin by the destruction complex",
                                    "C-type lectin receptors (CLRs)" ,                       
                                    "CLEC7A (Dectin-1) signaling" ,                          
                                    "Negative regulation of NOTCH4 signaling",               
                                    "Degradation of GLI1 by the proteasome" ,                
                                    "Signaling by WNT" ,                                     
                                    "Regulation of expression of SLITs and ROBOs",         
                                    "Metabolism of amino acids and derivatives" ,            
                                    "Axon guidance",                                         
                                    "Misspliced GSK3beta mutants stabilize beta-catenin",    
                                    "Nervous system development",                            
                                    "Beta-catenin phosphorylation cascade",                  
                                    "Metabolism of steroids"
                                    ),
              node_label="gene",
              categorySize="pvalue", colorEdge = F,
              order=TRUE, by="Count", foldChange=geneList.down,
              layout = "kk")+
  scale_color_gradient(#high='#deebf7', low='#2171b5'
    high="lightblue", low="blue3") +
  labs(edge="c",color = "log2 fold change",size="Gene count") #
p

description <- unlist(as.character(react_down.BS$Description))
View(description)
levels(react_down.BS$Description)
p <- cnetplot(react.d, showCategory=17,
              node_label="category",cex_label_category = 0.1,
              categorySize="pvalue", colorEdge = TRUE,
              order=F, foldChange=geneList,
              layout = "dh")+
  scale_color_gradient(high='#deebf7', low='#2171b5')+
  labs(color = "log2 fold change",size="gene counts")
p
#layout: kk,fr,graphopt,dh,mds,gem,lgl,drl...
##









#================================================================================================
# pathway analysis for DE genes

## create backgroud gene list
norm <- read.table(file="output/Norm_data_filtered.txt",row.names=1)
dim(norm)
head(norm)
View(norm)
#norm$entrez <- mapIds(org.Hs.eg.db, keys=rownames(norm), column="ENTREZID", keytype="SYMBOL", multiVals="first")
norm$entrez = mapIds(org.Hs.eg.db, keys=rownames(norm), column="ENTREZID", keytype="ALIAS", multiVals="first")
#norm1 <- norm %>%  mutate(entrez = coalesce(entrez,entrez2)) 
norm$sum <- rowSums(norm[,1:12])
norm <- norm[order(norm[,'entrez'],-norm[,'sum']),] # order by total counts high to low
dup1=norm %>% group_by(entrez) %>% dplyr::filter(n() > 1) 
head(dup1)
norm1 <- norm[!duplicated(norm$entrez),] # some essembl and symbol are the same gene, kept one with more total/basemean counts.
dup1=norm1 %>% group_by(entrez) %>% dplyr::filter(n() > 1) 
head(dup1)
norm2=norm1 %>% dplyr::filter(!entrez %in% c(NA)) 
head(norm2)
all_genes = unique(norm2$entrez) # all detected genes as the backgroud of ORA
head(all_genes)
str(all_genes)

## Over representation analyses (ORA)
res.merge <- read.table(file="output/Results_DE_merge.txt",row.names=1)
View(res.merge)
colnames(res.merge)

ora_gse<- function(x){
  res.sig <- res.merge %>% select(c("symbol",
                                    paste0("log2FoldChange_Protocol1.Protocol",x),
                                    paste0("padj_Protocol1.Protocol",x),
                                    "entrez")) %>% 
    na.omit() #select the comparison from the merged DE gene spreadsheet
  View(res.sig)
  
  #KEGG ORA
  KEGG <- enrichKEGG(gene         = res.sig$entrez,
                     organism     = 'hsa',
                     universe      = all_genes,
                     pvalueCutoff=0.05, pAdjustMethod="BH", 
                     qvalueCutoff=0.05)
  KEGG <- setReadable(KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  head(KEGG)
  KEGG.df <- as.data.frame(KEGG)
  View(KEGG.df)
  KEGG.df2 <- KEGG.df
  KEGG.df2$Description <- gsub(x = KEGG.df2$Description, pattern = "\\ ",replacement = "_") 
  KEGG.df2$Description <- gsub(x = KEGG.df2$Description, pattern = "\\,",replacement = ".")
  write.csv(KEGG.df, sep="\t",file=paste0("output/KEGG_protocol1_",x,".csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  write.csv(KEGG.df2, sep="\t",file=paste0("output/KEGG_protocol1_",x,"_read.csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  
  #Reactome ORA
  react <- enrichPathway(gene         = res.sig$entrez,
                         organism     = 'human',
                         universe      = all_genes,
                         minGSSize = 10,
                         maxGSSize = 500,
                         pvalueCutoff=0.05, pAdjustMethod="BH", 
                         qvalueCutoff=0.05)
  react <- setReadable(react, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  react.df=as.data.frame(react)
  View(react.df)
  react.df2 <- react.df
  react.df2$Description <- gsub(x = react.df2$Description, pattern = "\\ ",replacement = "_") 
  react.df2$Description <- gsub(x = react.df2$Description, pattern = "\\,",replacement = ".")
  write.csv(react.df, sep="\t",file=paste0("output/Reactome_protocol1_",x,".csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  write.csv(react.df, sep="\t",file=paste0("output/Reactome_protocol1_",x,"_read.csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  
  
  
  # GSEA gene list
  res <- read.table(file=paste0("output/Results_Protocol1_Protocol",x,".txt"),row.names=1)
  res$entrez = mapIds(org.Hs.eg.db, keys=rownames(res), column="ENTREZID", keytype="ALIAS", multiVals="first")
  res <- res[order(-res[,"baseMean"]),] # order by basemean (counts)
  res <- res[!duplicated(res$entrez),] %>% dplyr::filter(!entrez %in% c(NA))  # some essembl and symbol are the same gene, kept one with more total/basemean counts.
  dup1=res %>% group_by(entrez) %>% dplyr::filter(n() > 1) 
  head(dup1)
  geneList = res[,"log2FoldChange"]
  ## feature 2: named vector
  names(geneList) = as.character(res$entrez)
  ## feature 3: decreasing orde
  geneList= sort(geneList, decreasing = TRUE)
  head(geneList)
  
  # KEGG GSEA
  KEGG.gse <- gseKEGG(
    geneList,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    nPerm=1000000
  )
  KEGG.gse <- setReadable(KEGG.gse, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  KEGG.gse.df <- as.data.frame(KEGG.gse)
  View(KEGG.gse.df)
  KEGG.gse.df2 <- KEGG.gse.df
  KEGG.gse.df2$Description <- gsub(x = KEGG.gse.df2$Description, pattern = "\\ ",replacement = "_") 
  KEGG.gse.df2$Description <- gsub(x = KEGG.gse.df2$Description, pattern = "\\,",replacement = ".")
  write.csv(KEGG.gse.df, sep="\t",file=paste0("output/KEGG_GSEA_protocol1_",x,".csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  write.csv(KEGG.gse.df2, sep="\t",file=paste0("output/KEGG_GSEA_protocol1_",x,"_read.csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  
  #Reactome GSEA
  react.gse <- gsePathway(
    geneList,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    nPerm=1000000
  )
  react.gse <- setReadable(react.gse, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  react.gse.df <- as.data.frame(react.gse)
  View(react.gse.df)
  react.gse.df2 <- react.gse.df
  react.gse.df2$Description <- gsub(x = react.gse.df2$Description, pattern = "\\ ",replacement = "_") 
  react.gse.df2$Description <- gsub(x = react.gse.df2$Description, pattern = "\\,",replacement = ".")
  write.csv(react.gse.df, sep="\t",file=paste0("output/Reactome_GSEA_protocol1_",x,".csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  write.csv(react.gse.df2, sep="\t",file=paste0("output/Reactome_GSEA_protocol1_",x,"_read.csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  
  GO <- enrichGO(gene          = res.sig$entrez,
                 universe      = all_genes,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
  head(GO)
  GO.df <- as.data.frame(GO)
  View(GO.df)
  GO.df2 <- GO.df
  GO.df2$Description <- gsub(x = GO.df2$Description, pattern = "\\ ",replacement = "_") 
  GO.df2$Description <- gsub(x = GO.df2$Description, pattern = "\\,",replacement = ".")
  write.csv(GO.df, sep="\t",file=paste0("output/GO_protocol1_",x,".csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  write.csv(GO.df2, sep="\t",file=paste0("output/GO_protocol1_",x,"_read.csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  
}


for(i in 2:6) {
  ora_gse(i)
}

go<- function(x){
  res.sig <- res.merge %>% select(c("symbol",
                                    paste0("log2FoldChange_Protocol1.Protocol",x),
                                    paste0("padj_Protocol1.Protocol",x),
                                    "entrez")) %>% 
    na.omit() #select the comparison from the merged DE gene spreadsheet
  View(res.sig)
  
  # GSEA gene list
  res <- read.table(file=paste0("output/Results_Protocol1_Protocol",x,".txt"),row.names=1)
  res$entrez = mapIds(org.Hs.eg.db, keys=rownames(res), column="ENTREZID", keytype="ALIAS", multiVals="first")
  res <- res[order(-res[,"baseMean"]),] # order by basemean (counts)
  res <- res[!duplicated(res$entrez),] %>% dplyr::filter(!entrez %in% c(NA))  # some essembl and symbol are the same gene, kept one with more total/basemean counts.
  dup1=res %>% group_by(entrez) %>% dplyr::filter(n() > 1) 
  head(dup1)
  geneList = res[,"log2FoldChange"]
  ## feature 2: named vector
  names(geneList) = as.character(res$entrez)
  ## feature 3: decreasing orde
  geneList= sort(geneList, decreasing = TRUE)
  head(geneList)
  
  
  #Reactome GSEA
  react.gse <- gsePathway(
    geneList,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    nPerm=1000000
  )
  react.gse <- setReadable(react.gse, OrgDb = org.Hs.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
  react.gse.df <- as.data.frame(react.gse)
  View(react.gse.df)
  react.gse.df2 <- react.gse.df
  react.gse.df2$Description <- gsub(x = react.gse.df2$Description, pattern = "\\ ",replacement = "_") 
  react.gse.df2$Description <- gsub(x = react.gse.df2$Description, pattern = "\\,",replacement = ".")
  write.csv(react.gse.df, sep="\t",file=paste0("output/Reactome_GSEA_protocol1_",x,".csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  write.csv(react.gse.df2, sep="\t",file=paste0("output/Reactome_GSEA_protocol1_",x,"_read.csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  
  GO <- enrichGO(gene          = res.sig$entrez,
                 universe      = all_genes,
                 OrgDb         = org.Hs.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
  head(GO)
  GO.df <- as.data.frame(GO)
  View(GO.df)
  GO.df2 <- GO.df
  GO.df2$Description <- gsub(x = GO.df2$Description, pattern = "\\ ",replacement = "_") 
  GO.df2$Description <- gsub(x = GO.df2$Description, pattern = "\\,",replacement = ".")
  write.csv(GO.df, sep="\t",file=paste0("output/GO_protocol1_",x,".csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  write.csv(GO.df2, sep="\t",file=paste0("output/GO_protocol1_",x,"_read.csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  
  GO.gse <- enrichGO(gene     = geneList,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     readable      = TRUE,
                     minGSSize = 10,
                     maxGSSize = 500,
                     nPerm=1000000)
  head(GO)
  GO.df <- as.data.frame(GO)
  View(GO.df)
  GO.df2 <- GO.df
  GO.df2$Description <- gsub(x = GO.df2$Description, pattern = "\\ ",replacement = "_") 
  GO.df2$Description <- gsub(x = GO.df2$Description, pattern = "\\,",replacement = ".")
  write.csv(GO.df, sep="\t",file=paste0("output/GO_protocol1_",x,".csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  write.csv(GO.df2, sep="\t",file=paste0("output/GO_protocol1_",x,"_read.csv"), row.names=TRUE,col.names=NA,quote=FALSE)
  
}


for(i in 2:6) {
  ora_gse(i)
}

GO

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,)

GO

ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
KEGG <- enrichKEGG(gene         = Sig.up.BS$entrez,
                   organism     = 'mmu',
                   universe      = all_genes,
                   pvalueCutoff=0.005, pAdjustMethod="BH", 
                   qvalueCutoff=0.05)
KEGG_up.BS <- setReadable(KEGG_up.BS, OrgDb = org.Mm.eg.db, keyType="ENTREZID") # The geneID column is translated from EntrezID to symbol
head(KEGG_up.BS)
KEGG_up.BS.df=as.data.frame(KEGG_up.BS)
View(KEGG_up.BS.df)
write.csv(KEGG_up.BS.df, sep="\t",file="data/KEGG_up_BS200-0_noNA_BG_filter.csv", row.names=TRUE,col.names=NA,quote=FALSE)
#KEGG_up.BS.df$Description <- gsub(x = KEGG_up.BS.df$Description, pattern = "\\ ",replacement = "_") 
#KEGG_up.BS.df$Description <- gsub(x = KEGG_up.BS.df$Description, pattern = "\\,",replacement = ".") 
#write.table(KEGG_up.BS.df, sep="\t",file="data/KEGG_up_BS200-0_noNA_BG_filter.txt", row.names=TRUE,col.names=NA,quote=FALSE)

#-------------------------------------------------------------
# Correlation with S4 efficiency
#S4 <- read.xlsx2(file="input/S4_efficiency.xlsx", sheetIndex = 1)
#row.names(S4) <- S4$S1.Protocol
#S4 <- t(S4[,-1])
#View(S4)
# Shapiro-Wilk normality test
#shapiro.test(S4)
#library(ggpubr)
#ggqqplot(S4, ylab = "S4 efficiency")
#head(norm)
#correlation_function.p<- function(df){
#  r <- (1:ncol(df))
#  df[ , r] <- apply(df[ , r], 2,    # Specify own function within apply, character to numeric
#                    function(x) as.numeric(as.character(x)))
#  df.t = as.data.frame(t(df)) #convert row to column
#  S4 <- as.numeric(S4)
#  n <- ncol(df.t)
#  correlations <-lapply(1:n, function(x) (cor.test(S4,df.t[,x],method='pearson')))
#r value
#  R <- sapply(correlations,function(x) x[["estimate"]])
#p value
#  p<- sapply(correlations,function(x) x[["p.value"]])
#  corr<-data.frame(R=R,p=p,R2=R^2)
#  rownames(corr) <- colnames(df.t)
#  corr$adj.p <- p.adjust(corr$p, method = 'BH', n = length(corr$p)) # p.adjust.methods c("holm", "hochberg", "hommel", "bonferroni", "BH" (same as "fdr"), "BY", none")
#  corr$symbol <- row.names(corr)
#  corr$entrez = mapIds(org.Hs.eg.db, keys=rownames(corr), column="ENTREZID", keytype="ALIAS", multiVals="first")
#  corr$description = mapIds(org.Hs.eg.db, keys=rownames(corr), column="GENENAME", keytype="ALIAS", multiVals="first")
#  corr <<- corr[order(corr$R2,decreasing=T),]
#  head(corr)
#write.xlsx2(corr, file="output/correlated_gene_test.xlsx", sheetName = paste(deparse(substitute(df))),col.names = TRUE, row.names = TRUE, append = TRUE)
#write.xlsx2(corr, file=paste0("output/correlated_gene_pearson.xlsx"),col.names = TRUE, row.names = TRUE, append = TRUE)

#  corr.sig <<- corr %>% filter(adj.p<0.05)
#  corr.sig.ua <<- corr %>% filter(p<0.05)
#  df.order <<- df[,c(5,6,7,8,11,12,9,10,3,4,1,2)]
#  df.cor <- corr.sig %>% merge(df.order,by.x='symbol',by.y='row.names',all.x=TRUE)
#  df.cor <<- df.cor[order(df.cor$R,decreasing = T),]

#  df.cor.ua <- corr.sig.ua %>% merge(df.order,by.x='symbol',by.y='row.names',all.x=TRUE)
#  df.cor.ua <<- df.cor.ua[order(df.cor.ua$R,decreasing = T),]

#write.xlsx2(df.cor, file="output/correlated_gene_heatmap.xlsx",col.names = TRUE, row.names = TRUE, append = TRUE)
}

#correlation_function.p(norm)
#S4.sig <- corr.sig
#S4.sig.ua <- corr.sig.ua
#---------------------------------------------------------------------------------
# compared rank vs S4 efficiency correlation
#listInput.c <- list("Rank correlated gene (adj. p<0.05)"=rank.sig$symbol,
#                   "S4 efficiency correlated gene (adj. p<0.05)"=S4.sig$symbol)
#upset(fromList(listInput.c),nsets=10, nintersects = NA,mb.ratio = c(0.6, 0.4), 
#      text.scale = c(2,2, 2, 1.5, 2.5, 2), point.size = 6, line.size = 2,  order.by = c("degree"))

#listInput.c.ua <- list("Rank correlated gene (p<0.05)"=rank.sig.ua$symbol,
#                    "S4 efficiency correlated gene (p<0.05)"=S4.sig.ua$symbol)
#upset(fromList(listInput.c.ua),nsets=10, nintersects = NA,mb.ratio = c(0.6, 0.4), 
#      text.scale = c(2,2, 2, 1.5, 2.5, 2), point.size = 6, line.size = 2,  order.by = c("degree"))
#-----------------------------------------------------------------------------------
