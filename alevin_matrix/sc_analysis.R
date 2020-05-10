#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("biomaRt")
#install.packages("tidyverse")
#install.packages("ggplot2")
library(ggplot2)
library(biomaRt)
library(dplyr)
#install.packages('Seurat')
library(Seurat)

###Processing UMI Count Matrix and filterling low quality cells and low variance genes####

#Load count matrix
data <- read.csv('UMI_counts_matrix.csv', header = TRUE)
colnames(data)[1] <- "gene"

#subsetting gene column and removing '.' and any character after '.' from every ensembl gene ID
genes <- data[, 1]
ng_data <- data[, -1] #all columns except the gene column 
genes <- as.character(genes)
genes <- gsub("\\..*","",genes)

#checking for duplicate genes
sum(duplicated(genes)) 

#mapping ensembl gene ID to gene symbols 
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
test <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
              filter = "ensembl_gene_id", 
              values = genes,
              mart = mart)
matched <- match(genes, test$ensembl_gene_id)
genesyms <-  test$hgnc_symbol[matched]

data_1 <- cbind(genesyms, ng_data)
#data$gene <- genesyms

#Summing values of duplicated genes and removing genes with NA values 
data_1 <- na.omit(data_1)
data_1 <- data_1 %>% 
  group_by(genesyms) %>% 
  summarise_all(funs(sum))
#write.csv(data_1, file = "updated_count.csv")

#remove null values and empty rows with no mapping from dataframe 
rownames(data_1) <- data_1$genesyms
complete <- data_1[!(is.na(data_1$genesyms) | data_1$genesyms=="") ,]
rownames(complete) <- complete$genesyms

#Seurat object with non-normalized data 
pc_data <- CreateSeuratObject(counts = complete , project = "p5", min.cells = 3, min.features = 200)
pc_data
pc_data[["percent.mt"]] <- PercentageFeatureSet(pc_data, pattern = "^MT-")
head(pc_data@meta.data)

#Visualize QC metrics as a violin plot
VlnPlot(pc_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Visualize feature-feature relationships using FeatureScatter
plot1 <-FeatureScatter(pc_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pc_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#Filter cells that have unique feature counts over 2,500 or less than 200
#and cells with less than 10% mitochondrial genes
pc_data <- subset(pc_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

#Log normalize filtered expression measurements 
pc_data <- NormalizeData(pc_data, normalization.method = "LogNormalize", scale.factor = 10000)
head(pc_data@meta.data)

#Calculate a subset of features that exhibit high cell-to-cell variation in the dataset 
pc_data <- FindVariableFeatures(pc_data, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pc_data), 10)

# plot variable features with and without labels
plot3 <- VariableFeaturePlot(pc_data)
LabelPoints(plot = plot3, points = top10, repel = TRUE)

#Scaling data for PC analysis 
all.genes <- rownames(pc_data)
pc_data2 <- ScaleData(pc_data, features = all.genes)

#PCA
pc_data2 <- RunPCA(pc_data2, features = VariableFeatures(object = pc_data2))

print(pc_data2[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pc_data2, dims = 1:2, reduction = "pca")
DimPlot(pc_data2, reduction = "pca")

# Different methods to decide which PCs to include for further downstream analyses
DimHeatmap(pc_data2, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pc_data2, dims = 1:15, cells = 500, balanced = TRUE)
pc_data2 <- JackStraw(pc_data2, num.replicate = 100)
pc_data2 <- ScoreJackStraw(pc_data2, dims = 1:20)
JackStrawPlot(pc_data2, dims = 1:15)
ElbowPlot(pc_data2) #Elbow plot 

#Clustering 
pc_data2 <- FindNeighbors(pc_data2, dims = 1:10)
pc_data2 <- FindClusters(pc_data2, resolution = 0.7)
#Saving seruat object
#saveRDS(pc_data2, file = "seurat.rds")
####end####


####Identifying and clustering marker genes####

#RUN non linear dimensional reduction (TSNE)
pc_data2 <- RunTSNE(pc_data2, dims = 1:10, method= "FIt-SNE")
DimPlot(pc_data2, reduction= "tsne")


#Construction of barplot to show the number of cells in each cluster (#total #cells = 1881)
clusters <- Idents(pc_data2)
clusters_df <- as.data.frame(table(clusters))
clusters_df$percentage <- clusters_df$Freq/sum(clusters_df$Freq)*100
p<-ggplot(data=clusters_df, aes(x= clusters, y= Freq)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=Freq), vjust=1.6, color="white", size=3.5)+
  theme_minimal()+ 
  xlab("Clusters") + ylab("Number of Cells")
p

#Getting all marker genes for the cell clusters
marker_genes <- FindAllMarkers(pc_data2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#top_markers <-marker_genes %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#write.csv(marker_genes,  file ="markers.csv")

#Cluster ids retrieved from panglaodb for the top gene markers for the 11 clusters 
#Assigning cell type identity to clusters 
clusterids <- c("Beta cells", "Alpha cells", "Alpha cells", "Alpha cells", "Basal cells", " Pancreatic stellate cells", "Acinar cells", "Gamma cells","Ductal cells","Dendritic cells", "Endothelial cells")
names(clusterids) <- levels(pc_data2)
pc_data2 <- RenameIdents(pc_data2, clusterids)

#marker genes have logfc > 0.8 and adj pvalue < 0.05
top3 <- marker_genes %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)
DoHeatmap(pc_data2, features = top3$gene, angle = 90)
#tSNE for cell clusters
DimPlot(pc_data2, reduction= "tsne", label = TRUE, pt.size = 0.7)
####end####

####In-depth marker gene analysis####

#extracting each gene from each cluster 
cluster0 <- as.character(row.names(marker_genes)[marker_genes$cluster == 0])
cluster1 <- as.character(row.names(marker_genes)[marker_genes$cluster == 1])
cluster2 <- as.character(row.names(marker_genes)[marker_genes$cluster == 2])
cluster3 <- as.character(row.names(marker_genes)[marker_genes$cluster == 3])
cluster4 <- as.character(row.names(marker_genes)[marker_genes$cluster == 4])
cluster5 <- as.character(row.names(marker_genes)[marker_genes$cluster == 5])
cluster6 <- as.character(row.names(marker_genes)[marker_genes$cluster == 6])
cluster7 <- as.character(row.names(marker_genes)[marker_genes$cluster == 7])
cluster8 <- as.character(row.names(marker_genes)[marker_genes$cluster == 8])
cluster9 <- as.character(row.names(marker_genes)[marker_genes$cluster == 9])
cluster10 <-as.character(row.names(marker_genes)[marker_genes$cluster == 10])

#gene enrichment analysis for each gene in each cluster using enrichr  
#install.packages("enrichR")
library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("KEGG_2019_Human", "GO_Molecular_Function_2018")

enriched0 <- enrichr(c(cluster0), dbs)
enriched1 <- enrichr(c(cluster1), dbs)
enriched2 <- enrichr(c(cluster2), dbs)
enriched3 <- enrichr(c(cluster3), dbs)
enriched4 <- enrichr(c(cluster4), dbs)
enriched5 <- enrichr(c(cluster5), dbs)
enriched6 <- enrichr(c(cluster6), dbs)
enriched7 <- enrichr(c(cluster7), dbs)
enriched8 <- enrichr(c(cluster8), dbs)
enriched9 <- enrichr(c(cluster9), dbs)
enriched10 <- enrichr(c(cluster10), dbs)


kegg0 <- enriched0[["KEGG_2019_Human"]]
molecular0 <- enriched0[["GO_Molecular_Function_2018"]]

kegg1 <- enriched1[["KEGG_2019_Human"]]
molecular1 <- enriched1[["GO_Molecular_Function_2018"]]

kegg2 <- enriched2[["KEGG_2019_Human"]]
molecular2 <- enriched2[["GO_Molecular_Function_2018"]]

kegg3 <- enriched3[["KEGG_2019_Human"]]
molecular3 <- enriched3[["GO_Molecular_Function_2018"]]

kegg4 <- enriched4[["KEGG_2019_Human"]]
molecular4 <- enriched4[["GO_Molecular_Function_2018"]]

kegg5 <- enriched5[["KEGG_2019_Human"]]
molecular5 <- enriched5[["GO_Molecular_Function_2018"]]

kegg6 <- enriched6[["KEGG_2019_Human"]]
molecular6 <- enriched6[["GO_Molecular_Function_2018"]]

kegg7 <- enriched7[["KEGG_2019_Human"]]
molecular7 <- enriched7[["GO_Molecular_Function_2018"]]

kegg8 <- enriched8[["KEGG_2019_Human"]]
molecular8 <- enriched8[["GO_Molecular_Function_2018"]]

kegg9 <- enriched9[["KEGG_2019_Human"]]
molecular9 <- enriched9[["GO_Molecular_Function_2018"]]

kegg10 <- enriched10[["KEGG_2019_Human"]]
molecular10 <- enriched10[["GO_Molecular_Function_2018"]]
####end####


all <- rbind(kegg0[1:2,],kegg1[1,], kegg2[1:2,], kegg3[1,], kegg4[1:4,], kegg5[1:70,], kegg6[1:3,], kegg10[1:2,])
#write.csv(all,"pathways.csv")
