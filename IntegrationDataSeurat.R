library(dplyr)
library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(tidyverse)
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

samples.list <- list.dirs(getwd(),recursive = FALSE)
samples <- basename(samples.list)
remove(samples.list)

if (!file.exists("DataQuality")) {dir.create("DataQuality")}

###Input Data

for (i in samples) {
  assign("rawdata",Read10X(i))
  assign(i,CreateSeuratObject(rawdata,project = i,min.cells = 3, min.features = 500))
  remove("rawdata")
}

datasets <- lapply(samples, get)
names(datasets) <- samples

remove(list=samples)

for (i in seq_along(datasets)){
  datasets[[i]][["percent.mt"]] <- PercentageFeatureSet(datasets[[i]], pattern = "^MT-")
}

### Data quality check

for (i in seq_along(datasets)){
  
  pdf(file = paste("DataQuality/",samples[i],"_Vlnplot_beforeQC.pdf",sep = ""),title="Vlnplot_QC", width=20)
  VlnPlot(datasets[[i]], features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
  dev.off()
  
  pdf(file = paste("DataQuality/",samples[i],"_FeatureScatter_beforeQC.pdf",sep = ""),title="FeatureScatter_QC", width=20)
  plot1 <- FeatureScatter(datasets[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(datasets[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
  dev.off()
}
remove(list = c("plot1","plot2"))

##Quality Control 

f <- function(x) median(x) + 3*mad(x)

for (i in seq_along(datasets)){
  stats <- datasets[[i]]@meta.data[, c("nCount_RNA", "nFeature_RNA","percent.mt")]
  bounds <- apply(stats, 2,f)
  datasets[[i]] <- subset(datasets[[i]], subset = nFeature_RNA >200 &  
                            nCount_RNA < bounds[1] &
                            nFeature_RNA < bounds[2] & 
                            percent.mt < bounds[3])
  remove(list = c("stats","bounds"))
  }



### Data quality check after that cells are filtered out through QC metric

for (i in seq_along(datasets)){

  pdf(file = paste("DataQuality/",samples[i],"_Vlnplot_afterQC.pdf",sep = ""),title="Vlnplot_QC", width=20)
VlnPlot(datasets[[i]], features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
dev.off()

pdf(file = paste("DataQuality/",samples[i],"_FeatureScatter_afterQC.pdf",sep = ""),title="FeatureScatter_QC", width=20)
  plot1 <- FeatureScatter(datasets[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(datasets[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
dev.off()
}
remove(list = c("plot1","plot2"))

###Normalization

datasets <-lapply(datasets,\(i) SCTransform(i, vst.flavor = "v2", verbose = FALSE))


###select repeatedly variable features 

features <- SelectIntegrationFeatures(object.list = datasets,nfeatures = 5000)
datasets <- lapply(X =  datasets, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

### Integration

Anchors <- FindIntegrationAnchors(object.list = datasets, anchor.features = features, reduction = "rpca")
Combined <- IntegrateData(anchorset = Anchors)

#save (Run HPC)

save("DataSeurat.RData",list = c("datasets","Combined"))


