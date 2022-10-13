library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(harmony)
library(EnhancedVolcano)
set.seed(1234)

# Explore metadata
View(sce.big@meta.data)

# Add number of genes per UMI for each cell to metadata
sce.big$log10GenesPerUMI <- log10(sce.big$nFeature_RNA) / log10(sce.big$nCount_RNA)

# Compute percent mito ratio
sce.big$mitoRatio <- PercentageFeatureSet(object = sce.big, pattern = "^MT-")
sce.big$mitoRatio <- sce.big@meta.data$mitoRatio / 100

# Create metadata dataframe
metadata_rna <- sce.big@meta.data

# Add cell IDs to metadata
metadata_rna$cells <- rownames(metadata_rna)

# Rename columns
metadata_rna <- metadata_rna %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Create sample column
metadata_rna$sample <- NA
metadata_rna$sample[which(str_detect(metadata_rna$cells, "^DMSO_"))] <- "DMSO"
metadata_rna$sample[which(str_detect(metadata_rna$cells, "^ENZ48_"))] <- "ENZ48"
metadata_rna$sample[which(str_detect(metadata_rna$cells, "^RESA_"))] <- "RESA"
metadata_rna$sample[which(str_detect(metadata_rna$cells, "^RESB_"))] <- "RESB"

# Add metadata back to Seurat object
sce.big@meta.data <- metadata_rna

# Quality control
filtered_DMSO <- subset(x = sce.big, 
                        subset = (sample == "DMSO") &
                          (nUMI > 16000) & (nUMI < 50000) &
                          (nGene > 3000 ) & (nGene < 7000 ) &
                          (mitoRatio < 0.15))

filtered_ENZ48 <- subset(x = sce.big, 
                         subset= (sample == "ENZ48") &
                           (nUMI > 5000) & (nUMI < 25000) &
                           (nGene > 1500 ) & (nGene < 5000 ) &
                           (mitoRatio < 0.15))

filtered_RESA <- subset(x = sce.big, 
                        subset= (sample == "RESA") &
                          (nUMI > 5000) & (nUMI < 25000) &
                          (nGene > 1500 ) & (nGene < 5000 ) &
                          (mitoRatio < 0.17))

filtered_RESB <- subset(x = sce.big, 
                        subset= (sample == "RESB") &
                          (nUMI > 5000) & (nUMI < 25000) &
                          (nGene > 1500 ) & (nGene < 5000 ) &
                          (mitoRatio < 0.20))

# Merge all samples
filtered_rna <- merge(filtered_DMSO, 
                      y = c(filtered_ENZ48,filtered_RESA,filtered_RESB),  
                      project = "scRNA")

save(filtered_rna, file="/home/jli39/Code/filtered_rna.RData")

# Normalization and dimensionality reduction
filtered_rna <- NormalizeData(filtered_rna) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)

# Integration with harmony
rna.integrated <- RunHarmony(
  object = filtered_rna,
  group.by.vars = 'seq_folder'
)

# Produce clusters
rna.integrated <- RunUMAP(rna.integrated, reduction = "harmony", dims = 1:20)
rna.integrated <- FindNeighbors(rna.integrated, reduction = "harmony", dims = 1:20) 
rna.integrated <- FindClusters(object = rna.integrated,resolution = 0.9)
DimPlot(rna.integrated, split.by = "seq_folder", ncol = 4)

## one plot with shapes
library(cowplot);
rna_one<-DimPlot(rna.integrated, shape.by = "seq_folder", pt.size = .3) +  coord_fixed() +
  theme(legend.position="bottom");
zoom.in.purple<-DimPlot(rna.integrated, shape.by = "seq_folder", pt.size = 1.2) +
  xlim(c(-7, -5)) + ylim(c(-3, 0)) + theme(legend.position = "none") +
  coord_fixed() 

plot_grid(rna_one, '', '', zoom.in.purple, '', '', ncol=3); ## re-scale the panel for later assembly
dev.copy2pdf(file='merged rna-seq umap.pdf', family='ArialMT', height=8, width=12);

# Barplot
ptr <- table(Idents(rna.integrated), rna.integrated$seq_folder)
ptr <- as.data.frame(ptr)
ptr$Var1 <- as.character(ptr$Var1)
ptr$pro <- ifelse(ptr$Var2 == "DMSO", ptr$Freq/1782, 
                 ifelse(ptr$Var2 == "ENZ48", ptr$Freq/4315, 
                        ifelse(ptr$Var2 == "RESA", ptr$Freq/4569, ptr$Freq/4061)))

xx<-ggplot(ptr, aes(x = factor(Var1, levels = c(as.character(seq(0, 12)))), y = pro, fill = Var2)) +
  xlab("scRNA-seq cluster") +
  ylab("Proportion of cells in cluster") +
  geom_bar(position="dodge", stat="identity")
xx<-as_grob(xx);
plot_grid(xx, '', '', '', '', '', ncol=3);

dev.copy2pdf(file='RNA_barplot.pdf', family='ArialMT', paper='a4')

# Regroup clusters
rna.integrated$group <- ifelse(rna.integrated$seurat_clusters == 1 | rna.integrated$seurat_clusters == 4 | rna.integrated$seurat_clusters == 9, "initial cluster",
                              ifelse(rna.integrated$seurat_clusters == 0 | rna.integrated$seurat_clusters == 5 | rna.integrated$seurat_clusters == 8 | rna.integrated$seurat_clusters == 12, "ENZ-induced cluster",
                                     ifelse(rna.integrated$seurat_clusters == 10 | rna.integrated$seurat_clusters == 11 | rna.integrated$seurat_clusters == 3 | rna.integrated$seurat_clusters == 6 | rna.integrated$seurat_clusters == 7, "Persistent cluster","NNA")))

table(rna.integrated$group)
save(rna.integrated, file="/Users/jli39/Desktop/rna.integrated.RData")

# Find differentially expressed genes 
### initial cluster 1, 4, 9
### induced cluster 0, 5, 8, 12
rna.integrated$group2<-'None';
rna.integrated$group2[which(rna.integrated$seurat_clusters%in%c(1, 4, 9))]<-'Initial';
rna.integrated$group2[which(rna.integrated$seurat_clusters%in%c(3, 6, 7, 10, 11))]<-'Persistent';
rna.integrated$group2[which(rna.integrated$seurat_clusters%in%c(0, 5, 8, 12))]<-'Induced';

initial.cells<-colnames(rna.integrated)[which(rna.integrated$group2=='Initial')];
induced.cells<-colnames(rna.integrated)[which(rna.integrated$group2=='Induced')];
persistent.cells<-colnames(rna.integrated)[which(rna.integrated$group2=='Persistent')];

## filter fc first then subset then feed to findmarkers for significance calculation
#-------------------induced-initial-------------------
relog2fc.all.rna.initial_induced <-mclapply(rownames(rna.integrated), function(x){
  t1<-rna.integrated@assays$RNA@data[x, initial.cells];
  t2<-rna.integrated@assays$RNA@data[x, induced.cells];
  tt1<-sum(t1, na.rm = T);
  tt2<-sum(t2, na.rm = T);
  log2(tt2/tt1);
}, mc.cores = 20);

relog2fc.all.rna.t.initial_induced <- t(as.data.frame(relog2fc.all.rna.initial_induced))
row.names(relog2fc.all.rna.t.initial_induced) <- rownames(rna.integrated)

rna.markers1 <- FindMarkers(rna.integrated, ident.1 = c("0","5","8","12"), ident.2 = c("1","4","9"),
                            min.pct = 0, logfc.threshold=0);

deg.initial_induced.f <- merge(rna.markers1, relog2fc.all.rna.t.initial_induced, by="row.names", all.x=T)
rownames(deg.initial_induced.f) <- deg.initial_induced.f$Row.names
deg.initial_induced.f <- deg.initial_induced.f[, c(2:7) ]
colnames(deg.initial_induced.f)[6] <- "relog2fc"
deg.initial_induced.p <- dplyr::filter(deg.initial_induced.f, relog2fc != "Inf" & relog2fc != "-Inf")

rr1<-EnhancedVolcano(deg.initial_induced.p, lab = rownames(deg.initial_induced.p), 
                    x='relog2fc', y='p_val_adj');
plot_grid(rr1, '', '', '', ncol=2);
dev.copy2pdf(file='RNA-df.initial_induced.pdf', family='ArialMT', width = 8, height = 12);

ss.deg.initial_induced <- dplyr::filter(deg.initial_induced.p, p_val_adj < 0.05 & abs(relog2fc)>=log2(1.2))
save(ss.deg.initial_induced, file='ss.deg.initial_induced_0128.rda');

#-------------------persistent-induced-------------------
rna.markers2 <- FindMarkers(rna.integrated, ident.1 = c("0","5","8","12"), ident.2 = c("3","6","7","10","11"),
                            min.pct = 0, logfc.threshold=0);

deg.persistent_induced.f <- merge(rna.markers2, relog2fc.all.rna.t, by="row.names", all.x=T)
rownames(deg.persistent_induced.f) <- deg.persistent_induced.f$Row.names
deg.persistent_induced.f <- deg.persistent_induced.f[, c(2:7) ]
colnames(deg.persistent_induced.f)[6] <- "relog2fc"
deg.persistent_induced.p <- dplyr::filter(deg.persistent_induced.f, relog2fc != "Inf" & relog2fc != "-Inf")

rr2<-EnhancedVolcano(deg.persistent_induced.p, lab = rownames(deg.persistent_induced.p), 
                    x='relog2fc', y='p_val_adj', pCutoff=1e-6);
plot_grid(rr2, '', '', '', ncol=2);
dev.copy2pdf(file='RNA-df.persistent_induced.pdf', family='ArialMT', width = 8, height = 12);

ss.deg.persistent_induced <- dplyr::filter(deg.persistent_induced.p, p_val_adj < 1e-6)
save(ss.deg.persistent_induced, file='ss.deg.persistent_induced_20211202.rda');

#------------------persistent-initial-----------------
relog2fc.all.rna.initial_persistent <-mclapply(rownames(rna.integrated), function(x){
  t1<-rna.integrated@assays$RNA@data[x, initial.cells];
  t2<-rna.integrated@assays$RNA@data[x, persistent.cells];
  tt1<-sum(t1, na.rm = T);
  tt2<-sum(t2, na.rm = T);
  log2(tt2/tt1);
}, mc.cores = 8);

relog2fc.all.rna.t.initial_persistent <- t(as.data.frame(relog2fc.all.rna.initial_persistent))
row.names(relog2fc.all.rna.t.initial_persistent) <- rownames(rna.integrated)

rna.markers3 <- FindMarkers(rna.integrated, ident.1 = c("3","6","7","10","11"), ident.2 = c("1","4","9"),
                            min.pct = 0, logfc.threshold=0);

deg.persistent_initial.f <- merge(rna.markers3, relog2fc.all.rna.t.initial_persistent, by="row.names", all.x=T)
rownames(deg.persistent_initial.f) <- deg.persistent_initial.f$Row.names
deg.persistent_initial.f <- deg.persistent_initial.f[, c(2:7) ]
colnames(deg.persistent_initial.f)[6] <- "relog2fc"
deg.persistent_initial.p <- dplyr::filter(deg.persistent_initial.f, relog2fc != "Inf" & relog2fc != "-Inf")

rr3<-EnhancedVolcano(deg.persistent_initial.p, lab = rownames(deg.persistent_initial.p), 
                    x='relog2fc', y='p_val_adj', pCutoff=1e-6);
plot_grid(rr3, '', '', '', ncol=2);
dev.copy2pdf(file='RNA-df.persistent-initial.pdf', family='ArialMT', width = 8, height = 12);

ss.deg.persistent_initial <- dplyr::filter(deg.persistent_initial.p, p_val_adj < 1e-6 & abs(relog2fc)>=log2(1.2))
save(ss.deg.persistent_initial, file='ss.deg.persistent_initial_0128.rda');

# Add regulate direction
ss.deg.initial_induced$direction <- ifelse(ss.deg.initial_induced$avg_log2FC > 0, "up", "down")
ss.deg.persistent_induced$direction <- ifelse(ss.deg.persistent_induced$avg_log2FC > 0, "up", "down")
ss.deg.persistent_initial$direction <- ifelse(ss.deg.persistent_initial$avg_log2FC > 0, "up", "down")

write.csv(ss.deg.initial_induced,'ss.deg.initial_induced.csv')
write.csv(ss.deg.persistent_induced,'ss.deg.persistent_induced.csv')
write.csv(ss.deg.persistent_initial,'ss.deg.persistent_induced.csv')



