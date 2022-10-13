## package version issue
## install GenomeInfoDb from github for the newest version
### ### ### ### ### ### ### 
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager");
BiocManager::install(version = "3.14");

library(devtools);
install_git('https://github.com/Bioconductor/BiocManager.git', force = TRUE);
install_git('https://github.com/Bioconductor/GenomeInfoDb.git', force = TRUE);

### ### ### ### ### ### ### 
library(Signac);
library(Seurat);
library(ggplot2);
library(patchwork);
library(EnsDb.Hsapiens.v86);
library(dplyr);
library(rtracklayer);
library(hdf5r)
set.seed(1234);

setwd('G:/.shortcut-targets-by-id/1guiMrPeJVqjpkEwof0EhFkw3RxoiI6cZ/Jinze_scATAC/data/scATAC-seq');

# gene annotation
stardb<-import('genes.gtf.gz');
stardb$gene_biotype<-stardb$gene_type;
seqlevelsStyle(stardb) <- 'UCSC';
genome(stardb) <- "hg38";

# peaks in each dataset
counts_DMSO <- Read10X_h5(filename = "GSM5155451/filtered_peak_bc_matrix.h5");
metadata_DMSO <- read.csv(
  file = "GSM5155451/singlecell.csv",
  header = TRUE,
  row.names = 1
)

counts_ENZ48 <- Read10X_h5(filename = "GSM5155452/filtered_peak_bc_matrix.h5");
metadata_ENZ48 <- read.csv(
  file = "GSM5155452/singlecell.csv",
  header = TRUE,
  row.names = 1
)

counts_RESA <- Read10X_h5(filename = "GSM5155453/filtered_peak_bc_matrix.h5");
metadata_RESA <- read.csv(
  file = "GSM5155453/singlecell.csv",
  header = TRUE,
  row.names = 1
)

counts_RESB <- Read10X_h5(filename = "GSM5155454/filtered_peak_bc_matrix.h5");
metadata_RESB <- read.csv(
  file = "GSM5155454/singlecell.csv",
  header = TRUE,
  row.names = 1
)

# Create Objects
chrom_assay_DMSO <- CreateChromatinAssay(
  counts = counts_DMSO,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = 'GSM5155451/fragments.tsv.gz',
  annotation=stardb
)
pbmc_DMSO <- CreateSeuratObject(
  counts = chrom_assay_DMSO,
  assay = "peaks",
  project = 'ATAC',
  meta.data = metadata_DMSO
)

chrom_assay_ENZ48 <- CreateChromatinAssay(
  counts = counts_ENZ48,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = 'GSM5155452/fragments.tsv.gz',
  annotation=stardb
)
pbmc_ENZ48 <- CreateSeuratObject(
  counts = chrom_assay_ENZ48,
  assay = "peaks",
  project = 'ATAC',
  meta.data = metadata_ENZ48
)

chrom_assay_RESA <- CreateChromatinAssay(
  counts = counts_RESA,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = 'GSM5155453/fragments.tsv.gz',
  annotation=stardb
)
pbmc_RESA <- CreateSeuratObject(
  counts = chrom_assay_RESA,
  assay = "peaks",
  project = 'ATAC',
  meta.data = metadata_RESA
)

chrom_assay_RESB <- CreateChromatinAssay(
  counts = counts_RESB,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = 'GSM5155454/fragments.tsv.gz',
  annotation=stardb
)
pbmc_RESB <- CreateSeuratObject(
  counts = chrom_assay_RESB,
  assay = "peaks",
  project = 'ATAC',
  meta.data = metadata_RESB
)

# compute nucleosome signal score per cell
pbmc_DMSO <- NucleosomeSignal(object = pbmc_DMSO);
pbmc_ENZ48 <- NucleosomeSignal(object = pbmc_ENZ48);
pbmc_RESA <- NucleosomeSignal(object = pbmc_RESA);
pbmc_RESB <- NucleosomeSignal(object = pbmc_RESB);

# compute TSS enrichment score per cell
pbmc_DMSO <- TSSEnrichment(object = pbmc_DMSO, fast = T);
pbmc_ENZ48 <- TSSEnrichment(object = pbmc_ENZ48, fast = T);
pbmc_RESA <- TSSEnrichment(object = pbmc_RESA, fast = T);
pbmc_RESB <- TSSEnrichment(object = pbmc_RESB, fast = T);

# add blacklist ratio and fraction of reads in peaks
pbmc_DMSO$pct_reads_in_peaks <- pbmc_DMSO$peak_region_fragments / pbmc_DMSO$passed_filters * 100
pbmc_ENZ48$pct_reads_in_peaks <- pbmc_ENZ48$peak_region_fragments / pbmc_ENZ48$passed_filters * 100
pbmc_RESA$pct_reads_in_peaks <- pbmc_RESA$peak_region_fragments / pbmc_RESA$passed_filters * 100
pbmc_RESB$pct_reads_in_peaks <- pbmc_RESB$peak_region_fragments / pbmc_RESB$passed_filters * 100

# QC
pbmc_DMSO <- pbmc_DMSO[,pbmc_DMSO@meta.data$peak_region_fragments > 2000 & pbmc_DMSO@meta.data$peak_region_fragments < 20000
             & pbmc_DMSO@meta.data$pct_reads_in_peaks > 30 
             & pbmc_DMSO@meta.data$TSS.enrichment > 2
             & pbmc_DMSO@meta.data$nucleosome_signal < 9]


pbmc_ENZ48 <- pbmc_ENZ48[,pbmc_ENZ48@meta.data$peak_region_fragments > 1000 & pbmc_ENZ48@meta.data$peak_region_fragments < 20000
             & pbmc_ENZ48@meta.data$pct_reads_in_peaks > 30 
             & pbmc_ENZ48@meta.data$TSS.enrichment > 2
             & pbmc_ENZ48@meta.data$nucleosome_signal < 9]

pbmc_RESA <- pbmc_RESA[,pbmc_RESA@meta.data$peak_region_fragments > 2000 & pbmc_RESA@meta.data$peak_region_fragments < 20000
             & pbmc_RESA@meta.data$pct_reads_in_peaks > 40 
             & pbmc_RESA@meta.data$TSS.enrichment > 2
             & pbmc_RESA@meta.data$nucleosome_signal < 8]

pbmc_RESB <- pbmc_RESB[,pbmc_RESB@meta.data$peak_region_fragments > 2000 & pbmc_RESB@meta.data$peak_region_fragments < 25000
             & pbmc_RESB@meta.data$pct_reads_in_peaks > 30 
             & pbmc_RESB@meta.data$TSS.enrichment > 2
             & pbmc_RESB@meta.data$nucleosome_signal < 8]

# merge all 4 smaples
pbmc_DMSO$orig.ident <- "DMSO"
save(pbmc_DMSO, file="pbmc_DMSO.RData")

pbmc_ENZ48$orig.ident <- "ENZ48"
save(pbmc_ENZ48, file="pbmc_ENZ48.RData")

pbmc_RESA$orig.ident <- "RESA"
save(pbmc_RESA, file="pbmc_RESA.RData")

pbmc_RESB$orig.ident <- "RESB"
save(pbmc_RESB, file="pbmc_RESB.RData")

unintegrated <- merge(
  x = pbmc_DMSO,
  y = list(pbmc_ENZ48, pbmc_RESA, pbmc_RESB),
  add.cell.ids = c("DMSO", "ENZ48", "RESA", "RESB")
)

# Normalization and dimensionality reduction
unintegrated <- RunTFIDF(unintegrated)
unintegrated <- FindTopFeatures(unintegrated, min.cutoff = 50)
unintegrated <- RunSVD(unintegrated, n = 50, reduction.name = 'lsi', reduction.key = 'LSI_')

# Integration with Harmony
hm.integrated <- RunHarmony(
  object = unintegrated,
  group.by.vars = 'orig.ident',
  reduction = 'lsi',
  assay.use = 'peaks',  ## this is not running through, the assay is peaks right now
  project.dim = FALSE
)

hm.integrated <- RunUMAP(hm.integrated, reduction = "harmony", dims = 1:20)
hm.integrated <- FindNeighbors(hm.integrated, reduction = "harmony", dims = 1:20) 
hm.integrated <- FindClusters(object = hm.integrated, resolution = 0.9, algorithm = 3)
save(hm.integrated, file="Clustered_ATAC.RData");

DimPlot(hm.integrated, split.by = "orig.ident", ncol = 4);
dev.copy2pdf(file='atac_cluster.pdf', family='ArialMT', paper='a4');

## one plot with shapes
library(cowplot);
hui1<-DimPlot(hm.integrated, shape.by = "orig.ident", pt.size = .3) +  coord_fixed() +
  theme(legend.position="bottom");
zoom.in.purple<-DimPlot(hm.integrated, shape.by = "orig.ident", pt.size = 1.2) +
  xlim(c(-3.5, -1)) + ylim(c(-5.5, -3.5)) + theme(legend.position = "none") +
  coord_fixed() 

plot_grid(hui1, '', '', zoom.in.purple, '', '', ncol=3); ## re-scale the panel for later assembly
dev.copy2pdf(file='merged atac-seq umap.pdf', family='ArialMT', height=8, width=12);

# barplot
pt <- table(Idents(hm.integrated), hm.integrated$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
pt$pro <- ifelse(pt$Var2 == "DMSO", pt$Freq/3284, 
                 ifelse(pt$Var2 == "ENZ48", pt$Freq/3115, 
                        ifelse(pt$Var2 == "RESA", pt$Freq/3823, pt$Freq/3227)))

xx<-ggplot(pt, aes(x = factor(Var1, levels = c(as.character(seq(0, 12)))), y = pro, fill = Var2)) +
  xlab("scATAC-seq cluster") +
  ylab("Proportion of cells in cluster") +
  geom_bar(position="dodge", stat="identity")
xx<-as_grob(xx);
plot_grid(xx, '', '', '', '', '', ncol=3);
dev.copy2pdf(file='atac-seq barchart for cluster count.pdf', family='ArialMT', height=8, width=12)

# matching cluster why this way?
# hm.integrated$group <- ifelse(hm.integrated$seurat_clusters ==3 | hm.integrated$seurat_clusters ==6, "initial cluster",
#                              ifelse(hm.integrated$seurat_clusters ==10 | hm.integrated$seurat_clusters == 7 | hm.integrated$seurat_clusters == 8, "ENZ-induced cluster", 
#                                     ifelse(hm.integrated$seurat_clusters ==9 | hm.integrated$seurat_clusters ==4 | hm.integrated$seurat_clusters ==9 | hm.integrated$seurat_clusters ==11 | hm.integrated$seurat_clusters ==12,"Persistent cluster","NA")))

## persistent cluster 0 7
## induced cluster 2, 4, 5, 8, 9, 12
## initial cluster 1, 3
hm.integrated$group2<-'None';
hm.integrated$group2[which(hm.integrated$seurat_clusters%in%c(1, 3))]<-'Initial';
hm.integrated$group2[which(hm.integrated$seurat_clusters%in%c(0, 7))]<-'Persistent';
hm.integrated$group2[which(hm.integrated$seurat_clusters%in%c(2, 4))]<-'Induced';

table(hm.integrated$group2);
save(hm.integrated, file='hm.integrated_20211129.rda');

### go use server
## remove empty cells to refine average fold change 
setwd('/data2/jli39');
library(Seurat); library(parallel); library(Signac);

initial.cells<-colnames(hm.integrated)[which(hm.integrated$group2=='Initial')];
induced.cells<-colnames(hm.integrated)[which(hm.integrated$group2=='Induced')];
persistent.cells<-colnames(hm.integrated)[which(hm.integrated$group2=='Persistent')];

## filter fc first then subset then feed to findmarkers for significance calculation
##---------------initial-induced------------------
relog2fc.all.initial_induced <-mclapply(rownames(hm.integrated), function(x){
  t1<-hm.integrated@assays$peaks@data[x, initial.cells];
  t2<-hm.integrated@assays$peaks@data[x, induced.cells];
  tt1<-sum(t1, na.rm = T);
  tt2<-sum(t2, na.rm = T);
  log2(tt2/tt1);
}, mc.cores = 20);

relog2fc.all.t.initial_induced <- t(as.data.frame(relog2fc.all.initial_induced))
row.names(relog2fc.all.t.initial_induced) <- rownames(hm.integrated)

hm.log2.filter.intial_induced <- subset(hm.integrated, features=rownames(hm.integrated)[which(abs(as.numeric(relog2fc.all.initial_induced))>=1)]);

df.initial_induced <- FindMarkers(hm.log2.filter.intial_induced, ident.1 = c("2", "4"), ident.2 = c("1","3"), 
                                  min.pct = 0.001, logfc.threshold=0.05);

# relog2fcs<-mclapply(rownames(df.initial_induced), function(x){
#  t1<-hm.integrated@assays$peaks@data[x, initial.cells];
#  t2<-hm.integrated@assays$peaks@data[x, induced.cells];
#  tt1<-sum(t1, na.rm = T);
#  tt2<-sum(t2, na.rm = T);
#  log2(tt2/tt1);
# }, mc.cores = 20);

# df.initial_induced$relog2fc<-as.numeric(unlist(relog2fcs));

# length(which(df.initial_induced$p_val_adj<1e-6 & abs(df.initial_induced$avg_log2FC)>=1)); ## 3043 about 1% of all the peaks
# save(df.initial_induced, file='df.initial_induced_20211129.rda');

df.initial_induced.f <- merge(df.initial_induced, relog2fc.all.t.initial_induced, by="row.names", all.x=T)
rownames(df.initial_induced.f) <- df.initial_induced.f$Row.names
df.initial_induced.f <- df.initial_induced.f[, c(2:7) ]
colnames(df.initial_induced.f)[6] <- "relog2fc"
df.initial_induced.p <- dplyr::filter(df.initial_induced.f, relog2fc != "Inf" & relog2fc != "-Inf")

library(EnhancedVolcano);
hh1<-EnhancedVolcano(df.initial_induced.p, lab = rownames(df.initial_induced.p), 
                x='relog2fc', y='p_val_adj', pCutoff=1e-6);
plot_grid(hh1, '', '', '', ncol=2);
dev.copy2pdf(file='ATAC-df.initial_induced.pdf', family='ArialMT', width = 10, height = 7);

ss.initial_induced <- dplyr::filter(df.initial_induced.f, p_val_adj < 1e-6)
save(ss.initial_induced, file='ss.initial_induced_0128.rda');

#------------------induced-persistent-----------------

df.persistent_induced <- FindMarkers(hm.log2.filter, ident.1 = c("2", "4"), ident.2 = c("0", "7"), 
                                     min.pct = 0.001, logfc.threshold=0.05);

df.persistent_induced.f <- merge(df.persistent_induced, relog2fc.all.t, by="row.names", all.x=T)
rownames(df.persistent_induced.f) <- df.persistent_induced.f$Row.names
df.persistent_induced.f <- df.persistent_induced.f[, c(2:7) ]
colnames(df.persistent_induced.f)[6] <- "relog2fc"
df.persistent_induced.p <- dplyr::filter(df.persistent_induced.f, relog2fc != "Inf" & relog2fc != "-Inf")

hh2<-EnhancedVolcano(df.persistent_induced.p, lab = rownames(df.persistent_induced.p), 
                     x='relog2fc', y='p_val_adj', pCutoff=1e-6);
plot_grid(hh2, '', '', '', ncol=2);
dev.copy2pdf(file='ATAC-df.induced_persistent.pdf', family='ArialMT');

ss.persistent_induced <- dplyr::filter(df.persistent_induced.f, p_val_adj < 1e-6)
save(ss.persistent_induced, file='ss.persistent_induced_20211202.rda');

#------------------persistent-initial-----------------
relog2fc.all.initial_persistent <-mclapply(rownames(hm.integrated), function(x){
  t1<-hm.integrated@assays$peaks@data[x, initial.cells];
  t2<-hm.integrated@assays$peaks@data[x, persistent.cells];
  tt1<-sum(t1, na.rm = T);
  tt2<-sum(t2, na.rm = T);
  log2(tt2/tt1);
}, mc.cores = 8);

relog2fc.all.t.initial_persistent <- t(as.data.frame(relog2fc.all.initial_persistent))
row.names(relog2fc.all.t.initial_persistent) <- rownames(hm.integrated)

hm.log2.filter.intial_persistent <- subset(hm.integrated, features=rownames(hm.integrated)[which(abs(as.numeric(relog2fc.all.initial_induced))>=1)]);

df.persistent_initial <- FindMarkers(hm.log2.filter.intial_persistent, ident.1 = c("0", "7"), ident.2 = c("1","3"), 
                                  min.pct = 0.001, logfc.threshold=0.05);

df.persistent_initial.f <- merge(df.persistent_initial, relog2fc.all.t.initial_persistent, by="row.names", all.x=T)
rownames(df.persistent_initial.f) <- df.persistent_initial.f$Row.names
df.persistent_initial.f <- df.persistent_initial.f[, c(2:7) ]
colnames(df.persistent_initial.f)[6] <- "relog2fc"
df.persistent_initial.p <- dplyr::filter(df.persistent_initial.f, relog2fc != "Inf" & relog2fc != "-Inf")

hh3<-EnhancedVolcano(df.persistent_initial.p, lab = rownames(df.persistent_initial.p), 
                    x='relog2fc', y='p_val_adj', pCutoff=1e-6);
plot_grid(hh3, '', '', '', ncol=2);
dev.copy2pdf(file='ATAC-df.persistent-initial.pdf', family='ArialMT');

ss.persistent_initial <- dplyr::filter(df.persistent_initial.f, p_val_adj < 1e-6)
save(ss.persistent_initial, file='ss.persistent_initial_0128.rda');

# direction and annotation
tdf.initial_induced <- ss.initial_induced
tdf.initial_induced$direction <- ifelse(tdf.initial_induced$avg_log2FC > 0, "open", "closed")
tdf.initial_induced <- cbind(rownames(tdf.initial_induced), data.frame(tdf.initial_induced, row.names=NULL))
names(tdf.initial_induced)[names(tdf.initial_induced) == 'rownames(tdf.initial_induced)'] <- 'query_region'
clo_genes_initial_induced <- ClosestFeature(hm.integrated, regions = tdf.initial_induced$query_region)
fin.df.initial_induced <- merge(tdf.initial_induced, clo_genes_initial_induced, by="query_region")
write.csv(fin.df.initial_induced, file='fin.df.initial_induced_0128.csv')

tdf.persistent_induced <- ss.persistent_induced
tdf.persistent_induced$direction <- ifelse(tdf.persistent_induced$avg_log2FC > 0, "open", "closed")
tdf.persistent_induced <- cbind(rownames(tdf.persistent_induced), data.frame(tdf.persistent_induced, row.names=NULL))
names(tdf.persistent_induced)[names(tdf.persistent_induced) == 'rownames(tdf.persistent_induced)'] <- 'query_region'
clo_genes_persistent_induced <- ClosestFeature(hm.integrated, regions = tdf.persistent_induced$query_region)
fin.df.persistent_induced <- merge(tdf.persistent_induced, clo_genes_persistent_induced, by="query_region")
write.csv(fin.df.persistent_induced, file='fin.df.initial_induced_20211129.csv')

tdf.persistent_initial <- ss.persistent_initial.t
tdf.persistent_initial$direction <- ifelse(tdf.persistent_initial$avg_log2FC > 0, "open", "closed")
tdf.persistent_initial <- cbind(rownames(tdf.persistent_initial), data.frame(tdf.persistent_initial, row.names=NULL))
names(tdf.persistent_initial)[names(tdf.persistent_initial) == 'rownames(tdf.persistent_initial)'] <- 'query_region'
clo_genes_persistent_initial <- ClosestFeature(hm.integrated, regions = tdf.persistent_initial$query_region)
fin.df.persistent_initial <- merge(tdf.persistent_initial, clo_genes_persistent_initial, by="query_region")
write.csv(fin.df.persistent_initial, file='fin.df.persistent_initial_0128.csv')