#----overlap DEGs and DARs------
#-------------------------------------------------
#-----------initial-induced-----------------------
#-------------------------------------------------
overlap.1 <- as.data.frame(intersect(rownames(ss.deg.initial_induced),fin.df.initial_induced$gene_name))
colnames(overlap.1)[1] <- "gene_name"
initial_induced_overlap_rna <- ss.deg.initial_induced[rownames(ss.deg.initial_induced) %in% overlap.1$gene_name,]
names(initial_induced_overlap_rna)[1:7] <- c("rna_p_val","rna_avg_log2FC","rna_pct.1","rna_pct.2","rna_p_val_adj","rna_relog2fc","rna_direction")
initial_induced_overlap_rna <- cbind(rownames(initial_induced_overlap_rna), data.frame(initial_induced_overlap_rna, row.names=NULL))
colnames(initial_induced_overlap_rna)[1] <- "gene_name"
initial_induced_overlap_atac <- fin.df.initial_induced[fin.df.initial_induced$gene_name %in% overlap.1$gene_name,]
names(initial_induced_overlap_atac)[c(2:8)] <- c("atac_p_val","atac_avg_log2FC","atac_pct.1","atac_pct.2","atac_p_val_adj","atac_relog2fc","atac_direction")
initial_induced_merged <- merge(initial_induced_overlap_rna, initial_induced_overlap_atac, by="gene_name")
initial_induced_merged$accessibility_expression.direction <- ifelse(initial_induced_merged$atac_direction == "up" & initial_induced_merged$rna_direction == "up", "up-up",
                                                                    ifelse(initial_induced_merged$atac_direction == "up" & initial_induced_merged$rna_direction == "down", "up-down",
                                                                           ifelse(initial_induced_merged$atac_direction == "down" & initial_induced_merged$rna_direction == "up", "down-up", "down-down")))
write.csv(initial_induced_merged,'initial_induced_merged_0128.csv');
table(initial_induced_merged$accessibility_expression.direction)
library(cowplot)
library(ggplot2)
tmp1<-initial_induced_merged[, c('accessibility_expression.direction', 'gene_name')];
tmp1<-unique(tmp1);
jj<-ggplot(tmp1, aes(x=accessibility_expression.direction)) +
  geom_bar() +
  stat_count(geom='text', color='white', aes(label=..count..),
             position=position_stack(vjust=0.5))+xlab("Regulation Direction (Induced-Initial)")+ylab('Count')
plot_grid(jj, '', '', '', ncol=2);
dev.copy2pdf(file='Regulation Direction (Induced-Initial).pdf', family='ArialMT')

#?# pie(initial_induced_merged$gene_type[which(initial_induced_merged$accessibility_expression.direction=='down-down')])
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### functional analysis 
### ### ### 
### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
library(clusterProfiler);
library(org.Hs.eg.db);

initial_induced.degs.df<-tmp1;
initial_induced.name2id<-bitr(initial_induced.degs.df$gene_name, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db');
initial_induced.name2id<-initial_induced.name2id[-which(duplicated(initial_induced.name2id[, 1])), ]; rownames(initial_induced.name2id)<-initial_induced.name2id$ALIAS;

initial_induced.degs.df$id<-initial_induced.name2id[initial_induced.degs.df$gene_name, 2];
colnames(initial_induced.degs.df)[1]<-c('direction');
initial_induced.go.all<-compareCluster(id~direction, data=initial_induced.degs.df[which(initial_induced.degs.df$direction!='up-down'), ], minGSSize=15, fun='enrichGO', OrgDb='org.Hs.eg.db', ont = 'BP', qvalueCutoff=0.05);
initial_induced.kegg.all<-compareCluster(id~direction, data=initial_induced.degs.df[which(initial_induced.degs.df$direction!='up-down'), ], minGSSize=15, fun='enrichKEGG', qvalueCutoff=0.05);

all.gsea<-read.gmt('Analysis 1202/msigdb.v7.4.entrez.gmt.txt');
#### remove gene sets with out clear information in their names
## start with module_
pos<-grep('^MODULE_', all.gsea[, 1]);
all.gsea<-all.gsea[-pos, ];

initial_induced.degs.df<-initial_induced.degs.df[which(initial_induced.degs.df$direction!='up-down'), ];
initial_induced.deg.list<-tapply(initial_induced.degs.df$id, initial_induced.degs.df$direction, list);

initial_induced.all.gsea.lt<-lapply(initial_induced.deg.list, function(x)enricher(x, TERM2GENE = all.gsea, qvalueCutoff = .05));
save(initial_induced.all.gsea.lt, file='initial_induced.all.gsea.lt.rda');

library(ggnewscale); library(enrichplot); library(wordcloud); library(xlsx)

write.xlsx(initial_induced.all.gsea.lt$`down-down`, file="initial_induced.gsea.xlsx", sheetName="down-down", row.names=FALSE)
write.xlsx(initial_induced.all.gsea.lt$`down-up`, file="initial_induced.gsea.xlsx", sheetName="down-up", append=TRUE, row.names=FALSE)
write.xlsx(initial_induced.all.gsea.lt$`up-up`, file="initial_induced.gsea.xlsx", sheetName="up-up", append=TRUE, row.names=FALSE)

edox<-initial_induced.all.gsea.lt$`down-down`;
## add in gene names for readable information for top 50 only now
edox@result$geneName<-'';
edox@result$geneName[1:50]<-sapply(edox@result$geneID[1:50], function(x){
  t1<-strsplit(x, split='/')[[1]];
  t2<-bitr(t1, fromType = 'ENTREZID', toType = 'ALIAS', OrgDb = 'org.Hs.eg.db');
  t2<-unique(t2[, 2]);
  paste(t2, collapse ='; ');
  });

edox2 <- pairwise_termsim(edox);
p1 <- treeplot(edox2, offset=2, fontsize=3);
dev.copy2pdf(file='tree plot for up-up top50 in initial_induced.pdf', family='ArialMT',width = 14, height = 9);
## plot gene word cloud for each tree plot
p1.cluster<-p1$data; ## remove level 0
p1.cluster.list<-tapply(p1.cluster$label, p1.cluster$group, list)[-1];
p1.cluster.genes<-lapply(1:length(p1.cluster.list), function(y){
  x<-p1.cluster.list[[y]];
  t1<-edox@result[which(edox@result$ID%in%x), 'geneName'];
  t2<-sapply(t1,strsplit, split='; ');
  t3<-unlist(t2); t4<-sort(table(t3), decreasing = T);
  df<-data.frame(names(t4), as.numeric(t4), names(p1.cluster.list)[y]);
  colnames(df)<-c('word', 'freq', 'Cluster');
  return(df);
})

library(ggwordcloud);
loc.up2dn<-c(1,2,3,4,5); ### re-order the clusters based on the clustering from top to down
### plot word cloud for each group in the tree plot
tmp.word<-p1.cluster.genes[loc.up2dn];
# n=100; ### for an easy view
# tmp<-lapply(tmp, function(x){if(dim(x)[1]>=n){return(x[1:n, ])}else{return(x)}})
# tmp<-do.call('rbind', tmp);
## get a table
# set.seed(42);
# wc1<-ggplot(tmp, aes(label = word, size = freq, col=factor(sample.int(10, nrow(tmp), replace=T)))) +
#  geom_text_wordcloud_area(shape='circle') +
#  scale_size_area(max_size = 24) +
#  theme_minimal() + facet_wrap(~factor(Cluster, levels = names(p1.cluster.list)[loc.up2dn]), ncol = 1);
par(mfcol=c(length(loc.up2dn), 2), mar=rep(0, 4));
for(i in 1:length(loc.up2dn)){
  wordcloud(words = tmp.word[[i]]$word, 
            freq = tmp.word[[i]]$freq, min.freq = 1,     
            max.words=150, random.order=FALSE,             
            colors=brewer.pal(8, "Set2"));}
dev.copy2pdf(file='initial_induced.all.gsea.lt_down-down_wc.pdf', family='ArialMT');

# edox<-initial_induced.all.gsea.lt$`down-up`;
# edox2 <- pairwise_termsim(edox)
# p2 <- treeplot(edox2, offset=2, fontsize=3)

# edox<-initial_induced.all.gsea.lt$`up-up`;
# edox2 <- pairwise_termsim(edox)
# p3 <- treeplot(edox2, offset=2, fontsize=3)

# library(cowplot);
# plot_grid(p1, p2, p3, ncol=1, labels = c('A. Down-Down', 'B. Down-Up', 'C. Up-Up'));
# dev.copy2pdf(file='GSEA all top30 each class_initial_induced.pdf', family='ArialMT',width = 14, height = 26);

#-------------------------------------------------
#--------------------persistent-induced-----------
#-------------------------------------------------
overlap.2 <- as.data.frame(intersect(rownames(ss.deg.persistent_induced),fin.df.persistent_induced$gene_name))
colnames(overlap.2)[1] <- "gene_name"
persistent_induced_overlap_rna <- ss.deg.persistent_induced[rownames(ss.deg.persistent_induced) %in% overlap.2$gene_name,]
names(persistent_induced_overlap_rna)[1:7] <- c("rna_p_val","rna_avg_log2FC","rna_pct.1","rna_pct.2","rna_p_val_adj","rna_relog2fc","rna_direction")
persistent_induced_overlap_rna <- cbind(rownames(persistent_induced_overlap_rna), data.frame(persistent_induced_overlap_rna, row.names=NULL))
colnames(persistent_induced_overlap_rna)[1] <- "gene_name"
persistent_induced_overlap_atac <- fin.df.persistent_induced[fin.df.persistent_induced$gene_name %in% overlap.2$gene_name,]
names(persistent_induced_overlap_atac)[c(2:8)] <- c("atac_p_val","atac_avg_log2FC","atac_pct.1","atac_pct.2","atac_p_val_adj","atac_relog2fc","atac_direction")
persistent_induced_merged <- merge(persistent_induced_overlap_rna, persistent_induced_overlap_atac, by="gene_name")
persistent_induced_merged$accessibility_expression.direction <- ifelse(persistent_induced_merged$atac_direction == "up" & persistent_induced_merged$rna_direction == "up", "up-up",
                                                                    ifelse(persistent_induced_merged$atac_direction == "up" & persistent_induced_merged$rna_direction == "down", "up-down",
                                                                           ifelse(persistent_induced_merged$atac_direction == "down" & persistent_induced_merged$rna_direction == "up", "down-up", "down-down")))
write.csv(persistent_induced_merged,'persistent_induced_merged.csv')
table(persistent_induced_merged$accessibility_expression.direction)

tmp2<-persistent_induced_merged[, c('accessibility_expression.direction', 'gene_name')];
tmp2<-unique(tmp2);
jj<-ggplot(tmp2, aes(x=accessibility_expression.direction)) +
  geom_bar() +
  stat_count(geom='text', color='white', aes(label=..count..),
             position=position_stack(vjust=0.5))+xlab("Regulation Direction (Induced-Persistent)")+ylab('Count')
plot_grid(jj, '', '', '', ncol=2);
dev.copy2pdf(file='Regulation Direction (Induced-Persistent).pdf', family='ArialMT')

### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### functional analysis 
### ### ### 
### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
persistent_induced.degs.df<-tmp2;
persistent_induced.name2id<-bitr(persistent_induced.degs.df$gene_name, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db');
persistent_induced.name2id<-persistent_induced.name2id[-which(duplicated(persistent_induced.name2id[, 1])), ]; rownames(persistent_induced.name2id)<-persistent_induced.name2id$ALIAS;

persistent_induced.degs.df$id<-persistent_induced.name2id[persistent_induced.degs.df$gene_name, 2];
colnames(persistent_induced.degs.df)[1]<-c('direction');
persistent_induced.go.all<-compareCluster(id~direction, data=persistent_induced.degs.df[which(persistent_induced.degs.df$direction!='down-up'), ], minGSSize=15, fun='enrichGO', OrgDb='org.Hs.eg.db', ont = 'BP', qvalueCutoff=0.05);
persistent_induced.kegg.all<-compareCluster(id~direction, data=persistent_induced.degs.df[which(persistent_induced.degs.df$direction!='down-up'), ], minGSSize=15, fun='enrichKEGG', qvalueCutoff=0.05);

#### remove gene sets with out clear information in their names
## start with module_
persistent_induced.degs.df<-persistent_induced.degs.df[which(persistent_induced.degs.df$direction!='down-up'), ];
persistent_induced.deg.list<-tapply(persistent_induced.degs.df$id, persistent_induced.degs.df$direction, list);

persistent_induced.all.gsea.lt<-lapply(persistent_induced.deg.list, function(x)enricher(x, TERM2GENE = all.gsea, qvalueCutoff = .05));

library(xlsx)
write.xlsx(persistent_induced.all.gsea.lt$`down-down`, file="persistent_induced.gsea.xlsx", sheetName="down-down", row.names=FALSE)
write.xlsx(persistent_induced.all.gsea.lt$`down-up`, file="persistent_induced.gsea.xlsx", sheetName="down-up", append=TRUE, row.names=FALSE)

persistent_induced.edox<-persistent_induced.all.gsea.lt$`down-down`;
## add in gene names for readable information for top 50 only now
persistent_induced.edox@result$geneName<-'';
persistent_induced.edox@result$geneName[1:50]<-sapply(persistent_induced.edox@result$geneID[1:50], function(x){
  t1<-strsplit(x, split='/')[[1]];
  t2<-bitr(t1, fromType = 'ENTREZID', toType = 'ALIAS', OrgDb = 'org.Hs.eg.db');
  t2<-unique(t2[, 2]);
  paste(t2, collapse ='; ');
});

persistent_induced.edox2 <- pairwise_termsim(persistent_induced.edox);
persistent_induced.p1 <- treeplot(persistent_induced.edox2, offset=2, fontsize=3);
dev.copy2pdf(file='tree plot for down-down top50 in persistent_induced.pdf', family='ArialMT',width = 14, height = 6);

## plot gene word cloud for each tree plot
persistent_induced.p1.cluster<-persistent_induced.p1$data; ## remove level 0
persistent_induced.p1.cluster.list<-tapply(persistent_induced.p1.cluster$label, persistent_induced.p1.cluster$group, list)[-1];
persistent_induced.p1.cluster.genes<-lapply(persistent_induced.p1.cluster.list, function(x){
  t1<-persistent_induced.edox@result[which(persistent_induced.edox@result$ID%in%x), 'geneName'];
  t2<-sapply(t1,strsplit, split='; ');
  t3<-unlist(t2); t4<-sort(table(t3), decreasing = T);
  df<-data.frame(names(t4), as.numeric(t4));
  colnames(df)<-c('word', 'freq')
  return(df);
})

### plot word cloud for each group in the tree plot
pdf(file='Word Could for persistent-induced.pdf', family='ArialMT',width = 13, height = 13) 
for(i in persistent_induced.p1.cluster.genes) {
  wordcloud(words = i$word, 
            freq = i$freq, min.freq = 1,           
            max.words=50, random.order=FALSE, rot.per=0.35,            
            colors=brewer.pal(8, "Dark2"))
}
dev.off()

#-------------------------------------------------
#--------------------persistent-initial-----------
#-------------------------------------------------
overlap.3 <- as.data.frame(intersect(rownames(ss.deg.persistent_initial),fin.df.persistent_initial$gene_name))
colnames(overlap.3)[1] <- "gene_name"
persistent_initial_overlap_rna <- ss.deg.persistent_initial[rownames(ss.deg.persistent_initial) %in% overlap.3$gene_name,]
names(persistent_initial_overlap_rna)[1:7] <- c("rna_p_val","rna_avg_log2FC","rna_pct.1","rna_pct.2","rna_p_val_adj","rna_relog2fc","rna_direction")
persistent_initial_overlap_rna <- cbind(rownames(persistent_initial_overlap_rna), data.frame(persistent_initial_overlap_rna, row.names=NULL))
colnames(persistent_initial_overlap_rna)[1] <- "gene_name"
persistent_initial_overlap_atac <- fin.df.persistent_initial[fin.df.persistent_initial$gene_name %in% overlap.3$gene_name,]
names(persistent_initial_overlap_atac)[c(2:8)] <- c("atac_p_val","atac_avg_log2FC","atac_pct.1","atac_pct.2","atac_p_val_adj","atac_relog2fc","atac_direction")
persistent_initial_merged <- merge(persistent_initial_overlap_rna, persistent_initial_overlap_atac, by="gene_name")
persistent_initial_merged$accessibility_expression.direction <- ifelse(persistent_initial_merged$atac_direction == "up" & persistent_initial_merged$rna_direction == "up", "up-up",
                                                                       ifelse(persistent_initial_merged$atac_direction == "up" & persistent_initial_merged$rna_direction == "down", "up-down",
                                                                              ifelse(persistent_initial_merged$atac_direction == "down" & persistent_initial_merged$rna_direction == "up", "down-up", "down-down")))
write.csv(persistent_initial_merged,'persistent_initial_merged_0128.csv')
table(persistent_initial_merged$accessibility_expression.direction)

tmp3<-persistent_initial_merged[, c('accessibility_expression.direction', 'gene_name')];
tmp3<-unique(tmp3);
jj<-ggplot(tmp3, aes(x=accessibility_expression.direction)) +
  geom_bar() +
  stat_count(geom='text', color='white', aes(label=..count..),
             position=position_stack(vjust=0.5))+xlab("Regulation Direction (Persistent-Initial)")+ylab('Count')
plot_grid(jj, '', '', '', ncol=2);
dev.copy2pdf(file='Regulation Direction (Persistent-Initial).pdf', family='ArialMT')

### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### functional analysis 
### ### ### 
### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
persistent_initial.degs.df<-tmp3;
persistent_initial.name2id<-bitr(persistent_initial.degs.df$gene_name, fromType = 'ALIAS', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db');
persistent_initial.name2id<-persistent_initial.name2id[-which(duplicated(persistent_initial.name2id[, 1])), ]; rownames(persistent_initial.name2id)<-persistent_initial.name2id$ALIAS;

persistent_initial.degs.df$id<-persistent_initial.name2id[persistent_initial.degs.df$gene_name, 2];
colnames(persistent_initial.degs.df)[1]<-c('direction');
persistent_initial.go.all<-compareCluster(id~direction, data=persistent_initial.degs.df[which(persistent_initial.degs.df$direction!='down-down'), ], minGSSize=15, fun='enrichGO', OrgDb='org.Hs.eg.db', ont = 'BP', qvalueCutoff=0.05);
persistent_initial.kegg.all<-compareCluster(id~direction, data=persistent_initial.degs.df[which(persistent_initial.degs.df$direction!='down-down'), ], minGSSize=15, fun='enrichKEGG', qvalueCutoff=0.05);

#### remove gene sets with out clear information in their names
## start with module_
persistent_initial.degs.df<-persistent_initial.degs.df[which(persistent_initial.degs.df$direction!='down-down'), ];
persistent_initial.deg.list<-tapply(persistent_initial.degs.df$id, persistent_initial.degs.df$direction, list);

persistent_initial.all.gsea.lt<-lapply(persistent_initial.deg.list, function(x)enricher(x, TERM2GENE = all.gsea, qvalueCutoff = .05));
save(persistent_initial.all.gsea.lt, file='persistent_initial.all.gsea.lt.rda');

write.xlsx(persistent_initial.all.gsea.lt$`up-down`, file="persistent_initial.gsea.xlsx", sheetName="up-down", row.names=FALSE)
write.xlsx(persistent_initial.all.gsea.lt$`up-up`, file="persistent_initial.gsea.xlsx", sheetName="up-up", append=TRUE, row.names=FALSE)

persistent_initial.edox<-persistent_initial.all.gsea.lt$`up-up`;

## add in gene names for readable information for top 50 only now
persistent_initial.edox@result$geneName<-'';
persistent_initial.edox@result$geneName[1:50]<-sapply(persistent_initial.edox@result$geneID[1:50], function(x){
  t1<-strsplit(x, split='/')[[1]];
  t2<-bitr(t1, fromType = 'ENTREZID', toType = 'ALIAS', OrgDb = 'org.Hs.eg.db');
  t2<-unique(t2[, 2]);
  paste(t2, collapse ='; ');
});

edox2 <- pairwise_termsim(persistent_initial.edox);
p2 <- treeplot(edox2, offset=2, fontsize=3);
dev.copy2pdf(file='tree plot for up-up top50 in persistent_initial.pdf', family='ArialMT',width = 14, height = 9);

## plot gene word cloud for each tree plot
p2.cluster<-p2$data; ## remove level 0
p2.cluster.list<-tapply(p2.cluster$label, p2.cluster$group, list)[-1];
p2.cluster.genes<-lapply(1:length(p2.cluster.list), function(y){
  x<-p2.cluster.list[[y]];
  t1<-edox@result[which(edox@result$ID%in%x), 'geneName'];
  t2<-sapply(t1,strsplit, split='; ');
  t3<-unlist(t2); t4<-sort(table(t3), decreasing = T);
  df<-data.frame(names(t4), as.numeric(t4), names(p2.cluster.list)[y]);
  colnames(df)<-c('word', 'freq', 'Cluster');
  return(df);
})

loc.up2dn.2<-c(1, 4, 2, 5, 3); ### re-order the clusters based on the clustering from top to down
### plot word cloud for each group in the tree plot
tmp2<-p2.cluster.genes[loc.up2dn.2];
par(mfcol=c(length(loc.up2dn.2), 2), mar=rep(0, 4));
for(i in 1:length(loc.up2dn.2)){
  wordcloud(words = tmp2[[i]]$word, 
            freq = tmp2[[i]]$freq, min.freq = 1,     
            max.words=200, random.order=FALSE,             
            colors=brewer.pal(8, "Set2"));}
dev.copy2pdf(file='persistent_initial.all.gsea.lt_up-up_wc.pdf', family='ArialMT', height=13);

# --------------------------
# -------venn plot----------
# --------------------------
library(eulerr)
set.seed(1234)
d1 <- filter(ss.deg.initial_induced, direction == 'down')
d2 <- filter(fin.df.initial_induced, direction == 'down')
u1 <- filter(ss.deg.initial_induced, direction == 'up')
u2 <- filter(fin.df.initial_induced, direction == 'up')

dd <- list(rna.down = unique(rownames(d1)), atac.down = unique(d2$gene_name))
du <- list(rna.down = unique(rownames(d1)), atac.up = unique(u2$gene_name))
ud <- list(rna.up = unique(rownames(u1)), atac.down = unique(d2$gene_name))
uu <- list(rna.up = unique(rownames(u1)), atac.up = unique(u2$gene_name))

v1 <- plot(euler(dd, shape = "ellipse"), quantities = TRUE)
v2 <- plot(euler(du, shape = "ellipse"), quantities = TRUE)
v3 <- plot(euler(ud, shape = "ellipse"), quantities = TRUE)
v4 <- plot(euler(uu, shape = "ellipse"), quantities = TRUE)

library(ggpubr)
ggarrange(v1, v2, v3, v4, 
          labels = c("down-down", "down-up", "up-down","up-up"),
          ncol = 2, nrow = 2)
dev.copy2pdf(file='initial_induced.overlap.pdf', family='ArialMT');

#-----------------venn plot--------------------
d1 <- filter(ss.deg.persistent_initial, direction == 'down')
d2 <- filter(fin.df.persistent_initial, direction == 'down')
u1 <- filter(ss.deg.persistent_initial, direction == 'up')
u2 <- filter(fin.df.persistent_initial, direction == 'up')

dd <- list(rna.down = unique(rownames(d1)), atac.down = unique(d2$gene_name))
du <- list(rna.down = unique(rownames(d1)), atac.up = unique(u2$gene_name))
ud <- list(rna.up = unique(rownames(u1)), atac.down = unique(d2$gene_name))
uu <- list(rna.up = unique(rownames(u1)), atac.up = unique(u2$gene_name))

v1 <- plot(euler(dd, shape = "ellipse"), quantities = TRUE)
v2 <- plot(euler(du, shape = "ellipse"), quantities = TRUE)
v3 <- plot(euler(ud, shape = "ellipse"), quantities = TRUE)
v4 <- plot(euler(uu, shape = "ellipse"), quantities = TRUE)

ggarrange(v1, v2, v3, v4, 
          labels = c("down-down", "down-up", "up-down","up-up"),
          ncol = 2, nrow = 2)
dev.copy2pdf(file='persistent_induced.overlap.pdf', family='ArialMT')

#---------------upset plot----------------------
library(UpSetR)
library(dplyr)
initial_induced.direction <- initial_induced_merged %>% dplyr::select(c('gene_name', 'accessibility_expression.direction'))
initial_induced.up_up <- initial_induced.direction[initial_induced.direction$accessibility_expression.direction == "up-up", ]    
colnames(initial_induced.up_up)[2] <- "Induced vs. Initial (Up/Up)"
initial_induced.up_up <- unique(initial_induced.up_up, by="gene_name")
initial_induced.up_down <- initial_induced.direction[initial_induced.direction$accessibility_expression.direction == "up-down", ]    
colnames(initial_induced.up_down)[2] <- "Induced vs. Initial (Up/Down)"
initial_induced.up_down <- unique(initial_induced.up_down, by="gene_name")
initial_induced.down_up <- initial_induced.direction[initial_induced.direction$accessibility_expression.direction == "down-up", ]    
colnames(initial_induced.down_up)[2] <- "Induced vs. Initial (Down/Up)"
initial_induced.down_up <- unique(initial_induced.down_up, by="gene_name")
initial_induced.down_down <- initial_induced.direction[initial_induced.direction$accessibility_expression.direction == "down-down", ]    
colnames(initial_induced.down_down)[2] <- "Induced vs. Initial (Down/Down)"
initial_induced.down_down <- unique(initial_induced.down_down, by="gene_name")

persistent_initial.direction <- persistent_initial_merged %>% dplyr::select(c('gene_name', 'accessibility_expression.direction'))
persistent_initial.up_up <- persistent_initial.direction[persistent_initial.direction$accessibility_expression.direction == "up-up", ]    
colnames(persistent_initial.up_up)[2] <- "Persistent vs. Initial (Up/Up)"
persistent_initial.up_up <- unique(persistent_initial.up_up, by="gene_name")
persistent_initial.up_down <- persistent_initial.direction[persistent_initial.direction$accessibility_expression.direction == "up-down", ]    
colnames(persistent_initial.up_down)[2] <- "Persistent vs. Initial (Up/Down)"
persistent_initial.up_down <- unique(persistent_initial.up_down, by="gene_name")
persistent_initial.down_down <- persistent_initial.direction[persistent_initial.direction$accessibility_expression.direction == "down-down", ]    
colnames(persistent_initial.down_down)[2] <- "Persistent vs. Initial (Down/Down)"
persistent_initial.down_down <- unique(persistent_initial.down_down, by="gene_name")

library(tidyverse)
plot_list <- list(initial_induced.up_up, initial_induced.up_down, initial_induced.down_up, initial_induced.down_down,
                  persistent_initial.up_up, persistent_initial.up_down, persistent_initial.down_down)
all_plot <- plot_list %>% purrr::reduce(full_join, by='gene_name')

all_plot$`Induced vs. Initial (Up/Up)`[which(all_plot$`Induced vs. Initial (Up/Up)`=='up-up')]<-'TRUE'
all_plot$`Induced vs. Initial (Up/Down)`[which(all_plot$`Induced vs. Initial (Up/Down)`=="up-down")]<-'TRUE'
all_plot$`Induced vs. Initial (Down/Up)`[which(all_plot$`Induced vs. Initial (Down/Up)`=="down-up")]<-'TRUE'
all_plot$`Induced vs. Initial (Down/Down)`[which(all_plot$`Induced vs. Initial (Down/Down)`=="down-down")]<-'TRUE'
all_plot$`Persistent vs. Initial (Up/Up)`[which(all_plot$`Persistent vs. Initial (Up/Up)`=="up-up")]<-'TRUE'
all_plot$`Persistent vs. Initial (Up/Down)`[which(all_plot$`Persistent vs. Initial (Up/Down)`=="up-down")]<-'TRUE'
all_plot$`Persistent vs. Initial (Down/Down)`[which(all_plot$`Persistent vs. Initial (Down/Down)`=="down-down")]<-'TRUE'

all_plot[is.na(all_plot)==T] <- "FALSE"

up_down <- colnames(all_plot)[-1];
names(up_down) <- up_down

library(ComplexUpset)
h1<-upset(data = all_plot, intersect = up_down, 
          name="Overlapped Genes Frequency", 
          min_size = 0,
          width_ratio = 0.125) +
  labs(title = "")
plot_grid(h1, '', '', '', ncol=2);
dev.copy2pdf(file='UpSet for all gene sets.pdf', family='ArialMT', height=21, width=17);

######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### 
######### a very simply upset plot  ##############
plot_list <- list(initial_induced.up_up, initial_induced.up_down, initial_induced.down_up, initial_induced.down_down,
                  persistent_initial.up_up, persistent_initial.up_down, persistent_initial.down_down)
names(plot_list)<-as.character(unlist(lapply(plot_list, function(x)colnames(x)[2])))
ecm.venn.up<-lapply(plot_list, function(x)unique(x[, 1])); 
names(ecm.venn.up)<-names(plot_list);
ecm.venn.up.df<-matrix(0, ncol=length(plot_list), nrow=length(unique(unlist(ecm.venn.up))));
colnames(ecm.venn.up.df)<-names(ecm.venn.up); 
rownames(ecm.venn.up.df)<-unique(unlist(ecm.venn.up));
for(i in 1:length(ecm.venn.up)){ecm.venn.up.df[ecm.venn.up[[i]], names(ecm.venn.up)[i]]<-1;}
t1<-upset(as.data.frame(ecm.venn.up.df), nsets = length(ecm.venn.up), mb.ratio = c(0.5, 0.5),
          order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))
t1
dev.copy2pdf(file='UpSet for all gene sets.pdf', family='ArialMT', paper='a4');


