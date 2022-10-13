### ### ### ### ### ### ### ### ### ### ### 
library(ggplot2)
library(readxl)
library(survival)
library(survminer)
library(dplyr)
library(AnnotationDbi)
library(rlist)
library(scater)

initial_induced_merged <- read.csv("Analysis 1202/results 1202/initial_induced_merged.csv");
persistent_initial_merged <- read.csv("Analysis 1202/results 1202/persistent_initial_merged.csv");

all.test <- list(persistent_initial.down_down, persistent_initial.up_down, persistent_initial.up_up, 
                 initial_induced.up_up, initial_induced.up_down, initial_induced.down_up, initial_induced.down_down);
for (i in 1:7) {
  all.test[[i]] <- filter(all.test[[i]], all.test[[i]][,"gene_name"] %in% rownames(tcga.expr.sybl)==TRUE);
  all.test[[i]][, 2]<-colnames(all.test[[i]])[2];
  colnames(all.test[[i]])<-NULL; rownames(all.test[[i]])<-NULL
  }
all.test.df<-cbind(gene=c(all.test[[1]][, 1], all.test[[2]][, 1], all.test[[3]][, 1], all.test[[4]][, 1], 
                     all.test[[5]][, 1], all.test[[6]][, 1], all.test[[7]][, 1]),
                   class=c(all.test[[1]][, 2], all.test[[2]][, 2], all.test[[3]][, 2], all.test[[4]][, 2], 
                     all.test[[5]][, 2], all.test[[6]][, 2], all.test[[7]][, 2]));

prad<-tcga.samples$case_submitter_id[which(tcga.samples$project_id=='TCGA-PRAD')];
expr.tumor<-names(tcga.expr.sybl)[which(substring(names(tcga.expr.sybl), 14, 15)!='11' &
                                        substring(names(tcga.expr.sybl), 1, 12)%in%prad)]; ## remove non-tumor samples
pos<-which(duplicated(substring(expr.tumor, 1, 12)));

expr.tumor.sht<-substring(expr.tumor, 1, 12)[-pos];

clin.data<-read_xlsx('Analysis 1202/TCGA clinical survival information_published.xlsx');
clin.data<-as.data.frame(clin.data);
clin.data$Time<-clin.data$PFI.time; ## change different survival time OS; DSS; DFI; PFI
rownames(clin.data)<-clin.data$bcr_patient_barcode;
clin.data$status<-clin.data$PFI;   ## change different survival event OS; DSS; DFI; PFI
clin.data.PC <- filter(clin.data, type == 'PRAD')

library(parallel);
pvals<-NULL;
  for (i in 1:length(all.test.df[, 1])){
  x<-all.test.df[, 1][i];
  print(i);
  gene1 <- log2(as.matrix(tcga.expr.sybl)[x, expr.tumor]+1);
  gene1 <- gene1[-pos];
  
  gene.df1<- data.frame(expr = gene1, status=clin.data.PC[expr.tumor.sht, 'status'], 
                             time=clin.data.PC[expr.tumor.sht, 'Time'], Tumor=clin.data.PC[expr.tumor.sht, 'type'])
  
  gene.df1<- filter(gene.df1, is.na(Tumor)==F)
  gene.df1$time <- as.numeric(gene.df1$time);
  
  med<-median(gene.df1$expr, na.rm = T);
  if(med==0){
    pvals<-c(pvals, NA);
  }else{
  positive <- paste(x, "+", sep = '')
  negative <- paste(x, "-", sep = '')
  gene.updn1 <- ifelse(gene.df1$expr>med, positive, negative)
  names(gene.updn1) <- rownames(gene.df1)
  gene.df1$gene <- unlist2(gene.updn1)[rownames(gene.df1)]
  
    fit1 <- survfit(Surv(gene.df1$time/365, gene.df1$status) ~ gene.df1$gene, data = gene.df1)
    pvalue1<- surv_pvalue(fit1, data = gene.df1, method = "survdiff")$pval;
   pvals<-c(pvals, pvalue1);
       }
  }

all.test.df<-as.data.frame(all.test.df);
all.test.df$pval<-pvals;
all.test.df$FDR<-p.adjust(pvals);
all.test.df.sig <- filter(all.test.df, FDR <= 0.05);
save(all.test.df.sig, file='all.test.df.sig.rda');

sig.df<-all.test.df.sig;
sig.df[5, 2]<-paste(sig.df[5, 2], '\n', sig.df[7, 2], sep='');
sig.df<-sig.df[-7, ]; rownames(sig.df)<-sig.df$gene;

sig.df<-lapply(sig.df$gene, function(x){
  gene1 <- log2(as.matrix(tcga.expr.sybl)[x, expr.tumor]+1);
  gene1 <- gene1[-pos];
  
  gene.df1<- data.frame(expr = gene1, status=clin.data.PC[expr.tumor.sht, 'status'], 
                        time=clin.data.PC[expr.tumor.sht, 'Time'], Tumor=clin.data.PC[expr.tumor.sht, 'type'])
  
  gene.df1<- filter(gene.df1, is.na(Tumor)==F)
  gene.df1$time <- as.numeric(gene.df1$time);
  
    med<-median(gene.df1$expr, na.rm = T);
    positive <- paste('Median', "+", sep = '')
    negative <- paste('Median', "-", sep = '')
    gene.updn1 <- ifelse(gene.df1$expr>med, positive, negative)
    names(gene.updn1) <- rownames(gene.df1)
    gene.df1$gene <- unlist2(gene.updn1)[rownames(gene.df1)]
    
    gene.df1$name=x; gene.df1$class=sig.df[x, 'class'];
    gene.df1$pval<-sig.df[x, 'pval'];
    return(gene.df1);
})

sig.df<-do.call('rbind', sig.df);
fit2 <- survfit( Surv(time/365, status) ~ gene, data = sig.df )
sig.df$Class<-paste(sig.df$class, sig.df$name, 
format(sig.df$pval, digits = 4, scientific = T),
              sep='\n');
ggsurvplot_facet(fit2, sig.df, facet.by = c('Class'), 
                 palette = c("blue", "brown"), pval = F) 
dev.copy2pdf(file='survival.sig.pdf', family='ArialMT', height=7, width=7)

################### a more complicated way #####################
geneset1 <- list()
for (i in all.test[[1]]$gene_name) {
  geneset1[[i]] <- log2(as.matrix(tcga.expr.sybl)[i, expr.tumor]+1)
  geneset1[[i]] <- geneset1[[i]][-pos]
}

expr.tumor<-substring(expr.tumor, 1, 12)[-pos];

clin.data<-read_xlsx('TCGA clinical survival information_published.xlsx');
clin.data<-as.data.frame(clin.data);
clin.data$Time<-clin.data$DSS.time; ## change different survival time OS; DSS; DFI; PFI
rownames(clin.data)<-clin.data$bcr_patient_barcode;
clin.data$status<-clin.data$DSS;   ## change different survival event OS; DSS; DFI; PFI
clin.data.PC <- filter(clin.data, type == 'PRAD')

# ----------------persistent_initial.down_down------------
gene.df1 <- list()
for (i in names(geneset1)) {
  gene.df1[[i]]<- data.frame(expr = geneset1[[i]], status=clin.data.PC[expr.tumor, 'status'], 
                            time=clin.data.PC[expr.tumor, 'Time'], Tumor=clin.data.PC[expr.tumor, 'type'])
  gene.df1[[i]] <- filter(gene.df1[[i]],is.na(Tumor)==F)
  gene.df1[[i]]$time <- as.numeric(gene.df1[[i]]$time);
}

gene.updn1 <- list()
for (i in names(gene.df1)) {
  med<-median(gene.df1[[i]]$expr)
  positive <- paste(i, "+", sep = '')
  negative <- paste(i, "-", sep = '')
  gene.updn1[[i]] <- ifelse(gene.df1[[i]]$expr>med, positive, negative)
  names(gene.updn1[[i]]) <- rownames(gene.df1[[i]])
  gene.df1[[i]]$gene <- unlist2(gene.updn1[[i]])[rownames(gene.df1[[i]])]
}

fit1 <- list()
pvalue1 <- list()
significant1 <- list()
for (i in names(gene.df1)) {
  fit1[[i]] <- survfit(Surv(gene.df1[[i]]$time/365, gene.df1[[i]]$status) ~ gene.df1[[i]]$Tumor + gene.df1[[i]]$gene, data = gene.df1[[i]])
  pvalue1[[i]] <- surv_pvalue(fit1[[i]], data = gene.df1[[i]], method = "survdiff")
  significant1[[i]] <- filter(pvalue1[[i]], pval <= 0.05)
}

# ----------------persistent_initial.up_down--------------
expr.tumor<-names(tcga.expr.sybl)[which(substring(names(tcga.expr.sybl), 14, 15)!='11')]; ## remove non-tumor samples
geneset2 <- list()
for (i in all.test[[2]]$gene_name) {
  geneset2[[i]] <- as.matrix(tcga.expr.sybl)[i, expr.tumor]
  geneset2[[i]] <- geneset2[[i]][-pos]
}

expr.tumor<-substring(expr.tumor, 1, 12)[-pos];

gene.df2 <- list()
for (i in names(geneset2)) {
  gene.df2[[i]]<- data.frame(expr = geneset2[[i]], status=clin.data.PC[expr.tumor, 'status'], 
                             time=clin.data.PC[expr.tumor, 'Time'], Tumor=clin.data.PC[expr.tumor, 'type'])
  gene.df2[[i]] <- filter(gene.df2[[i]],is.na(Tumor)==F)
  gene.df2[[i]]$time <- as.numeric(gene.df2[[i]]$time);
}

gene.updn2 <- list()
for (i in names(gene.df2)) {
  med<-median(gene.df2[[i]]$expr)
  positive <- paste(i, "+", sep = '')
  negative <- paste(i, "-", sep = '')
  gene.updn2[[i]] <- ifelse(gene.df2[[i]]$expr>med, positive, negative)
  names(gene.updn2[[i]]) <- rownames(gene.df2[[i]])
  gene.df2[[i]]$gene <- unlist2(gene.updn2[[i]])[rownames(gene.df2[[i]])]
}

fit2 <- list()
pvalue2 <- list()
significant2 <- list()
for (i in names(gene.df2)) {
  fit2[[i]] <- survfit(Surv(gene.df2[[i]]$time/365, gene.df2[[i]]$status) ~ gene.df2[[i]]$Tumor + gene.df2[[i]]$gene, data = gene.df2[[i]])
  pvalue2[[i]] <- surv_pvalue(fit2[[i]], data = gene.df2[[i]], method = "survdiff")
  significant2[[i]] <- filter(pvalue2[[i]], pval <= 0.05)
}

gplot2 <- list()
for (i in c("B3GNTL1", "C22orf39", "DUSP8", "GOLGA2", "KCNMA1", "KLF9", "NLGN1", "RBPMS" ,"UBE2H")) {
  fit2[[i]] <- survfit(Surv(time/365, status) ~ Tumor + gene, data = gene.df2[[i]])
  gplot2[[i]] <- ggsurvplot(fit2[[i]], gene.df2[[i]], palette = "jco", pval = TRUE)
}
res2 <- arrange_ggsurvplots(gplot2, print = FALSE, title = NA, ncol = 3, nrow = 3)
ggsave("persistent_initial.up_down survival plots.pdf", res2, width = 17, height = 16)

# ----------------persistent_initial.up_up--------------
expr.tumor<-names(tcga.expr.sybl)[which(substring(names(tcga.expr.sybl), 14, 15)!='11')]; ## remove non-tumor samples
geneset3 <- list()
for (i in all.test[[3]]$gene_name) {
  geneset3[[i]] <- as.matrix(tcga.expr.sybl)[i, expr.tumor]
  geneset3[[i]] <- geneset3[[i]][-pos]
}

expr.tumor<-substring(expr.tumor, 1, 12)[-pos];

gene.df3 <- list()
for (i in names(geneset3)) {
  gene.df3[[i]]<- data.frame(expr = geneset3[[i]], status=clin.data.PC[expr.tumor, 'status'], 
                             time=clin.data.PC[expr.tumor, 'Time'], Tumor=clin.data.PC[expr.tumor, 'type'])
  gene.df3[[i]] <- filter(gene.df3[[i]],is.na(Tumor)==F)
  gene.df3[[i]]$time <- as.numeric(gene.df3[[i]]$time);
}

gene.updn3 <- list()
for (i in names(gene.df3)) {
  med<-median(gene.df3[[i]]$expr)
  positive <- paste(i, "+", sep = '')
  negative <- paste(i, "-", sep = '')
  gene.updn3[[i]] <- ifelse(gene.df3[[i]]$expr>med, positive, negative)
  names(gene.updn3[[i]]) <- rownames(gene.df3[[i]])
  gene.df3[[i]]$gene <- unlist2(gene.updn3[[i]])[rownames(gene.df3[[i]])]
}

fit3 <- list()
pvalue3 <- list()
significant3 <- list()
for (i in names(gene.df3)) {
  fit3[[i]] <- survfit(Surv(gene.df3[[i]]$time/365, gene.df3[[i]]$status) ~ gene.df3[[i]]$Tumor + gene.df3[[i]]$gene, data = gene.df3[[i]])
  pvalue3[[i]] <- surv_pvalue(fit3[[i]], data = gene.df3[[i]], method = "survdiff")
  significant3[[i]] <- filter(pvalue3[[i]], pval <= 0.05)
}

gplot3 <- list()
for (i in c("ACTL6A", "AIFM1", "ATAD3A", "CCDC18", "DHX15", "FASN", "FRYL", "GOLIM4" ,"GSTA2", "HS6ST1", 
            "IGF1", "MBNL2", "MKNK2", "MRPS30", "MYOZ1", "NUP205", "PPIF", "PPP2R5E", "QSOX1", "RUVBL1",
            "SDCCAG8", "SGK2", "SMYD2", "TRA2B", "U2AF2", "ZMYND11")) {
  fit3[[i]] <- survfit(Surv(time/365, status) ~ Tumor + gene, data = gene.df3[[i]])
  gplot3[[i]] <- ggsurvplot(fit3[[i]], gene.df3[[i]], palette = "jco", pval = TRUE)
}
res3 <- arrange_ggsurvplots(gplot3, print = FALSE, title = NA, ncol = 3, nrow = 3)
ggsave("persistent_initial.up_up survival plots.pdf", res3, width = 17, height = 16)

# ----------------induced_initial.up_up--------------
expr.tumor<-names(tcga.expr.sybl)[which(substring(names(tcga.expr.sybl), 14, 15)!='11')]; ## remove non-tumor samples
geneset4 <- list()
for (i in all.test[[4]]$gene_name) {
  geneset4[[i]] <- as.matrix(tcga.expr.sybl)[i, expr.tumor]
  geneset4[[i]] <- geneset4[[i]][-pos]
}

expr.tumor<-substring(expr.tumor, 1, 12)[-pos];

gene.df4 <- list()
for (i in names(geneset4)) {
  gene.df4[[i]]<- data.frame(expr = geneset4[[i]], status=clin.data.PC[expr.tumor, 'status'], 
                             time=clin.data.PC[expr.tumor, 'Time'], Tumor=clin.data.PC[expr.tumor, 'type'])
  gene.df4[[i]] <- filter(gene.df4[[i]],is.na(Tumor)==F)
  gene.df4[[i]]$time <- as.numeric(gene.df4[[i]]$time);
}

gene.updn4 <- list()
for (i in names(gene.df4)) {
  med<-median(gene.df4[[i]]$expr)
  positive <- paste(i, "+", sep = '')
  negative <- paste(i, "-", sep = '')
  gene.updn4[[i]] <- ifelse(gene.df4[[i]]$expr>med, positive, negative)
  names(gene.updn4[[i]]) <- rownames(gene.df4[[i]])
  gene.df4[[i]]$gene <- unlist2(gene.updn4[[i]])[rownames(gene.df4[[i]])]
}

fit4 <- list()
pvalue4 <- list()
significant4 <- list()
for (i in names(gene.df4)) {
  fit4[[i]] <- survfit(Surv(time/365, status) ~ Tumor + gene, data = gene.df4[[i]])
  pvalue4[[i]] <- surv_pvalue(fit4[[i]], data = gene.df4[[i]], method = "survdiff")
  significant4[[i]] <- filter(pvalue4[[i]], pval <= 0.05)
}

gplot4 <- list()
for (i in c("FOXQ1", "MBNL2", "MGST3", "MRPS30", "PPARGC1A", "PPIF", "QSOX1", "SGK2" ,"SMYD2", "ZMYND11")) {
  fit4[[i]] <- survfit(Surv(time/365, status) ~ Tumor + gene, data = gene.df4[[i]])
  gplot4[[i]] <- ggsurvplot(fit4[[i]], gene.df4[[i]], palette = "jco", pval = TRUE)
}
res4 <- arrange_ggsurvplots(gplot4, print = FALSE, title = NA, ncol = 3, nrow = 3)
ggsave("induced_initial.up_up survival plots.pdf", res4, width = 17, height = 16)

# ----------------induced_initial.up_down--------------
expr.tumor<-names(tcga.expr.sybl)[which(substring(names(tcga.expr.sybl), 14, 15)!='11')]; ## remove non-tumor samples
geneset5 <- list()
for (i in all.test[[5]]$gene_name) {
  geneset5[[i]] <- as.matrix(tcga.expr.sybl)[i, expr.tumor]
  geneset5[[i]] <- geneset5[[i]][-pos]
}

expr.tumor<-substring(expr.tumor, 1, 12)[-pos];

gene.df5 <- list()
for (i in names(geneset5)) {
  gene.df5[[i]]<- data.frame(expr = geneset5[[i]], status=clin.data.PC[expr.tumor, 'status'], 
                             time=clin.data.PC[expr.tumor, 'Time'], Tumor=clin.data.PC[expr.tumor, 'type'])
  gene.df5[[i]] <- filter(gene.df5[[i]],is.na(Tumor)==F)
  gene.df5[[i]]$time <- as.numeric(gene.df5[[i]]$time);
}

gene.updn5 <- list()
for (i in names(gene.df5)) {
  med<-median(gene.df5[[i]]$expr)
  positive <- paste(i, "+", sep = '')
  negative <- paste(i, "-", sep = '')
  gene.updn5[[i]] <- ifelse(gene.df5[[i]]$expr>med, positive, negative)
  names(gene.updn5[[i]]) <- rownames(gene.df5[[i]])
  gene.df5[[i]]$gene <- unlist2(gene.updn5[[i]])[rownames(gene.df5[[i]])]
}

fit5 <- list()
pvalue5 <- list()
significant5 <- list()
for (i in names(gene.df5)) {
  fit5[[i]] <- survfit(Surv(gene.df5[[i]]$time/365, gene.df5[[i]]$status) ~ gene.df5[[i]]$Tumor + gene.df5[[i]]$gene, data = gene.df5[[i]])
  pvalue5[[i]] <- surv_pvalue(fit5[[i]], data = gene.df5[[i]], method = "survdiff")
  significant5[[i]] <- filter(pvalue5[[i]], pval <= 0.05)
}

gplot5 <- list()
for (i in c("B3GNTL1", "NLGN1", "PECR")) {
  fit5[[i]] <- survfit(Surv(time/365, status) ~ Tumor + gene, data = gene.df5[[i]])
  gplot5[[i]] <- ggsurvplot(fit5[[i]], gene.df5[[i]], palette = "jco", pval = TRUE)
}
res5 <- arrange_ggsurvplots(gplot5, print = FALSE, title = NA, ncol = 3, nrow = 3)
ggsave("induced_initial.up_down survival plots.pdf", res5, width = 17, height = 16)

# ----------------induced_initial.down_up--------------
expr.tumor<-names(tcga.expr.sybl)[which(substring(names(tcga.expr.sybl), 14, 15)!='11')]; ## remove non-tumor samples
geneset6 <- list()
for (i in all.test[[6]]$gene_name) {
  geneset6[[i]] <- as.matrix(tcga.expr.sybl)[i, expr.tumor]
  geneset6[[i]] <- geneset6[[i]][-pos]
}

expr.tumor<-substring(expr.tumor, 1, 12)[-pos];

gene.df6 <- list()
for (i in names(geneset6)) {
  gene.df6[[i]]<- data.frame(expr = geneset6[[i]], status=clin.data.PC[expr.tumor, 'status'], 
                             time=clin.data.PC[expr.tumor, 'Time'], Tumor=clin.data.PC[expr.tumor, 'type'])
  gene.df6[[i]] <- filter(gene.df6[[i]],is.na(Tumor)==F)
  gene.df6[[i]]$time <- as.numeric(gene.df6[[i]]$time);
}

gene.updn6 <- list()
for (i in names(gene.df6)) {
  med<-median(gene.df6[[i]]$expr)
  positive <- paste(i, "+", sep = '')
  negative <- paste(i, "-", sep = '')
  gene.updn6[[i]] <- ifelse(gene.df6[[i]]$expr>med, positive, negative)
  names(gene.updn6[[i]]) <- rownames(gene.df6[[i]])
  gene.df6[[i]]$gene <- unlist2(gene.updn6[[i]])[rownames(gene.df6[[i]])]
}

fit6 <- list()
pvalue6 <- list()
significant6 <- list()
for (i in names(gene.df6)) {
  fit6[[i]] <- survfit(Surv(gene.df6[[i]]$time/365, gene.df6[[i]]$status) ~ gene.df6[[i]]$Tumor + gene.df6[[i]]$gene, data = gene.df6[[i]])
  pvalue6[[i]] <- surv_pvalue(fit6[[i]], data = gene.df6[[i]], method = "survdiff")
  significant6[[i]] <- filter(pvalue6[[i]], pval <= 0.05)
}

gplot6 <- list()
for (i in c("ARFGEF1", "CBLN2", "CYP4B1", "DHFR", "HIBCH", "PEBP4", "PPP4R2", "THRB")) {
  fit6[[i]] <- survfit(Surv(time/365, status) ~ Tumor + gene, data = gene.df6[[i]])
  gplot6[[i]] <- ggsurvplot(fit6[[i]], gene.df6[[i]], palette = "jco", pval = TRUE)
}
res6 <- arrange_ggsurvplots(gplot6, print = FALSE, title = NA, ncol = 3, nrow = 3)
ggsave("induced_initial.down_up survival plots.pdf", res6, width = 17, height = 16)

# ----------------induced_initial.down_down--------------
expr.tumor<-names(tcga.expr.sybl)[which(substring(names(tcga.expr.sybl), 14, 15)!='11')]; ## remove non-tumor samples
geneset7 <- list()
for (i in all.test[[7]]$gene_name) {
  geneset7[[i]] <- as.matrix(tcga.expr.sybl)[i, expr.tumor]
  geneset7[[i]] <- geneset7[[i]][-pos]
}

expr.tumor<-substring(expr.tumor, 1, 12)[-pos];

gene.df7 <- list()
for (i in names(geneset7)) {
  gene.df7[[i]]<- data.frame(expr = geneset7[[i]], status=clin.data.PC[expr.tumor, 'status'], 
                             time=clin.data.PC[expr.tumor, 'Time'], Tumor=clin.data.PC[expr.tumor, 'type'])
  gene.df7[[i]] <- filter(gene.df7[[i]],is.na(Tumor)==F)
  gene.df7[[i]]$time <- as.numeric(gene.df7[[i]]$time);
}

gene.updn7 <- list()
for (i in names(gene.df7)) {
  med<-median(gene.df7[[i]]$expr)
  positive <- paste(i, "+", sep = '')
  negative <- paste(i, "-", sep = '')
  gene.updn7[[i]] <- ifelse(gene.df7[[i]]$expr>med, positive, negative)
  names(gene.updn7[[i]]) <- rownames(gene.df7[[i]])
  gene.df7[[i]]$gene <- unlist2(gene.updn7[[i]])[rownames(gene.df7[[i]])]
}

gene.df7 <- list.remove(gene.df7, c('CNBD1','DDX53',"MAGEA10"))

fit7 <- list()
pvalue7 <- list()
significant7 <- list()
for (i in names(gene.df7)) {
  fit7[[i]] <- survfit(Surv(time/365, status) ~ Tumor + gene, data = gene.df7[[i]])
  pvalue7[[i]] <- surv_pvalue(fit7[[i]], data = gene.df7[[i]], method = "survdiff")
  significant7[[i]] <- filter(pvalue7[[i]], pval <= 0.05)
}

gplot7 <- list()
for (i in c("DUSP8", "EPHA3", "H2AFJ", "HERC5", "KCNMA1", "KLF9", "LRBA", "MAP1LC3A", "NLGN1", "OPRK1", 
            "PKD2L2", "PRKG2", "PRPF39", "RALYL", "SMO", "TP53I11", "ZFP28", "ZNF814")) {
  fit7[[i]] <- survfit(Surv(time/365, status) ~ Tumor + gene, data = gene.df7[[i]])
  gplot7[[i]] <- ggsurvplot(fit7[[i]], gene.df7[[i]], palette = "jco", pval = TRUE)
}
res7 <- arrange_ggsurvplots(gplot7, print = FALSE, title = NA, ncol = 3, nrow = 3)
ggsave("induced_initial.down_down survival plots.pdf", res6, width = 17, height = 16)

