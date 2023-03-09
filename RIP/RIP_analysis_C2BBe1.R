library(tximport)
library(DESeq2)
library(org.Hs.eg.db)
library(ggfortify)
library(dplyr)
library(ggrepel)
library(enrichplot)
library(clusterProfiler)

rm(list=ls())
setwd("/Users/gilbertgiri/Desktop/Dominguez_Lab/PCBP1/RIPs/C2BBe1/")

files <- list.files(pattern = "\\.genes.results$")
Counts<-list()

#Reading the files
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
txi.rsem$length[txi.rsem$length == 0] <- 1

##Design for DESeq2
sample<-sapply(files,function(x){strsplit(x,"[.]")[[1]][1]})
assay<-sapply(sample,function(x){strsplit(x,"_|-")[[1]][3]})
condition<-sapply(sample,function(x){strsplit(x,"_|-")[[1]][2]})
sampleInfo<-as.data.frame(cbind(sample,condition,assay))
sampleInfo

#Reading into DESEQ2
deseqdata <- DESeqDataSetFromTximport(txi.rsem, colData = sampleInfo, design = ~assay+condition+assay:condition)
deseqdata

# pre-filtering- Averaging at 5 per sample
keep <- rowSums(counts(deseqdata)) >= 60
deseqdata <- deseqdata[keep,]

#checking the factor levels
deseqdata$condition<-relevel(deseqdata$condition,ref = "gfp")
deseqdata$assay<-relevel(deseqdata$assay,ref="in")

###Running DESEQ
deseqdata<-estimateSizeFactors(deseqdata)
deseqdata<-estimateDispersions(deseqdata)
deseqdata<-nbinomLRT(deseqdata,reduced=~assay+condition)

###Looking at the variability between replicates
par(mfrow=c(1,1))
df <- as.data.frame(t(counts(deseqdata,normalize=T)))
df$Sample<-paste0(sampleInfo$condition,sampleInfo$assay)
pca_res <- prcomp(df[1:dim(df)[2]-1])
pca_res<-as.list(pca_res)
autoplot(pca_res,data=df,colour="Sample")

### Comparing fold change in mutant vs wildtype
resultsNames(deseqdata)
res<-results(deseqdata,contrast = list(c("assayflagip.conditionwt","assayflagip.conditionmut")))
res
summary(res)

#Plotting
res_df<-as.data.frame(res)
res_df$Genes<-rownames(res_df)


####Volcano Plot
res_df <- res_df %>% mutate(gene_type = case_when(log2FoldChange >= 1 & padj <= 0.05 ~ "up",
                                                  log2FoldChange <= -1 & padj <= 0.05 ~ "down",
                                                  TRUE ~ "ns"))


cols <- c("up" = "#FE0000", "down" = "#079BF5", "ns" = "grey") 
sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)


sig_up_genes <- res_df %>%
  dplyr::filter(Genes %in% c())



volcano<-res_df %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             fill = gene_type,    
             size = gene_type,
             alpha =gene_type)) + 
  geom_point(shape = 21,
             colour = "black") +
  geom_point(data = sig_up_genes,
             shape = 21,
             alpha = 0.05,
             size = 2, 
             fill = "firebrick", 
             colour = "black")+
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(log2(0.5), log2(2)),
             linetype = "dashed")+
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-20, 20, 5)),       
                     limits = c(-20, 20))+
  theme_classic()+
  geom_text_repel(data = sig_up_genes, # Add labels last to appear as the top layer  
                  aes(label = Genes),
                  force = 2,
                  nudge_y = 1)+
  theme(legend.position = "none",
        axis.title.y = element_text(size=10),
        axis.text.x=element_text(size=10,color="black"),
        axis.text.y=element_text(size=10,color="black")
  )

volcano

#GO Analysis
UP<-res_df[which(res_df$gene_type=="up"),]
Down<-res_df[which(res_df$gene_type=="down"),]

Go_up <-enrichGO(gene          = UP$Genes,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 universe     = res_df$Genes,
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05)


Go_down<-enrichGO(gene          = Down$Genes,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = "BP",
                  universe     = res_df$Genes,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05)


dotplot(Go_up, showCategory=25,orderBy="GeneRatio",font.size = 10)+
  theme(axis.title.y = element_text(size=10),
        axis.text.x=element_text(size=10,color="black"),
        axis.text.y=element_text(size=10,color="black")
  )


dotplot(Go_down, showCategory=25,orderBy="GeneRatio",font.size = 10)+
  theme(axis.title.y = element_text(size=10),
        axis.text.x=element_text(size=10,color="black"),
        axis.text.y=element_text(size=10,color="black")
  )

res_df$Gene_Symbol = mapIds(org.Hs.eg.db,
                              keys = res_df$Genes,
                              column = "SYMBOL",
                              keytype="ENSEMBL",
                              multiVals = "first")

write.table(res_df,"../RIP_DESeq_summary_C2BBe1.txt",quote=F,row.names = F,col.names = T,sep="\t")
