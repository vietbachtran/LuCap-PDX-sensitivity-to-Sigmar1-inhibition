##LuCap PDX transcriptomics data
rm(list=ls())
library(dplyr)
library(limma)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(annotate)
library(DESeq2)
library(edgeR)
library(enrichplot)

## Get RNAseq data from GEO
# Compile downloaded individual files from GEO
ls<-list.files("GSE199596_RAW/")

Exprs<-read.table(paste("GSE199596_RAW/",ls[1],sep=""), header=TRUE)
Exprs<-Exprs %>% dplyr::select(c(Entrez_ID,Gene_symbol,Frag_count))
colnames(Exprs)[3]<-sub("\\_.*","",ls[1])

for (i in ls[-1]) {
  temp<-read.table(paste("GSE199596_RAW/",i,sep=""), header=TRUE)
  temp<-temp %>% dplyr::select(c(Entrez_ID,Frag_count))
  colnames(temp)[2]<-sub("\\_.*","",i)
  Exprs <- merge(Exprs,temp, by='Entrez_ID')
}

# Get metadata
GEO <-"GSE199596"
gset <- getGEO(GEO, GSEMatrix =TRUE, AnnotGPL=FALSE)
pheno_1<-pData(phenoData(gset[[1]]))
pheno_2<-pData(phenoData(gset[[2]]))
pheno<-rbind(pheno_1,pheno_2)
pheno<-pheno %>% dplyr::select(c(title,geo_accession,source_name_ch1,`molecular phenotype:ch1`))

# Get samples treated with CT-1107
# Tumors 147CR, 96CR, 96CS were identified as unresponsive to the drug
# Tumors 23.1, 58, 77, 136 were identified as responsive to the drug
pheno<-filter(pheno,grepl(paste(c("147CR","58 ","136","96","77 ","23.1 "),collapse="|"), pheno$`source_name_ch1`))
pheno$treatment_1107<-ifelse(grepl(paste(c("147CR","96"),collapse="|"),pheno$source_name_ch1),"Unresponsive","Responsive")
Exprs<-Exprs[,c("Entrez_ID","Gene_symbol",pheno$geo_accession)]

write.table(Exprs, file= "Raw_read.tsv",row.names = FALSE, sep='\t')
write.table(pheno, file= "Metadata.tsv",row.names = FALSE, sep='\t')


## Differential Expression (DE) Analysis
# Import data
metadata<-read.table(file="Metadata.tsv",header=TRUE, sep='\t')
raw_count<-read.table("Raw_read.tsv", header=TRUE, sep='\t')

rownames(metadata)<- metadata$geo_accession
metadata$title=NULL

rownames(raw_count)<-raw_count$Entrez_ID
raw_count$Entrez_ID=NULL
raw_count$Gene_symbol=NULL

#Create a DESeqDataSet
all(colnames(raw_count) %in% rownames(metadata)) #Check row names and columns
all(colnames(raw_count) == rownames(metadata))

dds <- DESeqDataSetFromMatrix(countData = raw_count, 
                              colData = metadata, 
                              design = ~ treatment_1107)

dds$treatment_1107<-relevel(dds$treatment_1107, ref="Unresponsive")

#Remove lowly expressed genes
nrow(dds)
smallestGroupSize <- 7
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
nrow(dds)

#Normalize count with median of ratios method
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F, col.names=NA)

#Visualize the gene count normalization
boxplot(log10(counts(dds)), target="core", xlab="Sample", ylab="Raw gene counts (Log10 transformed)",las=2)
abline(h=median(log10(counts(dds))),col="blue")

boxplot(log10(normalized_counts), target="core",xlab="Sample", ylab="Normalized gene counts (Log10 transformed)",las=2)
abline(h=median(log10(normalized_counts)),col="blue")

#Variance stabilization
rld <- rlog(dds, blind = FALSE)

# Plot Principal Component Analysis (PCA)
plotPCA(rld, intgroup = "treatment_1107")

## DE Analysis
dds <- DESeq(dds)
DE_res <- results(dds)
DE_res<- DE_res %>% data.frame() %>% drop_na()

#Get Gene symbol
DE_res$SYMBOL<-AnnotationDbi::select(org.Hs.eg.db, keys=rownames(DE_res), columns='SYMBOL', keytype='ENTREZID')$SYMBOL
DE_res$Gene.ID<-rownames(DE_res)

#Volcano plot
library(RColorBrewer)
library(ggrepel)

Vol_df<-DE_res
Vol_df$DE<-"NO"
Vol_df$DE[Vol_df$padj<0.05 & Vol_df$log2FoldChange> log2(1.5)] <-"UP"
Vol_df$DE[Vol_df$padj<0.05 & Vol_df$log2FoldChange< -log2(1.5)] <-"DOWN"

#Plot volcanoe plot with significant genes
ggplot(data = Vol_df, aes(x = log2FoldChange, y = -log10(padj), col = DE)) +
  geom_hline(yintercept = -log10(0.05), col = "darkgray", linetype = 'dashed') +
  geom_point(size = 2,alpha = 0.5,shape = 16) +
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable
                     labels = c("Down regulated", "Not significant", "Up regulated"))+ # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  labs(color = NULL, x = expression("log"[2]*"FC"), y = expression("-log"[10]*"(adjusted P-value)"))

#Expression heatmap
library(ComplexHeatmap)

# Plot Expression Heatmap of significant genes (adjusted p-value < 0.05 and Fold Change > 1.5)
DE_genes<-DE_res %>% filter(padj<0.05 & abs(log2FoldChange)>log2(1.5))
TopDE <- counts(dds, normalized=TRUE)[DE_genes$Gene.ID,]

TopDE.z <- t(apply(TopDE, 1, scale)) #Normalize to Z-score

COLN<-as.character(dds$treatment_1107)
COLN<-ifelse(COLN=="Responsive",0,1)

#Color annotation
meta_anno<-metadata %>% dplyr::select(source_name_ch1,treatment_1107)
meta_anno$source_name_ch1<-substr(meta_anno$source_name_ch1,6,nchar(meta_anno$source_name_ch1)-6)
colnames(meta_anno)<-c("Tumor type","CT-1107 Treatment")
ha = HeatmapAnnotation(df= meta_anno, col = list(`CT-1107 Treatment`=c("Responsive"="green","Unresponsive"="red")))

png("DE_Heatmap.png",width = 3000, height = 3000,res = 300)
Heatmap(TopDE.z,cluster_rows = T,show_row_dend = FALSE, cluster_columns = T,show_column_dend = FALSE,show_column_names = FALSE,show_row_names = FALSE,column_title =NULL, column_split =COLN,cluster_column_slices = FALSE,
        bottom_annotation =ha, row_labels = rownames(TopDE.z),name = "Z-score", row_names_gp = gpar(fontsize = 3),use_raster=FALSE)
dev.off()


## Gene set enrichment analysis (GSEA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

#Rank genes by t stats
ordered_genes <- distinct(DE_res, Gene.ID, .keep_all = TRUE)$stat
names(ordered_genes) <- distinct(DE_res, Gene.ID, .keep_all = TRUE)$Gene.ID
ordered_genes <- sort(ordered_genes, decreasing = TRUE)

#GO pathways
gsea.go <- gseGO(geneList = ordered_genes, ont = "BP",OrgDb =org.Hs.eg.db,keyType= "ENTREZID",
                 pvalueCutoff = 0.05,pAdjustMethod = "BH",eps = 0)

enrichplot::dotplot(gsea.go, showCategory=30,orderBy = "NES",x=~NES) + ggtitle("GO pathways") +xlab("NES")

#KEGG pathways
gsea.kegg <- gseKEGG(geneList = ordered_genes,organism = 'hsa',keyType = "ncbi-geneid",
                     pvalueCutoff = 0.05,pAdjustMethod = "BH",eps=0)

enrichplot::dotplot(gsea.kegg, showCategory=30 ,x = ~NES) + ggtitle("KEGG pathways")+xlab("NES")

#Hallmark gene sets
msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::rename(ont = gs_name, gene = entrez_gene)%>% 
  distinct(ont,gene, .keep_all= TRUE)

msig_gsea <- GSEA(ordered_genes, TERM2GENE = msig_h,
                  eps= 0, pvalueCutoff = 1,pAdjustMethod = "BH",seed = 10)

enrichplot::dotplot(msig_gsea, showCategory=25,x = ~NES)+ggtitle("GSEA Hallmark Gene sets")+xlab("NES")

# Save enrichment plot
# Helper function
GSEA_enrichplot<-function(Pathway,GSEA_res){
  Pathway_title<-gsub(Pathway, pattern="_", replacement=" ") #Edit name
  Pathway_title<-gsub(Pathway_title, pattern="HALLMARK", replacement="")
  Path_index=which(GSEA_res@result$ID == Pathway)
  #Get pval and NES
  NES.label = paste("NES = ",signif(GSEA_res@result$NES[Path_index],digits=3))
  FDR.label = paste("p-value = ",signif(GSEA_res@result$p.adjust[Path_index],digits=1))
  label = paste(NES.label,FDR.label,sep = "\n")
  #Plot
  enrich_plot<-gseaplot2(GSEA_res, geneSetID = Pathway,title = Pathway_title,subplots = 1:2)
  enrich_plot[[1]]<-enrich_plot[[1]]+annotate("text",x=Inf, y = Inf, label= label,vjust=1, hjust=1, size=5)
  ggsave(paste0(Pathway,".png"), enrich_plot, width = 5, height = 5, dpi = 300)
}

GSEA_enrichplot("HALLMARK_DNA_REPAIR",msig_gsea)


## Gene set variation analysis (GSVA) Molecular signature scores
library(GSVA)
dds_norm <- vst(dds)
vst_df <- assay(dds_norm) %>% as.matrix()

#Create hallmark gene sets
hallmark_gene_sets <- msigdbr::msigdbr(species = "Homo sapiens",category = "H")
hallmarks_list <- split(hallmark_gene_sets$entrez_gene,hallmark_gene_sets$gs_name)

#Perform GSVA
gsvapar <- gsvaParam(vst_df,hallmarks_list,kcdf = "Gaussian",minSize = 5,maxSize = 1000,maxDiff = TRUE)
gsva_es <- gsva(gsvapar, verbose = FALSE)

# Clean up result data
gsva_es<-t(gsva_es) %>% as.data.frame()
gsva_es$geo_accession=rownames(gsva_es)
gsva_df<-merge(metadata,gsva_es,by='geo_accession')

# Plot GSVA score
# Helper function
Plot_GSVA <- function(Y_Value){
  Y_title <- gsub(Y_Value, pattern="HALLMARK_", replacement="")
  Y_title <- gsub(Y_title, pattern="_", replacement=" ")
  plot<-ggplot(gsva_df, aes(x = treatment_1107 , y = gsva_df[,Y_Value] , colour= treatment_1107)) +
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter() +
    xlab("Group")+
    ylab(Y_title)+
    guides(colour="none")+
    theme_minimal()
  ggsave(paste0(Y_Value,"_GSVA.png"), plot, width = 3, height = 3, dpi = 300)
}

Plot_GSVA("HALLMARK_DNA_REPAIR")


## Transcription factor activity inference
library(tibble)
library(tidyr)
library(stats)

# DeCoupleR package
library(decoupleR)

net <- decoupleR::get_collectri(organism='human', split_complexes=FALSE)

# Infer TF activity between treated vs untreated
DE_res<-DE_res[!duplicated(DE_res$SYMBOL),] %>% data.frame()
DE_res<-na.omit(DE_res)
rownames(DE_res)<-DE_res$SYMBOL

df_TF <- run_ulm(mat=DE_res[,'stat',drop=FALSE], net=net, .source='source', .target='target',.mor='mor', minsize = 5)

# Rank TF based on their score
f_contrast_acts <- df_TF %>% filter(p_value<0.05) %>%mutate(rnk = NA)
msk <- f_contrast_acts$score > 0
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
tfs <- f_contrast_acts %>%arrange(rnk) %>%pull(source)
f_contrast_acts <- f_contrast_acts %>%filter(source %in% tfs)
f_contrast_acts$Regulation<-ifelse(f_contrast_acts$score>0,"Up","Down")

# Plot top transcription factors
ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = Regulation), stat = "identity") +
  scale_fill_manual(values = c("#00AFBB","#bb0c00")) +
  guides(fill="none")+
  xlab("Transcription Factors")+
  ylab("Biological activity")

# Upstream regulation factor analysis
library(QuaternaryProd)

QP<-DE_res %>% select(Gene.ID,padj,log2FoldChange) %>% na.omit()
names(QP)<-c("entrez", "pvalue", "fc") 
QP$entrez<-as.character(QP$entrez)

QP_res <- RunCRE_HSAStringDB(QP, method = "Ternary", 
                                      fc.thresh = log2(1.5), pval.thresh = 0.05, 
                                      only.significant.pvalues = TRUE, 
                                      significance.level = 0.05, epsilon = 1e-16, 
                                      progressBar = TRUE, relations = NULL, entities = NULL)
QP_res<-QP_res %>% filter(symbol!="No-Symbol" & pvalue>0)

# Plot top 50 significant regulators
df<-head(QP_res,50)
df$RegulationCode=ifelse(df$regulation=="up",1,-1)

ggplot(df, aes(x =reorder(symbol, -log10(pvalue)*RegulationCode), y = -log10(pvalue)*RegulationCode)) +
  geom_bar(aes(fill = regulation), stat = "identity")+
  scale_fill_manual(values = c("#00AFBB","#bb0c00"),labels = c("Down Regulated","Up Regulated")) +
  labs(fill = NULL, x= "Regulators", y = expression("±log"[10]*"(P-value)"))