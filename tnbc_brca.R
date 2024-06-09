########################################################
######### TCGA BRCA TNBC ANALYSIS ######################

# ----- TNBC BRCA lncRNAs

# Loading required packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(ComplexHeatmap)
library(cowplot)
library(circlize)
library(RColorBrewer)


# Download GDC BRCA clinical data for TNBC samples
clinical_gdc <- read.delim("CollabGen/ClinicalDataGDC/nationwidechildrens.org_clinical_patient_brca.txt",
                           header = T, sep = "\t", dec = ".")
clinical_gdc <- clinical_gdc[3:1099,]

# Subset TNBC samples
TNBC_samples <- subset(clinical_gdc, er_status_by_ihc == "Negative" &
                         pr_status_by_ihc == "Negative" &
                         her2_status_by_ihc == "Negative")

tnbc_barcodes <- TNBC_samples$bcr_patient_barcode

# Download gene expression data from TCGAbiolinks
query = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor",
  barcode = tnbc_barcodes)

GDCdownload(query)
tcga_brca = GDCprepare(query)
dim(tcga_brca)

save(file = "CollabGen/tcga_BRCA.RData", tcga_brca)
load("CollabGen/tcga_BRCA.RData")

colnames(colData(tcga_brca))
gexp_tnbc <- assay(tcga_brca)
clin_tnbc <- colData(tcga_brca)

ens.row <- rownames(gexp_tnbc)
ens.row <- sub('\\.[0-9]*$', '', ens.row)
rownames(gexp_tnbc) <- ens.row

# Separating lncRNAs from the gene expression matrix from biomaRt
library(biomaRt)
mart <- useMart(biomart = "ensembl", 
                dataset = "hsapiens_gene_ensembl")
gene_names <- getBM(attributes = c("ensembl_gene_id",
                                   "hgnc_symbol","external_gene_name",
                                   "gene_biotype"), 
                    filters = "ensembl_gene_id", 
                    values = rownames(gexp_tnbc),
                    bmHeader = T, 
                    mart = mart)
lnc_names <- subset(gene_names, gene_names$`Gene type`=="lncRNA")
lnc_names$ens.symbol <- ifelse(lnc_names$`HGNC symbol`=="",
                               lnc_names$`Gene name`,
                               lnc_names$`HGNC symbol`)

lnc_names <- lnc_names[,c("Gene stable ID", "ens.symbol")]
lnc_names <- lnc_names[!duplicated(lnc_names$`Gene stable ID`),] #16775
rownames(lnc_names) <- lnc_names$`Gene stable ID`
colnames(lnc_names) <- c("Gene stable ID", "Gene name")

# Filtering gene expression matrix for lncRNAs
gexp_lnc_tnbc <- gexp_tnbc[lnc_names$`Gene stable ID`,] #16775

# Selecting primary tumor samples - 01A
patients_tnbc <- colnames(gexp_lnc_tnbc)
patients_tnbc <- patients_tnbc[grepl("01A", patients_tnbc)]
patients_tnbc <- substr(patients_tnbc, 1, 16)

# Checking if there are duplicated barcodes
patients_tnbc[duplicated(patients_tnbc)] #0
colnames(gexp_lnc_tnbc) <- substr(colnames(gexp_lnc_tnbc), 1, 16)
gexp_lnc_tnbc <- gexp_lnc_tnbc[,colnames(gexp_lnc_tnbc) %in% patients_tnbc] #114

rm(mart,gene_names, ens.row, tnbc_barcodes, TNBC_samples, clinical_gdc)

## Download Thorsson et al (2018) data set
repo_link <- "https://www.cell.com/cms/10.1016/j.immuni.2018.03.023/attachment/"
name.file <- "1b63d4bc-af31-4a23-99bb-7ca23c7b4e0a/mmc2.xlsx"
download.file(
  url= paste0(repo_link,name.file),
  destfile = "thorsson2018.xlsx")
thorsson <- data.frame(readxl::read_xlsx(path= "thorsson2018.xlsx"))

## Selecting BRCA patients
thorsson <- thorsson[thorsson$TCGA.Study=="BRCA",]
rownames(thorsson) <- thorsson$TCGA.Participant.Barcode

# Filtering Thorsson data set for primary samples
patients_tnbc <- substr(patients_tnbc,1,12)
thorsson <- thorsson[rownames(thorsson) %in% patients_tnbc,] #113

# Preparing Thorsson Data set

thorsson %>% group_by(Immune.Subtype) %>% count(n=n())
# C1: 35
# C2: 69
# C3: 6
# C4: 3

thorsson$Immune.Subtype <- factor(thorsson$Immune.Subtype,
                                  levels = c("C1","C2","C3","C4"))
thorsson <- thorsson[thorsson$TCGA.Subtype!="NA",]

save(thorsson, file = "thorsson_tnbc_brca_2703.RData")
load("thorsson_tnbc_brca_2703.RData")

# Filtering basal-like and normal-like samples
thorsson_basal <- subset(thorsson, TCGA.Subtype == "BRCA.Basal")

rm(thorsson,name.file, patients_tnbc, repo_link)

# filtering basal samples by clinical data
colnames(gexp_lnc_tnbc) <- substr(colnames(gexp_lnc_tnbc), 1, 12)
gexp_lnc_tnbc <- gexp_lnc_tnbc[,colnames(gexp_lnc_tnbc) %in% thorsson_basal$TCGA.Participant.Barcode]
# 89
thorsson_basal <- thorsson_basal[match(colnames(gexp_lnc_tnbc), rownames(thorsson_basal)),]
identical(thorsson_basal$TCGA.Participant.Barcode, colnames(gexp_lnc_tnbc))#TRUE

identical(colnames(gexp_lnc_tnbc), rownames(thorsson_basal)) #TRUE

rm(gexp_tnbc, clin_tnbc, tcga_brca)

#---------------------------------------------------------------------

# lncRNAs groups from table 2 (gene symbol)
lncRNAs_tnbc <- sort(c("OTUD6B-AS1", "LINC00578", "LINC01871",
                           "MAPT-IT1", "SLC26A4-AS1", "VPS9D1-AS1",
                           "PCAT18", "LINC01234", "SPATA41",
                           "LINC01215", "MIAT", "LINC01010", 
                           "LINC00668", "LINC02418"))

lnc_names <- lnc_names[rownames(lnc_names) %in% rownames(gexp_lnc_tnbc),]
lnc_names <- lnc_names[match(rownames(lnc_names), rownames(gexp_lnc_tnbc)),]
identical(rownames(gexp_lnc_tnbc), rownames(lnc_names)) #TRUE

lnc_names_tnbc <- lnc_names[lnc_names$`Gene name` %in% lncRNAs_tnbc,]
lnc_names_tnbc <- lnc_names_tnbc[,2, drop = F]

# Extracting lncs from gene expression matrix
dat_lnc_BRCA <- gexp_lnc_tnbc[rownames(gexp_lnc_tnbc) %in% rownames(lnc_names_tnbc),] 
# 14 89
identical(rownames(dat_lnc_BRCA), rownames(lnc_names_tnbc)) #TRUE
rownames(dat_lnc_BRCA) <- lnc_names_tnbc$`Gene name`

# Preparing dataframe for heatmap
unique(factor(thorsson_basal$TCGA.Subtype))
thorsson_basal <- thorsson_basal[thorsson_basal$TCGA.Subtype!="NA",]
dat_lnc_BRCA <- dat_lnc_BRCA[,colnames(dat_lnc_BRCA) %in% thorsson_basal$TCGA.Participant.Barcode]
# 14 102
thorsson_basal <- thorsson_basal[match(colnames(dat_lnc_BRCA), rownames(thorsson_basal)),]
identical(thorsson_basal$TCGA.Participant.Barcode, colnames(dat_lnc_BRCA)) #TRUE

# log2 + 1 
dat_lnc_BRCA <- log2(dat_lnc_BRCA + 1)

# Column-wise Z-score and setting maximum and minimum to +2 and -2 
dat_lnc_BRCA <- data.frame(t(dat_lnc_BRCA))
dat_lnc_BRCA <- data.frame(scale(dat_lnc_BRCA))
range(dat_lnc_BRCA) # -2.910468  6.668068
summary(dat_lnc_BRCA) # mean 0
dat_lnc_BRCA[dat_lnc_BRCA>2]<-2
dat_lnc_BRCA[dat_lnc_BRCA<(-2)]<-(-2)
range(dat_lnc_BRCA) # -2 2

# Generate top annotation dataframe
annot_lnc_BRCA <- thorsson_basal[,c("TCGA.Participant.Barcode", "TCGA.Subtype", "Immune.Subtype")]
annot_lnc_BRCA <- annot_lnc_BRCA[,-1]
annot_lnc_BRCA[,1:2] <- lapply(annot_lnc_BRCA[,1:2], as.factor)

identical(rownames(annot_lnc_BRCA), rownames(dat_lnc_BRCA)) #TRUE 

library(dplyr)
thorsson_basal %>% group_by(Immune.Subtype) %>% count(n=n())
which(thorsson_basal$Immune.Subtype == "C3") #86
sample_cols3 <- rownames(thorsson_basal)[86] # "TCGA-LL-A73Y"
which(rownames(dat_lnc_BRCA) == sample_cols3) #86
sample_cols3 <- dat_lnc_BRCA[86,]
cols3 <- rownames(sample_cols3)

# Clustering samples between Immune Subtypes
## Clustering C1 samples
cols <-hclust(dist(dat_lnc_BRCA[rownames(annot_lnc_BRCA)[annot_lnc_BRCA$Immune.Subtype=="C1"],],
                   method="euclidean"),
              method = "ward.D2")
cols <- cols$labels[rev(cols$order)]
## Clustering C2 samples
cols2 <-hclust(dist(dat_lnc_BRCA[rownames(annot_lnc_BRCA)[annot_lnc_BRCA$Immune.Subtype=="C2"],],
                    method="euclidean"),
               method = "ward.D2")
cols2 <- cols2$labels[rev(cols2$order)]
## Clustering C3 samples
# cols3 <-hclust(dist(dat_lnc_BRCA[rownames(annot_lnc_BRCA)[annot_lnc_BRCA$Immune.Subtype=="C3"],],
#                     method="euclidean"),
#                method = "ward.D2")
# cols3 <- cols3$labels[rev(cols3$order)]
## Clustering C4 samples
cols4 <-hclust(dist(dat_lnc_BRCA[rownames(annot_lnc_BRCA)[annot_lnc_BRCA$Immune.Subtype=="C4"],],
                    method="euclidean"),
               method = "ward.D2")
cols4 <- cols4$labels[rev(cols4$order)]

# Reorder data frames with clustered samples 
annot_lnc_BRCA<- annot_lnc_BRCA[c(cols,cols2, cols3, cols4),]
dat_lnc_BRCA <- dat_lnc_BRCA[rownames(annot_lnc_BRCA),]

#------------------------------------------------------
# Setting heatmap colors
x <- t(dat_lnc_BRCA)

cpall <- colorRampPalette(rev(brewer.pal(8, "RdBu")))(11)
bks <- as.numeric(x)
bks <- quantile(bks, probs = seq(0,1,0.1))
cpall <- circlize::colorRamp2(breaks = bks, cpall)


# Plotting Heatmap
h1<- Heatmap(matrix = t(dat_lnc_BRCA),
             name = "Log2 Gene\nexpression\nz-score",
             col = cpall,
             show_column_names = F, 
             show_column_dend = F,
             column_names_side="bottom",
             column_names_gp= gpar(col="white", fontsize=1),
             column_title_gp = gpar(fontsize=8), 
             column_title = c("C1\n(32)", "C2\n(53)",
                              "C3\n(1)", "C4\n(3)"),
             column_split = factor(annot_lnc_BRCA$Immune.Subtype,
                                   levels=c("C1", "C2",
                                            "C3", "C4" )),
             cluster_columns = F, 
             cluster_column_slices = F,
             cluster_rows = F, 
             row_title_rot = 0, 
             row_title_gp = gpar(fontsize=9),
             row_names_gp = gpar(fontsize=8),
             row_title_side ="left",
             height = unit(5, "cm"), width = unit(7, "cm"),
             raster_device = "png", use_raster = T, raster_quality = 6,
             heatmap_legend_param = list(
               labels_gp=gpar(fontsize=8),
               at=c(-2,0,2),
               title_gp=gpar(fontsize=8, fontface="bold"),
               direction = "horizontal",
               border=T,
               legend_width= unit(1.7,"cm")),
             top_annotation = columnAnnotation(
               df=annot_lnc_BRCA,
               simple_anno_size= unit(3, "mm"),
               height=unit(1,"cm"),
               annotation_name_gp = gpar(fontsize=8),
               annotation_name_side = "right",
               na_col="grey95", 
               annotation_legend_param = list(
                 # TCGA.Subtype=list(
                 #   labels_gp=gpar(fontsize=8),
                 #   title_gp=gpar(fontsize=8, fontface="bold"),
                 #   nrow=2, 
                 #   title="TCGA Subtype"),
                 Immune.Subtype=list(
                   labels_gp=gpar(fontsize=8),
                   title_gp=gpar(fontsize=8, fontface="bold"),
                   grid_height = unit(3, "mm"), grid_width = unit(3, "mm"),
                   nrow=2, 
                   title="Immune\nSubtype")),
               col=list(
                 Immune.Subtype=c("C1"="#FFFF99","C2"="#B15928","C3"="#B2DF8A",
                                  "C4"="#33A02C"),
                 TCGA.Subtype = c("BRCA.Basal" = "lightpink3"))))

g<- grid.grabExpr(
  draw(h1,
       heatmap_legend_side="bottom",
       annotation_legend_side="bottom",
       merge_legends=T),
  height = 5,
  width = 7)

plot_grid(g)

library(ggplot2)
ggsave(plot=g, device = "svg", width = 7, height = 5, 
       filename = "CollabGen/PlotsBRCA/heatmap_table2.svg", dpi=300)

# Clean environment
rm(annot_lnc_BRCA, annot_row_lnc, dat_lnc_BRCA, g, h1,
   cols, cols2, cols3, cols4, lnc_os_imm_inf, lnc_overall_surv,
   lncRNAs_BRCA, lncs, split, TNBC_samples, symbols, gexp_lnc, x, bks, cpall)

#--------------------------------------------

# lncRNAs from table 3 (ens)
lncRNAs_other_cancers <- c("ENSG00000266445", "ENSG00000231439",
                           "ENSG00000228315", "ENSG00000182057",
                           "ENSG00000152931", "ENSG00000265688",
                           "ENSG00000245937", "ENSG00000214049",
                           "ENSG00000225733", "ENSG00000253563",
                           "ENSG00000206337")

gexp_lnc_others <- gexp_lnc_tnbc[rownames(gexp_lnc_tnbc) %in% lncRNAs_other_cancers,]

# filtering samples basal and normal by Thorsson
gexp_lnc_others <- gexp_lnc_others[,colnames(gexp_lnc_others) %in% thorsson_basal$TCGA.Participant.Barcode]
thorsson_basal <- thorsson_basal[match(colnames(gexp_lnc_others), rownames(thorsson_basal)),]
identical(thorsson_basal$TCGA.Participant.Barcode, colnames(gexp_lnc_others))#TRUE

identical(colnames(gexp_lnc_others), rownames(thorsson_basal))

lnc_others_names <- lnc_names[rownames(lnc_names) %in% lncRNAs_other_cancers,]
lnc_others_names <- lnc_others_names[,2, drop = F]

# Extracting lncs from gene expression matrix
dat_lnc_others <- gexp_lnc_others[rownames(gexp_lnc_others) %in% rownames(lnc_others_names),] 
# 11 89
identical(rownames(dat_lnc_others), rownames(lnc_others_names)) #TRUE
rownames(dat_lnc_others) <- lnc_others_names$`Gene name`

# Preparing dataframe for heatmap
unique(factor(thorsson_basal$TCGA.Subtype))
thorsson_basal <- thorsson_basal[thorsson_basal$TCGA.Subtype!="NA",]
dat_lnc_others <- dat_lnc_others[,colnames(dat_lnc_others) %in% thorsson_basal$TCGA.Participant.Barcode]
# 11 102
thorsson_basal <- thorsson_basal[match(colnames(dat_lnc_others), rownames(thorsson_basal)),]
identical(thorsson_basal$TCGA.Participant.Barcode, colnames(dat_lnc_others))

# log2 + 1
dat_lnc_others <- log2(dat_lnc_others + 1)

# Column-wise Z-score and setting maximum and minimum to +2 and -2 
dat_lnc_others <- data.frame(t(dat_lnc_others))
dat_lnc_others <- data.frame(scale(dat_lnc_others))
range(dat_lnc_others) # -3.674187  5.295571
summary(dat_lnc_others) #mean 0
dat_lnc_others[dat_lnc_others>2]<-2
dat_lnc_others[dat_lnc_others<(-2)]<-(-2)
range(dat_lnc_others) # -2 2


# Generating top annotation dataframe
annot_others <- thorsson_basal[,c("TCGA.Participant.Barcode", "TCGA.Subtype", "Immune.Subtype")]
annot_others <- annot_others[,-1]
annot_others[,1:2] <- lapply(annot_others[,1:2], as.factor)


identical(rownames(annot_others), rownames(dat_lnc_others)) #TRUE 

thorsson_basal %>% group_by(Immune.Subtype) %>% count(n=n())
which(thorsson_basal$Immune.Subtype == "C3") #86
sample_cols3 <- rownames(thorsson_basal)[86] # "TCGA-LL-A73Y"
which(rownames(dat_lnc_others) == sample_cols3) #86
sample_cols3 <- dat_lnc_others[86,]
cols3 <- rownames(sample_cols3)

# Clustering samples between Immune Subtypes
## Clustering C1 samples
cols <-hclust(dist(dat_lnc_others[rownames(annot_others)[annot_others$Immune.Subtype=="C1"],],
                   method="euclidean"),
              method = "ward.D2")
cols <- cols$labels[rev(cols$order)]
## Clustering C2 samples
cols2 <-hclust(dist(dat_lnc_others[rownames(annot_others)[annot_others$Immune.Subtype=="C2"],],
                    method="euclidean"),
               method = "ward.D2")
cols2 <- cols2$labels[rev(cols2$order)]
# ## Clustering C3 samples
# cols3 <-hclust(dist(dat_lnc_others[rownames(annot_others)[annot_others$Immune.Subtype=="C3"],],
#                     method="euclidean"),
#                method = "ward.D2")
# cols3 <- cols3$labels[rev(cols3$order)]
## Clustering C4 samples
cols4 <-hclust(dist(dat_lnc_others[rownames(annot_others)[annot_others$Immune.Subtype=="C4"],],
                    method="euclidean"),
               method = "ward.D2")
cols4 <- cols4$labels[rev(cols4$order)]

# Reorder dataframes with clustered samples
annot_others<- annot_others[c(cols,cols2,cols3,cols4),]
dat_lnc_others <- dat_lnc_others[rownames(annot_others),]

# Setting heatmap colors
x <- t(dat_lnc_others)

cpall <- colorRampPalette(rev(brewer.pal(8, "RdBu")))(11)
bks <- as.numeric(x)
bks <- quantile(bks, probs = seq(0,1,0.1))
cpall <- circlize::colorRamp2(breaks = bks, cpall)


# Plotting Heatmap
h1<- Heatmap(matrix = t(dat_lnc_others),
             name = "Log2 Gene\nexpression\nz-score",
             col = cpall,
             show_column_names = F, 
             show_column_dend = F,
             column_names_side="bottom",
             column_names_gp= gpar(col="white", fontsize=1),
             column_title_gp = gpar(fontsize=8), 
             column_title = c("C1\n(32)", "C2\n(53)",
                              "C3\n(1)", "C4\n(3)"),
             column_split = factor(annot_others$Immune.Subtype,
                                   levels=c("C1", "C2", "C3", "C4")),
             cluster_columns = F, 
             cluster_column_slices = F,
             cluster_rows = F, 
             row_title_rot = 0, 
             row_title_gp = gpar(fontsize=9),
             row_names_gp = gpar(fontsize=8),
             row_title_side ="left",
             height = unit(5, "cm"), width = unit(7, "cm"),
             raster_device = "png", use_raster = T, raster_quality = 6,
             heatmap_legend_param = list(
               labels_gp=gpar(fontsize=8),
               at=c(-2,0,2),
               title_gp=gpar(fontsize=8, fontface="bold"),
               direction = "horizontal",
               border=T,
               legend_width= unit(1.7,"cm")),
             top_annotation = columnAnnotation(
               df=annot_others,
               simple_anno_size= unit(3, "mm"),
               height=unit(1,"cm"),
               annotation_name_gp = gpar(fontsize=8),
               annotation_name_side = "right",
               na_col="white", 
               annotation_legend_param = list(
                 # TCGA.Subtype=list(
                 #   labels_gp=gpar(fontsize=8),
                 #   title_gp=gpar(fontsize=8, fontface="bold"),
                 #   nrow=2, 
                 #   title="TCGA Subtype"),
                 Immune.Subtype=list(
                   labels_gp=gpar(fontsize=8),
                   title_gp=gpar(fontsize=8, fontface="bold"),
                   grid_height = unit(3, "mm"), grid_width = unit(3, "mm"),
                   nrow=2, 
                   title="Immune\nSubtype")),
               col=list(
                 Immune.Subtype=c("C1"="#FFFF99","C2"="#B15928","C3"="#B2DF8A",
                                  "C4"="#33A02C"),
                 TCGA.Subtype=c("BRCA.Basal"="lightpink3"))))

g<- grid.grabExpr(
  draw(h1,
       heatmap_legend_side="bottom",
       annotation_legend_side="bottom",
       merge_legends=T),
  height = 5,
  width = 7)

plot_grid(g)
ggsave(plot=g, device = "svg", width = 7, height = 5, 
       filename = "CollabGen/PlotsBRCA/heatmap_table3.svg", dpi=300)

rm(g, h1, lnc_names_tnbc, lnc_others_names, sample_cols3, x, bks,
   cols, cols2, cols3, cols4, lncRNAs_other_cancers, lncRNAs_tnbc)

#-------------------- DESEQ2 ANALYSIS

thorsson_c1_c2 <- subset(thorsson_basal, Immune.Subtype == "C1" | Immune.Subtype == "C2")
#85

# filtering basal samples by clinical data
colnames(gexp_lnc_tnbc) <- substr(colnames(gexp_lnc_tnbc), 1, 12)
gexp_lnc_tnbc <- gexp_lnc_tnbc[,colnames(gexp_lnc_tnbc) %in% thorsson_c1_c2$TCGA.Participant.Barcode]
# 89
thorsson_c1_c2 <- thorsson_c1_c2[match(colnames(gexp_lnc_tnbc), rownames(thorsson_c1_c2)),]
identical(thorsson_c1_c2$TCGA.Participant.Barcode, colnames(gexp_lnc_tnbc))#TRUE

identical(colnames(gexp_lnc_tnbc), rownames(thorsson_c1_c2)) #TRUE

rm(gexp_tnbc, clin_tnbc, tcga_brca)


dds <- DESeqDataSetFromMatrix(countData = gexp_lnc_tnbc, colData = thorsson_c1_c2, 
                              design = ~ Immune.Subtype)


dds <- dds[rowSums(counts(dds)) >= 10,]
dds <- DESeq(dds)
resultsNames(dds) # "Immune.Subtype_C2_vs_C1"
res <- results(dds, contrast = c("Immune.Subtype", "C1", "C2"))
res
res_padj <- res[order(res$padj),] 

summary(results(dds, alpha=0.05))
# LFC > 0: 219
# LFC < 0: 143
# OUTLIERS: 0
# LOW COUNTS: 4241

save(dds, file = "CollabGen/dds_TNBC_IS_3103.RData")
load("CollabGen/dds_TNBC_IS_3103.RData")

lnc_deseq_names <- lnc_names[rownames(lnc_names) %in% rownames(res),]
lnc_deseq_names <- lnc_deseq_names[match(rownames(lnc_deseq_names), rownames(res)),]
res$SYMBOL <- lnc_deseq_names$`Gene name`
res_ordered <- res[order(res$padj),]

write.csv(res_ordered, 
          file="CollabGen//basal_tnbc_IS_symbol_3103.csv")

vsd <- vst(dds, blind = F)

save(file = "CollabGen/vsd_basal_3103.RData", vsd)
load("CollabGen/vsd_basal_3103.RData")

res_sig_tnbc <- subset(res, padj < 0.05 & abs(log2FoldChange) > 2) #79 ens

mat1 <- assay(vsd)
idx <- rownames(res_sig_tnbc)
DEgenes <- mat1[idx,]

annotation1 <- as.data.frame(colData(vsd)[, c("TCGA.Participant.Barcode","Immune.Subtype")])
annotation1 <- annotation1[,-1,drop=F]
identical(rownames(annotation1), colnames(DEgenes))

#-- plot heatmap of top DE genes
ordered_sig_genes <- res_sig_tnbc[order(res_sig_tnbc$padj),]
id1 <- rownames(ordered_sig_genes)
id2 <- ordered_sig_genes$SYMBOL
topDE <- mat1[id1,]
rownames(topDE) <- id1
#rownames(topDE) <- id1 #id2
#top50DE <- head(topDE, n = 50)

# log2 + 1 
topDE_log2 <- log2(topDE)

# Column-wise Z-score and setting maximum and minimum to +2 and -2 
topDE_log2 <- data.frame(t(topDE_log2))
topDE_log2 <- data.frame(scale(topDE_log2))
range(topDE_log2) # -4.509081  3.950693
summary(topDE_log2) # mean 0
topDE_log2[topDE_log2>2]<-2
topDE_log2[topDE_log2<(-2)]<-(-2)
range(topDE_log2) # -2 2

annotation1 %>% group_by(Immune.Subtype) %>% count(n=n())
#C1 32
#C2 53


# Clustering samples between Immune Subtypes
## Clustering C1 samples
cols <-hclust(dist(topDE_log2[rownames(annotation1)[annotation1$Immune.Subtype=="C1"],],
                   method="euclidean"),
              method = "ward.D2")
cols <- cols$labels[rev(cols$order)]
## Clustering C2 samples
cols2 <-hclust(dist(topDE_log2[rownames(annotation1)[annotation1$Immune.Subtype=="C2"],],
                    method="euclidean"),
               method = "ward.D2")
cols2 <- cols2$labels[rev(cols2$order)]

# Reorder dataframes with clustered samples
annotation1<- annotation1[c(cols,cols2),, drop=F]
topDE_log2 <- topDE_log2[rownames(annotation1),]

# Clustering rows
mat_final <- topDE_log2
mat_final <- t(mat_final)
mat_final <- hclust(dist(mat_final, method = "euclidean"), method = "ward.D2")
mat_final <- mat_final$labels[rev(mat_final$order)]

topDE_log2 <- topDE_log2[,mat_final]

# Setting heatmap colors
x <- t(topDE_log2)

cpall <- colorRampPalette(rev(brewer.pal(8, "RdBu")))(11)
bks <- as.numeric(x)
bks <- quantile(bks, probs = seq(0,1,0.1))
cpall <- circlize::colorRamp2(breaks = bks, cpall)

mat_h <- t(topDE_log2)

write.csv(file = "top_DE_c1c2.csv", topDE_log2)

# Plotting Heatmap
h1<- Heatmap(matrix = mat_h[1:20,],
             name = "Log2 Gene\nexpression\nz-score",
             col = cpall,
             show_column_names = F,
             show_row_names = T,
             show_column_dend = F,
             show_row_dend = F,
             column_names_side="bottom",
             column_names_gp= gpar(col="white", fontsize=1),
             column_title_gp = gpar(fontsize=10), 
             column_title = c("C1\n(n=32)", "C2\n(n=53)"),
             column_split = factor(annotation1$Immune.Subtype,
                                   levels=c("C1", "C2")),
             cluster_columns = F, 
             cluster_column_slices = F,
             cluster_rows = F, 
             height = unit(5, "cm"), width = unit(7, "cm"),
             raster_device = "png", use_raster = T, raster_quality = 6,
             row_title_rot = 0, 
             row_title_gp = gpar(fontsize=9),
             row_names_gp = gpar(fontsize=7.5),
             row_title_side ="left",
             heatmap_legend_param = list(
               labels_gp=gpar(fontsize=8),
               at=c(-2,0,2),
               title_gp=gpar(fontsize=8, fontface="bold"),
               direction = "horizontal",
               border=T,
               legend_width= unit(1.7,"cm")),
             top_annotation = columnAnnotation(
               df=annotation1,
               simple_anno_size= unit(3, "mm"),
               height=unit(1,"cm"),
               annotation_name_gp = gpar(fontsize=8),
               annotation_name_side = "right",
               na_col="grey95", 
               annotation_legend_param = list(
                 Immune.Subtype=list(
                   labels_gp=gpar(fontsize=8),
                   title_gp=gpar(fontsize=8, fontface="bold"),
                   nrow=2, 
                   title="Immune Subtype")),
               col=list(
                 Immune.Subtype=c("C1"="#FFFF99","C2"="#B15928"))))    

g<- grid.grabExpr(
  draw(h1,
       heatmap_legend_side="bottom",
       annotation_legend_side="bottom",
       merge_legends=T))
  # height = 5,
  # width = 7)

plot_grid(g)

ggsave(plot=g, device = "svg", width = 7, height = 5, filename = "CollabGen/PlotsBRCA/heatmap_c1_c2_20.svg", dpi=300)
