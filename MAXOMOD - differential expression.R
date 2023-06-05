suppressMessages( library(ggplot2))
suppressMessages( library(dplyr))
suppressMessages( library(readxl))
suppressMessages( library(DESeq2))
suppressMessages( library(writexl))
suppressMessages( library(vsn))
suppressMessages( library(ggrepel))
suppressMessages( library(outliers))
suppressMessages( library(flashClust))
BiocManager::install("sva")
BiocManager::install("GOSim")
devtools::install_github('juanbot/CoExpNets')
suppressMessages( library(CoExpNets))

########### Download data ###########

# load all species 
load("~/Documents/HMGU/MAXOMOD/data to be used/exceRpt_smallRNAQuants_ReadCounts.RData")

# -> piRNA
counts_piRNA <- exprs.piRNA

# -> circRNA
counts_circRNA <- exprs.circRNA

# -> miRNA
counts_miRNA <- exprs.miRNA

# -> tRNA
counts_tRNA <- exprs.tRNA

# -> all classes of smallRNA
counts_smallRNA <- do.call("rbind",list(counts_piRNA,counts_tRNA,counts_miRNA,counts_circRNA))

# Clinical information 
clinical_info <- read_excel("~/Documents/HMGU/MAXOMOD/data to be used/Samples that passed QC smallRNA_clean cohort.xlsx")

########### Create column data ###########
col_data <- clinical_info %>%
  select(`SAMPLES MAXOMOD (Munich)`,condition,sex) %>%
  as.data.frame()
rownames(col_data) <- col_data[,1]
col_data <- col_data[,-1]
col_data_male <- col_data %>%
  filter(sex == "male")
col_data_female <- col_data %>%
  filter(sex == "female")

########### Make sure count matrix has same IDs and in the same order as column data ###########
counts_smallRNA <- counts_smallRNA[,rownames(col_data)]
counts_smallRNA_male <- counts_smallRNA[,rownames(col_data_male)]
counts_smallRNA_female <- counts_smallRNA[,rownames(col_data_female)]
counts_miRNA <- counts_miRNA[,rownames(col_data)]
counts_miRNA_male <- counts_miRNA[,rownames(col_data_male)]
counts_miRNA_female <- counts_miRNA[,rownames(col_data_female)]

########### DESeq2 Data Matrix ###########
# -> all smallRNA species
dds_smallRNA <- DESeqDataSetFromMatrix(countData = round(counts_smallRNA),
                                        colData = col_data,
                                        design = ~ condition)
dds_smallRNA_sex <- DESeqDataSetFromMatrix(countData = round(counts_smallRNA), # sex as covariate 
                                           colData = col_data,
                                           design = ~ condition + sex)
dds_smallRNA_female <- DESeqDataSetFromMatrix(countData = round(counts_smallRNA_female),
                                              colData = col_data_female,
                                              design = ~ condition)
dds_smallRNA_male <- DESeqDataSetFromMatrix(countData = round(counts_smallRNA_male),
                                              colData = col_data_male,
                                              design = ~ condition)

# -> miRNA
dds_miRNA <- DESeqDataSetFromMatrix(countData = round(counts_miRNA),
                                       colData = col_data,
                                       design = ~ condition)
dds_miRNA_sex <- DESeqDataSetFromMatrix(countData = round(counts_miRNA), # sex as covariate
                                           colData = col_data,
                                           design = ~ condition + sex)
dds_miRNA_female <- DESeqDataSetFromMatrix(countData = round(counts_miRNA_female),
                                              colData = col_data_female,
                                              design = ~ condition)
dds_miRNA_male <- DESeqDataSetFromMatrix(countData = round(counts_miRNA_male),
                                            colData = col_data_male,
                                            design = ~ condition)

########### Pre-filter: eliminate omics with single digit counts ###########
# -> smallRNA
count_keep <- rowSums(counts(dds_smallRNA)) >= 10
dds_smallRNA <- dds_smallRNA[count_keep]
dds_smallRNA_sex <- dds_smallRNA_sex[count_keep]
count_keep <- rowSums(counts(dds_smallRNA_female)) >= 10
dds_smallRNA_female <- dds_smallRNA_female[count_keep]
count_keep <- rowSums(counts(dds_smallRNA_male)) >= 10
dds_smallRNA_male <- dds_smallRNA_male[count_keep]

# -> miRNA
count_keep <- rowSums(counts(dds_miRNA)) >= 10
dds_miRNA <- dds_miRNA[count_keep]
dds_miRNA_sex <- dds_miRNA_sex[count_keep]
count_keep <- rowSums(counts(dds_miRNA_female)) >= 10
dds_miRNA_female <- dds_miRNA_female[count_keep]
count_keep <- rowSums(counts(dds_miRNA_male)) >= 10
dds_miRNA_male <- dds_miRNA_male[count_keep]

########### Filtering based on quantile info (TCGA package) ###########
# dds <- estimateSizeFactors(dds)
# mRNA.count.data.new.temp <- TCGAanalyze_Filtering(assay(dds), method = "quantile")
# mRNA.count.data.new <- DESeqDataSetFromMatrix(countData = mRNA.count.data.new.temp,
#                                               colData = col.count.data,
#                                               design = ~ condition)

########### Define factor levels ALS vs Ctrl ###########
# -> smallRNA
dds_smallRNA$condition <- factor(dds_smallRNA$condition, levels = c("ctrl","als"))
dds_smallRNA_sex$condition <- factor(dds_smallRNA_sex$condition, levels = c("ctrl","als"))
dds_smallRNA_female$condition <- factor(dds_smallRNA_female$condition, levels = c("ctrl","als"))
dds_smallRNA_male$condition <- factor(dds_smallRNA_male$condition, levels = c("ctrl","als"))

# -> miRNA
dds_miRNA$condition <- factor(dds_miRNA$condition, levels = c("ctrl","als"))
dds_miRNA_sex$condition <- factor(dds_miRNA_sex$condition, levels = c("ctrl","als"))
dds_miRNA_female$condition <- factor(dds_miRNA_female$condition, levels = c("ctrl","als"))
dds_miRNA_male$condition <- factor(dds_miRNA_male$condition, levels = c("ctrl","als"))

########### Count transformation of original count data for visualization purposes: PCA, variance dist, heatmaps, etc. ###########
# -> smallRNA
rlog_smallRNA <- rlog(dds_smallRNA)
rlog_smallRNA_female <- rlog(dds_smallRNA_female)
rlog_smallRNA_male <- rlog(dds_smallRNA_male)

# -> miRNA
rlog_miRNA <- rlog(dds_miRNA)
rlog_miRNA_female <- rlog(dds_miRNA_female)
rlog_miRNA_male <- rlog(dds_miRNA_male)

########### Visualization 1: Standard deviation distribution ###########
cidr <- getwd()
dir.create(file.path(cidr,"plots_output"),recursive = T)
# -> smallRNA
pdf("plots_output/Mean variance smallRNA.pdf")
meanSdPlot(assay(rlog_smallRNA), xlab = "Mean normalized smallRNA") 
dev.off()
pdf("plots_output/Mean variance smallRNA_female.pdf")
meanSdPlot(assay(rlog_smallRNA_female), xlab = "Mean normalized smallRNA") 
dev.off()
pdf("plots_output/Mean variance smallRNA_male.pdf")
meanSdPlot(assay(rlog_smallRNA_male), xlab = "Mean normalized smallRNA") 
dev.off()

# -> miRNA
pdf("plots_output/Mean variance miRNA.pdf")
meanSdPlot(assay(rlog_miRNA), xlab = "Mean normalized miRNA") 
dev.off()
pdf("plots_output/Mean variance miRNA_female.pdf")
meanSdPlot(assay(rlog_miRNA_female), xlab = "Mean normalized miRNA") 
dev.off()
pdf("plots_output/Mean variance miRNA_male.pdf")
meanSdPlot(assay(rlog_miRNA_male), xlab = "Mean normalized miRNA") 
dev.off()

########### Visualization 2: PCA ###########
plot_PCA <- function(rlog_data,title_plot){
  PCA_data <- plotPCA(rlog_data, returnData = T)
  percentVar <- round(100 * attr(PCA_data, "percentVar"))
  ggplot(PCA_data, aes(x=PC1,y=PC2,color = condition,size=2)) + 
    geom_point(alpha = 0.7) +  
    geom_text_repel(aes(label = name)) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    guides(size = FALSE) +
    scale_color_manual(values = c("#81a25d","#7E5DA2"),breaks = c("ctrl", "als"), labels = c("Control", "ALS")) + 
    ggtitle(title_plot) + theme_minimal() +
    theme(legend.position = "bottom",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          plot.title = element_text(size = 15)) 
}

# -> smallRNA
pdf("plots_output/PCA smallRNA.pdf")
plot_PCA(rlog_smallRNA,"PCA smallRNA")
dev.off()
pdf("plots_output/PCA smallRNA female.pdf")
plot_PCA(rlog_smallRNA_female,"PCA smallRNA Female")
dev.off()
pdf("plots_output/PCA smallRNA male.pdf")
plot_PCA(rlog_smallRNA_male,"PCA smallRNA Male")
dev.off()

# -> miRNA
pdf("plots_output/PCA miRNA.pdf")
plot_PCA(rlog_miRNA,"PCA miRNA")
dev.off()
pdf("plots_output/PCA miRNA female.pdf")
plot_PCA(rlog_miRNA_female,"PCA miRNA Female")
dev.off()
pdf("plots_output/PCA miRNA male.pdf")
plot_PCA(rlog_miRNA_male,"PCA miRNA Male")
dev.off()

########### Outlier detection ###########
# -> based on PCAs we have 3 potential outliers to investigate: G103, MM77 and CC55 in both smallRNA and miRNA data
outlier_detection <- function(rlog_data){
  data_outlier <- assay(rlog_data)
  outlist = NULL
  #Traspose data frame and get a distance matrix between pairs of samples
  sampdist = as.matrix(dist(t(data_outlier)))
  result = grubbs.test(apply(sampdist,1,sum))
  #Is the test significant at 95% confidence level?
  if(result$p.value < 0.05){
    cat("We should drop out",names(result$p.value),"\n")
    plot(flashClust(dist(CoExpNets::trasposeDataFrame(data_outlier,F)), method = "average"),cex=0.5,
         main=paste0("Detected outlier ",names(result$p.value)," in ROS/MAP cases + ctrls"))
  }else
    cat("No plausible outlier detected\n")
}

# -> smallRNA
outlier_detection(rlog_smallRNA) # G103 identified as outlier with 5% significance 
outlier_detection(rlog_smallRNA_male) # no outlier detected
outlier_detection(rlog_smallRNA_female) # no outlier detected

# -> miRNA
outlier_detection(rlog_miRNA) # G103 identified as outlier with 5% significance 
outlier_detection(rlog_miRNA_female) # no outlier detected
outlier_detection(rlog_miRNA_male) # no outlier detected

########### Differential expression analysis (with outlier) ###########
# -> smallRNA 
DE_smallRNA <- DESeq(dds_smallRNA)
DE_smallRNA_sex <- DESeq(dds_smallRNA_sex)
DE_smallRNA_female <- DESeq(dds_smallRNA_female)
DE_smallRNA_male <- DESeq(dds_smallRNA_male)

# -> miRNA
DE_miRNA <- DESeq(dds_miRNA)
DE_miRNA_sex <- DESeq(dds_miRNA_sex)
DE_miRNA_female <- DESeq(dds_miRNA_female)
DE_miRNA_male <- DESeq(dds_miRNA_male)

########### Differential expression analysis (without outlier G103) ###########
# -> smallRNA 
DE_smallRNA_noutlier <- DESeq(dds_smallRNA[,-which(colnames(dds_smallRNA) == "G103")])
DE_smallRNA_sex_noutlier <- DESeq(dds_smallRNA_sex[,-which(colnames(dds_smallRNA_sex) == "G103")])

# -> miRNA
DE_miRNA_noutlier <- DESeq(dds_miRNA[,-which(colnames(dds_miRNA) == "G103")])
DE_miRNA_sex_noutlier <- DESeq(dds_miRNA_sex[,-which(colnames(dds_miRNA_sex) == "G103")])

########### Retrieve ordered results from differential expression analysis (with outlier) ###########
# -> smallRNA
res_smallRNA <- results(DE_smallRNA)
res_smallRNA <- res_smallRNA[order(res_smallRNA$padj),]
res_smallRNA_sex <- results(DE_smallRNA_sex)
res_smallRNA_sex <- res_smallRNA_sex[order(res_smallRNA_sex$padj),]
res_smallRNA_female <- results(DE_smallRNA_female)
res_smallRNA_female <- res_smallRNA_female[order(res_smallRNA_female$padj),]
res_smallRNA_male <- results(DE_smallRNA_male)
res_smallRNA_male <- res_smallRNA_male[order(res_smallRNA_male$padj),]

# -> miRNA
res_miRNA <- results(DE_miRNA)
res_miRNA <- res_miRNA[order(res_miRNA$padj),]
res_miRNA_sex <- results(DE_miRNA_sex)
res_miRNA_sex <- res_miRNA_sex[order(res_miRNA_sex$padj),]
res_miRNA_female <- results(DE_miRNA_female)
res_miRNA_female <- res_miRNA_female[order(res_miRNA_female$padj),]
res_miRNA_male <- results(DE_miRNA_male)
res_miRNA_male <- res_miRNA_male[order(res_miRNA_male$padj),]

########### Retrieve ordered results from differential expression analysis (without outlier) ###########
# -> smallRNA
res_smallRNA_noutlier <- results(DE_smallRNA_noutlier)
res_smallRNA_noutlier <- res_smallRNA_noutlier[order(res_smallRNA_noutlier$padj),]
res_smallRNA_sex_noutlier <- results(DE_smallRNA_sex_noutlier)
res_smallRNA_sex_noutlier <- res_smallRNA_sex_noutlier[order(res_smallRNA_sex_noutlier$padj),]

# -> miRNA
res_miRNA_noutlier <- results(DE_miRNA_noutlier)
res_miRNA_noutlier <- res_miRNA_noutlier[order(res_miRNA_noutlier$padj),]
res_miRNA_sex_noutlier <- results(DE_miRNA_sex_noutlier)
res_miRNA_sex_noutlier <- res_miRNA_sex_noutlier[order(res_miRNA_sex_noutlier$padj),]

########### Save results ###########
dir.create(file.path(cidr,"data_output"),recursive = T)
# -> smallRNA
write_xlsx(data.frame(smallRNA = row.names(res_smallRNA),res_smallRNA),"data_output/DE_smallRNA.xlsx")
write_xlsx(data.frame(smallRNA = row.names(res_smallRNA_sex),res_smallRNA_sex),"data_output/DE_smallRNA_sex.xlsx")
write_xlsx(data.frame(smallRNA = row.names(res_smallRNA_noutlier),res_smallRNA_noutlier),"data_output/DE_smallRNA_noOutlier.xlsx")
write_xlsx(data.frame(smallRNA = row.names(res_smallRNA_sex_noutlier),res_smallRNA_sex_noutlier),"data_output/DE_smallRNA_sex_noOutlier.xlsx")
write_xlsx(data.frame(smallRNA = row.names(res_smallRNA_male),res_smallRNA_male),"data_output/DE_smallRNA_female.xlsx")
write_xlsx(data.frame(smallRNA = row.names(res_smallRNA_female),res_smallRNA_female),"data_output/DE_smallRNA_male.xlsx")

# -> miRNA
write_xlsx(data.frame(miRNA = row.names(res_miRNA),res_miRNA),"data_output/DE_miRNA.xlsx")
write_xlsx(data.frame(miRNA = row.names(res_miRNA_sex),res_miRNA_sex),"data_output/DE_miRNA_sex.xlsx")
write_xlsx(data.frame(miRNA = row.names(res_miRNA_noutlier),res_miRNA_noutlier),"data_output/DE_miRNA_noOutlier.xlsx")
write_xlsx(data.frame(miRNA = row.names(res_miRNA_sex_noutlier),res_miRNA_sex_noutlier),"data_output/DE_miRNA_sex_noOutlier.xlsx")
write_xlsx(data.frame(miRNA = row.names(res_miRNA_female),res_miRNA_female),"data_output/DE_miRNA_female.xlsx")
write_xlsx(data.frame(miRNA = row.names(res_miRNA_male),res_miRNA_male),"data_output/DE_miRNA_male.xlsx")

########### Visualization 3: Volcano plots ###########
volcano_plot <- function(data_res,alpha_sig,name_title){
  log2fc <- grep("log2FoldChange",colnames(data_res))
  padj <- grep("padj",colnames(data_res))
  df <- data.frame(x = data_res[,log2fc], 
                   y = -log10(data_res[,padj]),
                   name = rownames(data_res)) %>%
    mutate(name = lapply(strsplit(name,"\\|"), function(x) x[1]))
  names(df) <- c("x","y","name")
  df <- df %>%
    mutate(omic_type = case_when(x >= 0 & y >= (-log10(alpha_sig)) ~ "up",
                                 x <= (0) & y >= (-log10(alpha_sig)) ~ "down",
                                 TRUE ~ "ns")) 
  cols <- c("up" = "#d4552b", "down" = "#26b3ff", "ns" = "grey") 
  sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
  alphas <- c("up" = 0.7, "down" = 0.7, "ns" = 0.5)
  ggplot(data = df, aes(x,y)) + 
    geom_point(aes(colour = omic_type), 
               alpha = 0.5, 
               shape = 16,
               size = 3) + 
    geom_hline(yintercept = log10(10),
               linetype = "dashed") + 
    geom_vline(xintercept = 0,linetype = "dashed") +
    geom_point(data = filter(df, y >= (-log10(alpha_sig))),
               aes(colour = omic_type), 
               alpha = 0.5, 
               shape = 16,
               size = 4) + 
    #annotate(geom="text", x=-1.9, y= (-log10(alpha_sig)) + 0.15, label="FDR = 10%",size = 5) +
    geom_text_repel(data = filter(df, y >= (-log10(alpha_sig)) & y > 0),
                     aes(label = name),
                     force = 1,
                    hjust = 1,
                     nudge_x = - 0.3,
                    nudge_y = 0.1,
                    #direction = "x",
                     max.overlaps = 5,
                    segment.size = 0.2,
                     size = 4) +
    geom_text_repel(data = filter(df, y >= (-log10(alpha_sig)) & y < 0),
                    aes(label = name),
                    force = 1,
                    hjust = 0,
                    nudge_x = 0.3,
                    nudge_y = 0.1,
                    #direction = "y",
                    max.overlaps = 5,
                    size = 4) +
    scale_colour_manual(values = cols) + 
    scale_fill_manual(values = cols) + 
    scale_x_continuous(expand = c(0, 0), 
                       limits = c(-7, 8.5)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(-0.1, NA)) +
    labs(title = name_title,
         x = "log2(fold change)",
         y = expression(-log[10] ~ "(adjusted p-value)"),
         colour = "Differential \nExpression") +
    theme_classic() + # Select theme with a white background  
    theme(axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 15, hjust = 0.5),
          text = element_text(size = 14))
}

# -> smallRNA
pdf("plots_output/volcano_smallRNA.pdf",paper = "a4")
volcano_plot(as.data.frame(res_smallRNA),0.1,"Volcano plot smallRNA (ALS vs Ctrl)")
dev.off()
pdf("plots_output/volcano_smallRNA_noOutlier.pdf",paper = "a4")
volcano_plot(as.data.frame(res_smallRNA_noutlier),0.1,"Volcano plot smallRNA without G103 (ALS vs Ctrl)")
dev.off()
pdf("plots_output/volcano_smallRNA_sex.pdf",paper = "a4")
volcano_plot(as.data.frame(res_smallRNA_sex),0.1,"Volcano plot smallRNA + sex covariate (ALS vs Ctrl)")
dev.off()
pdf("plots_output/volcano_smallRNA_sex_noOutlier.pdf",paper = "a4")
volcano_plot(as.data.frame(res_smallRNA_sex_noutlier),0.1,"Volcano plot smallRNA + sex covariate without G103 (ALS vs Ctrl)")
dev.off()
pdf("plots_output/volcano_smallRNA_female.pdf",paper = "a4")
volcano_plot(as.data.frame(res_smallRNA_female),0.1,"Volcano plot smallRNA - Female only (ALS vs Ctrl)")
dev.off()
pdf("plots_output/volcano_smallRNA_male.pdf",paper = "a4")
volcano_plot(as.data.frame(res_smallRNA_male),0.1,"Volcano plot smallRNA - Male only (ALS vs Ctrl)")
dev.off()

# -> miRNA
pdf("plots_output/volcano_miRNA.pdf",paper = "a4")
volcano_plot(as.data.frame(res_miRNA),0.1,"Volcano plot miRNA (ALS vs Ctrl)")
dev.off()
pdf("plots_output/volcano_miRNA_noOutlier.pdf",paper = "a4")
volcano_plot(as.data.frame(res_miRNA_noutlier),0.1,"Volcano plot miRNA without G103 (ALS vs Ctrl)")
dev.off()
pdf("plots_output/volcano_miRNA_sex.pdf",paper = "a4")
volcano_plot(as.data.frame(res_miRNA_sex),0.1,"Volcano plot miRNA + sex covariate (ALS vs Ctrl)")
dev.off()
pdf("plots_output/volcano_miRNA_sex_noOutlier.pdf",paper = "a4")
volcano_plot(as.data.frame(res_miRNA_sex_noutlier),0.1,"Volcano plot miRNA + sex covariate without G103 (ALS vs Ctrl)")
dev.off()
pdf("plots_output/volcano_miRNA_female.pdf",paper = "a4")
volcano_plot(as.data.frame(res_miRNA_female),0.1,"Volcano plot miRNA - Female only (ALS vs Ctrl)")
dev.off()
pdf("plots_output/volcano_miRNA_male.pdf",paper = "a4")
volcano_plot(as.data.frame(res_miRNA_male),0.1,"Volcano plot miRNA - Male only (ALS vs Ctrl)")
dev.off()
