library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tximport)
library(grid)
library(tidyr)
set.seed(16)
setwd("/n/home00/mlaubstein/climatechangeRNA/Thryomanes_Discrete/")
system("rm *sig_DEGs.csv")
system("rm *all_genes.csv")
system("rm combined_output_Thryomanes.csv")
system("rm *Volcano.pdf")
system("rm BEWR_Specific_PCA_All_Tissues.pdf")
system("rm Summary_Plot.pdf")
system("rm Thryomanes_Figures.pdf")

load("/n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/climatechangeRNA/Thryomanes_output/kallisto/deseq2_qc/deseq2.dds.RData")
dds_og <- dds
rm("dds")
metadata <- read.csv('~/climatechangeRNA/Metadata_Climate_Change_RNA_Final.csv', header = TRUE, sep = ",", check.names = FALSE)

metadata$id <- paste0(metadata$`Catalog Number`, "_", metadata$Tissue)
metadata <- subset(metadata, metadata$Species == "Thryomanes bewickii")
metadata <- metadata %>% dplyr::select("id", "Tissue", "Treatment", "Sex", Individual = "Catalog Number")

rownames(metadata) <- metadata$id
metadata <- metadata[colnames(dds_og), ]
all(rownames(metadata) == colnames(dds_og))


####Tissue & Treatment Effects:
dds <- dds_og
colData(dds) <- DataFrame(metadata)
dds$Tissue <- factor(dds$Tissue)
dds$Treatment <- factor(dds$Treatment)
design(dds) <- ~ Treatment + Tissue
keep <- rowSums(counts(dds) >= 10) >= 3 # only keep data points where the gene count is greater than 10 in at least 3 samples 
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast = list("Treatment_Heat_vs_Control"), alpha = 0.05)
summary(res)
#perform variance stabilizing transformation:
vsd <- vst(dds, blind=FALSE)
pca_df <- plotPCA(vsd, intgroup=c("Treatment", "Tissue"), returnData = TRUE)
all_PCA_plot <- ggplot(pca_df)+
  geom_point(aes(x = PC1, y = PC2, shape = Tissue, fill = Treatment, color = Treatment), size = 3, stroke = 1)+
  scale_fill_manual("Treatment", values = c(Heat = "#E63946", Control = "#74C0E3"))+
  scale_color_manual("Treatment", values = c(Heat = "#E63946", Control = "#74C0E3"))+
  xlab(paste0("PC1 (", round(attr(pca_df, "percentVar")[1]*100, 2),"%)" ))+
  ylab(paste0("PC2 (", round(attr(pca_df, "percentVar")[2]*100, 2),"%)" ))+
  scale_shape_manual("Tissue", values = 21:25)+
  theme_minimal()+
  ggtitle("Figure 1. PCA of Gene Expression Across Treatments and Tissues")

all_PCA_plot_fig <- ggplot(pca_df)+
  geom_point(aes(x = PC1, y = PC2, shape = Tissue, fill = Treatment, color = Treatment), size = 3, stroke = 1)+
  xlab(paste0("PC1 (", round(attr(pca_df, "percentVar")[1]*100, 2),"%)" ))+
  ylab(paste0("PC2 (", round(attr(pca_df, "percentVar")[2]*100, 2),"%)" ))+
  scale_fill_manual("Treatment", values = c(Heat = "#E63946", Control = "#74C0E3"))+
  scale_color_manual("Treatment", values = c(Heat = "#E63946", Control = "#74C0E3"))+
  scale_shape_manual("Tissue", values = 21:25)+
  theme_minimal()
ggsave("BEWR_Specific_PCA_All_Tissues.pdf", all_PCA_plot_fig, width = 6, height = 4, units = "in")

####Treatment Effects Within Tissues
tissues <- c("heart", "liver", "muscle", "kidney", "brain")

for(i in 1:length(tissues)){
  message(paste0("Analyzing ", tissues[i], "..."))
  dds <- dds_og
  colData(dds) <- DataFrame(metadata)
  dds$Tissue <- factor(dds$Tissue)
  dds$Treatment <- factor(dds$Treatment)
  dds <- dds[, colData(dds)$Tissue == tissues[i]]
  design(dds) <- ~ Treatment
  keep <- rowSums(counts(dds) >= 10) >= 3 # only keep data points where the gene count is greater than 10 in at least 3 samples 
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds, contrast = list("Treatment_Heat_vs_Control"), alpha = 0.05)
  vsd <- vst(dds, blind = FALSE)
  
  pca_df <- DESeq2::plotPCA(vsd, intgroup=c("Treatment"), returnData = TRUE)
  pcaplot <- ggplot(pca_df)+
    geom_point(aes(x = PC1, y = PC2, color = Treatment), size = 3)+
    scale_color_manual("Treatment", values = c("#E63946", "#74C0E3"))+
    theme_minimal()+
    ggtitle(paste0("Figure ", i+1, "A. PCA of Gene Expression in ", tools::toTitleCase(tissues[i]), " Across Treatments"))
  assign(paste0(tissues[i], "_PCA_plot"), pcaplot)
  message(paste0("PCA plot stored in '", paste0(tissues[i], "_PCA_plot","'")))
  
  volcano <- as.data.frame(res)
  volcano$gene <- rownames(volcano)
  volcano$sig <- volcano$padj < 0.05 & abs(volcano$log2FoldChange) >= 1
  volcano <- na.omit(volcano)
  
  volcanoplot <- ggplot(volcano, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(color=sig), alpha=0.7, size=2, na.rm=TRUE)+
    scale_color_manual(values=c("grey60","red"), labels=NULL, guide="none") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="gray70") +
    geom_vline(xintercept=c(-1,1), linetype="dashed", color="gray70") +
    xlab(expression(log[2] ~ "Fold Change"))+
    ylab(expression(-log[10]*P[adj.]))+
    xlim(-10,10)+
    ylim(0,6)+
    theme_minimal()+
    ggtitle(paste0("Figure ", i+1, "B. Volcano Plot of Heated vs. Control in ", tools::toTitleCase(tissues[i])))
  assign(paste0(tissues[i], "_volcano_plot"), volcanoplot)
  message(paste0("Volcano plot stored in '", paste0(tissues[i], "_volcano_plot","'")))
  
  sig_DEG_df <- subset(volcano, abs(volcano$log2FoldChange) >= 1 & volcano$padj <= 0.05)
  
  write.csv(sig_DEG_df, file = paste0(tissues[i], "_sig_DEGs.csv"))
  message(paste0("Significant DEGs file stored in '", paste0(tissues[i], "_sig_DEGs.csv","'")))
  if(nrow(sig_DEG_df) == 1){
    message(paste0("1 significant DEG was identified in ", tissues[i]))
  }
  if(nrow(sig_DEG_df) > 1){
    message(paste0(nrow(sig_DEG_df), " DEGs were identified in ", tissues[i]))
  }
  if(nrow(sig_DEG_df) == 0){
    message(paste0("No significant DEGs were identified in ", tissues[i]))
  }
  
  all_genes_df <- volcano
  all_genes_df$sig <- NULL
  write.csv(all_genes_df, file = paste0(tissues[i], "_all_genes.csv"))
  
  message("\n")
}

all_plots <- list(
  all_PCA = all_PCA_plot,
  heart_PCA = heart_PCA_plot,
  heart_volcano = heart_volcano_plot,
  
  liver_PCA = liver_PCA_plot,
  liver_volcano = liver_volcano_plot,
  
  muscle_PCA = muscle_PCA_plot,
  muscle_volcano = muscle_volcano_plot,
  
  kidney_PCA = kidney_PCA_plot,
  kidney_volcano = kidney_volcano_plot,
  
  brain_PCA = brain_PCA_plot,
  brain_volcano = brain_volcano_plot
)

pdf("Thryomanes_Figures.pdf", width = 8.5, height = 11)

for (p in all_plots) {
  
  grid::grid.newpage()
  
  pushViewport(
    grid::viewport(
      width = unit(8, "in"),
      height = unit(6, "in"),
      x = 0.5,
      y = 0.5,
      just = "center"
    )
  )
  
  if (inherits(p, "ggplot")) {
    print(p, newpage = FALSE)
  } else {
    grid::grid.draw(p$gtable)
  }
  
  popViewport()
}

dev.off()

brain_volcano_plot_fig <- brain_volcano_plot + ggtitle(NULL)
ggsave("BEWR_Brain_Volcano.pdf", brain_volcano_plot_fig, width = 6, height = 4, units = "in")

heart_volcano_plot_fig <- heart_volcano_plot + ggtitle(NULL)
ggsave("BEWR_Heart_Volcano.pdf", heart_volcano_plot_fig, width = 6, height = 4, units = "in")

liver_volcano_plot_fig <- liver_volcano_plot + ggtitle(NULL)
ggsave("BEWR_Liver_Volcano.pdf", liver_volcano_plot_fig, width = 6, height = 4, units = "in")

muscle_volcano_plot_fig <- muscle_volcano_plot + ggtitle(NULL)
ggsave("BEWR_Muscle_Volcano.pdf", muscle_volcano_plot_fig, width = 6, height = 4, units = "in")

kidney_volcano_plot_fig <- kidney_volcano_plot + ggtitle(NULL)
ggsave("BEWR_Kidney_Volcano.pdf", kidney_volcano_plot_fig, width = 6, height = 4, units = "in")

summary_plot_df <- data.frame(
  Tissue = c("Heart", "Liver", "Muscle", "Kidney", "Brain"),
  nDEGs = c(
    nrow(read.csv("heart_sig_DEGs.csv")),
    nrow(read.csv("liver_sig_DEGs.csv")),
    nrow(read.csv("muscle_sig_DEGs.csv")),
    nrow(read.csv("kidney_sig_DEGs.csv")),
    nrow(read.csv("brain_sig_DEGs.csv"))
  ),
  nUP = c(
    sum(read.csv("heart_sig_DEGs.csv")$log2FoldChange > 0),
    sum(read.csv("liver_sig_DEGs.csv")$log2FoldChange > 0),
    sum(read.csv("muscle_sig_DEGs.csv")$log2FoldChange > 0),
    sum(read.csv("kidney_sig_DEGs.csv")$log2FoldChange > 0),
    sum(read.csv("brain_sig_DEGs.csv")$log2FoldChange > 0)
  ),
  nDOWN = c(
    sum(read.csv("heart_sig_DEGs.csv")$log2FoldChange < 0),
    sum(read.csv("liver_sig_DEGs.csv")$log2FoldChange < 0),
    sum(read.csv("muscle_sig_DEGs.csv")$log2FoldChange < 0),
    sum(read.csv("kidney_sig_DEGs.csv")$log2FoldChange < 0),
    sum(read.csv("brain_sig_DEGs.csv")$log2FoldChange < 0)
  )
)

summary_plot_df_long <- summary_plot_df |>
  dplyr::select(Tissue, nUP, nDOWN) |>
  pivot_longer(cols = c(nUP, nDOWN),
               names_to = "UPorDOWN",
               values_to = "Count")


summary_plot <- ggplot(summary_plot_df_long, aes(x = Tissue, y = Count, fill = UPorDOWN))+
  geom_col()+
  scale_fill_manual(NULL, values = c(nUP = "#E63946", nDOWN = "#74C0E3"),
                    labels = c(nUP = "Upregulated", nDOWN = "Downregulated"))+
  ylab("# Significant DEGs")+
  theme_minimal()

ggsave("Summary_Plot.pdf", summary_plot, width = 6, height = 4, units = "in")

heart_genes <- read.csv("heart_all_genes.csv")
heart_genes$tissue <- "heart"

liver_genes <- read.csv("liver_all_genes.csv")
liver_genes$tissue <- "liver"

muscle_genes <- read.csv("muscle_all_genes.csv")
muscle_genes$tissue <- "muscle"

kidney_genes <- read.csv("kidney_all_genes.csv")
kidney_genes$tissue <- "kidney"

brain_genes <- read.csv("brain_all_genes.csv")
brain_genes$tissue <- "brain"

combined_output <- rbind(heart_genes, liver_genes, muscle_genes, kidney_genes, brain_genes)
write.csv(combined_output, file = "combined_output_Thryomanes.csv")

