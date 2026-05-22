#This script runs the cross-species differential expression analysis on all...
#...single copy orthologous genes between Bewick's Wren and Ash-throated flycatcher...
#... treating treatment versus control as a discrete, binary independent variable.
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tximport)
library(grid)
library(stringr)
setwd("~/climatechangeRNA/Combined_Discrete/")

system("rm *sig_DEGs.csv")
system("rm *all_genes.csv")
system("rm *pdf")

load("/n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/climatechangeRNA/Thryomanes_output/kallisto/deseq2_qc/deseq2.dds.RData")
BEWR_dds <- dds
rm(dds)
load("/n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/climatechangeRNA/Myiarchus_output//kallisto/deseq2_qc/deseq2.dds.RData")
ATFL_dds <- dds
rm(dds)
metadata <- read.csv('~/climatechangeRNA/Metadata_Climate_Change_RNA_Final.csv', header = TRUE, sep = ",", check.names = FALSE)

orthologs <- read.csv("/n/holylfs06/LABS/edwards_lab/Lab/maxlaubstein/climatechangeRNA/TOGA_orthology/taeGut.myiCay.orthology/orthology_classification.tsv", sep = "\t")
one2one <- subset(orthologs, orthologs$orthology_class == "one2one")
one2one_genes <- unique(one2one$t_gene)

sum(one2one_genes %in% rownames(BEWR_dds))
sum(one2one_genes %in% rownames(ATFL_dds))

common_genes <- Reduce(intersect, list(one2one_genes, rownames(BEWR_dds), rownames(ATFL_dds)))

BEWR_dds <- BEWR_dds[common_genes,]
ATFL_dds <- ATFL_dds[common_genes,]

message(paste("There are", length(common_genes), "shared one-to-one orthologs"))

combined_dds <- cbind(ATFL_dds, BEWR_dds)

colData(combined_dds) <- DataFrame(metadata)
combined_dds$Tissue <- factor(combined_dds$Tissue)
combined_dds$Treatment <- factor(combined_dds$Treatment)
combined_dds$Species <- factor(combined_dds$Species)

keep <- rowSums(counts(combined_dds) >= 10) >= 3 # only keep data points where the gene count is greater than 10 in at least 3 samples 
combined_dds <- combined_dds[keep,]
design(combined_dds) <- ~ Treatment + Tissue
combined_dds <- DESeq(combined_dds)
resultsNames(combined_dds)
res <- results(combined_dds, contrast = list("Treatment_Heat_vs_Control"), alpha = 0.05)
summary(res)
#perform variance stabilizing transformation:
vsd <- vst(combined_dds, blind=FALSE)

pca_df <- plotPCA(vsd, intgroup=c("Species", "Tissue"), returnData = TRUE)
all_PCA_plot <- ggplot(pca_df)+
  geom_point(aes(x = PC1, y = PC2, shape = Tissue, fill = Species, color = Species), size = 3, stroke = 1)+
  scale_fill_manual("Species", values = c(`Thryomanes bewickii` = "#8B5A2B", `Myiarchus cinerascens` = "#E1C85A"))+
  scale_color_manual("Species", values = c(`Thryomanes bewickii` = "#5B3A1C", `Myiarchus cinerascens` = "#7F7F7F"))+
  scale_shape_manual("Tissue", values = 21:25)+
  theme_minimal()+
  ggtitle("Figure 1. PCA of Gene Expression Across Species and Tissues")
all_PCA_plot_fig <- ggplot(pca_df)+
  geom_point(aes(x = PC1, y = PC2, shape = Tissue, fill = Species, color = Species), size = 3, stroke = 1)+
  scale_fill_manual("Species", values = c(`Thryomanes bewickii` = "#8B5A2B", `Myiarchus cinerascens` = "#E1C85A"))+
  scale_color_manual("Species", values = c(`Thryomanes bewickii` = "#5B3A1C", `Myiarchus cinerascens` = "#7F7F7F"))+
  scale_shape_manual("Tissue", values = 21:25, labels = function(x) tools::toTitleCase(tolower(x)))+
  xlab(paste0("PC1 (", round(attr(pca_df, "percentVar")[1]*100, 2),"%)" ))+
  ylab(paste0("PC2 (", round(attr(pca_df, "percentVar")[2]*100, 2),"%)" ))+
  theme_minimal()
ggsave("Cross_Species_PCA_All_Tissues.pdf", all_PCA_plot_fig, width = 6, height = 4, units = "in")



####Treatment Effects Within Tissues
tissues <- c("heart", "liver", "muscle", "kidney", "brain")

for(i in 1:length(tissues)){
  message(paste0("Analyzing ", tissues[i], "..."))
  dds <- combined_dds
  colData(dds) <- DataFrame(metadata)
  dds$Tissue <- factor(dds$Tissue)
  dds$Treatment <- factor(dds$Treatment)
  dds$Species <- factor(dds$Species)
  dds <- dds[, colData(dds)$Tissue == tissues[i]]
  design(dds) <- ~ Treatment + Species
  keep <- rowSums(counts(dds) >= 10) >= 3 # only keep data points where the gene count is greater than 10 in at least 3 samples 
  dds <- dds[keep,]
  dds <- DESeq(dds)
  res <- results(dds, contrast = list("Treatment_Heat_vs_Control"), alpha = 0.05)
  vsd <- vst(dds, blind = FALSE)
  
  pca_df <- DESeq2::plotPCA(vsd, intgroup=c("Treatment", "Species"), returnData = TRUE)
  pcaplot <- ggplot(pca_df)+
    geom_point(aes(x = PC1, y = PC2, color = Treatment, shape = Species), size = 3)+
    scale_color_manual("Treatment", values = c("#E63946", "#74C0E3"))+
    theme_minimal()+
    ggtitle(paste0("Figure ", i+1, "A. PCA of Gene Expression in ", tools::toTitleCase(tissues[i]), " In Both Species Across Treatments"))+
    xlab(paste0("PC1 (", round(attr(pca_df, "percentVar")[1]*100, 2),"%)" ))+
    ylab(paste0("PC2 (", round(attr(pca_df, "percentVar")[2]*100, 2),"%)" ))
  assign(paste0(tissues[i], "_PCA_plot"), pcaplot)
  message(paste0("PCA plot stored in '", paste0(tissues[i], "_Combined_PCA_plot","'")))
  
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
    theme_minimal()+
    ggtitle(paste0("Figure ", i+1, "B. Cross-Species Volcano Plot of Heated vs. Control in ", tools::toTitleCase(tissues[i])))
  assign(paste0(tissues[i], "_volcano_plot"), volcanoplot)
  message(paste0("Volcano plot stored in '", paste0(tissues[i], "_Combined_volcano_plot","'")))
  
  DEG_df <- subset(volcano, abs(volcano$log2FoldChange) >= 1 & volcano$padj < 0.05)
  
  write.csv(DEG_df, file = paste0(tissues[i], "_sig_DEGs.csv"))
  message(paste0("Significant DEGs file stored in '", paste0(tissues[i], "_Combined_sig_DEGs.csv","'")))
  if(nrow(DEG_df) == 1){
    message(paste0("1 significant DEG was identified in ", tissues[i]))
  }
  if(nrow(DEG_df) > 1){
    message(paste0(nrow(DEG_df), " DEGs were identified in ", tissues[i]))
  }
  if(nrow(DEG_df) == 0){
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

pdf("Cross_Species_Figures.pdf", width = 8.5, height = 11)

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
ggsave("Combined_Brain_Volcano.pdf", brain_volcano_plot_fig, width = 6, height = 4, units = "in")

heart_volcano_plot_fig <- heart_volcano_plot + ggtitle(NULL)
ggsave("Combined_Heart_Volcano.pdf", heart_volcano_plot_fig, width = 6, height = 4, units = "in")

liver_volcano_plot_fig <- liver_volcano_plot + ggtitle(NULL)
ggsave("Combined_Liver_Volcano.pdf", liver_volcano_plot_fig, width = 6, height = 4, units = "in")

muscle_volcano_plot_fig <- muscle_volcano_plot + ggtitle(NULL)
ggsave("Combined_Muscle_Volcano.pdf", muscle_volcano_plot_fig, width = 6, height = 4, units = "in")

kidney_volcano_plot_fig <- kidney_volcano_plot + ggtitle(NULL)
ggsave("Combined_Kidney_Volcano.pdf", kidney_volcano_plot_fig, width = 6, height = 4, units = "in")





