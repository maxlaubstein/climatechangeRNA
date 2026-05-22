#Venn diagram
setwd("/n/home00/mlaubstein/climatechangeRNA/")
set.seed(16)
library(ggVennDiagram)
library(ggplot2)
library(ggforce)

system("rm *GSEA_rect_venn.pdf Combined_Rect_Venn_fig.pdf *GSEA_circle_venn.pdf")

load("Myiarchus_Discrete/GSEA_Myiarchus/ATFL_gsea_output.RData")
ATFL_gsea <- gsea_all_outputs
rm("gsea_all_outputs")

load("Thryomanes_Discrete/GSEA_Thryomanes/BEWR_gsea_output.RData")
BEWR_gsea <- gsea_all_outputs
rm("gsea_all_outputs")

for(tissue in c("heart", "liver", "muscle", "kidney", "brain")){
  tissue = tissue
  BEWR_df <- BEWR_gsea[[tissue]]
  ATFL_df <- ATFL_gsea[[tissue]]
  both <- intersect(unique(BEWR_df$ID), unique(ATFL_df$ID))
  BEWR_only <- setdiff(unique(BEWR_df$ID), both) 
  ATFL_only <- setdiff(unique(ATFL_df$ID), both) 
  n_both <- length(both)
  n_BEWR_only <- length(BEWR_only)
  n_ATFL_only <- length(ATFL_only)
  plot <- ggplot()+
    annotate("rect", xmin = 0, xmax = n_BEWR_only+n_both, ymin = 0, ymax = 2, fill = '#8B5A2B', alpha = 0.5)+ #left rectange for BEWR
    annotate("rect", xmin = n_BEWR_only, xmax = n_BEWR_only+n_both+n_ATFL_only, ymin = 0, ymax = 2, fill = '#E1C85A', alpha = 0.5)+ #right rectangle for ATFL
    annotate("text", x = 1, y = 1.6, label = n_BEWR_only, hjust = 0)+
    scale_x_continuous(limits = c(0, 200))+
    theme_void()
  assign(paste0(tissue,"_rect_venn"), plot)
  ggsave(filename=paste0(tissue,"_GSEA_rect_venn.pdf"), plot, width = 6, height = 0.5, units = "in")
  message(tissue)
  cat(c(n_BEWR_only, n_both, n_ATFL_only, "\n"))
  cat(dQuote(unique(subset(rbind(BEWR_df, ATFL_df), rbind(BEWR_df, ATFL_df)$ID %in% both)$Description)), sep = ", ")
  print(" ")
  if(tissue == "brain"){
    axis <- ggplot()+
      theme_minimal()+
      xlab("# GO Terms")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      scale_x_continuous(limits = c(0, 200))
    ggsave(filename="axis.pdf", axis, width = 6, height = 0.4, units = "in")
    
  }
}

pdf_files <- c("brain_GSEA_rect_venn.pdf",
               "kidney_GSEA_rect_venn.pdf",
               "liver_GSEA_rect_venn.pdf",
               "muscle_GSEA_rect_venn.pdf",
               "heart_GSEA_rect_venn.pdf")
tex <- c(
  "\\documentclass{article}",
  "\\usepackage[paperwidth=7in,paperheight=8in,",
  "            left=0in,right=0in,top=0in,bottom=0in]{geometry}",
  "\\usepackage{graphicx}",
  "\\begin{document}",
  "\\thispagestyle{empty}",
  "\\noindent",
  sprintf("\\includegraphics[width=7in]{%s}\\\\[0.075in]", pdf_files[1]),
  sprintf("\\includegraphics[width=7in]{%s}\\\\[0.075in]", pdf_files[2]),
  sprintf("\\includegraphics[width=7in]{%s}\\\\[0.075in]", pdf_files[3]),
  sprintf("\\includegraphics[width=7in]{%s}\\\\[0.075in]", pdf_files[4]),
  sprintf("\\includegraphics[width=7in]{%s}\\\\", pdf_files[5]),
  sprintf("\\includegraphics[width=7in]{%s}", "axis.pdf"),
  "\\end{document}"
)
writeLines(tex, "Combined_Rect_Venn_fig.tex")
tinytex::latexmk("Combined_Rect_Venn_fig.tex")
system("rm Combined_Rect_Venn_fig.tex axis.pdf")



##########################################
#Regular Circle Venn Diagrams
##########################################

for(tissue in c("heart", "liver", "muscle", "kidney", "brain")){
  tissue = tissue
  BEWR_df <- BEWR_gsea[[tissue]]
  ATFL_df <- ATFL_gsea[[tissue]]
  both <- intersect(unique(BEWR_df$ID), unique(ATFL_df$ID))
  BEWR_only <- setdiff(unique(BEWR_df$ID), both) 
  ATFL_only <- setdiff(unique(ATFL_df$ID), both) 
  n_both <- length(both)
  n_BEWR_only <- length(BEWR_only)
  n_ATFL_only <- length(ATFL_only)
  
  
  if(n_both > 0){
    plot <- ggplot()+
      geom_circle(aes(x0 = 0, y0 = 0, r = 1), fill = '#8B5A2B', alpha = 0.7)+
      geom_circle(aes(x0 = 1.25, y0 = 0, r = 1), fill = '#E1C85A', alpha = 0.7)+
      annotate("text", x = -0.2, y = 0, label = n_BEWR_only)+
      annotate("text", x = 0.625, y = 0, label = n_both)+
      annotate("text", x = 1.45, y = 0, label = n_ATFL_only)+
      theme_void()+
      coord_equal()
  }
  
  if(n_both == 0){
    plot <- ggplot()+
      geom_circle(aes(x0 = 0, y0 = 0, r = 1), fill = '#8B5A2B', alpha = 0.7)+
      geom_circle(aes(x0 = 2.2, y0 = 0, r = 1), fill = '#E1C85A', alpha = 0.7)+
      annotate("text", x = 0, y = 0, label = n_BEWR_only)+
      annotate("text", x = 2.2, y = 0, label = n_ATFL_only)+
      theme_void()+
      coord_equal()
  }
  
  
  
  assign(paste0(tissue,"_circle_venn"), plot)
  ggsave(filename=paste0(tissue,"_GSEA_circle_venn.pdf"), plot = plot, width = 3.5, height = 2, units = "in")
  message(tissue)
  cat(c(n_BEWR_only, n_both, n_ATFL_only, "\n"))
  cat(dQuote(unique(subset(rbind(BEWR_df, ATFL_df), rbind(BEWR_df, ATFL_df)$ID %in% both)$Description)), sep = ", ")
  print(" ")
}

pdf_files <- c(
  "brain_GSEA_circle_venn.pdf",
  "kidney_GSEA_circle_venn.pdf",
  "liver_GSEA_circle_venn.pdf",
  "muscle_GSEA_circle_venn.pdf",
  "heart_GSEA_circle_venn.pdf"
)

tex <- c(
  "\\documentclass{article}",
  "\\usepackage[paperwidth=7in,paperheight=8in,",
  "            left=0in,right=0in,top=0in,bottom=0in]{geometry}",
  "\\usepackage{graphicx}",
  "\\begin{document}",
  "\\thispagestyle{empty}",
  "\\noindent",
  
  sprintf("\\includegraphics[width=3.5in]{%s}", pdf_files[1]),
  sprintf("\\includegraphics[width=3.5in]{%s}\\\\[0.075in]", pdf_files[2]),
  
  sprintf("\\includegraphics[width=3.5in]{%s}", pdf_files[3]),
  sprintf("\\includegraphics[width=3.5in]{%s}\\\\[0.075in]", pdf_files[4]),
  
  sprintf(
    "\\centerline{\\includegraphics[width=3.5in]{%s}}",
    pdf_files[5]
  ),
  
  "\\end{document}"
)

writeLines(tex, "Combined_Circle_Venn_fig.tex")

tinytex::latexmk("Combined_Circle_Venn_fig.tex")

system("rm Combined_Circle_Venn_fig.tex")
