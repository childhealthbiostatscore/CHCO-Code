###Making Multiomic VSG Figures for Long 



################################################################################
# FIGURE 1 ASSEMBLY - VSG vs SMT Manuscript
# 
# Your file location: 
# C:\Users\netio\Documents\UofW\Projects\Long_Paper\VSG_analysis-main\VSG_analysis-main\output\Figure1
################################################################################

library(magick)
library(pdftools)
library(cowplot)
library(ggplot2)
library(patchwork)

# Set your base path
base_path <- "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"

# Helper function
pdf_to_plot <- function(path, density = 300) {
  if (!file.exists(path)) {
    stop(paste("File not found:", path))
  }
  img <- image_read_pdf(path, density = density)
  ggdraw() + draw_image(img)
}

################################################################################
# FIGURE 1: TRANSCRIPTIONAL LANDSCAPE
################################################################################
#
# Available files in your Figure1 folder:
# - BMI_lineage_genes.pdf
# - BMI_lineage_genes_scatter.pdf
# - DiD_PrePost_UpSet.pdf
# - GSEA_Gene_Heatmaps_DiD.pdf
# - GSEA_Gene_Heatmaps_STD.pdf
# - GSEA_Gene_Heatmaps_VSG.pdf
# - MixedModel_DEG_Heatmaps.pdf
# - std_PrePost_UpSet.pdf
# - VSG_PrePost_UpSet.pdf
#
# Based on what you showed me earlier:
# - Panel A: BMI_lineage_genes_scatter.pdf (the FDR scatter plot)
# - Panel B: MixedModel_DEG_Heatmaps.pdf (the large functional heatmap)
# - Panel C: DiD_PrePost_UpSet.pdf (gene overlap)
# - Panel D: GSEA_Gene_Heatmaps_DiD.pdf (pathway enrichment)
#
# Alternative panels for supplemental:
# - std_PrePost_UpSet.pdf
# - VSG_PrePost_UpSet.pdf
# - GSEA_Gene_Heatmaps_STD.pdf
# - GSEA_Gene_Heatmaps_VSG.pdf
################################################################################

make_figure1 <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    
    # Main panels - UPDATE THESE if my guesses are wrong
    panel_A = "BMI_lineage_genes_scatter.pdf",  # FDR scatter VSG vs Standard
    
    panel_B = "MixedModel_DEG_Heatmaps.pdf",     # Large functional heatmap
    panel_C = "DiD_PrePost_UpSet.pdf",           # UpSet plot
    panel_D = "GSEA_Gene_Heatmaps_DiD.pdf",      # DiD pathway enrichment
    
    output = "Figure1_Transcriptional_Landscape.pdf",
    width = 16,
    height = 16,
    density = 300
) {
  
  # Build full paths
  fig1_path <- file.path(base_path, "Figure1")
  
  path_A <- file.path(fig1_path, panel_A)
  path_B <- file.path(fig1_path, panel_B)
  path_C <- file.path(fig1_path, panel_C)
  path_D <- file.path(fig1_path, panel_D)
  
  # Check files exist
  for (p in c(path_A, path_B, path_C, path_D)) {
    if (!file.exists(p)) {
      warning(paste("File not found:", p))
    }
  }
  
  message("Loading panels...")
  A <- pdf_to_plot(path_A, density)
  B <- pdf_to_plot(path_B, density)
  C <- pdf_to_plot(path_C, density)
  D <- pdf_to_plot(path_D, density)
  
  message("Assembling figure...")
  
  # Layout Option 1: 2x2 grid
  # fig1 <- (A | B) / (C | D) +
  #   plot_annotation(tag_levels = 'A') &
  #   theme(plot.tag = element_text(size = 24, face = "bold"))
  
  # Layout Option 2: A small, B large on top; C and D on bottom
  fig1 <- (A + B + plot_layout(widths = c(1, 1.5))) /
    (C + D + plot_layout(widths = c(1, 1.5))) +
    plot_layout(heights = c(1, 1.2)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 24, face = "bold"))
  
  # Save
  output_path <- file.path(base_path, output)
  ggsave(output_path, fig1, width = width, height = height, dpi = 300)
  
  message("✓ Figure 1 saved: ", output_path)
  return(fig1)
}

# Alternative layouts you can try:

# Layout with 3 rows (if heatmap is very tall)
make_figure1_v2 <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"
) {
  
  fig1_path <- file.path(base_path, "Figure1")
  
  A <- pdf_to_plot(file.path(fig1_path, "BMI_lineage_genes_scatter.pdf"))
  B <- pdf_to_plot(file.path(fig1_path, "MixedModel_DEG_Heatmaps.pdf"))
  C <- pdf_to_plot(file.path(fig1_path, "DiD_PrePost_UpSet.pdf"))
  D <- pdf_to_plot(file.path(fig1_path, "GSEA_Gene_Heatmaps_DiD.pdf"))
  
  # A on top left, B takes more space
  # C below spanning width
  # D at bottom
  
  fig1 <- (A | B) / C / D +
    plot_layout(heights = c(1, 0.6, 1.2)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 24, face = "bold"))
  
  output_path <- file.path(base_path, "Figure1_v2.pdf")
  ggsave(output_path, fig1, width = 16, height = 18, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig1)
}

# Supplemental Figure S1: All UpSet plots
make_figure_s1_upset <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"
) {
  
  fig1_path <- file.path(base_path, "Figure1")
  
  A <- pdf_to_plot(file.path(fig1_path, "VSG_PrePost_UpSet.pdf"))
  B <- pdf_to_plot(file.path(fig1_path, "std_PrePost_UpSet.pdf"))
  C <- pdf_to_plot(file.path(fig1_path, "DiD_PrePost_UpSet.pdf"))
  
  fig_s1 <- A / B / C +
    plot_annotation(
      title = "Supplemental Figure S1: Gene Overlap Analysis",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  output_path <- file.path(base_path, "FigureS1_UpSet_plots.pdf")
  ggsave(output_path, fig_s1, width = 12, height = 14, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig_s1)
}

# Supplemental Figure S2: Individual pathway enrichments
make_figure_s2_pathways <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"
) {
  
  fig1_path <- file.path(base_path, "Figure1")
  
  A <- pdf_to_plot(file.path(fig1_path, "GSEA_Gene_Heatmaps_VSG.pdf"))
  B <- pdf_to_plot(file.path(fig1_path, "GSEA_Gene_Heatmaps_STD.pdf"))
  
  fig_s2 <- A / B +
    plot_annotation(
      title = "Supplemental Figure S2: Pathway Enrichment by Treatment",
      subtitle = "A) VSG Pre vs Post   B) Standard Medical Therapy Pre vs Post",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  output_path <- file.path(base_path, "FigureS2_Pathways.pdf")
  ggsave(output_path, fig_s2, width = 14, height = 16, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig_s2)
}

################################################################################
# RUN FIGURE 1
################################################################################

# Run this to generate Figure 1:
fig1 <- make_figure1()

# Or try the alternative layout:
# fig1_v2 <- make_figure1_v2()

# Generate supplemental figures:
# fig_s1 <- make_figure_s1_upset()
# fig_s2 <- make_figure_s2_pathways()

################################################################################
# TROUBLESHOOTING
################################################################################

# If you get an error, first check your files are there:
# list.files("C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output/Figure1")

# If panels look wrong, you might need to swap the file assignments.
# Tell me which PDF corresponds to which panel (from our earlier discussion):
# - PDF 1 (FDR scatter): Which file is this?
# - PDF 7 (functional heatmap): Which file is this?
# - PDF 9 (ranked genes): Which file is this?
# - PDF 4 (DiD pathways): Which file is this?











###### Figure 2

################################################################################
# FIGURE 2 ASSEMBLY - Cell-Cell Communication & Trajectories
# 
# Your file location: 
# C:\Users\netio\Documents\UofW\Projects\Long_Paper\VSG_analysis-main\VSG_analysis-main\output\Figure2
################################################################################

library(magick)
library(pdftools)
library(cowplot)
library(ggplot2)
library(patchwork)

# Set your base path
base_path <- "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"

# Helper function
pdf_to_plot <- function(path, density = 300) {
  if (!file.exists(path)) {
    stop(paste("File not found:", path))
  }
  img <- image_read_pdf(path, density = density)
  ggdraw() + draw_image(img)
}

################################################################################
# YOUR FIGURE 2 FILES:
################################################################################
#
# TRAJECTORY FILES:
# - EC_trajectory_median_pseudotime.pdf   (PDF 16 - boxplot)
# - EC_trajectory_ridge.pdf               (PDF 17 - density)
# - EC_trajectory_sling.pdf               (slingshot visualization)
# - PT_trajectory_median_pseudotime.pdf
# - PT_trajectory_ridge.pdf
# - PT_trajectory_sling.pdf
# - TAL_trajectory_median_pseudotime.pdf
# - TAL_trajectory_ridge.pdf
# - TAL_trajectory_sling.pdf
#
# CELLCHAT - VSG:
# - vsg_pre_fibremod.pdf      (PDF 14 - Pre-VSG fibroblast remodeling)
# - vsg_post_fibremod.pdf     (PDF 23 - Post-VSG fibroblast remodeling)
# - vsg_pre_fibtoendo.pdf     (fibroblast to endothelial)
# - vsg_post_fibtoendo.pdf
# - vsg_pre_immtoepi.pdf      (immune to epithelial)
# - vsg_post_immtoepi.pdf     (PDF 11 - Post-VSG immune to epithelial)
# - vsg_pre_immtoimm.pdf      (immune to immune)
# - vsg_post_immtoimm.pdf     (PDF 12 - Post-VSG immune to immune)
# - vsg_pre_tls.pdf           (PDF 15 - Pre-VSG TLS signature)
# - vsg_post_tls.pdf          (PDF 13 - Post-VSG TLS signature)
#
# CELLCHAT - SMT:
# - SMT_pre_fibremod.pdf
# - SMT_post_fibremod.pdf     (PDF 19 - Post-SMT fibroblast remodeling)
# - SMT_pre_fibtoendo.pdf
# - SMT_post_fibtoendo.pdf    (PDF 20 - Post-SMT fibroblast to EC - large!)
# - SMT_pre_immtoepi.pdf
# - SMT_post_immtoepi.pdf     (PDF 21 - Post-SMT immune to epithelial)
# - SMT_pre_immtoimm.pdf
# - SMT_post_immtoimm.pdf
# - SMT_pre_tls.pdf
#
# MILO:
# - Milo_DA_combined_plot.pdf
# - Milo_DA_heatmap_celltype_vs_test.pdf
################################################################################

################################################################################
# MAIN FIGURE 2: Communication & Trajectory Overview
################################################################################

make_figure2 <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    output = "Figure2_CellCell_Communication.pdf",
    width = 18,
    
    height = 20,
    density = 300
) {
  
  fig2_path <- file.path(base_path, "Figure2")
  
  message("Loading panels...")
  
  # Row A: EC Trajectory
  A1 <- pdf_to_plot(file.path(fig2_path, "EC_trajectory_median_pseudotime.pdf"), density)
  A2 <- pdf_to_plot(file.path(fig2_path, "EC_trajectory_ridge.pdf"), density)
  
  # Row B: Fibroblast Remodeling (Pre-VSG | Post-VSG | Post-SMT)
  B1 <- pdf_to_plot(file.path(fig2_path, "vsg_pre_fibremod.pdf"), density)
  B2 <- pdf_to_plot(file.path(fig2_path, "vsg_post_fibremod.pdf"), density)
  B3 <- pdf_to_plot(file.path(fig2_path, "SMT_post_fibremod.pdf"), density)
  
  # Row C: Immune to Epithelial (Post-VSG | Post-SMT)
  C1 <- pdf_to_plot(file.path(fig2_path, "vsg_post_immtoepi.pdf"), density)
  C2 <- pdf_to_plot(file.path(fig2_path, "SMT_post_immtoepi.pdf"), density)
  
  # Row D: TLS Signatures (Pre-VSG | Post-VSG)
  D1 <- pdf_to_plot(file.path(fig2_path, "vsg_pre_tls.pdf"), density)
  D2 <- pdf_to_plot(file.path(fig2_path, "vsg_post_tls.pdf"), density)
  
  message("Assembling figure...")
  
  # Layout:
  # Row A: EC trajectory (boxplot | ridge density)
  # Row B: Fibroblast remodeling comparison (3 panels)
  # Row C: Immune-epithelial communication (2 panels)
  # Row D: TLS signatures (2 panels)
  
  row_A <- (A1 | A2) + plot_layout(widths = c(1, 2.5))
  row_B <- (B1 | B2 | B3) + plot_layout(widths = c(1, 0.7, 1))
  row_C <- (C1 | C2) + plot_layout(widths = c(1, 1.2))
  row_D <- (D1 | D2) + plot_layout(widths = c(1.5, 1))
  
  fig2 <- row_A / row_B / row_C / row_D +
    plot_layout(heights = c(0.8, 1.2, 1.2, 0.8)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 24, face = "bold"))
  
  # Save
  output_path <- file.path(base_path, output)
  ggsave(output_path, fig2, width = width, height = height, dpi = 300)
  
  message("✓ Figure 2 saved: ", output_path)
  return(fig2)
}

################################################################################
# ALTERNATIVE: Figure 2 with all trajectories (PT, TAL, EC)
################################################################################

make_figure2_all_trajectories <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    output = "Figure2_with_all_trajectories.pdf",
    width = 18,
    height = 22
) {
  
  fig2_path <- file.path(base_path, "Figure2")
  
  # Trajectories
  EC_box <- pdf_to_plot(file.path(fig2_path, "EC_trajectory_median_pseudotime.pdf"))
  EC_ridge <- pdf_to_plot(file.path(fig2_path, "EC_trajectory_ridge.pdf"))
  PT_box <- pdf_to_plot(file.path(fig2_path, "PT_trajectory_median_pseudotime.pdf"))
  TAL_box <- pdf_to_plot(file.path(fig2_path, "TAL_trajectory_median_pseudotime.pdf"))
  
  # Fibroblast
  fib_pre_vsg <- pdf_to_plot(file.path(fig2_path, "vsg_pre_fibremod.pdf"))
  fib_post_vsg <- pdf_to_plot(file.path(fig2_path, "vsg_post_fibremod.pdf"))
  fib_post_smt <- pdf_to_plot(file.path(fig2_path, "SMT_post_fibremod.pdf"))
  
  # Immune
  imm_post_vsg <- pdf_to_plot(file.path(fig2_path, "vsg_post_immtoepi.pdf"))
  imm_post_smt <- pdf_to_plot(file.path(fig2_path, "SMT_post_immtoepi.pdf"))
  
  # Layout
  row1 <- (PT_box | TAL_box | EC_box) + plot_layout(widths = c(1, 1, 1))
  row2 <- EC_ridge
  row3 <- (fib_pre_vsg | fib_post_vsg | fib_post_smt)
  row4 <- (imm_post_vsg | imm_post_smt)
  
  fig2 <- row1 / row2 / row3 / row4 +
    plot_layout(heights = c(0.7, 1, 1.2, 1.2)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 24, face = "bold"))
  
  output_path <- file.path(base_path, output)
  ggsave(output_path, fig2, width = width, height = height, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig2)
}

################################################################################
# ALTERNATIVE: Figure 2 with Milo
################################################################################

make_figure2_with_milo <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    output = "Figure2_with_Milo.pdf",
    width = 18,
    height = 22
) {
  
  fig2_path <- file.path(base_path, "Figure2")
  
  # EC Trajectory
  A1 <- pdf_to_plot(file.path(fig2_path, "EC_trajectory_median_pseudotime.pdf"))
  A2 <- pdf_to_plot(file.path(fig2_path, "EC_trajectory_ridge.pdf"))
  
  # Milo
  B <- pdf_to_plot(file.path(fig2_path, "Milo_DA_combined_plot.pdf"))
  
  # Fibroblast
  C1 <- pdf_to_plot(file.path(fig2_path, "vsg_pre_fibremod.pdf"))
  C2 <- pdf_to_plot(file.path(fig2_path, "vsg_post_fibremod.pdf"))
  C3 <- pdf_to_plot(file.path(fig2_path, "SMT_post_fibremod.pdf"))
  
  # Immune
  D1 <- pdf_to_plot(file.path(fig2_path, "vsg_post_immtoepi.pdf"))
  D2 <- pdf_to_plot(file.path(fig2_path, "SMT_post_immtoepi.pdf"))
  
  fig2 <- ((A1 | A2) + plot_layout(widths = c(1, 2))) /
    B /
    (C1 | C2 | C3) /
    (D1 | D2) +
    plot_layout(heights = c(0.8, 1, 1.2, 1.2)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 24, face = "bold"))
  
  output_path <- file.path(base_path, output)
  ggsave(output_path, fig2, width = width, height = height, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig2)
}

################################################################################
# SUPPLEMENTAL FIGURES
################################################################################

# Supp: All trajectory plots
make_supp_trajectories <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"
) {
  
  fig2_path <- file.path(base_path, "Figure2")
  
  # PT
  PT_sling <- pdf_to_plot(file.path(fig2_path, "PT_trajectory_sling.pdf"))
  PT_box <- pdf_to_plot(file.path(fig2_path, "PT_trajectory_median_pseudotime.pdf"))
  PT_ridge <- pdf_to_plot(file.path(fig2_path, "PT_trajectory_ridge.pdf"))
  
  # TAL
  TAL_sling <- pdf_to_plot(file.path(fig2_path, "TAL_trajectory_sling.pdf"))
  TAL_box <- pdf_to_plot(file.path(fig2_path, "TAL_trajectory_median_pseudotime.pdf"))
  TAL_ridge <- pdf_to_plot(file.path(fig2_path, "TAL_trajectory_ridge.pdf"))
  
  # EC
  EC_sling <- pdf_to_plot(file.path(fig2_path, "EC_trajectory_sling.pdf"))
  EC_box <- pdf_to_plot(file.path(fig2_path, "EC_trajectory_median_pseudotime.pdf"))
  EC_ridge <- pdf_to_plot(file.path(fig2_path, "EC_trajectory_ridge.pdf"))
  
  # 3x3 grid: rows = PT/TAL/EC, cols = sling/boxplot/ridge
  fig <- (PT_sling | PT_box | PT_ridge) /
    (TAL_sling | TAL_box | TAL_ridge) /
    (EC_sling | EC_box | EC_ridge) +
    plot_annotation(
      title = "Supplemental: Pseudotemporal Trajectory Analysis",
      subtitle = "Columns: Slingshot embedding | Median pseudotime | Ridge density",
      tag_levels = list(c("A", "", "", "B", "", "", "C", "", ""))
    ) &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  output_path <- file.path(base_path, "FigureS_Trajectories.pdf")
  ggsave(output_path, fig, width = 18, height = 16, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

# Supp: Complete CellChat comparison (2x2 for each interaction type)
make_supp_cellchat_fibroblast <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"
) {
  
  fig2_path <- file.path(base_path, "Figure2")
  
  # 2x2: Pre/Post × VSG/SMT
  vsg_pre <- pdf_to_plot(file.path(fig2_path, "vsg_pre_fibremod.pdf"))
  vsg_post <- pdf_to_plot(file.path(fig2_path, "vsg_post_fibremod.pdf"))
  smt_pre <- pdf_to_plot(file.path(fig2_path, "SMT_pre_fibremod.pdf"))
  smt_post <- pdf_to_plot(file.path(fig2_path, "SMT_post_fibremod.pdf"))
  
  fig <- (vsg_pre | vsg_post) / (smt_pre | smt_post) +
    plot_annotation(
      title = "Supplemental: Fibroblast/Perivascular Remodeling",
      subtitle = "Top: VSG (Pre | Post)   Bottom: SMT (Pre | Post)",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  output_path <- file.path(base_path, "FigureS_Fibroblast_Remodeling.pdf")
  ggsave(output_path, fig, width = 14, height = 14, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

make_supp_cellchat_immune <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"
) {
  
  fig2_path <- file.path(base_path, "Figure2")
  
  # Immune to Epithelial - 2x2
  vsg_pre <- pdf_to_plot(file.path(fig2_path, "vsg_pre_immtoepi.pdf"))
  vsg_post <- pdf_to_plot(file.path(fig2_path, "vsg_post_immtoepi.pdf"))
  smt_pre <- pdf_to_plot(file.path(fig2_path, "SMT_pre_immtoepi.pdf"))
  smt_post <- pdf_to_plot(file.path(fig2_path, "SMT_post_immtoepi.pdf"))
  
  fig <- (vsg_pre | vsg_post) / (smt_pre | smt_post) +
    plot_annotation(
      title = "Supplemental: Immune to Epithelial/Endothelial Communication",
      subtitle = "Top: VSG (Pre | Post)   Bottom: SMT (Pre | Post)",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  output_path <- file.path(base_path, "FigureS_Immune_to_Epithelial.pdf")
  ggsave(output_path, fig, width = 16, height = 14, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

make_supp_cellchat_tls <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"
) {
  
  fig2_path <- file.path(base_path, "Figure2")
  
  # TLS - VSG Pre/Post + SMT Pre (no SMT post TLS file)
  vsg_pre <- pdf_to_plot(file.path(fig2_path, "vsg_pre_tls.pdf"))
  vsg_post <- pdf_to_plot(file.path(fig2_path, "vsg_post_tls.pdf"))
  smt_pre <- pdf_to_plot(file.path(fig2_path, "SMT_pre_tls.pdf"))
  
  fig <- (vsg_pre | vsg_post) / (smt_pre | plot_spacer()) +
    plot_annotation(
      title = "Supplemental: TLS-like Signatures",
      subtitle = "Top: VSG (Pre | Post)   Bottom: SMT Pre",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  output_path <- file.path(base_path, "FigureS_TLS_Signatures.pdf")
  ggsave(output_path, fig, width = 16, height = 12, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

################################################################################
# RUN FIGURE 2
################################################################################

# Main Figure 2
fig2 <- make_figure2()

# Alternative versions:
# fig2_traj <- make_figure2_all_trajectories()
# fig2_milo <- make_figure2_with_milo()

# Supplemental figures:
# supp_traj <- make_supp_trajectories()
# supp_fib <- make_supp_cellchat_fibroblast()
# supp_imm <- make_supp_cellchat_immune()
# supp_tls <- make_supp_cellchat_tls()





#########Figure 3

################################################################################
# FIGURE 3 ASSEMBLY - hdWGCNA Modules & Transcription Factors
# 
# Your file location: 
# C:\Users\netio\Documents\UofW\Projects\Long_Paper\VSG_analysis-main\VSG_analysis-main\output\Figure3\scRNA
################################################################################

library(magick)
library(pdftools)
library(cowplot)
library(ggplot2)
library(patchwork)

# Set your base path
base_path <- "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"

# Helper function
pdf_to_plot <- function(path, density = 300) {
  if (!file.exists(path)) {
    stop(paste("File not found:", path))
  }
  img <- image_read_pdf(path, density = density)
  ggdraw() + draw_image(img)
}

################################################################################
# YOUR FIGURE 3 FILES (in scRNA subfolder):
################################################################################
#
# DENDROGRAMS:
# - Dendrogram_hdWGCNA_EC.pdf      (PDF 24)
# - Dendrogram_hdWGCNA_Immune.pdf  (PDF 25)
# - Dendrogram_hdWGCNA_PT.pdf      (PDF 26)
# - Dendrogram_hdWGCNA_TAL.pdf     (PDF 27)
#
# DME HEATMAPS (Module × Subtype × Condition comparisons):
# - DME_Heatmap_EC.pdf             (PDF 28)
# - DME_Heatmap_Immune.pdf         (PDF 29)
# - DME_Heatmap_PT.pdf             (PDF 30)
# - DME_Heatmap_TAL.pdf            (PDF 31)
#
# MODULE-TRAIT CORRELATIONS:
# - module_trait_correlation_EC.pdf      (PDF 43)
# - module_trait_correlation_Immune.pdf  (PDF 44)
# - module_trait_correlation_PT.pdf      (PDF 45)
# - module_trait_correlation_TAL.pdf     (PDF 46)
#
# TF REGULON - BAR PLOTS:
# - EC_RegulonBarPlots_topTFs.pdf        (PDF 32)
# - Immune_RegulonBarPlots_topTFs.pdf    (PDF 51)
# - PT_RegulonBarPlots_topTFs.pdf        (PDF 49)
# - TAL_RegulonBarPlots_topTFs.pdf       (PDF 63)
#
# TF REGULON - NETWORKS:
# - EC_TF_Network_Pos_Neg_Regulons.pdf   (PDF 33)
# - Immune_TF_Network_Pos_Neg_Regulons.pdf (PDF 52)
# - PT_TF_Network_Pos_Neg_Regulons.pdf   (PDF 50)
# - TAL_TF_Network_Pos_Neg_Regulons.pdf  (PDF 64)
#
# TF REGULON - VOLCANO PLOTS:
# - EC_TF_volcano_plot.pdf               (PDF 34)
# - Immune_TF_volcano_plot.pdf           (PDF 53)
# - PT_TF_volcano_plot.pdf               (PDF 58)
# - TAL_TF_volcano_plot.pdf              (PDF 65)
#
# PATHWAY ENRICHMENT:
# - Enrichr_EC.pdf                       (PDF 35)
# - Enrichr_Immune.pdf                   (PDF 36)
# - Enrichr_PT.pdf                       (PDF 37)
# - Enrichr_TAL.pdf                      (PDF 38)
#
# KME PLOTS (hub genes):
# - KMEs_plot_EC.pdf                     (PDF 54)
# - KMEs_plot_Immune.pdf                 (PDF 55)
# - KMEs_plot_PT.pdf                     (PDF 56)
# - KMEs_plot_TAL.pdf                    (PDF 57)
#
# RADAR PLOTS (module composition):
# - Radar_Plot_EC.pdf                    (PDF 59)
# - Radar_Plot_Immune.pdf                (PDF 60)
# - Radar_Plot_PT.pdf                    (PDF 61)
# - Radar_Plot_TAL.pdf                   (PDF 62)
#
# PATIENT-LEVEL:
# - patient_pca_trajectory.pdf           (PDF 47)
# - PCA_sensitivity_HC_T2D_PostVSG.pdf   (PDF 48)
# - lineage_delta_mahalanobis_pc5.pdf    (PDF 39)
# - lineage_delta_mahalanobis_pc10.pdf   (PDF 40)
# - lineage_delta_mahalanobis_pc20.pdf   (PDF 41)
# - lineage_delta_mahalanobis_pc30.pdf   (PDF 42)
################################################################################

################################################################################
# MAIN FIGURE 3: hdWGCNA Overview
################################################################################

make_figure3 <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    output = "Figure3_hdWGCNA.pdf",
    width = 20,
    height = 22,
    density = 300
) {
  
  fig3_path <- file.path(base_path, "Figure3", "scRNA")
  
  message("Loading panels...")
  
  # Panel A: Patient PCA trajectory
  A <- pdf_to_plot(file.path(fig3_path, "patient_pca_trajectory.pdf"), density)
  
  # Panels B-E: DME Heatmaps (module comparisons)
  B <- pdf_to_plot(file.path(fig3_path, "DME_Heatmap_PT.pdf"), density)
  C <- pdf_to_plot(file.path(fig3_path, "DME_Heatmap_TAL.pdf"), density)
  D <- pdf_to_plot(file.path(fig3_path, "DME_Heatmap_EC.pdf"), density)
  E <- pdf_to_plot(file.path(fig3_path, "DME_Heatmap_Immune.pdf"), density)
  
  # Panels F-I: Module-trait correlations
  F <- pdf_to_plot(file.path(fig3_path, "module_trait_correlation_PT.pdf"), density)
  G <- pdf_to_plot(file.path(fig3_path, "module_trait_correlation_TAL.pdf"), density)
  H <- pdf_to_plot(file.path(fig3_path, "module_trait_correlation_EC.pdf"), density)
  I <- pdf_to_plot(file.path(fig3_path, "module_trait_correlation_Immune.pdf"), density)
  
  message("Assembling figure...")
  
  # Layout:
  # Row 1: A (PCA trajectory - full width)
  # Row 2: B | C (PT and TAL DME heatmaps)
  # Row 3: D | E (EC and Immune DME heatmaps)
  # Row 4: F | G | H | I (all trait correlations)
  
  fig3 <- A /
    (B | C) /
    (D | E) /
    (F | G | H | I) +
    plot_layout(heights = c(0.8, 1.3, 1.5, 0.7)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 24, face = "bold"))
  
  output_path <- file.path(base_path, output)
  ggsave(output_path, fig3, width = width, height = height, dpi = 300)
  
  message("✓ Figure 3 saved: ", output_path)
  return(fig3)
}

################################################################################
# FIGURE 3 VERSION 2: With TF Networks
################################################################################

make_figure3_with_TFs <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    output = "Figure3_with_TF_networks.pdf",
    width = 18,
    height = 26,
    density = 300
) {
  
  fig3_path <- file.path(base_path, "Figure3", "scRNA")
  
  message("Loading panels...")
  
  # Patient trajectory
  A <- pdf_to_plot(file.path(fig3_path, "patient_pca_trajectory.pdf"), density)
  
  # DME Heatmaps
  B <- pdf_to_plot(file.path(fig3_path, "DME_Heatmap_PT.pdf"), density)
  C <- pdf_to_plot(file.path(fig3_path, "DME_Heatmap_TAL.pdf"), density)
  D <- pdf_to_plot(file.path(fig3_path, "DME_Heatmap_EC.pdf"), density)
  E <- pdf_to_plot(file.path(fig3_path, "DME_Heatmap_Immune.pdf"), density)
  
  # Module-trait correlations (select 2)
  F <- pdf_to_plot(file.path(fig3_path, "module_trait_correlation_PT.pdf"), density)
  G <- pdf_to_plot(file.path(fig3_path, "module_trait_correlation_Immune.pdf"), density)
  
  # TF Networks (select 2)
  H <- pdf_to_plot(file.path(fig3_path, "PT_TF_Network_Pos_Neg_Regulons.pdf"), density)
  I <- pdf_to_plot(file.path(fig3_path, "Immune_TF_Network_Pos_Neg_Regulons.pdf"), density)
  
  message("Assembling figure...")
  
  fig3 <- A /
    (B | C) /
    (D | E) /
    (F | G) /
    (H | I) +
    plot_layout(heights = c(0.7, 1.2, 1.4, 0.6, 1.2)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 24, face = "bold"))
  
  output_path <- file.path(base_path, output)
  ggsave(output_path, fig3, width = width, height = height, dpi = 300)
  
  message("✓ Figure 3 saved: ", output_path)
  return(fig3)
}

################################################################################
# FIGURE 3 VERSION 3: Compact (for journals with space limits)
################################################################################

make_figure3_compact <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    output = "Figure3_compact.pdf",
    width = 16,
    height = 16,
    density = 300
) {
  
  fig3_path <- file.path(base_path, "Figure3", "scRNA")
  
  # Minimal set
  A <- pdf_to_plot(file.path(fig3_path, "patient_pca_trajectory.pdf"), density)
  B <- pdf_to_plot(file.path(fig3_path, "DME_Heatmap_PT.pdf"), density)
  C <- pdf_to_plot(file.path(fig3_path, "DME_Heatmap_Immune.pdf"), density)
  D <- pdf_to_plot(file.path(fig3_path, "module_trait_correlation_PT.pdf"), density)
  E <- pdf_to_plot(file.path(fig3_path, "module_trait_correlation_Immune.pdf"), density)
  
  fig3 <- A /
    (B | C) /
    (D | E) +
    plot_layout(heights = c(0.8, 1.4, 0.7)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 24, face = "bold"))
  
  output_path <- file.path(base_path, output)
  ggsave(output_path, fig3, width = width, height = height, dpi = 300)
  
  message("✓ Figure 3 (compact) saved: ", output_path)
  return(fig3)
}

################################################################################
# SUPPLEMENTAL FIGURES FOR FIGURE 3
################################################################################

# Supp: All Dendrograms
make_supp_dendrograms <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"
) {
  
  fig3_path <- file.path(base_path, "Figure3", "scRNA")
  
  A <- pdf_to_plot(file.path(fig3_path, "Dendrogram_hdWGCNA_PT.pdf"))
  B <- pdf_to_plot(file.path(fig3_path, "Dendrogram_hdWGCNA_TAL.pdf"))
  C <- pdf_to_plot(file.path(fig3_path, "Dendrogram_hdWGCNA_EC.pdf"))
  D <- pdf_to_plot(file.path(fig3_path, "Dendrogram_hdWGCNA_Immune.pdf"))
  
  fig <- (A | B) / (C | D) +
    plot_annotation(
      title = "Supplemental: hdWGCNA Module Dendrograms",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  output_path <- file.path(base_path, "FigureS_Dendrograms.pdf")
  ggsave(output_path, fig, width = 14, height = 10, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

# Supp: All Pathway Enrichments
make_supp_enrichr <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"
) {
  
  fig3_path <- file.path(base_path, "Figure3", "scRNA")
  
  A <- pdf_to_plot(file.path(fig3_path, "Enrichr_PT.pdf"))
  B <- pdf_to_plot(file.path(fig3_path, "Enrichr_TAL.pdf"))
  C <- pdf_to_plot(file.path(fig3_path, "Enrichr_EC.pdf"))
  D <- pdf_to_plot(file.path(fig3_path, "Enrichr_Immune.pdf"))
  
  fig <- A / B / C / D +
    plot_annotation(
      title = "Supplemental: Module Pathway Enrichment",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  output_path <- file.path(base_path, "FigureS_Enrichr.pdf")
  ggsave(output_path, fig, width = 14, height = 20, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

# Supp: All TF Volcano Plots
make_supp_tf_volcanos <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"
) {
  
  fig3_path <- file.path(base_path, "Figure3", "scRNA")
  
  A <- pdf_to_plot(file.path(fig3_path, "PT_TF_volcano_plot.pdf"))
  B <- pdf_to_plot(file.path(fig3_path, "TAL_TF_volcano_plot.pdf"))
  C <- pdf_to_plot(file.path(fig3_path, "EC_TF_volcano_plot.pdf"))
  D <- pdf_to_plot(file.path(fig3_path, "Immune_TF_volcano_plot.pdf"))
  
  fig <- (A | B) / (C | D) +
    plot_annotation(
      title = "Supplemental: Transcription Factor Regulon Changes",
      subtitle = "Volcano plots showing Pre→Post changes (SMT vs VSG)",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  output_path <- file.path(base_path, "FigureS_TF_Volcanos.pdf")
  ggsave(output_path, fig, width = 16, height = 14, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

# Supp: All TF Networks
make_supp_tf_networks <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"
) {
  
  fig3_path <- file.path(base_path, "Figure3", "scRNA")
  
  A <- pdf_to_plot(file.path(fig3_path, "PT_TF_Network_Pos_Neg_Regulons.pdf"))
  B <- pdf_to_plot(file.path(fig3_path, "TAL_TF_Network_Pos_Neg_Regulons.pdf"))
  C <- pdf_to_plot(file.path(fig3_path, "EC_TF_Network_Pos_Neg_Regulons.pdf"))
  D <- pdf_to_plot(file.path(fig3_path, "Immune_TF_Network_Pos_Neg_Regulons.pdf"))
  
  fig <- (A | B) / (C | D) +
    plot_annotation(
      title = "Supplemental: Transcription Factor Networks",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  output_path <- file.path(base_path, "FigureS_TF_Networks.pdf")
  ggsave(output_path, fig, width = 18, height = 16, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

# Supp: Mahalanobis Distance (all PC versions)
make_supp_mahalanobis <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"
) {
  
  fig3_path <- file.path(base_path, "Figure3", "scRNA")
  
  A <- pdf_to_plot(file.path(fig3_path, "lineage_delta_mahalanobis_pc5.pdf"))
  B <- pdf_to_plot(file.path(fig3_path, "lineage_delta_mahalanobis_pc10.pdf"))
  C <- pdf_to_plot(file.path(fig3_path, "lineage_delta_mahalanobis_pc20.pdf"))
  D <- pdf_to_plot(file.path(fig3_path, "lineage_delta_mahalanobis_pc30.pdf"))
  
  fig <- (A | B) / (C | D) +
    plot_annotation(
      title = "Supplemental: Mahalanobis Distance (Post - Pre)",
      subtitle = "Sensitivity analysis across PC dimensions (5, 10, 20, 30)",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  output_path <- file.path(base_path, "FigureS_Mahalanobis.pdf")
  ggsave(output_path, fig, width = 14, height = 12, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

# Supp: KME Hub Genes
make_supp_kme <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"
) {
  
  fig3_path <- file.path(base_path, "Figure3", "scRNA")
  
  A <- pdf_to_plot(file.path(fig3_path, "KMEs_plot_PT.pdf"))
  B <- pdf_to_plot(file.path(fig3_path, "KMEs_plot_TAL.pdf"))
  C <- pdf_to_plot(file.path(fig3_path, "KMEs_plot_EC.pdf"))
  D <- pdf_to_plot(file.path(fig3_path, "KMEs_plot_Immune.pdf"))
  
  fig <- A / B / C / D +
    plot_annotation(
      title = "Supplemental: Module Hub Genes (kME)",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  output_path <- file.path(base_path, "FigureS_KME_HubGenes.pdf")
  ggsave(output_path, fig, width = 14, height = 20, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

# Supp: Radar Plots (module composition)
make_supp_radar <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"
) {
  
  fig3_path <- file.path(base_path, "Figure3", "scRNA")
  
  A <- pdf_to_plot(file.path(fig3_path, "Radar_Plot_PT.pdf"))
  B <- pdf_to_plot(file.path(fig3_path, "Radar_Plot_TAL.pdf"))
  C <- pdf_to_plot(file.path(fig3_path, "Radar_Plot_EC.pdf"))
  D <- pdf_to_plot(file.path(fig3_path, "Radar_Plot_Immune.pdf"))
  
  fig <- (A | B) / (C | D) +
    plot_annotation(
      title = "Supplemental: Module Composition by Cell Subtype",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  output_path <- file.path(base_path, "FigureS_Radar_ModuleComposition.pdf")
  ggsave(output_path, fig, width = 14, height = 12, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

################################################################################
# RUN FIGURE 3
################################################################################

# Main Figure 3
fig3 <- make_figure3()

# Alternative versions:
# fig3_tf <- make_figure3_with_TFs()
# fig3_compact <- make_figure3_compact()

# Supplemental figures:
# supp_dendro <- make_supp_dendrograms()
# supp_enrichr <- make_supp_enrichr()
# supp_volcano <- make_supp_tf_volcanos()
# supp_network <- make_supp_tf_networks()
# supp_mahal <- make_supp_mahalanobis()
# supp_kme <- make_supp_kme()
# supp_radar <- make_supp_radar()




####################Figure 4

################################################################################
# FIGURE 4 ASSEMBLY - Spatial Transcriptomics
# 
# Your file location: 
# C:\Users\netio\Documents\UofW\Projects\Long_Paper\VSG_analysis-main\VSG_analysis-main\output\Figure4
################################################################################

library(magick)
library(pdftools)
library(cowplot)
library(ggplot2)
library(patchwork)

# Set your base path
base_path <- "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"

# Helper function
pdf_to_plot <- function(path, density = 300) {
  if (!file.exists(path)) {
    warning(paste("File not found:", path))
    return(ggplot() + theme_void() + 
             annotate("text", x = 0.5, y = 0.5, label = "Missing", size = 4))
  }
  img <- image_read_pdf(path, density = density)
  ggdraw() + draw_image(img)
}

################################################################################
# YOUR FIGURE 4 FILES:
################################################################################
#
# ROOT FILES:
# - module_score_correlation_heatmap.pdf  (PDF 66 - signature-trait correlations)
# - TLS_epi_heatmap.pdf                   (PDF 67 - co-localization matrix)
# - TLS_heatmap.pdf                       (PDF 68 - another co-localization)
#
# SUBFOLDERS with individual sample data:
#
# /Deconvolution/ - Cell type pie charts per sample
#   VSG patients with Pre/Post pairs:
#   - IT_07 (RH-59-T)_baseline_pie.pdf
#   - IT_07 (RH-59-T)_12_months_post_surgery_pie.pdf
#   - IT_08 (RH-60-T)_baseline_pie.pdf
#   - IT_08 (RH-60-T)_12_months_post_surgery_pie.pdf
#   - IT_11_baseline_pie.pdf
#   - IT_11_12_months_post_surgery_pie.pdf
#   - IT_12_baseline_pie.pdf
#   - IT_12_12_months_post_surgery_pie.pdf
#   - IT_14_baseline_pie.pdf
#   - IT_14_12_months_post_surgery_pie.pdf
#
#   SMT/Other patients (baseline only):
#   - RH-23-T_baseline_pie.pdf
#   - RH-50-T_baseline_pie.pdf
#   - RH-63-T_baseline_pie.pdf
#   - CRC-* files (healthy controls?)
#
# /Module_score/ - Immune signature spatial maps per sample
#   (same sample naming as Deconvolution)
#
# /TLS_score/ - TLS signature spatial maps per sample
#   (same sample naming as Deconvolution)
################################################################################

################################################################################
# MAIN FIGURE 4: Spatial Overview
################################################################################

make_figure4 <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    output = "Figure4_Spatial.pdf",
    width = 18,
    height = 16,
    density = 300
) {
  
  fig4_path <- file.path(base_path, "Figure4")
  deconv_path <- file.path(fig4_path, "Deconvolution")
  module_path <- file.path(fig4_path, "Module_score")
  tls_path <- file.path(fig4_path, "TLS_score")
  
  message("Loading panels...")
  
  # Panel A-B: Representative deconvolution (Pre vs Post VSG)
  # Using IT_11 as example (has both timepoints)
  A <- pdf_to_plot(file.path(deconv_path, "IT_11_baseline_pie.pdf"), density)
  B <- pdf_to_plot(file.path(deconv_path, "IT_11_12_months_post_surgery_pie.pdf"), density)
  
  # Panel C-D: Module scores (same patient)
  C <- pdf_to_plot(file.path(module_path, "IT_11_baseline_module_scores.pdf"), density)
  D <- pdf_to_plot(file.path(module_path, "IT_11_12_months_post_surgery_module_scores.pdf"), density)
  
  # Panel E-F: TLS scores (same patient)
  E <- pdf_to_plot(file.path(tls_path, "IT_11_baseline_TLS_scores.pdf"), density)
  F <- pdf_to_plot(file.path(tls_path, "IT_11_12_months_post_surgery_TLS_scores.pdf"), density)
  
  # Panel G: Signature-trait correlations
  G <- pdf_to_plot(file.path(fig4_path, "module_score_correlation_heatmap.pdf"), density)
  
  # Panel H: Co-localization
  H <- pdf_to_plot(file.path(fig4_path, "TLS_epi_heatmap.pdf"), density)
  
  message("Assembling figure...")
  
  # Layout:
  # Row 1: A | B (Deconvolution Pre vs Post)
  # Row 2: C | D (Module scores Pre vs Post)
  # Row 3: E | F (TLS Pre vs Post)
  # Row 4: G | H (Correlations | Co-localization)
  
  fig4 <- (A | B) / (C | D) / (E | F) / (G | H) +
    plot_layout(heights = c(1, 1, 1, 0.8)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 24, face = "bold"))
  
  output_path <- file.path(base_path, output)
  ggsave(output_path, fig4, width = width, height = height, dpi = 300)
  
  message("✓ Figure 4 saved: ", output_path)
  return(fig4)
}

################################################################################
# FIGURE 4 VERSION 2: Multiple patients comparison
################################################################################

make_figure4_multi_patient <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    output = "Figure4_multi_patient.pdf",
    width = 20,
    height = 18,
    density = 300
) {
  
  fig4_path <- file.path(base_path, "Figure4")
  deconv_path <- file.path(fig4_path, "Deconvolution")
  module_path <- file.path(fig4_path, "Module_score")
  
  message("Loading panels...")
  
  # Row 1: IT_07 Pre vs Post
  A1 <- pdf_to_plot(file.path(deconv_path, "IT_07 (RH-59-T)_baseline_pie.pdf"), density)
  A2 <- pdf_to_plot(file.path(deconv_path, "IT_07 (RH-59-T)_12_months_post_surgery_pie.pdf"), density)
  
  # Row 2: IT_11 Pre vs Post
  B1 <- pdf_to_plot(file.path(deconv_path, "IT_11_baseline_pie.pdf"), density)
  B2 <- pdf_to_plot(file.path(deconv_path, "IT_11_12_months_post_surgery_pie.pdf"), density)
  
  # Row 3: IT_12 Pre vs Post
  C1 <- pdf_to_plot(file.path(deconv_path, "IT_12_baseline_pie.pdf"), density)
  C2 <- pdf_to_plot(file.path(deconv_path, "IT_12_12_months_post_surgery_pie.pdf"), density)
  
  # Row 4: Correlation heatmap
  D <- pdf_to_plot(file.path(fig4_path, "module_score_correlation_heatmap.pdf"), density)
  
  message("Assembling figure...")
  
  # Layout with row labels
  row1 <- (A1 | A2) + plot_annotation(subtitle = "Patient IT_07")
  row2 <- (B1 | B2) + plot_annotation(subtitle = "Patient IT_11")
  row3 <- (C1 | C2) + plot_annotation(subtitle = "Patient IT_12")
  
  fig4 <- (A1 | A2) / (B1 | B2) / (C1 | C2) / D +
    plot_layout(heights = c(1, 1, 1, 0.8)) +
    plot_annotation(
      title = "Spatial Deconvolution: Pre-VSG (left) vs Post-VSG (right)",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 24, face = "bold"))
  
  output_path <- file.path(base_path, output)
  ggsave(output_path, fig4, width = width, height = height, dpi = 300)
  
  message("✓ Figure 4 saved: ", output_path)
  return(fig4)
}

################################################################################
# FIGURE 4 VERSION 3: Compact with key findings
################################################################################

make_figure4_compact <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    output = "Figure4_compact.pdf",
    width = 16,
    height = 14,
    density = 300
) {
  
  fig4_path <- file.path(base_path, "Figure4")
  deconv_path <- file.path(fig4_path, "Deconvolution")
  module_path <- file.path(fig4_path, "Module_score")
  tls_path <- file.path(fig4_path, "TLS_score")
  
  # Select one representative patient
  A <- pdf_to_plot(file.path(deconv_path, "IT_11_baseline_pie.pdf"), density)
  B <- pdf_to_plot(file.path(deconv_path, "IT_11_12_months_post_surgery_pie.pdf"), density)
  C <- pdf_to_plot(file.path(module_path, "IT_11_baseline_module_scores.pdf"), density)
  D <- pdf_to_plot(file.path(module_path, "IT_11_12_months_post_surgery_module_scores.pdf"), density)
  E <- pdf_to_plot(file.path(fig4_path, "module_score_correlation_heatmap.pdf"), density)
  F <- pdf_to_plot(file.path(fig4_path, "TLS_epi_heatmap.pdf"), density)
  
  fig4 <- (A | B) / (C | D) / (E | F) +
    plot_layout(heights = c(1, 1.2, 0.8)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 24, face = "bold"))
  
  output_path <- file.path(base_path, output)
  ggsave(output_path, fig4, width = width, height = height, dpi = 300)
  
  message("✓ Figure 4 (compact) saved: ", output_path)
  return(fig4)
}

################################################################################
# SUPPLEMENTAL: All samples deconvolution
################################################################################

make_supp_deconv_vsg <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    density = 250
) {
  
  deconv_path <- file.path(base_path, "Figure4", "Deconvolution")
  
  # VSG patients with Pre/Post pairs
  IT07_pre <- pdf_to_plot(file.path(deconv_path, "IT_07 (RH-59-T)_baseline_pie.pdf"), density)
  IT07_post <- pdf_to_plot(file.path(deconv_path, "IT_07 (RH-59-T)_12_months_post_surgery_pie.pdf"), density)
  
  IT08_pre <- pdf_to_plot(file.path(deconv_path, "IT_08 (RH-60-T)_baseline_pie.pdf"), density)
  IT08_post <- pdf_to_plot(file.path(deconv_path, "IT_08 (RH-60-T)_12_months_post_surgery_pie.pdf"), density)
  
  IT11_pre <- pdf_to_plot(file.path(deconv_path, "IT_11_baseline_pie.pdf"), density)
  IT11_post <- pdf_to_plot(file.path(deconv_path, "IT_11_12_months_post_surgery_pie.pdf"), density)
  
  IT12_pre <- pdf_to_plot(file.path(deconv_path, "IT_12_baseline_pie.pdf"), density)
  IT12_post <- pdf_to_plot(file.path(deconv_path, "IT_12_12_months_post_surgery_pie.pdf"), density)
  
  IT14_pre <- pdf_to_plot(file.path(deconv_path, "IT_14_baseline_pie.pdf"), density)
  IT14_post <- pdf_to_plot(file.path(deconv_path, "IT_14_12_months_post_surgery_pie.pdf"), density)
  
  # 5 rows × 2 columns (Pre | Post)
  fig <- (IT07_pre | IT07_post) /
    (IT08_pre | IT08_post) /
    (IT11_pre | IT11_post) /
    (IT12_pre | IT12_post) /
    (IT14_pre | IT14_post) +
    plot_annotation(
      title = "Supplemental: Spatial Deconvolution - All VSG Patients",
      subtitle = "Left: Baseline | Right: 12 months post-surgery",
      tag_levels = list(c("IT_07", "", "IT_08", "", "IT_11", "", "IT_12", "", "IT_14", ""))
    )
  
  output_path <- file.path(base_path, "FigureS_Deconv_VSG_all.pdf")
  ggsave(output_path, fig, width = 14, height = 20, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

make_supp_deconv_controls <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    density = 250
) {
  
  deconv_path <- file.path(base_path, "Figure4", "Deconvolution")
  
  # Healthy controls (CRC samples)
  CRC02 <- pdf_to_plot(file.path(deconv_path, "CRC-02_baseline_pie.pdf"), density)
  CRC03 <- pdf_to_plot(file.path(deconv_path, "CRC-03_baseline_pie.pdf"), density)
  CRC10 <- pdf_to_plot(file.path(deconv_path, "CRC-10_baseline_pie.pdf"), density)
  CRC13 <- pdf_to_plot(file.path(deconv_path, "CRC-13_baseline_pie.pdf"), density)
  CRC14 <- pdf_to_plot(file.path(deconv_path, "CRC-14_baseline_pie.pdf"), density)
  CRC39 <- pdf_to_plot(file.path(deconv_path, "CRC-39_baseline_pie.pdf"), density)
  CRC40 <- pdf_to_plot(file.path(deconv_path, "CRC-40_baseline_pie.pdf"), density)
  
  # SMT patients
  RH23 <- pdf_to_plot(file.path(deconv_path, "RH-23-T_baseline_pie.pdf"), density)
  RH50 <- pdf_to_plot(file.path(deconv_path, "RH-50-T_baseline_pie.pdf"), density)
  RH63 <- pdf_to_plot(file.path(deconv_path, "RH-63-T_baseline_pie.pdf"), density)
  
  # 3×3 grid of controls + SMT
  fig <- (CRC02 | CRC03 | CRC10) /
    (CRC13 | CRC14 | CRC39) /
    (CRC40 | RH23 | RH50) +
    plot_annotation(
      title = "Supplemental: Spatial Deconvolution - Controls & SMT",
      subtitle = "Top two rows: Healthy controls (CRC) | Bottom row: SMT patients",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 16, face = "bold"))
  
  output_path <- file.path(base_path, "FigureS_Deconv_Controls.pdf")
  ggsave(output_path, fig, width = 16, height = 14, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

################################################################################
# SUPPLEMENTAL: Module scores for all VSG patients
################################################################################

make_supp_module_scores_vsg <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    density = 250
) {
  
  module_path <- file.path(base_path, "Figure4", "Module_score")
  
  # VSG patients Pre/Post
  IT07_pre <- pdf_to_plot(file.path(module_path, "IT_07 (RH-59-T)_baseline_module_scores.pdf"), density)
  IT07_post <- pdf_to_plot(file.path(module_path, "IT_07 (RH-59-T)_12_months_post_surgery_module_scores.pdf"), density)
  
  IT11_pre <- pdf_to_plot(file.path(module_path, "IT_11_baseline_module_scores.pdf"), density)
  IT11_post <- pdf_to_plot(file.path(module_path, "IT_11_12_months_post_surgery_module_scores.pdf"), density)
  
  IT12_pre <- pdf_to_plot(file.path(module_path, "IT_12_baseline_module_scores.pdf"), density)
  IT12_post <- pdf_to_plot(file.path(module_path, "IT_12_12_months_post_surgery_module_scores.pdf"), density)
  
  fig <- (IT07_pre | IT07_post) /
    (IT11_pre | IT11_post) /
    (IT12_pre | IT12_post) +
    plot_annotation(
      title = "Supplemental: Spatial Module Scores - VSG Patients",
      subtitle = "Left: Baseline | Right: 12 months post-surgery"
    )
  
  output_path <- file.path(base_path, "FigureS_ModuleScores_VSG.pdf")
  ggsave(output_path, fig, width = 16, height = 18, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

################################################################################
# SUPPLEMENTAL: TLS scores for all VSG patients
################################################################################

make_supp_tls_scores_vsg <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    density = 250
) {
  
  tls_path <- file.path(base_path, "Figure4", "TLS_score")
  
  # VSG patients Pre/Post
  IT07_pre <- pdf_to_plot(file.path(tls_path, "IT_07 (RH-59-T)_baseline_TLS_scores.pdf"), density)
  IT07_post <- pdf_to_plot(file.path(tls_path, "IT_07 (RH-59-T)_12_months_post_surgery_TLS_scores.pdf"), density)
  
  IT11_pre <- pdf_to_plot(file.path(tls_path, "IT_11_baseline_TLS_scores.pdf"), density)
  IT11_post <- pdf_to_plot(file.path(tls_path, "IT_11_12_months_post_surgery_TLS_scores.pdf"), density)
  
  IT12_pre <- pdf_to_plot(file.path(tls_path, "IT_12_baseline_TLS_scores.pdf"), density)
  IT12_post <- pdf_to_plot(file.path(tls_path, "IT_12_12_months_post_surgery_TLS_scores.pdf"), density)
  
  fig <- (IT07_pre | IT07_post) /
    (IT11_pre | IT11_post) /
    (IT12_pre | IT12_post) +
    plot_annotation(
      title = "Supplemental: Spatial TLS Scores - VSG Patients",
      subtitle = "Left: Baseline | Right: 12 months post-surgery"
    )
  
  output_path <- file.path(base_path, "FigureS_TLS_VSG.pdf")
  ggsave(output_path, fig, width = 14, height = 16, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

################################################################################
# RUN FIGURE 4
################################################################################

# Main Figure 4
fig4 <- make_figure4()

# Alternative versions:
# fig4_multi <- make_figure4_multi_patient()
# fig4_compact <- make_figure4_compact()

# Supplemental figures:
# supp_deconv_vsg <- make_supp_deconv_vsg()
# supp_deconv_ctrl <- make_supp_deconv_controls()
# supp_module_vsg <- make_supp_module_scores_vsg()
# supp_tls_vsg <- make_supp_tls_scores_vsg()






###Figure 5 
################################################################################
# FIGURE 5 ASSEMBLY - Clinical Outcomes & Risk Stratification
# 
# Your file location: 
# C:\Users\netio\Documents\UofW\Projects\Long_Paper\VSG_analysis-main\VSG_analysis-main\output\Figure5
################################################################################

library(magick)
library(pdftools)
library(cowplot)
library(ggplot2)
library(patchwork)

# Set your base path
base_path <- "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"

# Helper function
pdf_to_plot <- function(path, density = 300) {
  if (!file.exists(path)) {
    stop(paste("File not found:", path))
  }
  img <- image_read_pdf(path, density = density)
  ggdraw() + draw_image(img)
}

################################################################################
# YOUR FIGURE 5 FILES:
################################################################################
#
# GLYCEMIC CONTROL (Loss of glycemic control outcome):
# - GLYC_dendrogram_risk.pdf   (PDF 75 - hierarchical clustering)
# - GLYC_risk_group_KM.pdf     (PDF 76 - Kaplan-Meier, HR=1.26)
# - coef_GLYC_all_genes.csv    (coefficients data)
#
# MICROALBUMINURIA (MIC - Severe/Moderate Albuminuria outcome):
# - MIC_dendrogram_risk.pdf    (PDF 77 - hierarchical clustering)
# - MIC_risk_group_KM.pdf      (PDF 78 - Kaplan-Meier, HR=1.74)
# - coef_MIC_all_genes.csv     (coefficients data)
################################################################################

################################################################################
# MAIN FIGURE 5: Clinical Outcomes
################################################################################

make_figure5 <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    output = "Figure5_Clinical_Outcomes.pdf",
    width = 14,
    height = 10,
    density = 300
) {
  
  fig5_path <- file.path(base_path, "Figure5")
  
  message("Loading panels...")
  
  # Panel A: Glycemic control dendrogram
  A <- pdf_to_plot(file.path(fig5_path, "GLYC_dendrogram_risk.pdf"), density)
  
  # Panel B: Glycemic control KM curve
  B <- pdf_to_plot(file.path(fig5_path, "GLYC_risk_group_KM.pdf"), density)
  
  # Panel C: Albuminuria dendrogram
  C <- pdf_to_plot(file.path(fig5_path, "MIC_dendrogram_risk.pdf"), density)
  
  # Panel D: Albuminuria KM curve
  D <- pdf_to_plot(file.path(fig5_path, "MIC_risk_group_KM.pdf"), density)
  
  message("Assembling figure...")
  
  # Layout:
  # Row 1: A (GLYC dendro) | B (GLYC KM)
  # Row 2: C (MIC dendro)  | D (MIC KM)
  
  fig5 <- (A | B) / (C | D) +
    plot_layout(widths = c(1.2, 1)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 24, face = "bold"))
  
  output_path <- file.path(base_path, output)
  ggsave(output_path, fig5, width = width, height = height, dpi = 300)
  
  message("✓ Figure 5 saved: ", output_path)
  return(fig5)
}

################################################################################
# FIGURE 5 VERSION 2: KM curves only (if dendrograms go to supplemental)
################################################################################

make_figure5_km_only <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    output = "Figure5_KM_only.pdf",
    width = 12,
    height = 5,
    density = 300
) {
  
  fig5_path <- file.path(base_path, "Figure5")
  
  A <- pdf_to_plot(file.path(fig5_path, "GLYC_risk_group_KM.pdf"), density)
  B <- pdf_to_plot(file.path(fig5_path, "MIC_risk_group_KM.pdf"), density)
  
  fig5 <- (A | B) +
    plot_annotation(
      title = "Risk Stratification Predicts Clinical Outcomes",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 24, face = "bold"))
  
  output_path <- file.path(base_path, output)
  ggsave(output_path, fig5, width = width, height = height, dpi = 300)
  
  message("✓ Figure 5 (KM only) saved: ", output_path)
  return(fig5)
}

################################################################################
# FIGURE 5 VERSION 3: Vertical layout
################################################################################

make_figure5_vertical <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output",
    output = "Figure5_vertical.pdf",
    width = 10,
    height = 14,
    density = 300
) {
  
  fig5_path <- file.path(base_path, "Figure5")
  
  A <- pdf_to_plot(file.path(fig5_path, "GLYC_dendrogram_risk.pdf"), density)
  B <- pdf_to_plot(file.path(fig5_path, "GLYC_risk_group_KM.pdf"), density)
  C <- pdf_to_plot(file.path(fig5_path, "MIC_dendrogram_risk.pdf"), density)
  D <- pdf_to_plot(file.path(fig5_path, "MIC_risk_group_KM.pdf"), density)
  
  # Vertical stacking
  fig5 <- A / B / C / D +
    plot_layout(heights = c(0.8, 1, 0.8, 1)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(size = 24, face = "bold"))
  
  output_path <- file.path(base_path, output)
  ggsave(output_path, fig5, width = width, height = height, dpi = 300)
  
  message("✓ Figure 5 (vertical) saved: ", output_path)
  return(fig5)
}

################################################################################
# SUPPLEMENTAL: Dendrograms only
################################################################################

make_supp_dendrograms_outcomes <- function(
    base_path = "C:/Users/netio/Documents/UofW/Projects/Long_Paper/VSG_analysis-main/VSG_analysis-main/output"
) {
  
  fig5_path <- file.path(base_path, "Figure5")
  
  A <- pdf_to_plot(file.path(fig5_path, "GLYC_dendrogram_risk.pdf"))
  B <- pdf_to_plot(file.path(fig5_path, "MIC_dendrogram_risk.pdf"))
  
  fig <- A / B +
    plot_annotation(
      title = "Supplemental: Risk Group Clustering",
      subtitle = "A) Loss of Glycemic Control | B) Moderate/Severe Albuminuria",
      tag_levels = 'A'
    ) &
    theme(plot.tag = element_text(size = 20, face = "bold"))
  
  output_path <- file.path(base_path, "FigureS_Risk_Dendrograms.pdf")
  ggsave(output_path, fig, width = 12, height = 10, dpi = 300)
  message("✓ Saved: ", output_path)
  return(fig)
}

################################################################################
# RUN FIGURE 5
################################################################################

# Main Figure 5
fig5 <- make_figure5()

# Alternative versions:
# fig5_km <- make_figure5_km_only()
# fig5_vert <- make_figure5_vertical()

# Supplemental:
# supp_dendro <- make_supp_dendrograms_outcomes()




















