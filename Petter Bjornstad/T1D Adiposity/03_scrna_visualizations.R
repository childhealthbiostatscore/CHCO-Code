library(tidyverse)
library(arsenal)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(Hmisc)
library(readxl)
library(purrr)
library(aws.s3)
library(jsonlite)
library(RColorBrewer)
library(ggtrace)

user <- Sys.info()[["user"]]

if (user == "choiyej") {
  root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive/"
  git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad/"
  keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
} else if (user == "yejichoi") {
  root_path <- "/mmfs1/gscratch/togo/yejichoi/"
  git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
  keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "pylell") {
  root_path <- "/Users/pylell/Library/CloudStorage/OneDrive-SharedLibraries-UW/Bjornstad/Biostatistics Core Shared Drive"
  git_path <- "/Users/pylell/Documents/GitHub/CHCO-Code/Petter Bjornstad"
} else if (user == "hhampson") {
  root_path <- "/Volumes/Peds Endo"
} else {
  stop("Unknown user: please specify root path for this user.")
}

source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_scRNA_functions.R"))
source(file.path(git_path, "Renal HEIRitage/RH_RH2_IMPROVE_functions.R"))

## Create an S3 client

Sys.setenv(
  "AWS_ACCESS_KEY_ID" = keys$MY_ACCESS_KEY,
  "AWS_SECRET_ACCESS_KEY" = keys$MY_SECRET_KEY,
  "AWS_DEFAULT_REGION" = "",
  "AWS_REGION" = "",
  "AWS_S3_ENDPOINT" = "s3.kopah.uw.edu"
)

pb90_attempt_subset <- s3readRDS(
          object = "data_clean/t1d_hc_scrna_w_clinical.rds", 
          bucket = "t1d.adiposity",
          region = "")

my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(43)
# Read Seurat object and run UMAP, obtain cell numbers, etc
pb90_attempt_subset_meta <- pb90_attempt_subset@meta.data %>%
  cbind(pb90_attempt_subset@reductions$umap.harmony@cell.embeddings)

s3saveRDS(pb90_attempt_subset_meta, object = "data_clean/pb90_attempt_subset_meta.rds", bucket = "t1d.adiposity", region = "", multipart = T)

# start from here if file already saved
pb90_attempt_subset_meta <- s3readRDS(object = "data_clean/pb90_attempt_subset_meta.rds", bucket = "t1d.adiposity", region = "")

centers <- pb90_attempt_subset_meta %>%
  group_by(KPMP_celltype) %>%
  summarise(x = median(umapharmony_1),
            y = median(umapharmony_2), .groups = "drop")

umap_p <- ggplot(pb90_attempt_subset_meta,
                 aes(x = umapharmony_1, y = umapharmony_2, color = KPMP_celltype)) +
  geom_point(alpha = 0.2, shape = 16) +
  geom_text_repel(data = centers,
                  aes(x = x, y = y, label = KPMP_celltype),
                  inherit.aes = FALSE,
                  size = 4, fontface = "bold", color = "black") +
  scale_color_manual(values = my_colors) +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(hjust = 0),   # move x label to left
        axis.title.y = element_text(hjust = 0))  + # push y label down to bottom) 
  labs(x = "UMAP1", y = "UMAP2")

# Save to a temporary connection
tmp <- tempfile(fileext = ".png")
ggsave(tmp, plot = umap_p, width = 10, height = 8, dpi = 300)

# Upload directly
put_object(
  file = tmp,
  object = "results/figures/t1d_hc_umap_kpmpcelltype.png",
  bucket = "t1d.adiposity",
  region = ""
)


# UMAP by source
umap_p_source <- ggplot(pb90_attempt_subset_meta,
                 aes(x = umapharmony_1, y = umapharmony_2, color = source)) +
  geom_point(alpha = 0.05, shape = 16) +
  scale_color_manual(values = c("attempt" = "#c1121f", 
                                "pb90" = "#62b6cb")) +
  theme_minimal() + 
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(hjust = 0),   # move x label to left
        axis.title.y = element_text(hjust = 0))  + # push y label down to bottom) 
  labs(x = "UMAP1", y = "UMAP2",
       color = "Source") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

# Save to a temporary connection
tmp <- tempfile(fileext = ".png")
ggsave(tmp, plot = umap_p_source, width = 10, height = 8, dpi = 300)

# Upload directly
put_object(
  file = tmp,
  object = "results/figures/t1d_hc_umap_source.png",
  bucket = "t1d.adiposity",
  region = ""
)


cell_counts <- pb90_attempt_subset_meta %>%
  dplyr::count(KPMP_celltype_general) %>%
  arrange(desc(n))

# Create bar plot
cell_count_bar <- ggplot(cell_counts, aes(x = reorder(KPMP_celltype_general, n), y = n)) +
  geom_bar(stat = "identity", fill = "#97a97c") +
  geom_text(aes(label = n), hjust = -0.1, size = 4) +
  coord_flip() +
  labs(x = NULL, 
       y = "Count") +
  theme(text = element_text(size = 15),
        legend.position = "none",
        panel.grid = element_blank()) +
  theme_transparent +
  scale_y_continuous(expand = expansion(add = 5000)) 

# Save to a temporary connection
tmp <- tempfile(fileext = ".png")
ggsave(tmp, plot = cell_count_bar, width = 10, height = 8, dpi = 300)

# Upload directly
put_object(
  file = tmp,
  object = "results/figures/t1d_hc_kpmp_cell_count_bar.png",
  bucket = "t1d.adiposity",
  region = ""
)

# cell counts per group
cell_counts_group <- pb90_attempt_subset_meta %>%
  group_by(group) %>%
  dplyr::count(KPMP_celltype_general) %>%
  arrange(desc(n)) %>%
  dplyr::mutate(group = case_when(group == "Lean Control" ~ "HC",
                                  group == "Type 1 Diabetes" ~ "T1D"))

ggplot(cell_counts_group, aes(x = reorder(KPMP_celltype_general, n), y = n, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("T1D" = "#f27059",
                               "HC" = "#81b29a")) +
  coord_flip() +
  labs(x = NULL, 
       y = "Count",
       fill = "Group") +
  theme(text = element_text(size = 15),
        panel.grid = element_blank()) +
  theme_transparent +
  scale_y_continuous(expand = expansion(add = 5000)) 

# Save to a temporary connection
tmp <- tempfile(fileext = ".png")
ggsave(tmp, width = 10, height = 8, dpi = 300)

# Upload directly
put_object(
  file = tmp,
  object = "results/figures/t1d_hc_grp_kpmp_cell_count_bar.png",
  bucket = "t1d.adiposity",
  region = ""
)

# cell counts per obesity group in T1D
cell_counts_dxa_ob <- pb90_attempt_subset_meta %>%
  filter(group == "Type 1 Diabetes" & !is.na(dxa_obesity)) %>%
  group_by(dxa_obesity) %>%
  dplyr::count(KPMP_celltype_general) %>%
  arrange(desc(n))

ggplot(cell_counts_dxa_ob, aes(x = reorder(KPMP_celltype_general, n), y = n, fill = dxa_obesity)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Normal" = "#86ba90",
                               "Overweight" = "#f4845f",
                               "Obese" = "#f25c54")) +
  coord_flip() +
  labs(x = NULL, 
       y = "Count",
       fill = "Group\nby DXA") +
  theme(text = element_text(size = 15),
        panel.grid = element_blank()) +
  theme_transparent +
  scale_y_continuous(expand = expansion(add = 5000)) 

# Save to a temporary connection
tmp <- tempfile(fileext = ".png")
ggsave(tmp, width = 10, height = 8, dpi = 300)

# Upload directly
put_object(
  file = tmp,
  object = "results/figures/t1d_dxa_ob_kpmp_cell_count_bar.png",
  bucket = "t1d.adiposity",
  region = ""
)

# cell counts per obesity group in T1D
cell_counts_bmi_ob <- pb90_attempt_subset_meta %>%
  filter(group == "Type 1 Diabetes" & !is.na(bmi_obesity)) %>%
  group_by(bmi_obesity) %>%
  dplyr::count(KPMP_celltype_general) %>%
  arrange(desc(n))

ggplot(cell_counts_bmi_ob, aes(x = reorder(KPMP_celltype_general, n), y = n, fill = bmi_obesity)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Normal" = "#86ba90",
                               "Overweight" = "#f4845f",
                               "Obese" = "#f25c54")) +
  coord_flip() +
  labs(x = NULL, 
       y = "Count",
       fill = "Group\nby BMI %ile") +
  theme(text = element_text(size = 15),
        panel.grid = element_blank()) +
  theme_transparent +
  scale_y_continuous(expand = expansion(add = 5000)) 

# Save to a temporary connection
tmp <- tempfile(fileext = ".png")
ggsave(tmp, width = 10, height = 8, dpi = 300)

# Upload directly
put_object(
  file = tmp,
  object = "results/figures/t1d_bmi_ob_kpmp_cell_count_bar.png",
  bucket = "t1d.adiposity",
  region = ""
)

