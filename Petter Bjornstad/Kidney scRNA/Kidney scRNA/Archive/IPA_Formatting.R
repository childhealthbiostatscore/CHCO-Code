library(ggrepel)
# library(qpcR)
library(ggpubr)

# test_ipa <- readxl::read_xls("/Users/choiyej/Dropbox/PANTHER/IPA results/tanner_stage_12_345_rh_de_table_export.xls", skip = 1)
test_ipa <- readxl::read_xls("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/AST_pathways.xls", skip = 1)

colnames(test_ipa) <- c("pathway", "neglog_p", "ratio", "zscore", "molecules")



# Volcano plot
sig_cutoff = -log(0.05, base = 10)

test_ipa <- test_ipa %>%
  dplyr::mutate(direction = case_when(zscore < 0 & neglog_p > sig_cutoff ~ "Negative", 
                                      zscore > 0 & neglog_p > sig_cutoff ~ "Positive",
                                      zscore == 0 & neglog_p > sig_cutoff ~ "Unknown",
                                      T ~ "NS")) %>%
  group_by(direction) %>%
  dplyr::mutate(count = row_number(),
                label = case_when(direction != "NS" ~ paste0(str_sub(direction, 1,1), count))) %>%
  ungroup()

test_ipa$direction <- factor(test_ipa$direction, levels = c("Negative", "Positive", "Unknown", "NS"))

ipa.plot <- test_ipa %>% filter(neglog_p > 0) %>%
  ggplot(aes(x = zscore, y = neglog_p)) + 
  geom_point(aes(color = direction), alpha = 0.5, size = 2) +
  geom_hline(aes(yintercept = sig_cutoff), color = "red", linetype = "dashed") +
  labs(x = "Z-Score",
       y = "-log(p-value)",
       color = "Direction") + 
  scale_color_manual(values = c("#bf0603", "#457b9d",
                                "#f6bd60", "#dad7cd")) +
  ggrepel::geom_label_repel(aes(label = label, color = direction),
                            label.size = 0.15,
                            max.overlaps = 100, 
                            min.segment.length = 0.01) +
  theme_bw()

# Table
ipa.table.neg <- test_ipa %>%
  dplyr::filter(direction == "Negative") %>%
  dplyr::select(label, pathway) 
ipa.table.pos <- test_ipa %>%
  dplyr::filter(direction == "Positive") %>%
  dplyr::select(label, pathway) 
ipa.table.unk <- test_ipa %>%
  dplyr::filter(direction == "Unknown") %>%
  dplyr::select(label, pathway) 

cbind.fill <- function(...) {
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function(x) {
    filled_matrix <- ifelse(is.na(x), "", x)
    rbind(filled_matrix, matrix("", n - nrow(filled_matrix), ncol(filled_matrix)))
  }))
}

ipa.table.comb <- cbind.fill(ipa.table.neg,ipa.table.pos,ipa.table.unk)
gg.ipa.table <- ggtexttable(ipa.table.comb, rows = NULL,
                            cols = c("Label", "Pathway", "Label", "Pathway","Label", "Pathway"),
                            theme = ttheme("blank",
                                           tbody.style = tbody_style(fill = "white", 
                                                                     size = 8.5, 
                                                                     hjust = 0,
                                                                     x = 0.1))) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 0.8)


pdf(fs::path(dir.results,"PT_T2D_sglt2i_glp1_vs_glp1.pdf"),width=15,height=10)
ggarrange(ipa.plot, gg.ipa.table,
          ncol = 2, heights = c(1,1), widths = c(0.5,1),
          legend = "top")
dev.off()

library(readxl)
library(ggplot2)
library(dplyr)
library(forcats)
library(ggpubr)

color_table <- tibble(
  up_down = c(-1, 0, 1),
  color = c("indianred", "gray", "steelblue")
)

ipa_plot <- function(data){
  data <- data %>% mutate(up_down = case_when(
    abs(`z-score`) < 0.0001 | is.na(`z-score`) ~ 0,
    `z-score` <= -0.0001 ~ -1,
    `z-score` >= 0.0001 ~ 1,
    .default = 0
  ))
  data$up_down <- factor(data$up_down, levels = unique(data$up_down))
  p <- ggplot(data, aes(reorder(`Ingenuity Canonical Pathways`, `-log(p-value)`), `-log(p-value)`, fill = up_down))+
    geom_col() +
    geom_hline(yintercept = -log(0.05), linetype = "dashed", color = "black") +
    theme_classic()+
    # scale_y_continuous(limits = c(0, 6)) +
    ylab("-log(p-value)") +
    xlab("Pathway") + 
    # theme(legend.position="none") +
    scale_fill_manual(values = c("-1" = "steelblue", "0" = "grey", "1" = "indianred")) 
  p <- p + coord_flip()
}

# alb <- read_xls("/Volumes/PEDS/PEDS/RI Biostatistics Core/Shared/Shared Projects/Laura/Peds Endo/Petter Bjornstad/TODAY subaward/Results/Linear and Cox models/IPA/Output from IPA/albuminuria.xls",
#                 skip = 1)
ast <- readxl::read_xls("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/AST_pathways.xls", skip = 1)

ast <- ast %>% arrange(desc("-log(p-value)")) 
ast_keep <- ast[1:60,]

ast_plot <- ipa_plot(ast_keep)
ast_plot 
pdf("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/AST_pathways.pdf",width=10,height=8)
plot(ast_plot)
dev.off()

ast <- readxl::read_xls("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/AST_pathways_02_pval.xls", skip = 1)

ast <- ast %>% arrange(desc("-log(p-value)")) 
ast_keep <- ast[1:60,]

ast_plot <- ipa_plot(ast_keep)
ast_plot 
pdf("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/AST_pathways_02_pval.pdf",width=10,height=8)
plot(ast_plot)
dev.off()

alt <- readxl::read_xls("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/ALT_pathways.xls", skip = 1)

alt <- alt %>% arrange(desc("-log(p-value)")) 
alt_keep <- alt[1:60,]

alt_plot <- ipa_plot(alt_keep)
alt_plot 
pdf("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/ALT_pathways.pdf",width=10,height=8)
plot(alt_plot)
dev.off()

alt <- readxl::read_xls("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/alt_pathways_02_pval.xls", skip = 1)

alt <- alt %>% arrange(desc("-log(p-value)")) 
alt_keep <- alt[1:60,]

alt_plot <- ipa_plot(alt_keep)
alt_plot 
pdf("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/alt_pathways_02_pval.pdf",width=10,height=8)
plot(alt_plot)
dev.off()
