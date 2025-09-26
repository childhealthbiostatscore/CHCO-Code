#Libraries
library(ggrepel)
# library(qpcR)
library(ggpubr)
library(readxl)
library(ggplot2)
library(dplyr)
library(forcats)
library(ggpubr)

#Directories
user <- "hhampson" 
dir.dat <- c(paste0("/Users/",user,"/Library/CloudStorage/OneDrive-SharedLibraries-UW/Laura Pyle - Biostatistics Core Shared Drive/Liver project/IPA Results"))


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
ast <- readxl::read_xls(fs::path(dir.dat,"AST_IPA_Results.xls"), skip = 1)

ast <- ast %>% arrange(desc("-log(p-value)")) 
ast_keep <- ast[1:60,]

ast_plot <- ipa_plot(ast_keep)
ast_plot 
pdf("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/AST_pathways_TeenLabs.pdf",width=10,height=8)
plot(ast_plot)
dev.off()

ast <- readxl::read_xls("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/ALT_TeenLabs.xls", skip = 1)

ast <- ast %>% arrange(desc("-log(p-value)")) 
ast_keep <- ast[1:60,]

ast_plot <- ipa_plot(ast_keep)
ast_plot 
pdf("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/ALT_pathways_TeenLabs.pdf",width=10,height=8)
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

#Ratio
ratio <- readxl::read_xls("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/ALT_AST_Ratio_pathways.xls", skip = 1)

ratio <- ratio %>% arrange(desc("-log(p-value)")) 
ratio_keep <- ratio[1:60,]

ratio_plot <- ipa_plot(ratio_keep)
ratio_plot 
pdf("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/Ratio_pathways.pdf",width=10,height=8)
plot(ratio_plot)
dev.off()

ratio <- readxl::read_xls("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/ALT_AST_Ratio_pathways_02_pval.xls", skip = 1)

ratio <- ratio %>% arrange(desc("-log(p-value)")) 
ratio_keep <- ratio[1:60,]

ratio_plot <- ipa_plot(ratio_keep)
ratio_plot 
pdf("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Liver project/IPA Results/ratio_pathways_02_pval.pdf",width=10,height=8)
plot(ratio_plot)
dev.off()
