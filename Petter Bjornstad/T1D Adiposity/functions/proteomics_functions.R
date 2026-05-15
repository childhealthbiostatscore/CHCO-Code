


########################
# specify user for paths
########################

user <- Sys.info()[["user"]]
  
if (user == "choiyej") { # local version
    root_path <- "/Users/choiyej/Library/CloudStorage/OneDrive-UW/Bjornstad/Biostatistics Core Shared Drive"
    git_path <- "/Users/choiyej/GitHub/CHCO-Code/Petter Bjornstad"
    keys <- fromJSON("/Users/choiyej/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Bjornstad Pyle Lab/keys.json")
} else if (user == "rameshsh") { # hyak version
    root_path <- ""
    git_path <- "/mmfs1/gscratch/togo/rameshsh/CHCO-Code/Petter Bjornstad"
} else if (user == "yejichoi") { # hyak version
    root_path <- "/mmfs1/gscratch/togo/yejichoi/"
    git_path <- "/mmfs1/gscratch/togo/yejichoi/CHCO-Code/Petter Bjornstad"
    keys <- fromJSON("/mmfs1/home/yejichoi/keys.json")
} else if (user == "leidholt") {
    root_path <- c("/mmfs1/gscratch/togo/leidholt/")
    git_path <- "/mmfs1/gscratch/togo/leidholt/CHCO-Code/"
    keys <- fromJSON("/mmfs1/home/leidholt/keys.json")
} else if (user == "sleidholt") {
    root_path <- c("/Users/Shared/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/")
    dir.results <- c("/Users/Shared/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/T1D Adiposity/Results/Compiled Results/Proteomics Results")
    git_path <- "/Users/sleidholt/Documents/GitHub/CHCO-Code/Petter Bjornstad/"
    keys <- fromJSON("/Users/Shared/OneDrive - UW/Personal_storage/keys.json")
} else if (user == "savanahleidholt") {
    root_path <- c("/Users/savanahleidholt/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/")
    dir.results <- c("/Users/savanahleidholt/OneDrive - UW/Laura Pyle's files - Biostatistics Core Shared Drive/T1D Adiposity/Results/Compiled Results/Proteomics Results")
    git_path <- "/Users/savanahleidholt/Desktop/CHCO-Code/Petter Bjornstad"
    keys <- fromJSON("/Users/savanahleidholt/OneDrive - UW/Personal_storage/keys.json")
} else {
    stop("Unknown user: please specify root path for this user.")
}


#---------------------------------------------------------
#Function for plotting clinical data of interest
#---------------------------------------------------------
plot_clinical_boxplot <- function(data, outcome_var, outcome_label) {
  plot_df <- data %>%
    dplyr::select(group_bmi, age, sex, study, dplyr::all_of(outcome_var)) %>%
    dplyr::mutate(
      y = as.numeric(.data[[outcome_var]]),
      group_bmi = factor(group_bmi),
      sex = factor(sex),
      study = factor(study)
    ) %>%
    tidyr::drop_na(y, group_bmi)
  
  ggplot2::ggplot(plot_df,
    ggplot2::aes(
      x = group_bmi,
      y = y,
      fill = group_bmi
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.65,
      alpha = 0.75,
      outlier.shape = NA
    ) +
    ggplot2::geom_jitter(
      width = 0.15,
      size = 2.5,
      alpha = 0.75
    ) +
    ggplot2::labs(
      title = outcome_label,
      x = NULL,
      y = outcome_label
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "LC_Normal" = "cyan4",
        "LC_Overweight_Obese" = "darkred",
        "T1D_Normal" = "gold",
        "T1D_Overweight" = "plum3",
        "T1D_Obese" = "black"
      )
    ) +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      panel.grid = ggplot2::element_blank()
    )
}

#-------------------------------------
#PCA plotting function for group_bmi
#-------------------------------------
make_pca <- function(mat, meta_df, title_str, add_ellipses = TRUE) {
  stopifnot(identical(rownames(mat), meta_df$record_id))
  
  pca <- prcomp(mat, center = TRUE, scale. = TRUE)
  vexp <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)
  
  pca_df <- tibble(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    group_bmi = meta_df$group_bmi,
    study = meta_df$study
  )
  
  p <- ggplot(
    pca_df,
    aes(
      x = PC1,
      y = PC2,
      color = group_bmi
    )
  ) +
    {
      if (add_ellipses) {
        stat_ellipse(
          aes(group = group_bmi, color = group_bmi),
          type = "norm",
          level = 0.68,
          linewidth = 1,
          linetype = "solid",
          alpha = 0.8
        )
      }
    } +
    geom_point(aes(shape = study), size = 3, alpha = 0.85) +
    scale_color_manual(
      values = c(
        LC_Normal  = "cyan4",
        LC_Overweight_Obese = "darkred",
        T1D_Normal = "gold",
        T1D_Overweight = "plum3",
        T1D_Obese = "black"
      ),
      name = "Disease-adiposity group"
    ) +
    labs(
      title = title_str,
      x = paste0("PC1 (", vexp[1], "%)"),
      y = paste0("PC2 (", vexp[2], "%)"),
      color = "Disease-adiposity group",
      shape = "Study"
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  return(p)
}

