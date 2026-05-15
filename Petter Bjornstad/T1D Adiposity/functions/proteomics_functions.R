---
title: "proteomics_functions"
author: "Ye Ji Choi & Savanah Leidholt"
format: html
---
  
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