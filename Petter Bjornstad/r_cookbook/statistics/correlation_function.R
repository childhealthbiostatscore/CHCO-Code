library(dplyr)
library(Hmisc)
library(kableExtra)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(tidyselect)
library(utils)
library(purrr)
library(RColorBrewer)
library(labelled)
library(corrplot)
library(gtools)

correlation_table_minimal <- function(data, relevant_vars, n_cols, cor_method = "pearson", dict = NULL, save_path = NULL) {
  # Filter relevant variables for the correlation matrix
  dat_correlation <- data %>%
    dplyr::select(all_of(relevant_vars))
  
  # Compute the correlation matrix
  if (n_cols > 1) {
    M_table <- cor(y = dat_correlation[, 1:n_cols],
                   x = dat_correlation[, (n_cols + 1):ncol(dat_correlation)], 
                   use = "pairwise.complete.obs",
                   method = cor_method) %>%
      as.data.frame() %>%
      dplyr::mutate_all(round, digits = 3)
  }
  if (n_cols == 1) {
    M_table <- cor(y = dat_correlation[1],
                   x = dat_correlation[, (n_cols + 1):ncol(dat_correlation)], 
                   use = "pairwise.complete.obs",
                   method = cor_method) %>%
      as.data.frame() %>%
      dplyr::mutate_all(round, digits = 3)
  }
  
  # Compute p-values 
  res2 <- Hmisc::rcorr(as.matrix(dat_correlation),
                       type = cor_method)
  corr_pval <- as.data.frame(res2$P) %>%
    dplyr::select(all_of(relevant_vars[1:n_cols])) %>%
    base::sapply(scales::pvalue) 
  corr_pval <- as.data.frame(corr_pval) %>% 
    dplyr::select(all_of(relevant_vars[1:n_cols]))
  colnames(corr_pval) <- rep("p_value", ncol(corr_pval))
  
  # Add p-values to matrix table
  corr_pval = corr_pval[(n_cols + 1):nrow(corr_pval),]
  M_table_p <- dplyr::bind_cols(M_table, corr_pval) %>%
    dplyr::relocate(purrr::map(1:ncol(M_table), ~ c(.x, .x + ncol(M_table))) %>% unlist())
  
  # Label table if a dictionary of labels is provided
  if (!is.null(dict)) {
    dict_corr <- dict
    dict_corr <- dict_corr[dplyr::intersect(names(M_table_p), names(dict_corr))]
    dict_corr[setdiff(names(M_table_p), names(dict_corr))] <- "p-value"
    Hmisc::label(M_table_p) = dict_corr[match(names(M_table_p), names(dict_corr))]
    row.names(M_table_p) <- dict[match(row.names(M_table_p), names(dict))]
  }

  # Set column labels to the labels of the column variables
  colnames(M_table_p) <- base::sapply(M_table_p, function(x) Hmisc::label(x))
  
  # Save file if path is provided
  if (!is.null(save_path)) {
    write.csv2(M_table_p, save_path)
  }
  
  # Format the correlation matrix table using kableExtra functions
  formatted_table <- knitr::kable(M_table_p,
                                  format = "html",
                                  table.attr = "style='width:100%;'",
                                  caption = "Correlation",
                                  format.args = list(nsmall = 3)) %>%
    kableExtra::kable_styling(bootstrap_options = c("hover"),
                              full_width = TRUE,
                              fixed_thead = TRUE) 
  
  return(formatted_table)
}

correlation_table_colored <- function(data, relevant_vars, n_cols, cor_method = "pearson", dict = NULL, save_path = NULL) {
  # Filter relevant variables for the correlation matrix
  dat_correlation <- data %>%
    dplyr::select(all_of(relevant_vars))
  
  # Compute the correlation matrix
  M_table <- cor(y = dat_correlation[, 1:n_cols],
                 x = dat_correlation[, (n_cols + 1):ncol(dat_correlation)], 
                 use = "pairwise.complete.obs",
                 method = cor_method) %>%
    as.data.frame() %>%
    dplyr::mutate_all(round, digits = 3)
  
  # Compute p-values 
  res2 <- Hmisc::rcorr(as.matrix(dat_correlation),
                       type = cor_method)
  corr_pval <- as.data.frame(res2$P) %>%
    dplyr::select(all_of(relevant_vars[1:n_cols])) %>%
    base::sapply(scales::pvalue) 
  corr_pval <- as.data.frame(corr_pval) %>% 
    dplyr::select(all_of(relevant_vars[1:n_cols]))
  colnames(corr_pval) <- rep("p_value", ncol(corr_pval))
  
  # Add p-values to matrix table
  corr_pval = corr_pval[(n_cols + 1):nrow(corr_pval),]
  M_table_p <- dplyr::bind_cols(M_table, corr_pval) %>%
    dplyr::relocate(purrr::map(1:ncol(M_table), ~ c(.x, .x + ncol(M_table))) %>% unlist())
  
  # Label table if a dictionary of labels is provided
  if (!is.null(dict)) {
    dict_corr <- dict
    dict_corr <- dict_corr[dplyr::intersect(names(M_table_p), names(dict_corr))]
    dict_corr[setdiff(names(M_table_p), names(dict_corr))] <- "p-value"
    Hmisc::label(M_table_p) = dict_corr[match(names(M_table_p), names(dict_corr))]
    row.names(M_table_p) <- dict[match(row.names(M_table_p), names(dict))]
  }
  
  # Set column labels to the labels of the column variables
  colnames(M_table_p) <- base::sapply(M_table_p, function(x) Hmisc::label(x))
  unlab <- remove_labels(M_table_p)
  
  # Save file if path is provided
  if (!is.null(save_path)) {
    write.csv2(M_table_p, save_path)
  }
  
  # Format the correlation matrix table using kableExtra functions
  colored_table <- M_table_p %>% 
    map2(seq_along(M_table_p), function(x, col_index) {
      if (col_index %% 2 == 1 && is.numeric(x)) {
        
        right_column <- M_table_p[[col_index + 1]]
        
        cell_spec(sprintf("%.3f", as.numeric(x)), 
                  bold = case_when((!is.na(as.numeric(right_column)) & as.numeric(right_column)<=0.05) |
                                     is.na(as.numeric(right_column)) ~ T, 
                                   T~ F),    
                  underline = case_when((!is.na(as.numeric(right_column)) & as.numeric(right_column)<=0.05) |
                                          is.na(as.numeric(right_column)) ~ T, 
                                        T~ F))
      } else {
        sprintf("%.3f", as.numeric(x))  # Keep the column as it is
      }
    })
  
  colored_table <- as.data.frame(colored_table) %>%
    dplyr::select(-starts_with("p.value"))
  M_table_p <- M_table_p %>% dplyr::select(-starts_with("p-value"))
  rownames(colored_table) <- rownames(M_table_p)
  
  formatted_table <- kbl(colored_table,
                         format = "html",
                         table.attr = "style='width:100%;'",
                         escape = F,
                         row.names = T,
                         col.names = colnames(M_table_p),
                         align = "c", 
                         digits = 3,
                         format.args = list(nsmall = 3))
  
  return(formatted_table)
}


scatter_plot_correlation <- function(data, x, y, group, cor_method = "pearson", group_colors = brewer.pal(9, "Set1")) {
  ggplot(data, aes(x, y)) +
  labs(x = Hmisc::label(x),
       y = Hmisc::label(y),
       color  = Hmisc::label(group)) +
  geom_point(aes(color = group)) +
  scale_colour_manual(values = group_colors) +
  sm_statCorr(color = '#5A5A5A', corr_method = cor_method,
              linetype = 'dashed') +
  sm_corr_theme()
}

correlation_p_value_matrix <- function(data, relevant_vars, n_cols, cor_method = "pearson") {
  # Filter relevant variables for the correlation matrix
  dat_correlation <- data %>%
    dplyr::select(all_of(relevant_vars))
  
  # Compute p-values 
  res2 <- Hmisc::rcorr(as.matrix(dat_correlation),
                       type = cor_method)
  corr_pval <- as.data.frame(res2$P) %>%
    dplyr::select(all_of(relevant_vars[1:n_cols]))
  corr_pval = corr_pval[(n_cols + 1):nrow(corr_pval),]
  corr_pval[is.na(corr_pval)] <- 1
  return(as.matrix(corr_pval))
}

corr_plot_modified <- function(
    data, X, Y, cor_method = "pearson", adj_var = NULL, dict = dict,
    method = "color", insig = "pch", coef_col = NULL,
    pch = 4, pch.col = "black",
    pch.cex = 0, tl.cex = .8, cl.cex = 0.8, number.cex = 1,
    coef_text_cex = 0.9,   # NEW: size for printed coefficients via graphics::text
    star_text_cex = 0.9,   # NEW: size for p-value stars via graphics::text
    col = colorRampPalette(c("#0096c7", "white", "#d62828"))(200)
) {
  n_cols <- length(Y)
  M <- cor(y = subset(data, select = Y),
           x = subset(data, select = X),
           use = "pairwise.complete.obs",
           method = cor_method)
  
  colnames(M) <- as.list(dict[match(colnames(M), names(dict))])
  rownames(M) <- as.list(dict[match(rownames(M), names(dict))])
  
  if (!is.null(adj_var)) {
    x_vars <- rep(X, times = length(Y))
    x_coord <- rep(seq(1, length(Y)), each = length(X))
    y_vars <- rep(Y, each = length(X))
    y_coord <- rep(seq(length(X), 1), times = length(Y))
    lm_extracted <- data.frame(yName = character(0),
                               xName = character(0),
                               x = numeric(0),
                               y = numeric(0),
                               adj_var = character(0),
                               adj_x_coef = numeric(0),
                               adj_x_pval = numeric(0))
    for (i in 1:length(x_vars)) {
      if (x_vars[i] != y_vars[i]) {
        lm_formula <- as.formula(paste0(y_vars[i], "~", x_vars[i], "+", adj_var))
        lm_model <- lm(lm_formula, data = data)
        lm_summary <- summary(lm_model)
        lm_coef <- lm_summary$coefficient[x_vars[i], "Estimate"]
        lm_p_val <- lm_summary$coefficient[x_vars[i], "Pr(>|t|)"]
        lm_extracted <- rbind(lm_extracted,
                              data.frame(
                                yName = y_vars[i],
                                xName = x_vars[i],
                                x = x_coord[i],
                                y = y_coord[i],
                                adj_var = adj_var,
                                adj_x_coef = lm_coef,
                                adj_x_pval = lm_p_val))
      }
    }
  }
  
  correlation_p_value <- correlation_p_value_matrix(
    data, relevant_vars = c(Y, X), n_cols = n_cols, cor_method = cor_method
  )
  
  corrplot(M,
           p.mat = correlation_p_value,
           method = method,
           tl.col = "black",
           tl.srt = 45,
           tl.cex = tl.cex,         # variable label size
           insig = insig,
           addCoef.col = coef_col,  # if non-NULL, corrplot prints numbers; size via number.cex
           addgrid.col = 'white',
           number.digits = 2,
           number.cex = number.cex, # on-cell number size (if using addCoef.col)
           pch.col = pch.col,
           pch = pch,
           pch.cex = pch.cex,
           col = col,
           cl.cex = cl.cex          # legend tick label size
  )$corrPos -> p1
  
  p1_sub  <- subset(p1, p.value <= 0.05)
  p1_sub2 <- subset(p1, p.value <= 0.05 & abs(corr) >= 0.7)
  
  # Custom overlays: control size with cex =
  if (nrow(p1) > 0) {
    graphics::text(p1$x, p1$y, sprintf("%.2f", p1$corr),
                   adj = c(0.5, 0), cex = coef_text_cex)
  }
  if (nrow(p1_sub) > 0) {
    graphics::text(p1_sub$x, p1_sub$y, stars.pval(p1_sub$p.value),
                   adj = c(0.5, 2), cex = star_text_cex)
  }
  if (nrow(p1_sub2) > 0) {
    graphics::text(p1_sub2$x, p1_sub2$y, sprintf("%.2f", p1_sub2$corr),
                   col = "white", adj = c(0.5, 0), cex = coef_text_cex)
    graphics::text(p1_sub2$x, p1_sub2$y, stars.pval(p1_sub2$p.value),
                   col = "white", adj = c(0.5, 2), cex = star_text_cex)
  }
  if (!is.null(adj_var)) {
    lm_extracted <- subset(lm_extracted, adj_x_pval <= 0.05)
    if (nrow(lm_extracted) > 0){
      graphics::rect(xleft = lm_extracted$x - 0.45, 
                     ybottom = lm_extracted$y - 0.45,
                     xright = lm_extracted$x + 0.45,
                     ytop = lm_extracted$y + 0.45)
    }
  }
}

corr_plot_ver2 <- function(data, X, Y, cor_method = "pearson", adj_var = NULL, dict = dict,
                               method = "color", insig = "pch", coef_col = NULL,
                               pch = 4, pch.col = "black", pch.cex = 0) {
  n_cols = length(Y)
  M <- cor(y = subset(data, select = Y),
           x = subset(data, select = X),
           use = "pairwise.complete.obs",
           method = cor_method)
  colnames(M) = as.list(dict[match(colnames(M), names(dict))])
  rownames(M) = as.list(dict[match(rownames(M), names(dict))])
  
  if (!is.null(adj_var)){
    x_vars <- rep(X, times = length(Y))
    x_coord <- rep(seq(1,length(Y)), each = length(X))
    y_vars <- rep(Y, each = length(X))
    y_coord <- rep(seq(length(X),1), times = length(Y))
    lm_extracted <- data.frame(yName = character(0),
                               xName = character(0),
                               x = numeric(0),
                               y = numeric(0),
                               adj_var = character(0),
                               adj_x_coef = numeric(0),
                               adj_x_pval = numeric(0))
    for (i in 1:length(x_vars)) {
      if(x_vars[i] != y_vars[i]) {
        lm_formula <- as.formula(paste0(y_vars[i], "~", x_vars[i], "+", adj_var))
        lm_model <- lm(lm_formula, data = data)
        lm_summary <- summary(lm_model)
        lm_coef <- lm_summary$coefficient[x_vars[i], "Estimate"]
        lm_p_val <- lm_summary$coefficient[x_vars[i], "Pr(>|t|)"]
        lm_extracted <- rbind(lm_extracted,
                              data.frame(
                                yName = y_vars[i],
                                xName = x_vars[i],
                                x = x_coord[i],
                                y = y_coord[i],
                                adj_var = adj_var,
                                adj_x_coef = lm_coef,
                                adj_x_pval = lm_p_val))
      }
    }
  }
  
  correlation_p_value <- correlation_p_value_matrix(data, relevant_vars = c(Y, X), n_cols = n_cols, cor_method = cor_method)
  
  corrplot(M,
           p.mat = correlation_p_value,
           method = method,
           tl.col = "black",
           tl.srt = 45,
           tl.cex = .8,
           insig = insig,
           addCoef.col = coef_col,
           addgrid.col = 'lightgray',
           number.digits = 2,
           pch.col = pch.col, 
           pch = pch,
           pch.cex = pch.cex)$corrPos -> p1
  p1_sub <- subset(p1, p.value <= 0.05)
  p1_sub2 <- subset(p1, p.value <= 0.05 & abs(corr) >= 0.7)
  
  if (nrow(p1_sub) > 0) {
    graphics::points(p1_sub$x, p1_sub$y+0.15, pch = 8)
  }
  if (nrow(p1_sub2) > 0) {
    graphics::points(p1_sub2$x, p1_sub2$y+0.15, pch = 8, col = "white")
  }  
  if (!is.null(adj_var)) {
    lm_combined <- left_join(p1, lm_extracted, by = c("x","y"))
    lm_combined <- subset(lm_combined, adj_x_pval <= 0.05)
    lm_subset <- subset(lm_combined, adj_x_pval <= 0.05 & abs(corr) >= 0.7)
    if (nrow(lm_combined) > 0){
      graphics::points(lm_combined$x, lm_combined$y-0.15, pch = 16)
    }
    if (nrow(lm_subset) > 0){
      graphics::points(lm_subset$x, lm_subset$y-0.15, pch = 16, col = "white")
    }
    
  }
}


pcor.plot <- function(data, X, Y, Z, cor_method = "pearson", dict = dict,
                             method = "color", insig = "pch", coef_col = NULL,
                             pch = 4, pch.col = "black", pch.cex = 0) {

  n_cols = length(Y)
  M <- data.frame(matrix(ncol = n_cols, nrow = length(X)))
  colnames(M) <- Y
  rownames(M) <- X
  P <- data.frame(matrix(ncol = n_cols, nrow = length(X)))
  colnames(P) <- Y
  rownames(P) <- X
  X_rep <- rep(X, each = length(Y))
  Y_rep <- rep(Y, times = length(X))
  dict_col <- dict[match(colnames(M), names(dict))]
  dict_row <- dict[match(rownames(M), names(dict))]
  for (i in 1:length(X_rep)) {
    if (X_rep[i] != Y_rep[i] && Y_rep[i] != paste0(X_rep[i], "_row")) {
      data_non_missing <- eval(parse(text = 
                   paste0("subset(data, !is.na(", X_rep[i], ") & !is.na(", Y_rep[i], ") & !is.na(", Z, "))")))
      pcor <- pcor.test(y = subset(data_non_missing, select = Y_rep[i]),
                x = subset(data_non_missing, select = X_rep[i]),
                z = subset(data_non_missing, select = Z),
                method = cor_method)  

      M[X_rep[i], Y_rep[i]] <- pcor$estimate
      P[X_rep[i], Y_rep[i]] <- pcor$p.value
    }
    else {
      M[X_rep[i], Y_rep[i]] <- 1.0
      P[X_rep[i], Y_rep[i]] <- 0
    }
  }

    P <- as.matrix(P)
    colnames(M) = as.list(dict_col)
    rownames(M) = as.list(dict_row)
    
    corrplot(as.matrix(M),
             is.corr = T,
             p.mat = P,
             method = method,
             tl.col = "black",
             tl.srt = 45,
             tl.cex = .8,
             insig = insig,
             addCoef.col = coef_col,
             addgrid.col = 'lightgray',
             number.digits = 2,
             pch.col = pch.col, 
             pch = pch,
             pch.cex = pch.cex)$corrPos -> p1
    p1_sub <- subset(p1, p.value <= 0.05)
    p1_sub2 <- subset(p1, p.value <= 0.05 & abs(corr) >= 0.7)
    
    if (nrow(p1) > 0) {
      graphics::text(p1$x, p1$y, sprintf("%.2f", p1$corr), adj = c(0.5, 0))
    }
    if (nrow(p1_sub) > 0) {
      graphics::text(p1_sub$x, p1_sub$y, stars.pval(p1_sub$p.value), adj = c(0.5, 2))
    }
    if (nrow(p1_sub2) > 0) {
      graphics::text(p1_sub2$x, p1_sub2$y, sprintf("%.2f", p1_sub2$corr), col = "white", adj = c(0.5, 0))
      graphics::text(p1_sub2$x, p1_sub2$y, stars.pval(p1_sub2$p.value), col = "white", adj = c(0.5, 2))
    }
}

  