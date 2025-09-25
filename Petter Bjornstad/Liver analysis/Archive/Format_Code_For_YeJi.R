#Run zlm
model_formula <- ~ group
zlm_results <- zlm(formula = model_formula, sca = sca_gene_set)
summary_zlm <- summary(zlm_results, doLRT = "groupType 2 Diabetes")
summary_dt <- summary_zlm$datatable

#Format results for barchart
fcHurdle <- merge(summary_dt[contrast=='groupType 2 Diabetes' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summary_dt[contrast=='groupType 2 Diabetes' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
m <- fcHurdle[,c("primerid","coef","fdr")]

top_pos <- as.data.frame(m) %>%
  filter(fdr<0.05) %>%
  filter(coef>0) %>% 
  dplyr::rename(Gene=primerid)
top_pos <- top_pos[c("Gene","coef","fdr")]
rownames(top_pos) <- NULL
top_pos$Direction <- "Positive"

top_neg <- as.data.frame(m) %>%
  filter(fdr<0.05) %>%
  filter(coef<0) %>% 
  dplyr::rename(Gene=primerid)
top_neg$Direction <- "Negative"

top_genes <- bind_rows(top_pos, top_neg) 
top_genes <- top_genes %>%
  mutate(Gene = factor(Gene, levels = Gene[order(Direction, coef)])) 

ggplot(top_genes, aes(x = coef, y = Gene, fill = Direction)) +
  geom_bar(stat = "identity") +
  # labs(x = "logFC", y = "Gene", title = "Overall Association of GLP-1s (Yes vs. No) and Gene Expression, Among Type 2") +
  scale_fill_manual(values = c("yellow", "purple"), labels = c("Lower Expression", "Higher Expression")) +
  theme_minimal()+
  theme(element_text(family="Times"))+
  theme(legend.position="none")