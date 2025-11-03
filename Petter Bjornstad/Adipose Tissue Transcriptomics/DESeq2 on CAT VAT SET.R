setwd('G:/My Drive/lab files/Bryan Bergman/cohort 11-8-24 CAT VAT SAT')
library(DESeq2)
library(reshape2)
library(apeglm)
library(biomaRt)
library(factoextra)
library(FactoMineR)

library(enrichR)
library(ggplot2)
library(org.Hs.eg.db)
library("pathview")
library(europepmc)
library(ggupset)
library(clusterProfiler)
library(enrichplot)
library(cowplot)
library(DOSE)
library(ggplot2)
pre_symp = read.csv('full melted kallisto counts and tpm.csv')
#coldata = as.data.frame(table(pre_symp$sample))
#head(pre_symp)
#write.csv(coldata, file = 'sample_table.csv', row.names = F)
#go_terms = read.delim('E:/My Drive/Datasets/Human/genome files/uniprot-human-genes and goterms mapping.tab')
coldata = read.csv('sample_table.csv')

coldata$group = paste0(coldata$Tissue, '_', coldata$Sample)
table(coldata$group)

#SAT_post SAT_solo VAT_post VAT_solo 
#10        9       13       14 
coldata$genotype = as.factor(coldata$group)
row.names(coldata) = coldata$SampleID
head(pre_symp)
pre_symp$est_counts = as.integer(pre_symp$est_counts)
pre_symp$individualID = pre_symp$sample
pre_symp = pre_symp[!is.na(pre_symp$gene_symbol),]
pre_symp = pre_symp[!pre_symp$gene_symbol=='',]

#We cast into counts
cnts_mat = dcast(pre_symp, gene_symbol ~ sample, value.var = 'est_counts', fun.aggregate = mean, na.rm=T )
row.names(cnts_mat) = cnts_mat$gene_symbol
cnts_mat$gene_symbol = NULL
cnts_mat[1:15,1:15]

#start right here if you use a counts matrix
cnts_mat1 = cnts_mat[rowSums(cnts_mat) > 10,]


iris.pca <- PCA(t(cnts_mat1), graph = FALSE)
#coldata
#sample count genotype
#A1     A1 31715        A
#A2     A2 31715        A
#A4     A4 31715        A
#B1     B1 31715        B
#B2     B2 31715        B

trait_cols = unique(coldata$genotype)
names(trait_cols) = MetBrewer::met.brewer('Moreau', length(trait_cols))
#new_trait_t = 
# col.ind = coldata$genotype, # color by groups
#palette = MetBrewer::met.brewer('Lakota', length(trait_cols)),

#addEllipses = TRUE, # Concentration ellipses
pdf(file = 'PCA - all samples.pdf')
g2 = fviz_pca_ind(iris.pca, col.ind = coldata$genotype) + ggtitle(paste0('PCA plot grouping by pre vs post - ALL samples'))
print(g2)
dev.off()

coldata1 = coldata
iris.pca <- PCA(t(cnts_mat1), graph = FALSE)


trait_cols = unique(coldata1$genotype)
names(trait_cols) = MetBrewer::met.brewer('Moreau', length(trait_cols))
new_trait_t = colnames(cnts_mat1)
names(new_trait_t) = coldata1$genotype[match(new_trait_t, coldata1$SampleID)]


pdf(file = paste0('PCA - full colored.pdf'))
g2 = fviz_pca_ind(iris.pca, col.ind = names(new_trait_t), palette = MetBrewer::met.brewer('Lakota', length(trait_cols)), addEllipses=T) + ggtitle(paste0('PCA '))
print(g2)
dev.off()





####################Run DE

coldata$condition = coldata$group
table(coldata$condition)
numerator_group = 'SAT_T2D' 
denominator_group = 'SAT_healthy'
coldata = coldata[!coldata$SampleID=='77_SAT_22_S22',]
go_terms = read.delim('G:/My Drive/Datasets/Human/genome files/uniprot-human-genes and goterms mapping.tab')
run_full_DE = function(numerator_group, denominator_group){
  coldata1 = coldata[coldata$condition %in% numerator_group | coldata$condition %in% denominator_group,]
  
  cnts_mat1 = cnts_mat[rowMeans(cnts_mat) > 10,]
  cnts_mat1 = cnts_mat1[,colnames(cnts_mat1) %in% coldata1$SampleID]
  cnts_mat1 = dplyr::mutate_all(cnts_mat1, function(x) as.integer(as.character(x)))
  coldata1$genotype = factor(coldata1$condition, levels=c(numerator_group, denominator_group))
  iris.pca <- PCA(t(cnts_mat1), graph = FALSE)
  trait_cols = unique(coldata1$condition)
  names(trait_cols) = MetBrewer::met.brewer('Moreau', length(trait_cols))
  new_trait_t = colnames(cnts_mat1)
  names(new_trait_t) = coldata1$genotype[match(new_trait_t, coldata1$SampleID)]
  
  
  pdf(file = paste0(numerator_group,  ' over ', denominator_group, ' PCA.pdf'))
  g2 = fviz_pca_ind(iris.pca, col.ind = names(new_trait_t), palette = MetBrewer::met.brewer('Lakota', length(trait_cols)), addEllipses=T) + ggtitle(paste0('PCA '))
  print(g2)
  dev.off()
  
  dds <- DESeqDataSetFromMatrix(countData = cnts_mat1,
                                colData = coldata1,
                                design = ~ genotype)
  
  
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("genotype",  numerator_group, denominator_group))
  results_table =as.data.frame(res)
  res1 = results_table[order(results_table$pvalue, decreasing = F),]
  head(res1)
  test1 = cnts_mat1[grepl('RPS4Y1', row.names(cnts_mat1)),]
  colnames(test1) = coldata$group[match(colnames(test1), coldata$SampleID)]
  res1$GOBP = go_terms$Gene.ontology..biological.process.[match(row.names(res1), go_terms$Gene.names...primary..)]
  res1$GOCC = go_terms$Gene.ontology..cellular.component.[match(row.names(res1), go_terms$Gene.names...primary..)]
  res1$GOMF = go_terms$Gene.ontology..molecular.function.[match(row.names(res1), go_terms$Gene.names...primary..)]
  
  write.csv(res1, file =  paste0( 'DESeq2 Full results ', numerator_group,  ' over ', denominator_group, '.csv'))
  
  
  #need to play around to assessing proper thresholds.  
  res1 = na.omit(res1)
  res1$gene_symbol = row.names(res1)
  label_key = res1$gene_symbol[1:20]
  
  res1$label2 = ifelse(res1$gene_symbol %in% label_key, paste0(res1$gene_symbol), '')
  table(res1$label2)
  res1$label_col1 = ifelse(res1$log2FoldChange>0, 'firebrick3', 'darkorchid')
  res1$label_col2 = ifelse(res1$padj<0.1, paste0(res1$label_col1), 'gray74')
  
  library(ggrepel)
  #Number of genes which will be labelled
  #Volcano plot
  #change labels here too
  pdf(file = paste0( 'Volcano plot ', numerator_group,  ' over ', denominator_group, '.pdf'))
  volc1 = ggplot(res1, aes(x=log2FoldChange, y=-log10(padj))) + theme_classic() +
    geom_point(aes(x=log2FoldChange, y=-log10(padj)), color=res1$label_col2) +
    geom_label_repel(aes(x=log2FoldChange, y=-log10(padj), label = res1$label2), color = res1$label_col2, size = 2, label.size=NA, box.padding = 0.8, point.padding = 0.5, max.overlaps = Inf, alpha = .6, segment.color = 'grey50')  +   ggtitle(paste0('DESeq2 Full results ', numerator_group,  ' over ', denominator_group)) + ylab('-log10(padj)')
  print(volc1)
  dev.off()
  
  heat_genes = res1$gene_symbol[1:20]
  
  htmap = cnts_mat1[row.names(cnts_mat1) %in% heat_genes,]
  colnames(htmap) = coldata1$condition[match(colnames(htmap), coldata1$sample)]
  ht1 = scale(t(htmap))
  pdf(file = paste0('heatmap of top genes ', numerator_group,  ' over ', denominator_group, '.pdf'))
  ggt1 = pheatmap::pheatmap(ht1, color=colorRampPalette(c("dodgerblue4", "white", "darkorange"))(100))
  print(ggt1)
  dev.off()
  deg_set1 = res1[res1$pvalue<0.1,]
  fc_dko = deg_set1$log2FoldChange
  library(clusterProfiler)
  library(enrichR)
  library(cowplot)
  names(fc_dko) <- deg_set1$gene_symbol
  fc_dko <- sort(fc_dko, decreasing = TRUE)
  fc_dko = fc_dko[!duplicated(names(fc_dko))]
  organism = "org.Hs.eg.db"
  gse <-gseGO(
    geneList=fc_dko,
    ont = "ALL",
    OrgDb= organism,
    keyType = "SYMBOL",
    exponent = 1,
    minGSSize = 2,
    maxGSSize = 500,
    eps = 0,
    pvalueCutoff = 1,
    pAdjustMethod = "BH") 
  
  gse_results = as.data.frame(gse@result)
  write.csv(gse_results, file = paste0('full pathway enrichment file for ', numerator_group,  ' over ', denominator_group, '.csv'), row.names = F)
  ## gsego https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/gseGO
  pdf(file =paste0('GSEA output ', numerator_group,  ' over ', denominator_group, '.pdf'))
  g3 = dotplot(gse, showCategory=15, split=".sign", color="pvalue") + facet_grid(.~.sign) + ggtitle(paste0( 'GSEA terms')) + theme(axis.text.y=element_text(size=4))
  print(g3)
  dev.off()
  
  pdf(file = paste0('GSEA network graph ', numerator_group,  ' over ', denominator_group, '.pdf'))
  x2<- pairwise_termsim(gse)
  g4 = emapplot(x2, showCategory = 15, color = "pvalue", cex_label_category=1)+ ggtitle(paste0('GSEA Network graph'))
  print(g4)
  dev.off()
  
  edox2 <- pairwise_termsim(gse)
  p1 <- treeplot(edox2)
  p2 <- treeplot(edox2, hclust_method = "average")
  gg42 = aplot::plot_list(p1, p2, tag_levels='A')
  pdf(file = paste0('pathway tree graph ', numerator_group,  ' over ', denominator_group, '.pdf'))
  print(gg42)
  dev.off()
  pdf(file = paste0('Ridgeplot ', numerator_group,  ' over ', denominator_group, '.pdf')) 
  g1 = ridgeplot(gse, showCategory = 20) + labs(x = "logFC") + ggtitle('Distribution of fold changes')+ theme(axis.text.y=element_text(size=5))
  print(g1)
  dev.off()
  #
  gse@result$Description[1:20]
  ##cnetplot gene concept networks
  pdf(file = paste0('cnet gene concept network ', numerator_group,  ' over ', denominator_group, '.pdf'), width = 11, height = 8)
  g22 = cnetplot(gse, foldChange=fc_dko, 
                 showCategory = gse@result$Description[1:3], 
                 circular = TRUE, colorEdge = TRUE) +
    ggtitle('Gene-Concept Network from GSEA-GO')
  print(g22)
  dev.off()
}
dev.off()

setwd('G:/My Drive/lab files/Bryan Bergman/cohort 11-8-24 CAT VAT SAT/DESeq2 outputs')
table(coldata$condition)
run_full_DE('CAT_T2D', 'CAT_healthy' )

run_full_DE('SAT_T2D', 'SAT_healthy' )
run_full_DE('VAT_T2D', 'VAT_healthy' )



###compare de
setwd('G:/My Drive/lab files/Bryan Bergman/cohort 11-8-24 CAT VAT SAT/DESeq2 outputs/')
sat_de = read.csv('./SAT/DESeq2 Full results SAT_T2D over SAT_healthy.csv')
vat_de = read.csv('./VAT/DESeq2 Full results VAT_T2D over VAT_healthy.csv')
cat_de = read.csv('./CAT/DESeq2 Full results CAT_T2D over CAT_healthy.csv')

library(ggVennDiagram)
pvalthr = 0.05
gene_set = list(SAT_DE = sat_de$X[sat_de$pvalue<pvalthr & sat_de$log2FoldChange>0], VAT_DE = vat_de$X[vat_de$pvalue<pvalthr & vat_de$log2FoldChange>0], CAT_DE = cat_de$X[cat_de$pvalue<pvalthr & cat_de$log2FoldChange>0])
pdf(file = paste0('intersection of T2D increased DEGs pval less ', pvalthr,'.pdf'))
ggVennDiagram::ggVennDiagram(gene_set)
dev.off()

gene_set = list(SAT_DE = sat_de$X[sat_de$pvalue<pvalthr & sat_de$log2FoldChange<0], VAT_DE = vat_de$X[vat_de$pvalue<pvalthr & vat_de$log2FoldChange<0], CAT_DE = cat_de$X[cat_de$pvalue<pvalthr & cat_de$log2FoldChange<0])
pdf(file = paste0('intersection of healthy increased DEGs pval less ', pvalthr,'.pdf'))
ggVennDiagram::ggVennDiagram(gene_set)
dev.off()


library(enrichR)

gg1 = sat_de$X[sat_de$pvalue<pvalthr]
gg2 = vat_de$X[vat_de$pvalue<pvalthr]
gg2 = gg2[!gg2 %in% gg1]
gg1 = gg1[!gg1 %in% gg2]
dbs = listEnrichrDbs()
dbs = dbs[order(dbs$geneCoverage, decreasing = T),]
setEnrichrSite("Enrichr")
dbs1 <- c("TRANSFAC_and_JASPAR_PWMs", "GO_Biological_Process_2023", "Reactome_2022", "Diabetes_Perturbations_GEO_2022")

#first for positive
enriched <- enrichr(gg1, dbs1)
g1 = plotEnrich(enriched[[1]], showTerms = 10, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(paste0(names(enriched[1]))) + theme(axis.text=element_text(size=5),
                                                                                                                                             axis.title=element_text(size=6,face="bold")) + theme(plot.title = element_text(size = 7, face = "bold"))
g2 = plotEnrich(enriched[[2]], showTerms = 10, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(paste0(names(enriched[2])))+ theme(axis.text=element_text(size=5),
                                                                                                                                            axis.title=element_text(size=6,face="bold"))+ theme(plot.title = element_text(size = 7, face = "bold"))

g3 = plotEnrich(enriched[[3]], showTerms = 10, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(paste0(names(enriched[3])))+ theme(axis.text=element_text(size=5),
                                                                                                                                            axis.title=element_text(size=6,face="bold"))+ theme(plot.title = element_text(size = 7, face = "bold"))

g4 = plotEnrich(enriched[[4]], showTerms = 10, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(paste0(names(enriched[4])))+ theme(axis.text=element_text(size=5),
                                                                                                                                            axis.title=element_text(size=6,face="bold"))+ theme(plot.title = element_text(size = 7, face = "bold"))

pdf(file = paste0('SAT ONLY pathways.pdf'))
gg2 = gridExtra::grid.arrange(g1, g2, g3, g4)
print(gg2)
dev.off()


sat_de$vatLFC = vat_de$log2FoldChange[match(sat_de$X, vat_de$X)]
sat_de$vatPvalue = vat_de$pvalue[match(sat_de$X, vat_de$X)]
sat_de$fcDELTA = sat_de$vatLFC-sat_de$log2FoldChange  
ss1 = sat_de[order(sat_de$fcDELTA, decreasing = T),]
ss1 = ss1[ss1$pvalue<0.05,]
ss1 = ss1[ss1$vatPvalue<0.05,]

new_df = as.data.frame(cbind(ss1$X, ss1$log2FoldChange, ss1$pvalue, ss1$vatLFC, ss1$vatPvalue, ss1$fcDELTA))
colnames(new_df) = c('gene_symbol', 'SAT_FC', 'SAT_Pvalue', 'VAT_FC', 'VAT_Pvalue', 'FC_delta')
heatdf1 = new_df[1:15,]
heatdf2 = new_df[54:69,]
heatdf = as.data.frame(rbind(heatdf1, heatdf2))

row.names(heatdf) = heatdf$gene_symbol
#heatdf = heatdf[!heatdf$gene_symbol=='DAGLA',]
heatdf$gene_symbol=NULL
heatdf$SAT_Pvalue=NULL
heatdf$VAT_Pvalue=NULL
heatdf$FC_delta=NULL
heatdf$SAT_FC = scale(as.numeric(heatdf$SAT_FC))
heatdf$VAT_FC = scale(as.numeric(heatdf$VAT_FC))

pdf(file = paste0('heatmap of top genes by FC pvalue filtered.pdf'))
ggt1 = pheatmap::pheatmap(heatdf, color=colorRampPalette(c("darkorchid", "white", "orange"))(100))
print(ggt1)
dev.off()

sec_proteins = read.delim('G:/My Drive/Datasets/Human/genome files/human secreted proteins.tab')
new_df$secreted = ifelse(new_df$gene_symbol %in% sec_proteins$Gene.names...primary.., 'Secreted', 'Not')
table(new_df$secreted)
write.csv(new_df, file = 'delta FC of DEGs VAT - SAT pvalue filtered.csv', row.names = F)
head(new_df)


###plot specific genes







genesetplot = c( 'COL3A1', 'COL6A2', 'LAMB3', 'EMILIN1', 'SDC4', 'MMP2', 'C1QTNF1', 'CFB', 'CSF1', 'CXCL8','TSLP', 'VEGFA')
sat_de$condition = paste0('SAT')
vat_de$condition = paste0('VAT')

full_de = as.data.frame(rbind(sat_de, vat_de))

bar_genes = full_de[full_de$X %in% genesetplot,]
bar_genes$X = factor(bar_genes$X, levels = c( 'COL3A1', 'COL6A2', 'LAMB3', 'EMILIN1', 'SDC4', 'MMP2', 'C1QTNF1', 'CFB', 'CSF1', 'CXCL8','TSLP', 'VEGFA'))
colnames(bar_genes)
pdf(file = 'LFC for select genes final set.pdf')
ggplot(bar_genes, aes(x=X, y=log2FoldChange, fill=condition)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width=.2,
                position=position_dodge(.9))  + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +scale_fill_manual(values = c('gray95', 'gray25')) + xlab('') + ggtitle('Log2FC for key genes')
dev.off()


###
###plot specific genes

genesetplot = c( 'COL3A1', 'COL6A2', 'LAMB3', 'EMILIN1', 'SDC4', 'MMP2', 'C1QTNF1', 'CFB','C5AR1', 'CSF1', 'CXCL8','TSLP', 'VEGFA')
sat_de$condition = paste0('SAT')
vat_de$condition = paste0('VAT')

full_de = as.data.frame(rbind(sat_de, vat_de))

bar_genes = full_de[full_de$X %in% genesetplot,]
bar_genes$X = factor(bar_genes$X, levels = c( 'COL3A1', 'COL6A2', 'LAMB3', 'EMILIN1', 'SDC4', 'MMP2', 'C1QTNF1', 'CFB', 'C5AR1', 'CSF1', 'CXCL8','TSLP', 'VEGFA'))
colnames(bar_genes)
table(bar_genes$X)
pdf(file = 'LFC for select genes final set WITH C5AR1.pdf')
ggplot(bar_genes, aes(x=X, y=log2FoldChange, fill=condition)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width=.2,
                position=position_dodge(.9))  + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +scale_fill_manual(values = c('gray95', 'gray25')) + xlab('') + ggtitle('Log2FC for key genes')
dev.off()














genesetplot = c( 'C1QTNF1', 'C1R', 'FSTL3', 'SEMA3C', 'SDC4', 'LGALS3BP', 'LAMC1', 'MMP2')
sat_de$condition = paste0('SAT')
vat_de$condition = paste0('VAT')

full_de = as.data.frame(rbind(sat_de, vat_de))

bar_genes = full_de[full_de$X %in% genesetplot,]
bar_genes$X = factor(bar_genes$X, levels =c( 'C1QTNF1', 'C1R', 'FSTL3', 'SEMA3C', 'SDC4', 'LGALS3BP', 'LAMC1', 'MMP2'))
colnames(bar_genes)
pdf(file = 'LFC for select genes revised set.pdf')
ggplot(bar_genes, aes(x=X, y=log2FoldChange, fill=condition)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE), width=.2,
                position=position_dodge(.9))  + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +scale_fill_manual(values = c('gray95', 'gray25')) + xlab('') + ggtitle('Log2FC for key genes')
dev.off()



#select pathwyas and plot - you can ignore
pathways = read.delim('G:/My Drive/Datasets/Rat/uniprotkb_rat_AND_model_organism_10116_2024_07_19.tsv')

colnames(pathways)
pathway_ofIGOBP = 'defense response to other organism' 

pathplot1 = function(pathway_ofIGOBP){
  nn1 = pathways[grepl(pathway_ofIGOBP, pathways$Gene.Ontology..biological.process.),]
  filtered_cnts = cnts_mat[row.names(cnts_mat) %in% nn1$Gene.Names..primary.,]
  
  plot_cnts = reshape2::melt(as.matrix(filtered_cnts))
  colnames(plot_cnts) = c('gene_symbol', 'sampleID', 'TPM')
  plot_cnts$group = sampel_table$Condition[match(plot_cnts$sampleID, sampel_table$RNAseq.Sample.Name)]
  plot_cnts$group = gsub('*', 'STAR', fixed=T, plot_cnts$group)
  
  library(ggpubr)
  my_comparisons <- list( c('M2-Norm-1M NRVM + 300K BMDM', 'M2-Norm-1M NRVM + 0 BMDM' ), c('M2-Hyp-1M NRVM + 300K BMDM', 'M2-Hyp-1M NRVM + 0 BMDM'), c('M2-Hyp-1M NRVM + 300K BMDM', 'M0-Hyp-1M NRVM + 0 BMDM'), c('M0-Hyp-1M NRVM + 0 BMDM', 'M0-Norm-1M NRVM + 0 BMDM') , c('M2-Hyp-1M NRVM + 0 BMDM', 'M2STAR-Hyp-1M NRVM + 0 BMDM') )
  
  gg33 = ggboxplot(plot_cnts, x = "group", y = "TPM",
                   color = "group",
                   add = "jitter") +
    stat_compare_means(comparisons = my_comparisons)+   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position="none")+ ggtitle(paste0(pathway_ofIGOBP, 'expression'))
  pdf(file = paste0(pathway_ofIGOBP, 'expression across conidtions.pdf'))
  
  print(gg33)
  dev.off()
}
setwd("G:/My Drive/lab files/anya_grosberg/macrphage_cardiomyocyte/RNA-seq/DE_seq2 outputs")
pathplot1('defense response to other organism' )
pathplot1('extracellular space' )

pathplot1('DNA damage response' )
pathplot1('lymphocyte activation' )

pathplot1('response to estradiol' )

pathplot1('ion transport' )
