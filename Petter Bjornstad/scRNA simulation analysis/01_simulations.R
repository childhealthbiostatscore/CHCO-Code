set.seed(1)

# Core
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(splatter)
  library(Matrix)
  library(nebula)
  library(Seurat)
  library(scater)
})

# Helper: quick checks
check_vec <- function(x, nm) {
  cat(nm, "class:", class(x), " length:", length(x), "\n")
  cat("  NA fraction:", mean(is.na(x)), "\n")
}

# --- simulation settings ---
n_genes <- 2000
n_subjects <- 20
cells_per_subject <- 6000
n_cells <- n_subjects * cells_per_subject

# Splatter params
params <- newSplatParams()  # default-ish params; later you can estimate from real data
params <- setParams(params,
                    nGenes = n_genes,
                    batchCells = rep(cells_per_subject, n_subjects),  # each batch = one subject
                    group.prob = c(0.5, 0.5),        # within each subject, half pre half post (approx)
                    de.prob = 0
)

sce <- splatSimulate(params,
                     method = "groups",
                     verbose = T)

sce <- logNormCounts(sce)
sce <- runPCA(sce)
plotPCA(sce, colour_by = "Group")

# "Batch" in splatter colData is a factor; we use it as subject id
meta <- as.data.frame(colData(sce))
meta$subject <- as.integer(factor(meta$Batch))
meta$time <- factor(meta$Group, levels = levels(meta$Group), labels = c("pre", "post"))

table(meta$time)
table(meta$subject)[1:5]

de_prop <- 0.10
lfc <- log(1.25)

n_de <- ceiling(de_prop * n_genes)
de_genes <- sample(seq_len(n_genes), n_de)

post_cells <- which(meta$time == "post")

counts <- assay(sce, "counts")
# Multiply counts for DE genes in post cells
counts[de_genes, post_cells] <-
  round(counts[de_genes, post_cells] * exp(lfc))

# run nebula
pred <- model.matrix(~ time, data = meta)

lib <- Matrix::colSums(counts)
lib[lib <= 0] <- 1
offset <- log(lib)

fit <- nebula(
  count  = counts,
  id     = meta$subject,
  pred   = pred,
  offset = offset,
  method = "LN",
  ncore  = 2,
  verbose = FALSE
)

sm <- as.data.frame(fit$summary)

coef_name <- colnames(pred)[grepl("^time", colnames(pred))]

logfc_col <- paste0("logFC_", coef_name)
p_col     <- paste0("p_", coef_name)

res <- data.frame(
  gene = rownames(sm),
  est_lfc = sm[[logfc_col]],
  pval = sm[[p_col]]
)

summary(res$pval)
mean(is.na(res$pval))
sum(res$pval < 0.05, na.rm = TRUE)

# Add groups A/B and inject a group × time effect

# Assign half subjects to A, half to B
subjects <- sort(unique(meta$subject))
n_subjects <- length(subjects)

group_by_subject <- rep(c("A", "B"), length.out = n_subjects)
group_by_subject <- sample(group_by_subject)  # randomize
names(group_by_subject) <- subjects

meta$group <- factor(group_by_subject[as.character(meta$subject)], levels = c("A","B"))

table(meta$group)
table(meta$group, meta$time)[, ]

# Re-simulate a clean dataset with same params
sce2 <- splatSimulate(params, verbose = FALSE)
counts2 <- assay(sce2, "counts")
meta2 <- as.data.frame(colData(sce2))

meta2$subject <- as.integer(factor(meta2$Batch))

# Same pre/post assignment scheme
meta2$time <- unlist(
  lapply(split(seq_len(nrow(meta2)), meta2$subject), function(idx) {
    n <- length(idx)
    sample(rep(c("pre", "post"), length.out = n))
  })
)
meta2$time <- factor(meta2$time, levels = c("pre","post"))

# Same group assignment scheme (reuse mapping so only counts changed)
meta2$group <- factor(group_by_subject[as.character(meta2$subject)], levels = c("A","B"))

table(meta2$group, meta2$time)

de_prop <- 0.10
lfc_int <- log(1.25)

n_genes <- nrow(counts2)
n_de <- ceiling(de_prop * n_genes)
de_genes <- sample(seq_len(n_genes), n_de)

B_post_cells <- which(meta2$group == "B" & meta2$time == "post")

counts2[de_genes, B_post_cells] <-
  round(counts2[de_genes, B_post_cells] * exp(lfc_int))

# sanity check: the effect should appear in B post vs B pre, but not in A
B_pre_cells <- which(meta2$group == "B" & meta2$time == "pre")
A_post_cells <- which(meta2$group == "A" & meta2$time == "post")
A_pre_cells  <- which(meta2$group == "A" & meta2$time == "pre")

ratio_B <- rowMeans(counts2[de_genes, B_post_cells]) / rowMeans(counts2[de_genes, B_pre_cells])
ratio_A <- rowMeans(counts2[de_genes, A_post_cells]) / rowMeans(counts2[de_genes, A_pre_cells])

summary(ratio_B)
summary(ratio_A)

pred2 <- model.matrix(~ group * time, data = meta2)

lib2 <- Matrix::colSums(counts2)
lib2[lib2 <= 0] <- 1
offset2 <- log(lib2)

fit2 <- nebula(
  count  = counts2,
  id     = meta2$subject,
  pred   = pred2,
  offset = offset2,
  method = "LN",
  ncore  = 2,
  verbose = FALSE
)

sm2 <- as.data.frame(fit2$summary)

# Interaction coefficient name should be "groupB:timepost"
coef_int <- "groupB:timepost"

res2 <- data.frame(
  gene = rownames(sm2),
  est_lfc = sm2[[paste0("logFC_", coef_int)]],
  pval = sm2[[paste0("p_", coef_int)]]
)

summary(res2$pval)
mean(is.na(res2$pval))
sum(res2$pval < 0.05, na.rm = TRUE)

truth <- rep(FALSE, n_genes)
truth[de_genes] <- TRUE

power <- mean(res2$pval[truth] < 0.05)
type1 <- mean(res2$pval[!truth] < 0.05)

c(power = power, type1 = type1)