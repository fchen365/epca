## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, 
                      tidy = TRUE, 
                      tidy.opts = list(comment = FALSE))
library(epca)
library(Matrix)
library(tidyverse)

## ----simu---------------------------------------------------------------------
## simulate a rank-5 data matrix with some additive Gaussian noise 
n <- 300
p <- 50
k <- 5 ## rank
Z <- shrinkage(svd(matrix(runif(n * k), n, k))$u, gamma = sqrt(n))
B <- diag(5) * 3
Y <- shrinkage(svd(matrix(runif(p * k), p, k))$u, gamma = sqrt(p))
E <- matrix(rnorm(n * p, sd = .01), n, p)
X <- scale(Z %*% B %*% t(Y) + E)

## ----sca----------------------------------------------------------------------
## perform sparse PCA
s.sca <- sca(X, k = 5)
s.sca

## ----sma----------------------------------------------------------------------
## perform sparse matrix approximation
s.sma <- sma(X, k = 5)
s.sma

## ----3pc----------------------------------------------------------------------
## find 3 sparse PCs
data("pitprops", package = "epca")
s.sca <- sca(pitprops, k = 3, gamma = 4.5)
print(s.sca, verbose = TRUE)

## ----6pc----------------------------------------------------------------------
## find 6 sparse PCs
s.sca <- sca(pitprops, 6, gamma = 6)
print(s.sca, verbose = TRUE)

## ----import results, echo=FALSE-----------------------------------------------
load("scrnaseq.rda")

## ----import scRNA-seq data, eval=FALSE----------------------------------------
#  # library(scRNAseq)
#  dat <- BaronPancreasData('human')
#  # dim(dat) ## 20125  8569
#  gene.select <- !!apply(counts(dat), 1, sd) ## remove non-variance gene
#  label.select <- colData(dat) %>%
#    data.frame() %>%
#    dplyr::count(label) %>%
#    filter(n > 100)
#  #   label                  n
#  # 1 acinar               958
#  # 2 activated_stellate   284
#  # 3 alpha               2326
#  # 4 beta                2525
#  # 5 delta                601
#  # 6 ductal              1077
#  # 7 endothelial          252
#  # 8 gamma                255
#  # 9 quiescent_stellate   173
#  dat1 <- dat[gene.select, colData(dat)$label %in% label.select$label]

## ----extract count matrix, eval=FALSE-----------------------------------------
#  count <- counts(dat1)
#  # dim(count) ## 17499  8451
#  # length(count@i) / length(count) ## %(nnz)
#  ## 10.80605% non-zeros

## ----extract cell label, eval=FALSE-------------------------------------------
#  label <- setNames(factor(dat1$label), colnames(dat1))

## ----apply sca to scRNA-seq, eval=FALSE---------------------------------------
#  scar <- sca(t(count), k = 9, gamma = 12,
#               center = F, scale = F,
#               epsilon = 1e-3)

## ----number of non-zeros------------------------------------------------------
n.gene <- apply(!!scar$loadings, 2, sum)
n.gene

## ----plot, fig.width=6, fig.height = 6, fig.cap="Scores of sparse gene principal components (PCs) stratified by cell types."----
scar$scores %>%
  reshape2::melt(varnames = c("cell", "PC"), 
                 value.name = "scores") %>% 
  mutate(PC = factor(PC), label = label[cell]) %>%
  ggplot(aes(PC, scores / 1000, fill = PC)) +
  geom_boxplot(color = "grey30", outlier.shape = NA, 
               show.legend = FALSE) + 
  labs(x = "gene PC", y = bquote("scores ("~10^3~")")) + 
  scale_x_discrete(labels = 1:9) + 
  facet_wrap(~ label, nrow = 3) + 
  scale_fill_brewer(palette = "Set3") +
  theme_classic() 

