library(bigsnpr)

celiac <- snp_attach("../paper1-packages/backingfiles/celiacQC.rds")
y <- celiac$fam$affection - 1
G <- celiac$genotypes
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos
NCORES <- nb_cores()

big_counts(G, ind.col = 1:10) # OK

#### GWAS ####
gwas1 <- big_univLogReg(G, y, ncores = NCORES)

library(ggplot2)
snp_qq(gwas1) + ylim(NA, 10)
# ggsave("figures/qqplot1.png", scale = 1/75, width = 670, height = 520)

svd <- snp_autoSVD(G, CHR, POS, ncores = NCORES)
gwas2 <- big_univLogReg(G, y, covar.train = svd$u, ncores = NCORES)
snp_qq(gwas2) + ylim(NA, 10)
# ggsave("figures/qqplot2.png", scale = 1/75, width = 670, height = 520)

#### SVD ####
svd1 <- big_randomSVD(G, big_scale(), ncores = NCORES)
# Get population from external files
pop.files <- list.files(path = "../thesis-celiac/Dubois2010_data/",
                        pattern = "cluster_*",
                        full.names = TRUE)
pop <- snp_getSampleInfos(celiac, pop.files)[[1]]
pop.names <- c("Netherlands", "Italy", "UK1", "UK2", "Finland")
plot(svd1, type = "scores") +
  aes(color = pop.names[pop]) +
  labs(color = "Population")
# ggsave("figures/PC-1-2.png", scale = 1/75, width = 670, height = 520)
plot(svd1, type = "scores", scores = 3:4) +
  aes(color = pop.names[pop]) +
  labs(color = "Population")
# ggsave("figures/PC-3-4.png", scale = 1/75, width = 670, height = 520)

## With clumping
ind.keep <- snp_clumping(G, CHR, ncores = NCORES)
svd2 <- big_randomSVD(G, big_scale(), ind.col = ind.keep, ncores = NCORES)
plot(svd2, type = "scores", scores = 3:4) +
  aes(color = pop.names[pop]) +
  labs(color = "Population")
# ggsave("figures/PC2-3-4.png", scale = 1/75, width = 670, height = 520)

theone <- ind.keep[which.max(abs(svd2$v[, 4]))]
CHR[theone]
plot(svd2, type = "scores", scores = 3:4) +
  aes(color = as.factor(G[, theone])) +
  labs(color = "Genotype")
# ggsave("figures/PC3-3-4.png", scale = 1/75, width = 670, height = 520)

plot(svd2, type = "loadings", loadings = 1:6, coeff = 0.8)
# ggsave("figures/load3-3-4.png", scale = 1/75, width = 900, height = 700)

# When removing long-range LD regions
ind.keep2 <- snp_clumping(G, CHR, exclude = snp_indLRLDR(CHR, POS), 
                          ncores = NCORES)
svd3 <- big_randomSVD(G, big_scale(), ind.col = ind.keep2, ncores = NCORES)
plot(svd3, type = "scores", scores = 3:4) +
  aes(color = pop.names[pop]) +
  labs(color = "Population")
# ggsave("figures/PC4-3-4.png", scale = 1/75, width = 670, height = 520)

#### PRS ####
ind.train <- sort(sample(nrow(G), size = 12e3))
ind.test <- setdiff(rows_along(G), ind.train)

gwas.train <- big_univLogReg(G, y01.train = y[ind.train], 
                             ind.train = ind.train,
                             covar.train = svd$u[ind.train, ],
                             ncores = NCORES)

ind.keep <- snp_clumping(G, CHR, ind.row = ind.train, 
                         S = abs(gwas.train$score),
                         ncores = NCORES)
summary(lpS.keep <- -predict(gwas.train)[ind.keep])
thrs <- c(seq(0, 10, by = 0.5), seq(15, 200, by = 5), seq(220, 540, by = 20))
nb.pred <- sapply(thrs, function(thr) sum(lpS.keep > thr))
prs <- snp_PRS(G, betas.keep = gwas.train$estim[ind.keep],
               ind.test = ind.test,
               ind.keep = ind.keep,
               lpS.keep = lpS.keep,
               thr.list = thrs)

aucs <- apply(prs, 2, AUC, target = y[ind.test])
library(ggplot2)
bigstatsr:::MY_THEME(qplot(nb.pred, aucs)) +
  geom_line() +
  scale_x_log10() +
  labs(x = "Number of predictors", y = "AUC")
# ggsave("figures/AUC-PRS.png", scale = 1/75, width = 600, height = 540)


splog <- big_spLogReg(G, y01.train = y[ind.train], 
                      ind.train = ind.train,
                      covar.train = svd$u[ind.train, ],
                      alpha = 0.5)
nb.pred2 <- splog$p - splog$rejections
preds <- predict(splog, G, ind.row = ind.test,
                 covar.row = svd$u[ind.test, ])
aucs2 <- apply(preds, 2, AUC, target = y[ind.test])
bigstatsr:::MY_THEME(qplot(nb.pred2[-1], aucs2[-1])) +
  # geom_line() +
  scale_x_log10() + 
  scale_y_continuous(breaks = seq(0.8, 0.9, by = 0.02)) +
  labs(x = "Number of predictors", y = "AUC")
# ggsave("figures/AUC-spLog.png", scale = 1/75, width = 800, height = 560)

cmsa <- big_CMSA(big_spLogReg, AUC, G, 
                 y.train = y[ind.train],
                 ind.train = ind.train, 
                 covar.train = svd$u[ind.train, ],
                 K = 12, ncores = NCORES,
                 alpha = 0.5)
preds2 <- predict(cmsa, G, ind.row = ind.test,
                  covar.row = svd$u[ind.test, ])
last_plot() +
  geom_hline(yintercept = AUC(preds2, y[ind.test]), col = "red")
# ggsave("figures/AUC-spLog-CMSA.png", scale = 1/75, width = 800, height = 560)

# library(dplyr)
# data.frame(
#   Population = pop.names[pop[ind.test]], 
#   pred1 = prs[, which.max(aucs)],
#   pred2 = preds2, 
#   true = y[ind.test]
# ) %>%
#   group_by(Population) %>%
#   summarise(
#     AUC_PRS = AUC(pred1, true),
#     AUC_CMSA = AUC(pred2, true)
#   )
  

