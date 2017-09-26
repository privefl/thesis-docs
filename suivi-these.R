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
