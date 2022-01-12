# computed with snp_cor() (do not hesitate to be stringent on alpha)
corr <- readRDS(url("https://www.dropbox.com/s/65u96jf7y32j2mj/spMat.rds?raw=1"))

library(tidyverse)

## Transform sparse representation into (i,j,x) triplets
corrT <- as(corr, "dgTMatrix")
upper <- (corrT@i <= corrT@j)
df <- data.frame(
  i = corrT@i[upper], 
  j = corrT@j[upper],
  r2 = corrT@x[upper]
) %>%
  mutate(y = (j - i) / 2) %>% 
  filter(220 < (i + y), (i + y) < 295)

ggplot(df) +
  theme_bw(15) +
  geom_point(aes(i + y, y, color = r2)) +
  coord_fixed() + 
  scale_colour_gradient(low = "white", high = "black") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Position of variant on the genome",
       y = NULL, color = expression(r^2)) +
  scale_alpha(guide = 'none')
# ggsave("figures/fig-LD.png", width = 8, height = 4)
