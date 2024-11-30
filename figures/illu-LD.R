corr <- readRDS(url("https://www.dropbox.com/s/65u96jf7y32j2mj/spMat.rds?raw=1"))
corrT <- as(corr, "dgTMatrix")
# THR_R2 <- 0.02
ind <- 81:116
library(magrittr)
df2 <- corrT[ind, ind] %>%
  { dplyr::filter(data.frame(i = .@i + ind[1], j = .@j + ind[1], r2 = .@x^2),
                  abs(i - j) <= 10) }

geom_square <- function(min, max) {
  geom_segment(data = data.frame(i = c(min, min, max, max), 
                                 j = c(min, min, max, max)),
               aes(xend = c(min, max, min, max), 
                   yend = c(max, min, max, min)),
               linetype = 3)
}

library(ggplot2)
ggplot(df2, aes(i, j)) +
  theme_minimal(16) +
  geom_raster(aes(fill = r2), alpha = 0.6) +
  scale_fill_viridis_c(direction = -1, breaks = c(0, 0.02, 0.1, 0.3, 0.5, 0.8, 1)) +
  coord_equal() +
  scale_x_continuous(breaks = 1:1000 - 0.5, minor_breaks = NULL) +
  scale_y_reverse(breaks = 1:1000 - 0.5, minor_breaks = NULL) +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank()) +
  labs(x = NULL, y = NULL, fill = expression(r^2)) +
  geom_square(80.5, 84.5) +
  geom_square(84.5, 95.5) +
  geom_square(95.5, 106.5) +
  geom_square(106.5, 110.5) +
  geom_square(110.5, 116.5) +
  geom_segment(aes(xend = c(117.5, 107.5), yend = c(107.5, 117.5)),
               data = data.frame(i = c(89.5, 79.5), j = c(79.5, 89.5)),
               linetype = 2, color = "red") +
  guides(fill = guide_colorbar(barheight = 15, ticks.linewidth = 2))

ggsave("figures/illu-LD.png", width = 7, height = 6)
