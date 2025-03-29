library(ggplot2)
library(dplyr)
res <- c("N.W. Europe\n(other individuals)" = 100, "S. Europe" = 86, "Middle East" = 72,
         "India" = 65, "E. Asia" = 49, "W. Africa" = 18) %>% 
  tibble::enframe() %>% 
  mutate(name = ordered(name, levels = name))
ggplot(res) +
  geom_col(aes(name, value / 100), color = "black", fill = "chartreuse4", alpha = 0.3,
           width = 0.7) +
  theme_bw(18) +
  # geom_text(aes(x = 5, y = 0.95, label = "Source: F. Privé et al, AJHG (2022)")) +
  geom_hline(yintercept = 1, linetype = 3) +
  labs(x = "PRS testing population (Pop)", 
       y = expression(frac(R[Pop]^2, R["NW Eur"]^2)),
       subtitle = "PRS training population: N.W. Europe") +
  theme(axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(angle = 0, vjust = 0.5, margin = margin(r = 15)),
        panel.grid.major.x = element_blank()) +
  geom_text(aes(name, 0.075, label = paste0(value, " %")), size = 5) +
  scale_y_continuous(limits = c(0, 1.01), expand = c(0, 0), breaks = 0:10 / 10)

            