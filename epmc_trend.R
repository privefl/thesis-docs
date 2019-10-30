library(europepmc)
library(cowplot)
library(tidyverse)

ot_trend <- europepmc::epmc_hits_trend(query = "polygenic risk scores", 
                                       period = 2000:2018)
# Standard plot

ot_trend %>% 
  ggplot(aes(year, query_hits / all_hits)) + 
  geom_point() + 
  geom_line()

# Nicer plot

ot_trend %>%
  ggplot(aes(x = factor(year), y = (query_hits / all_hits))) +
  geom_col(fill = "#56B4E9", width = 0.6, alpha = 0.9) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(breaks = seq(2000, 2018, by = 2)) +
  theme_minimal_hgrid(12) +
  labs(x = "Year", y = "Proportion of all published articles") +
  ggtitle("Interest in 'Polygenic Risk Scores' research since 2000",
          "(query of 'polygenic risk scores' in PubMed Central)")
