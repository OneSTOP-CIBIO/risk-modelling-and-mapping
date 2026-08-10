

library(tidyverse)


perf <- read.csv("./data/projects/onestop_v02/Overview_model_performance.csv")
colnames(perf)

perf_long <- pivot_longer(perf, cols = Specificity:Kappa, names_to = "perf_ind", values_to = "val")

g <- ggplot(perf_long, aes(y = val, x = perf_ind)) + 
  geom_boxplot(fill="lightgrey",alpha=0.7) +
  geom_jitter(alpha = 0.2, size = 1.7, color = "#228B22", width = 0.07) +
  xlab("Performance index") + 
  ylab("Value") +  
  theme_bw()

plot(g)


perf_long %>% 
  #filter(perf_ind == "AUC") %>% 
  group_by(perf_ind) %>% 
  summarise(avg = mean(val, na.rm=TRUE),
            std = sd(val, na.rm=TRUE))


length(unique(perf_long$acceptedScientificName))
