cat("\014")
rm(list=ls())
library(dplyr)
library(ggplot2)

data(CO2)

hist(CO2$uptake, prob=T)
rug(CO2$uptake)
density(CO2$uptake) %>%
  lines(col = "red", lwd = 2)

ggplot(CO2, aes(x = uptake)) +
  geom_histogram(aes(y = ..density..), bins = 13, color = "black", fill = "grey", alpha = 0.5) +
  geom_density(col = "red", size = 1) + 
  geom_rug()

hist(CO2$uptake, prob=T, breaks = seq(min(CO2$uptake), max(CO2$uptake),length.out = 15))
rug(CO2$uptake)
density(CO2$uptake) %>%
  lines(col = "red", lwd = 2)
density(CO2$uptake, adjust = 0.5, cut = 0) %>%
  lines(col = "blue", lwd = 2)
density(CO2$uptake, adjust = 2, cut = 0) %>%
  lines(col = "green", lwd = 2)


hist(CO2$uptake, prob=T, breaks = seq(min(CO2$uptake), max(CO2$uptake),length.out = 15))
rug(CO2$uptake)
density(CO2$uptake) %>%
  lines(col = "black", lwd = 1, lty = 2)
density(CO2$uptake, bw = "SJ", cut = 0) %>% ## Default kernel is "gaussian") 
  lines(col = "red", lwd = 2)
density(CO2$uptake, kernel = "epanechnikov", bw = "SJ", cut = 0) %>%
  lines(col = "blue", lwd = 2)
density(CO2$uptake, kernel = "rectangular", bw = "SJ", cut = 0) %>%
  lines(col = "blue", lwd = 2)
