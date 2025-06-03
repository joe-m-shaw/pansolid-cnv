# PanSolid CNV deletion service quality

library(tidyverse)

stdev_live <- read_csv(file = paste0(
  config::get("data_folderpath"),
  "live_service/collated/",
  "stdev_live.csv"))

percent_138x_live <- read_csv(file = paste0(
  config::get("data_folderpath"),
  "live_service/collated/",
  "percent_138x_live.csv"))

ggplot(stdev_live, aes(x = worksheet, y = stdev_noise)) +
  geom_boxplot(outliers = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "") +
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.5)) +
  geom_hline(yintercept = 0.7, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed")
  
ggplot(percent_138x_live, aes(x = worksheet, y = percent_138x)) +
  geom_boxplot(outliers = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "")
