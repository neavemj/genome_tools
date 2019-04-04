#!/usr/bin/env Rscript

library(ggplot2)
require(scales)

#tiling <- read.table("wssv_250k.kmer_a5.tiling", sep="\t", header=T)

file_handle <- commandArgs(TRUE)
tiling_handle <- file_handle[1]
tiling <- read.table(tiling_handle, sep="\t", header=T)


p <- ggplot(tiling, aes(x=start, xend=end, y=1, yend=1, color=wrapped)) +
  geom_segment() + 
  geom_point() +
  geom_point(aes(x=end)) +
  #geom_point() +
  #geom_point(aes(x=end, color='blue')) +
  scale_x_continuous(labels = comma) +
  theme_minimal() +
  theme(axis.text.y=element_blank(), strip.text.y=element_text(angle=180)) +
  ylab("") +
  xlab("reference position (bp)") +
  facet_wrap(~sample, ncol=1, strip.position='left')

ggsave(paste(tiling_handle, ".pdf", sep=""), p)

print("plotting complete")