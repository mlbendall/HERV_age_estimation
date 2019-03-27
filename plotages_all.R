#! /usr/bin/env Rscript
library(tidyverse)
library(cowplot)
library(gridExtra)

dn <- 'k2p.dist'
distfiles <- Sys.glob("*/distances.tsv")

d.est <- lapply(distfiles, function(f){
    read.table(f, sep='\t', header=T, na.strings='NA', stringsAsFactors=F)
    }) %>% 
    bind_rows() %>%
    select(locus, subfamily, dn) %>% 
    na.omit()

byfam <- d.est %>%
    group_by(subfamily) %>%
    summarize(count=n(), age=median(.data[[dn]])) %>%
    filter(count>2) %>%
    arrange(age)

z <- d.est %>%
    filter(subfamily %in% byfam$subfamily) %>%
    mutate(subfamily=factor(subfamily, levels=rev(byfam$subfamily)))

pdf('all_age_distribution.pdf', width=4, height=40)
z %>%
    ggplot(aes_string('subfamily', dn)) + 
            geom_violin() + 
            # geom_jitter(width=0.1,height=0.0,alpha=0.4) +
            scale_y_reverse() + coord_flip()
dev.off()

