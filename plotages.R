#! /usr/bin/env Rscript
library(tidyverse)
library(cowplot)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)

dn <- 'k2p.dist'
file <- args[1]
d.est <- read.table(file, sep='\t', header=T, na.strings='NA')

comparedists <- function(dist.df, xpre, ypre, distname) {
    xname <- paste(xpre, distname, sep='.')
    yname <- paste(ypre, distname, sep='.')
    dist.df %>%
        select(locus, subfamily, category, xname, yname) %>%
        na.omit() %>%
        ggplot(aes_string(x=xname, y=yname)) +
        geom_point() +
        geom_abline() + geom_smooth(method='lm',formula=y~x) + coord_fixed()
}

agedist <- function(dist.df, distname, ypre='') {
    if(ypre != '') {
        yname <- paste(ypre, distname, sep='.')
    } else {
        yname <- distname
    }
    dist.df %>%
        select(locus, subfamily, yname) %>%
        filter(!grepl(',', subfamily)) %>%
        na.omit() %>%
        mutate(subfamily=factor(subfamily)) %>%
        ggplot(aes_string('subfamily', yname)) + 
            geom_violin() + 
            geom_jitter(width=0.1,height=0.0,alpha=0.4) +
            scale_y_reverse() + coord_flip()
}


p1 <- d.est %>% agedist(dn)
# p2 <- d.est %>% comparedists('con5', 'con3', dn)
# p3 <- d.est %>% comparedists('ltr', 'con5', dn)
# p4 <- d.est %>% comparedists('ltr', 'con3', dn)

pdf(args[2], paper='USr', width=11, height=8.5)
p1 # grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()

