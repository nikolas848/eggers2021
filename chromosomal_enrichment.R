rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################## enrichment ########################
######################################################

library(dplyr)
library(tidyverse)
library(matrixStats)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(RColorBrewer)


##### change this:
#peaks to read in
my_peaks <- list.files(pattern="*.bed$")
# chromosomes to consider
my_chr <- c("chrX","chr2L","chr2R","chr3L","chr3R","chr4")
### determine color for chromosomes
chr.col <- function(feature) {
  return(switch(feature,
                "chr2L" = "#1F78B4",
                "chrY" = "#f5e618",
                "chr3L" = "#E31A1C",
                "chr3R" = "#FF7F00",
                "chr4" = "#B15928",
                "chrX" = "#6A3D9A",
                "chr2R" = "#33A02C"))
}

##read in data into tidyformat
my_names <- c(gsub("_all","",gsub("merged.","",(gsub(".bed","",my_peaks)))))

for(i in seq_along(my_peaks)){
  my_name <- paste("chip",i,sep="")
  my_chip <- read_tsv(my_peaks[i], col_names = c("chr","start","stop"),comment = "#")
  my_chip$sam <- my_names[i]
  assign(my_name, my_chip)
}

chiplist <- lapply(ls(pattern="^chip"), get)
tidychip<- as_tibble(bind_rows(chiplist))

random <- as.data.frame(seqlengths(keepSeqlevels(TxDb.Dmelanogaster.UCSC.dm6.ensGene, my_chr, pruning.mode = "coarse"))/1e6)

random$sam <- c("1 Genome")
random$chr <- rownames(random)
random1 <- random[,c(2,3,1)]
colnames(random1) <- c("sam","chr","hits")

tidychip %>% 
  group_by(sam,chr) %>% 
  filter(chr %in% my_chr) %>% 
  summarise(hits=length(chr)) -> chiphits

plotdata <- as_tibble(rbind.data.frame(random1,chiphits))
plotdata$chr <- factor(plotdata$chr, levels=rev(my_chr))

cols <- sapply(my_chr, chr.col)

theme_bars1 <- theme_grey(base_size = 10, base_family = "") %+replace% 
  theme(
    axis.title = element_text(size = rel(1.5)), 
    axis.text = element_text(size = rel(1)), 
    axis.ticks.x = element_line(colour = "black", size= rel(0.7)), 
    axis.ticks.y = element_blank(), 
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    strip.background = element_blank(),
    strip.text =  element_text(size = rel(1)))

p <- ggplot(plotdata, aes(x="", y=hits, fill=chr)) + geom_bar(width=1, stat="identity") +
  facet_grid(~sam, scales = "free", space='free') +
  scale_x_discrete(expand = c(0, 0.5)) +
  coord_cartesian(expand = F)
p <- p + facet_wrap(~sam, scales ="free",strip.position = "left", nrow=length(my_peaks)+1, ncol=1)
p <- p + labs(title = "Peak distribution",x="",y="")
p1 <- p + scale_fill_manual(values=cols) + theme_bars1 
p1 <- p1 + coord_flip() 
p1 <- p1 + theme(strip.text.y.left = element_text(size=rel(1),angle=0))

pdf("1enrichment_X2.pdf",height=length(my_peaks)+1,width=5)
p1
dev.off()
##########

