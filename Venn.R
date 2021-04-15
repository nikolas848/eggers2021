rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(rtracklayer)
library(GenomicRanges)
library(RColorBrewer)
library(GenomicRanges)
library(rtracklayer)

#####################################################################################################################
#####################################################################################################################

my_chromosomes <- c("chr2L","chr2R","chr3L","chr3R","chr4","chrX")
my_colors <- brewer.pal(7,"Set1")
#####################################################################################################################
#####################################################################################################################

### load new peaks ###


my_peak_files <- c(list.files(pattern = "*.bed$"))

#my_peak_files <- my_peak_files[c(1,2,3)]

i=1
for(i in seq_along(my_peak_files)){
  
  my_name <- paste(gsub("../","",my_peak_files[i]),".ranges", sep="")
  
  my_bed <- import.bed(my_peak_files[i])
  
  #strand(my_bed) <- "*"
  
  seqlevels(my_bed, pruning.mode="coarse") <- my_chromosomes
  
  my_bed <- keepSeqlevels(my_bed, my_chromosomes, pruning.mode = "coarse")
  
  assign(my_name, my_bed)
}

my_peak_file_names <- ls(pattern = "\\.ranges$")


#####################################################################################################################
#####################################################################################################################
### pooled peaks ###

my_pooled_peaks <- GRanges()

for(i in seq_along(my_peak_file_names)){
  
  my_pooled_peaks <- append(my_pooled_peaks, get(my_peak_file_names[i]))
  
}

my_pooled_peaks <- GenomicRanges::reduce(my_pooled_peaks,ignore.strand=T)

my_pooled_peaks$peak_id <- paste("peak", seq_along(my_pooled_peaks), sep="_")


for(i in seq_along(my_peak_file_names)){
  
  my_overlaps <- !(is.na(findOverlaps(maxgap = 100,ignore.strand=T,my_pooled_peaks, get(my_peak_file_names[i]), 
                                      select="arbitrary")))
  
  mcols(my_pooled_peaks) <- cbind(mcols(my_pooled_peaks), my_overlaps)
  
  colnames(mcols(my_pooled_peaks))[i+1] <- gsub("EX*","",gsub(".factor.*.bed.ranges", "", my_peak_file_names[i]))
  
}

my_overlaps_df <- as.data.frame(mcols(my_pooled_peaks))
my_overlaps_df[,-1] <- data.matrix(my_overlaps_df[,-1])


my_sample_ids <- 2:ncol(my_overlaps_df)

my_overlaps_list <- sapply(my_sample_ids, function(x){my_overlaps_df$peak_id[my_overlaps_df[,x] == 1]})
names(my_overlaps_list) <- gsub(".2.peaks","",gsub(".bed.ranges","",ls(pattern = "\\.ranges$")))
#names(my_overlaps_list) <- c("Msl2","Msl2 + Clamp","Msl2 vivo")

#####################################################################################################################
#####################################################################################################################



pdf(paste("venn.test",names(my_overlaps_list)[1],names(my_overlaps_list)[2],names(my_overlaps_list)[3],"pdf",sep="."), width = 6, height = 6)

library(Vennerable)
library(grid)
library(gridExtra)

my_overlaps_venn <-  Venn(my_overlaps_list)
my_overlaps_venn_plot <- compute.Venn(my_overlaps_venn,doWeights=T)

SetLabels <- VennGetSetLabels(my_overlaps_venn_plot)
####red
#SetLabels[SetLabels$Label==names(my_overlaps_list)[1],"y"] <- 75
####blue
#SetLabels[SetLabels$Label==names(my_overlaps_list)[2],"x"] <- 20
####green
#SetLabels[SetLabels$Label==names(my_overlaps_list)[3],"x"] <- -30

my_overlaps_venn_plot <- VennSetSetLabels(my_overlaps_venn_plot,SetLabels)
SetFaceLabels <- VennGetFaceLabels(my_overlaps_venn_plot)

#SetFaceLabels$x[1] <- SetFaceLabels$x[1] -15
#SetFaceLabels$y[1] <- SetFaceLabels$y[1] +1 
# 
#SetFaceLabels$x[3] <- SetFaceLabels$x[3] -5
# 
#SetFaceLabels$x[4] <- SetFaceLabels$x[4] +10
#SetFaceLabels$y[4] <- SetFaceLabels$y[4] +10
# 
#SetFaceLabels$x[6] <- SetFaceLabels$x[6] +2
#SetFaceLabels$y[6] <- SetFaceLabels$y[6] +2
#
#SetFaceLabels$x[7] <- SetFaceLabels$x[7] +5

my_overlaps_venn_plot <- Vennerable:::VennSetFaceLabels(my_overlaps_venn_plot,SetFaceLabels)
gp <- VennThemes(my_overlaps_venn_plot)


#gp$Set$Set1$col <- my_colors[6]
#gp$Set$Set2$col <- my_colors[2]
#gp$Set$Set3$col <- my_colors[7]

#gp$SetText$Set1$col <- my_colors[6]
#gp$SetText$Set2$col <- my_colors[2]
#gp$SetText$Set3$col <- my_colors[7]

plot(my_overlaps_venn_plot, show = list(Faces = FALSE), gp = gp
)
#grid.text("Comparison", x = 0.5, y = 0.95, gp = gpar(cex=1.75))

dev.off()

