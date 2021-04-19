rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(ggplot2)
library(ggpmisc)
library(tidyverse)
library(GGally)
library(IRanges)

#####################################################################################################################
coverageWindowsStranded <- function(centers, window.size = 2000, coverage) {
  
  #  cl <- makeCluster(getOption("cl.cores", 8))
  #  clusterExport(cl, list("centers","coverage","window.size") , envir=environment())
  
  centers <- centers[centers$chr %in% names(coverage),]
  
  #  result <- parSapply(cl, names(cov), function(x) {
  result <- sapply(names(coverage), function(x) {
    my.cov <- coverage[[x]]
    my.peaks <- centers[centers$chr==x,]
    mw.views <- Views(my.cov, start=my.peaks$center-ceiling(window.size/2), width=window.size+1)
    ## remove out-of bounds views
    flt <- start(mw.views)>0 & end(mw.views) < length(my.cov)
    mw.views <- mw.views[flt,]
    my.peaks <- my.peaks[flt,]
    if (length(mw.views) > 0) {
      mat <- as.matrix(mw.views)
      colnames(mat) <- seq(from=(0-ceiling(window.size/2)), to=0+ceiling(window.size/2))
      rownames(mat) <- rownames(my.peaks)
      return(mat)
    } else {
      return(NULL)
    }
  })
  #  stopCluster(cl)
  mat <- Reduce(rbind, result)
  centers <- centers[rownames(centers) %in% rownames(mat),]
  na.omit(match(rownames(centers), rownames(mat)) )-> o
  centers <- centers[o,]
  #mat[centers$strand=="-",] <- t(apply(mat[centers$strand=="-",],1,rev))
  mat
}
#####################################################################################################################

my_files <- list.files(path = ".",pattern = "normalized.*all.rda",full.names = T)
my_files <- my_files[grepl("chip_cm",my_files)==T]
my_centers <- ("centers.msl2_chip_cm_sites.rda") ### my_centers is 0


####### calculate m_area ########
i=1
f=1
my_areas <- list.files(pattern = "^corr")  ### my_means is still empty

parallel::mclapply(seq_along(my_files), mc.cores = 8, FUN = function(i){
  for(f in seq_along(my_centers)){
    my_name <- paste("corr", my_files[i],sep=".")
    my_name <- gsub("../normalized","", my_name)
    my_name <- gsub(".all.rda","", my_name)
    my_name <- gsub("centers.","",my_name)
    
    norm.cover <- get(load(my_files[i]))
    center<- get(load(my_centers[f]))
    my_window = 2000
    m_area <- coverageWindowsStranded(center,my_window,norm.cover)
    
    m_area_mean <- rowMeans(m_area[,950:1050])
    
    assign(my_name, m_area_mean)
    
    save(list = my_name, file = my_name)
    
    print(paste(my_name,"created"))
    
  }
})

rm(list=ls())
########################################
my_set2 <- list.files(pattern="area*")
samples <- c("alone","neg","cm")
i=2
for(i in seq_along(samples)){
  
my_set <- my_set2[grepl(samples[i],my_set2)==T]

for(t in seq_along(my_set)){
  load(my_set[t])
}
my_titles <- gsub("area.","",my_set)
my_df <-data.frame(mget(my_set))
colnames(my_df) <- my_titles
my_df
myrho1 <- cor.test(my_df[,1],my_df[,2], method="spearman")
myrho2 <- cor.test(my_df[,2],my_df[,3], method="spearman")
myrho3 <- cor.test(my_df[,1],my_df[,3], method="spearman")
title <- gsub("rho = ","",paste0("103_85",myrho1[4],"85_96",myrho2[4],"103_96",myrho3[4]))

myplot <- ggpairs(my_df,axisLabels = "show")
myplot <- myplot+labs(title=title)
ggsave(paste0(samples[i],".correlation.pdf"),plot = myplot)

}

dev.off()
