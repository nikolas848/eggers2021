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
my_set <- list.files(pattern="^corr*")
my_set <- my_set[grepl("sum",my_set)==FALSE]
for(t in seq_along(my_set)){
  load(my_set[t])
}
  my_titles <- gsub("area.","",my_set)
  my_df <-data.frame(mget(my_set))
  colnames(my_df) <- my_titles
  my_df
  #my_df <- as_tibble(my_df)
  #my_df <- log10(my_df)
  

  myplot <- ggpairs(log2(my_df),axisLabels = "internal",)
  myplot
  
  ggsave(paste0("a.correlation.png"),myplot, width=5,height=5)


library(corrplot)
library(RColorBrewer)

correlation_matrix <- cor(my_df)
x <- corrplot(correlation_matrix, method = "square", type = "upper", tl.col = "black", order = "hclust", col = brewer.pal(n = 5, name = "RdYlBu"))
x






library(bedr)


samplenames <- c("_alone.bed","_neg.bed","_cm.bed")
i=1
for(i in seq_along(samplenames)){
  
files <- list.files(pattern=samplenames[i])
file_list=file.info(files)
file_list=file_list[order(file_list$size),]
file_list=head(file_list,1)
x <- rownames(file_list)
samplesize <- length(import.bed(x))
f=1
for(f in seq_along(files)){
  name <- files[f]
beds <- import.bed(files[f])
export.bed(sample(beds,size = samplesize),gsub(".bed",".sample2.bed",name))
}
}



head(
  x <- import.bed("85_msl2_chip_cm.sample2.bed")
length(x)
