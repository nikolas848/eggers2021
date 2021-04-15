rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(LSD)
library(RColorBrewer)
library(ComplexHeatmap)
library(genefilter)
library(ggplot2)
library(zoo)
library(gridExtra)
library(circlize)
library(IRanges)
library(ShortRead)
library(rtracklayer)
########################

##############plot data######################################

rm(list=ls())
my_window = 2000

################################

myHeatmap <- function(mat, column_title, name, col1, clustering_method_rows = "complete",limit,i){
  Heatmap(mat,
          col = col1,
          column_title = column_title,
          name = name,
          use_raster = T,
          gap = 0,
          show_column_names = FALSE,
          show_row_names = FALSE,
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          clustering_method_rows = clustering_method_rows,
          column_title_gp = gpar(fontsize = 12),
          column_title_rot = 0,
          row_title = paste("n =",length(BoverA)),
          row_title_rot = 90,
          row_title_gp = gpar(fontsize = 12),
          heatmap_legend_param = list(legend_width = unit(2.5,"cm"),
                                      title = "scale",
                                      direction = "horizontal",
                                      title_gp = gpar(fontsize = 12),
                                      labels_gp = gpar(fontsize = 9)),
          show_heatmap_legend = TRUE,
          width = 300,
          na_col = "black",
          top_annotation=HeatmapAnnotation(sum = anno_lines(colMeans(mat, na.rm = FALSE, dims = 1),height = unit(50,"points"),
                                                            ylim=c(min(limit),max(limit)),
                                                            axis=if(i==1){T}else{F},
                                                            border = F,gp = gpar(lwd=3,col = col1(max(mat)))),
                                                            show_annotation_name = F),
          bottom_annotation=columnAnnotation(foo = anno_text(x=c(rep("",ncol(mat)/200*40),"-1 kb",rep("",ncol(mat)/200*67),"0",rep("",ncol(mat)/200*85),"+ 1 kb",rep("",ncol(mat)/200*5)),
                                                             rot=0,
                                                             gp =gpar(fontsize = 9)))
  )
}

bin.matrix <- function(m, bin.size) {
  bm <- c()
  for (i in 1:round((ncol(m)/bin.size))) {
    er <- (i*bin.size)
    sr <- er-(bin.size-1)
    bm <- cbind(bm,rowMeans(m[,sr:er]))
  }
  bm
}
color_functionG <- function(mat){
  x = as.vector(as.matrix(mat))
  colfun <- colorRamp2(breaks = seq(1, quantile(x, 0.9), length = 5),
                       colors = brewer.pal(5,"BuPu"))
  return(colfun)
}
color_functionM <- function(mat){
  x = as.vector(as.matrix(mat))
  colfun <- colorRamp2(breaks = seq(1, quantile(x, 0.9), length = 4),
                       colors = brewer.pal(4,"YlGnBu"))
  return(colfun)
}
color_functionC <- function(mat){
  x = as.vector(as.matrix(mat))
  colfun <- colorRamp2(breaks = seq(1, quantile(x, 0.9), length = 4),
                       colors = brewer.pal(4,"YlGn"))
  return(colfun)
}
color_functionN <- function(mat){
  x = as.vector(as.matrix(mat))
  colfun <- colorRamp2(breaks = seq(1, quantile(x, 0.9), length = 8),
                       colors = brewer.pal(8,"YlOrBr"))
  return(colfun)
}

color_functionMN <- function(mat){
  x = as.vector(as.matrix(mat))
  colfun <- colorRamp2(breaks = seq(5, 40, length = 3),
                       colors = c("blue","white","red"))
  return(colfun)
}

################################heatmap
##
my_search<- gsub(".bed","",gsub("centers.","",list.files(path = "../../centers2",pattern=".bed")))
#my_search <- my_search[c(2,3,5,6,7,9,11,12,14)]
#my_search <- c("clamp_dip_sites","clamp_vivo_sites","Phaser_sites","has_sites","pionx_sites","GAF_vitro_GCM","SuHW_sites")
my_search <- c("clamp_dip_sites")

for(f in seq_along(my_search)){
  #my_meansN <- list.files(path="../../chip_nurf/R_both",pattern= "^area",full.names = T)
  #my_meansN <- my_meansN[grepl(paste0("all.",my_search[f]), my_meansN) == TRUE][c(2,1)]
  
   #my_meansM <- list.files(path="../../chip_msl2_25/1_R_new_normalized_to neg",pattern= "^area",full.names = T)
   #my_meansM <- my_meansM[grepl(paste0("all.",my_search[f]), my_meansM) == TRUE]
   my_meansM2 <- list.files(path="../../chip_msl2_gaf",pattern= "^area",full.names = T)
   my_meansM2 <- my_meansM2[grepl(paste0("all.",my_search[f]), my_meansM2) == TRUE]
   my_meansM <- c(my_meansM2)
  
  #my_meansC <- list.files(path="../../chip_clamp",pattern= "^area",full.names = T)
  #my_meansC <- my_meansC[grepl(paste0("all.",my_search[f]), my_meansC) == TRUE][]
  
   my_meansG <- list.files(path="../../chip_gaf",pattern= "^area",full.names = T)
   my_meansG <- my_meansG[grepl(paste0("all.",my_search[f]), my_meansG) == TRUE]
   #my_meansG <- my_meansG[c(4)]
  
  #my_mnase <- list.files(path="../../mnase_clamp",pattern= "^area",full.names = T)
  #my_mnase <- my_mnase[grepl(paste0("sum.*.",my_search[f]), my_mnase) == TRUE][]

  ###### sort heatmaps by
  my_sorter <- my_meansM[1]
  ######
  my_titlesN <- c("Nurf cntrl","Nurf +C +M")
  my_titlesM <- c("Msl2 +CM","Msl2 +GCM1","Msl2 +GCM2","Msl2 +GCM3","Msl2 +GM")
  my_titlesG <- c("Gaf +G","Gaf +GCM1","Gaf +GCM2","Gaf +GCM3")
  my_titlesC <- c("Clamp alone","Clamp +M")
  my_titlesMN <- c("MNase cntrl","MNase +C")

  rowMeans(get(load(my_sorter))[,950:1050])-> NFR # 100 bp binding
  m <- NFR[is.finite(NFR)]
  
  m.sort <- sort(m, decreasing = T, na.last = T)
  BoverA <- names(m.sort)


if(exists("my_meansN")){
  mylimit=NULL
  for(i in seq_along(my_meansN)){
  maxsort <- c(min(colMeans(get(load(my_meansN[i])))),max(colMeans(get(load(my_meansN[i])))))
  mylimit <- c(mylimit,maxsort)
  }
  for(i in seq_along(my_meansN)){
    mysort <- get(load(my_meansN[2]))
    my_name <- paste("htN",i,sep="")
    mat.rng <- get(load(my_meansN[i]))
    match(BoverA, rownames(mat.rng), nomatch = NA) -> sorted.number # gene sorting
    mat.rng <- na.omit(mat.rng[sorted.number,])
    mat.rng <- bin.matrix(mat.rng,10)
    my_ht <- myHeatmap(mat=mat.rng, col =  color_functionN(mysort), column_title = my_titlesN[i],limit=mylimit,i=i)
    my_ht@matrix[,101] <- NA
    my_ht@matrix[,100] <- NA
    assign(my_name, my_ht)

  }}else{next} 
  
if(exists("my_meansC")){  
  mylimit=NULL
  for(i in seq_along(my_meansC)){
    maxsort <- c(min(colMeans(get(load(my_meansC[i])))),max(colMeans(get(load(my_meansC[i])))))
    mylimit <- c(mylimit,maxsort)
  }
for(i in seq_along(my_meansC)){
  mysort <- get(load(my_meansC[1]))
  my_name <- paste("htC",i,sep="")
  mat.rng <- get(load(my_meansC[i]))
  match(BoverA, rownames(mat.rng), nomatch = NA) -> sorted.number # gene sorting
  mat.rng <- na.omit(mat.rng[sorted.number,])
  mat.rng <- bin.matrix(mat.rng,10)
  my_ht <- myHeatmap(mat.rng, col =  color_functionC(mysort), column_title = my_titlesC[i],limit=mylimit,i=i)
  my_ht@matrix[,101] <- NA
  my_ht@matrix[,100] <- NA
  assign(my_name, my_ht)

}}else{next} 
  
  if(exists("my_meansG")){  
  mylimit=NULL
  for(i in seq_along(my_meansG)){
    maxsort <- c(min(colMeans(get(load(my_meansG[i])))),max(colMeans(get(load(my_meansG[i])))))
    mylimit <- c(mylimit,maxsort)
  }
  for(i in seq_along(my_meansG)){
    mysort <- get(load(my_meansG[1]))
    my_name <- paste("htG",i,sep="")
    mat.rng <- get(load(my_meansG[i]))
    match(BoverA, rownames(mat.rng), nomatch = NA) -> sorted.number # gene sorting
    mat.rng <- na.omit(mat.rng[sorted.number,])
    mat.rng <- bin.matrix(mat.rng,10)
    my_ht <- myHeatmap(mat.rng, col =  color_functionG(mysort), column_title = my_titlesG[i],limit=mylimit,i=i)
    my_ht@matrix[,101] <- NA
    my_ht@matrix[,100] <- NA
    assign(my_name, my_ht)

  }}else{next} 
  
  if(exists("my_meansM")){  
  mylimit=NULL
  for(i in seq_along(my_meansM)){
    maxsort <- c(min(colMeans(get(load(my_meansM[i])))),max(colMeans(get(load(my_meansM[i])))))
    mylimit <- c(mylimit,maxsort)
  }
  for(i in seq_along(my_meansM)){
    mysort <- get(load(my_meansM[1]))
    my_name <- paste("htM",i,sep="")
    mat.rng <- get(load(my_meansM[i]))
    match(BoverA, rownames(mat.rng), nomatch = NA) -> sorted.number # gene sorting
    mat.rng <- na.omit(mat.rng[sorted.number,])
    mat.rng <- bin.matrix(mat.rng,10)
    my_ht <- myHeatmap(mat.rng, col =  color_functionM(mysort), column_title = my_titlesM[i],limit=mylimit,i=i)
    my_ht@matrix[,101] <- NA
    my_ht@matrix[,100] <- NA
    assign(my_name, my_ht)

  }}else{next} 
  
  if(exists("my_mnase")){  
  mylimit=NULL 
  for(i in seq_along(my_mnase)){
    maxsort <- c(min(colMeans(get(load(my_mnase[i])))),max(colMeans(get(load(my_mnase[i])))))
    mylimit <- c(mylimit,maxsort)
  }
  for(i in seq_along(my_mnase)){
    mysort <- get(load(my_mnase[2]))
    my_name <- paste("htMN",i,sep="")
    mat.rng <- get(load(my_mnase[i]))
    match(BoverA, rownames(mat.rng), nomatch = NA) -> sorted.number # gene sorting
    mat.rng <- na.omit(mat.rng[sorted.number,])
    mat.rng <- bin.matrix(mat.rng,10)
    my_ht <- myHeatmap(mat.rng, col = color_functionMN(mysort), column_title = my_titlesMN[i],limit=mylimit,i=i)
    my_ht@matrix[,101] <- NA
    my_ht@matrix[,100] <- NA
    assign(my_name, my_ht)

  }}else{next} 
rm(my_ht)  
  #my_list <- (my_list[c(3,2,1,4,5)])
  df <- as.integer(grepl("chrX",BoverA))
  df[df=="1"]<-"Chr X"
  df[df=="0"]<-"Autosome"
  col_letters = c("Chr X" = "red","Autosome"="White")
  hanno1 = Heatmap(df,
                   cluster_rows = F,
                   col=col_letters,
                   column_title = "Chr X",
                   column_title_rot = 90,
                   width=50,
                   gap = 0,
                   show_column_names = FALSE,
                   show_row_names = FALSE,
                   cluster_columns = FALSE,
                   show_heatmap_legend = F,
                   heatmap_legend_param = list(title_gp = gpar(fontsize = 12,font=2),
                                              labels_gp = gpar(fontsize = 12)))
  
  df2 <- data.frame(chr=sapply(strsplit(BoverA," "),"[", 1),start=as.integer(sapply(strsplit(BoverA," "),"[", 2)))
  df2$stop <- df2$start+200
  df2$start <- df2$start-200
  df2 <- GRanges(df2)
  HAS <- import.bed("../../centers2/has_sites.bed")
  df2 <- overlapsAny(df2,HAS)
  df2[df2==TRUE]<-"HAS"
  df2[df2==FALSE]<-"other"
  col_letters = c("HAS" = "Blue","other"="White")
  hanno2 = Heatmap(df2,
                    name="",
                    cluster_rows = F,
                    col=col_letters,
                    column_title = "HAS",
                    column_title_rot = 90,
                    width=50,
                    gap=0,
                    show_column_names = FALSE,
                    show_row_names = FALSE,
                    cluster_columns = FALSE,
                    show_heatmap_legend = F,
                    heatmap_legend_param = list(title_gp = gpar(fontsize = 12,font=2),
                                                labels_gp = gpar(fontsize = 12)))
  
  
  
###resorting heatmaps
  my_list <- ls(pattern="ht")
  #my_list <- my_list[c(1,2,5,6,3,4)]
  #my_list <- c("htC1","htC2","htMN1","htMN2","htN2","htN1","htM1","htM2","htM3","htG1")
  
my_list1=NULL
  for(q in seq_along(my_list)){
    x <- get(my_list[q])
    my_list1 <- my_list1+x
  }

my_list2 <- my_list1+hanno1+hanno2
  
   pdf(paste0("heatmap.raster.300wid",my_search[f],".pdf"),compress = F,height=4,width=length(my_list2)+1,)
   draw(my_list2,
        #padding = unit(c(1, 1, 1, 1), "cm"),
        column_title = paste(gsub("sites sites","sites",gsub("_"," ",my_search[f]),"sites")),
        column_title_gp = gpar(fontsize = 12))

   dev.off()
   
}

  