rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(DNAshapeR)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(RColorBrewer)
library(stringr)

my_files <- list.files(".",pattern="*.tsv$")

i=1



for(i in seq_along(my_files)){

myname <- gsub(".tsv","",my_files[i])
myseq <-  read.delim(my_files[i])
myseq <-head(myseq, -3)

myseq %>% 
  distinct(sequence_name , .keep_all = T) -> my_seq


myseq2 <- data.frame(chr=sapply((str_split(my_seq$sequence_name,pattern = ":")), getElement, 1))
ranges <- sapply((str_split(my_seq$sequence_name,pattern = ":")), getElement, 2)   
myseq2$start <- as.double(sapply((str_split(ranges,pattern = "-")), getElement, 1))
myseq2$stop <- as.double(sapply((str_split(ranges,pattern = "-")), getElement, 2))

myseq2$stop <- c(myseq2$start+my_seq$stop)
myseq2$start <- c(myseq2$start+my_seq$start)
myseq2$strand <- my_seq$strand
myseq2
seqranges <- GRanges(myseq2)
seqranges <- resize(seqranges,width=27,fix="center")


getFasta(seqranges,BSgenome.Dmelanogaster.UCSC.dm6,width = width(seqranges),filename = paste0(myname,".fasta"))
         
pred <- getShape(paste0(myname,".fasta"))
assign(paste0(myname,"_pred"),pred)
rm(pred)
}

  my_shapes <- c("MGW","HelT","ProT","Roll","EP")
 
  
  my_pred <- ls(pattern="_pred")
  #my_pred <- my_pred[c(1,2,3)]
  
  i=1
  for(i in seq_along(my_shapes)){
    fullplotdf=NULL
    fullmediansdf =NULL
    if(i %in% c(1,2,5)){
      next
    }else{
    for(f in seq_along(my_pred)){
      myname <- gsub("_pred","",my_pred[f])
      plotdf <- as.data.frame(get(my_pred[f])[[i]])
      plotdf <- gather(plotdf,key = "position",value = "score")
      plotdf$position <- gsub("V","",plotdf$position)
      plotdf$position <- factor(plotdf$position,levels = c(1:max(as.double(plotdf$position))))
      plotdf$name <- myname
      plotdf %>% 
        group_by(position,name) %>% 
        summarise(median = mean(score)) ->   plotmedians
      fullplotdf <- bind_rows(fullplotdf,plotdf)
      fullmediansdf <- bind_rows(fullmediansdf,plotmedians)
    } 
    fullplotdf$name <- factor(fullplotdf$name, levels = c("GAF","MSL2","MSL2+C","MSL2+GC","HAS"))
    fullmediansdf$name <- factor(fullmediansdf$name, levels = c("GAF","MSL2","MSL2+C","MSL2+GC","HAS"))
    
    title <- names(my_shapes[i])
    my_colors1 <- c(brewer.pal(5,"BuPu")[5],brewer.pal(9,"GnBu")[c(7,8,9)],"darkgreen")
    #barplot(c(1:19),col=my_colors1)
    if(i %in% c(3)){
    xlimits <- c(min(as.double(na.omit(fullmediansdf)$position))-1,max(as.double(na.omit(fullmediansdf)$position)))
    }else{
    xlimits <- c(min(as.double(na.omit(fullmediansdf)$position)),max(as.double(na.omit(fullmediansdf)$position)))
    }
    ylimits <- c(quantile(na.omit(fullplotdf)$score,0.02),quantile(na.omit(fullplotdf)$score,0.98))
    
    p <- ggplot(fullplotdf,aes(x=position,y=score,fill=name,col=name))+
      #geom_boxplot(na.rm = T,width=0.9,lwd=0.3,outlier.size = 0.3,outlier.alpha = 0.1,alpha=0.6)+
      labs(title = title)+
      theme_classic()+
      scale_color_manual(values=my_colors1,aesthetics = c("color","fill"))+
      geom_path(na.rm=T,data=fullmediansdf,mapping=aes(x=position,y=median,col=name,group=1),lineend = "round",size=2)+
      theme(legend.position = "none",text = element_text(size=20))+
      coord_cartesian(xlim=xlimits,ylim=ylimits)
    p
    #stat_compare_means(method="wilcox.test",aes(label = ..p.signif..),hide.ns = T)

    ggsave(paste0("1_mean",my_shapes[i],".pdf"), height = 6, width = 6)
    
  } 
  } 
 
  
  
  
  
  my_pred <- ls(pattern="_pred")
  my_pred <- my_pred[c(1,2,3)]

  for(i in seq_along(my_shapes)){
    fullplotdf=NULL
    fullmediansdf =NULL
    if(i %in% c(1,2,5)){
      next
    }else{
      for(f in seq_along(my_pred)){
        myname <- gsub("_pred","",my_pred[f])
        plotdf <- as.data.frame(get(my_pred[f])[[i]])
        plotdf <- gather(plotdf,key = "position",value = "score")
        plotdf$position <- gsub("V","",plotdf$position)
        plotdf$position <- factor(plotdf$position,levels = c(1:max(as.double(plotdf$position))))
        plotdf$name <- myname
        plotdf %>% 
          group_by(position,name) %>% 
          summarise(median = mean(score)) ->   plotmedians
        fullplotdf <- bind_rows(fullplotdf,plotdf)
        fullmediansdf <- bind_rows(fullmediansdf,plotmedians)
      } 
      fullplotdf$name <- factor(fullplotdf$name, levels = c("GAF","MSL2","MSL2+C","MSL2+GC","HAS"))
      fullmediansdf$name <- factor(fullmediansdf$name, levels = c("GAF","MSL2","MSL2+C","MSL2+GC","HAS"))
      
      title <- names(my_shapes[i])
      my_colors1 <- c(brewer.pal(5,"BuPu")[5],brewer.pal(9,"GnBu")[c(7)],"Red4")
      if(i %in% c(3)){
        xlimits <- c(min(as.double(na.omit(fullmediansdf)$position))-1,max(as.double(na.omit(fullmediansdf)$position)))
      }else{
        xlimits <- c(min(as.double(na.omit(fullmediansdf)$position)),max(as.double(na.omit(fullmediansdf)$position)))
      }
      ylimits <- c(quantile(na.omit(fullplotdf)$score,0.02),quantile(na.omit(fullplotdf)$score,0.98))     
      p <- ggplot(fullplotdf,aes(x=position,y=score,fill=name,col=name))+
        #geom_boxplot(na.rm = T,width=0.9,lwd=0.3,outlier.size = 0.3,outlier.alpha = 0.1,alpha=0.6)+
        labs(title = title)+
        theme_classic()+
        scale_color_manual(values=my_colors1,aesthetics = c("color","fill"))+
        geom_path(na.rm=T,data=fullmediansdf,mapping=aes(x=position,y=median,col=name,group=1),lineend = "round",size=2)+
        theme(legend.position = "none",text = element_text(size=20))+
        coord_cartesian(xlim=xlimits,ylim=ylimits)
      ggsave(paste0("2_mean",my_shapes[i],".pdf"), height = 6, width = 6)
      
    } 
  } 
  
  
  my_pred <- ls(pattern="_pred")
  my_pred <- my_pred[c(1,2,4)]
  
  for(i in seq_along(my_shapes)){
    fullplotdf=NULL
    fullmediansdf =NULL
    if(i %in% c(1,2,5)){
      next
    }else{
      for(f in seq_along(my_pred)){
        myname <- gsub("_pred","",my_pred[f])
        plotdf <- as.data.frame(get(my_pred[f])[[i]])
        plotdf <- gather(plotdf,key = "position",value = "score")
        plotdf$position <- gsub("V","",plotdf$position)
        plotdf$position <- factor(plotdf$position,levels = c(1:max(as.double(plotdf$position))))
        plotdf$name <- myname
        plotdf %>% 
          group_by(position,name) %>% 
          summarise(median = mean(score)) ->   plotmedians
        fullplotdf <- bind_rows(fullplotdf,plotdf)
        fullmediansdf <- bind_rows(fullmediansdf,plotmedians)
      } 
      fullplotdf$name <- factor(fullplotdf$name, levels = c("GAF","MSL2","MSL2+C","MSL2+GC","HAS"))
      fullmediansdf$name <- factor(fullmediansdf$name, levels = c("GAF","MSL2","MSL2+C","MSL2+GC","HAS"))
      
      title <- names(my_shapes[i])
      my_colors1 <- c(brewer.pal(5,"BuPu")[5],brewer.pal(9,"GnBu")[c(8)],"Red4")
      if(i %in% c(3)){
        xlimits <- c(min(as.double(na.omit(fullmediansdf)$position))-1,max(as.double(na.omit(fullmediansdf)$position)))
      }else{
        xlimits <- c(min(as.double(na.omit(fullmediansdf)$position)),max(as.double(na.omit(fullmediansdf)$position)))
      }
      ylimits <- c(quantile(na.omit(fullplotdf)$score,0.02),quantile(na.omit(fullplotdf)$score,0.98))
      p <- ggplot(fullplotdf,aes(x=position,y=score,fill=name,col=name))+
        #geom_boxplot(na.rm = T,width=0.9,lwd=0.3,outlier.size = 0.3,outlier.alpha = 0.1,alpha=0.6)+
        labs(title = title)+
        theme_classic()+
        scale_color_manual(values=my_colors1,aesthetics = c("color","fill"))+
        geom_path(na.rm=T,data=fullmediansdf,mapping=aes(x=position,y=median,col=name,group=1),lineend = "round",size=2)+
        theme(legend.position = "none",text = element_text(size=20))+
        coord_cartesian(xlim=xlimits,ylim=ylimits)
      ggsave(paste0("3_mean",my_shapes[i],".pdf"), height = 6, width = 6)
      
    } 
  } 
  
  
  
  my_pred <- ls(pattern="_pred")
  my_pred <- my_pred[c(1,2,5)]
  
  for(i in seq_along(my_shapes)){
    fullplotdf=NULL
    fullmediansdf =NULL
    if(i %in% c(1,2,5)){
      next
    }else{
      for(f in seq_along(my_pred)){
        myname <- gsub("_pred","",my_pred[f])
        plotdf <- as.data.frame(get(my_pred[f])[[i]])
        plotdf <- gather(plotdf,key = "position",value = "score")
        plotdf$position <- gsub("V","",plotdf$position)
        plotdf$position <- factor(plotdf$position,levels = c(1:max(as.double(plotdf$position))))
        plotdf$name <- myname
        plotdf %>% 
          group_by(position,name) %>% 
          summarise(median = mean(score)) ->   plotmedians
        fullplotdf <- bind_rows(fullplotdf,plotdf)
        fullmediansdf <- bind_rows(fullmediansdf,plotmedians)
      } 
      fullplotdf$name <- factor(fullplotdf$name, levels = c("GAF","MSL2","MSL2+C","MSL2+GC","HAS"))
      fullmediansdf$name <- factor(fullmediansdf$name, levels = c("GAF","MSL2","MSL2+C","MSL2+GC","HAS"))
      
      title <- names(my_shapes[i])
      my_colors1 <- c(brewer.pal(5,"BuPu")[5],brewer.pal(9,"GnBu")[c(9)],"Red4")
      if(i %in% c(3)){
        xlimits <- c(min(as.double(na.omit(fullmediansdf)$position))-1,max(as.double(na.omit(fullmediansdf)$position)))
      }else{
        xlimits <- c(min(as.double(na.omit(fullmediansdf)$position)),max(as.double(na.omit(fullmediansdf)$position)))
      }
      ylimits <- c(quantile(na.omit(fullplotdf)$score,0.02),quantile(na.omit(fullplotdf)$score,0.98))
      p <- ggplot(fullplotdf,aes(x=position,y=score,fill=name,col=name))+
        #geom_boxplot(na.rm = T,width=0.9,lwd=0.3,outlier.size = 0.3,outlier.alpha = 0.1,alpha=0.6)+
        labs(title = title)+
        theme_classic()+
        scale_color_manual(values=my_colors1,aesthetics = c("color","fill"))+
        geom_path(na.rm=T,data=fullmediansdf,mapping=aes(x=position,y=median,col=name,group=1),lineend = "round",size=2)+
        theme(legend.position = "none",text = element_text(size=20))+
        coord_cartesian(xlim=xlimits,ylim=ylimits)
      ggsave(paste0("4_mean",my_shapes[i],".pdf"), height = 6, width = 6)
      
    } 
  } 
  