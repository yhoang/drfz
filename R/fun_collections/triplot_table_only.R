#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019

### uncomment if not created yet
# fcs = new.env()



######## CREATE TRIPLOT TABLE this$tab3D w/o PLOTTING ##############
##### @param
# data              data matrix to analyse
# stainTable        table from DB, information of column names and cutoffs
# feat.X            name of feature X
# feat.Y            name of feature Y
# feat.Z1           name of feature Z1
# calc              calculation method
# col1              monochrome color palette for frequency of Z1
# col2              monochrome color palette for frequency of Z2
# binsize           bin size
# mincells          minimum amount of cells in a bin
# plot.range        plot range, first x axis, second y axis, default= c(2,12,2,12)
fcs$bintriploT_table <- function(
  data,
  stainTable,
  feat.X, 
  feat.Y, 
  feat.Z1, 
  calc, 
  col1, 
  col2, 
  binsize, 
  mincells, 
  plot.range=c(1,12,1,12)) 
{ 
  this=fcs 
  
  if (FALSE) {
    feat.X = "PD-1"
    feat.Y = "IFNg"
    feat.Z1 = "IL21"
    calc = "MFI+"
    binsize = 0.2
    mincells = 10
  }
  
  ### remove all signs and write anything with capitals
  feat.X.clean = gsub("[^[:alnum:]]","",toupper(feat.X))
  feat.Y.clean = gsub("[^[:alnum:]]","",toupper(feat.Y))
  feat.Z1.clean = gsub("[^[:alnum:]]","",toupper(feat.Z1))
  
  ### axes range
  xmin.val = plot.range[1]
  xmax.val = plot.range[2]
  ymin.val = plot.range[3]
  ymax.val = plot.range[4]
  
  ############ LOAD DATA AND CUTOFFS
  cutoffs = as.numeric(stainTable[which(stainTable$file_ID==file.idx),5])
  
  features = colnames(data) 
  
  ### remove all signs and write anything with capitals
  features.clean = gsub("[^[:alnum:]]","",make.unique(unlist(lapply(features,function(x) {
    len = length(strsplit(x,"[.]")[[1]])
    y = toupper(strsplit(x,"[.]")[[1]][1])
    paste(y, collapse=".")
  }))))
  
  idx.X = which(features.clean == feat.X.clean)
  idx.Y = which(features.clean == feat.Y.clean)
  idx.Z1 = which(features.clean == feat.Z1.clean)
  
  if ( is.na(idx.X) | is.na(idx.Y) | is.na(idx.Z1) ) stop(sprintf("Either feat.X (%s) or feat.Y (%s) or feat.Z1 (%s) not in file index %s",feat.X, feat.Y, feat.Z1, i))
  
  feat.Z1 = paste(sprintf("%s(%s)",feat.Z1,cutoffs[idx.Z1]))
  
  
  ### cut data if calc method is MFI+
  if ( calc == "MFI+" ) data = data[which(data[,idx.Z1]>cutoffs[idx.Z1]),]
  
  ### construct bin table with number of cells per bin
  fX = cut(data[,idx.X],breaks=seq(xmin.val,xmax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
  fY = cut(data[,idx.Y],breaks=seq(ymin.val,ymax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
  fXY=as.factor(paste(fX,fY))
  
  tab = table(fX,fY)
  colnames(tab)=seq(ymin.val,ymax.val-binsize,by=binsize)
  rownames(tab)=seq(xmin.val,xmax.val-binsize,by=binsize)
  
  this$tab = tab
  
  printf("length(bins with tab>=mincells)=%s",length(which(tab>=mincells)))
  
  if ( grepl(calc,"MFI") ) {
    ########## CALCULATE MFI
    my.lengths = aggregate(data[,idx.Z1],by=list(fXY),length)
    idx.len = which(my.lengths$x >= mincells)
    
    ### start with MFI calculation of Z1
    my.calc = aggregate(data[,idx.Z1],by=list(fXY),mean)
    
    min.range.Z1 = floor(min(my.calc[idx.len,'x'])*10)/10
    max.range.Z1 = ceiling(max(my.calc[idx.len,'x'])*10)/10
    
    # get steps for Z1
    step=round(diff(range(max.range.Z1,min.range.Z1))/10,2) 
    steps.Z1=seq(min.range.Z1,max.range.Z1,by=step)
    
    # bin color factor Z1
    my.calc.fac.Z1=cut(my.calc$x,breaks=steps.Z1,labels=1:10,include.lowest=TRUE)
    names(my.calc.fac.Z1) = my.calc$x
    
    ### combine all frequencies in one table
    my.calc = cbind(my.calc,fac.Z1=as.numeric(my.calc.fac.Z1))
    my.calc = cbind(my.calc,ncells=my.lengths$x)
    
  } else if ( calc == "freq" ) {
    ########## CALCULATE FREQUENCIES
    ### frequency of feature Z1
    my.calc = aggregate(data[,idx.Z1],by=list(fXY),function(x) {
      y= round( 100 * length(which(x >= cutoffs[idx.Z1])) / length(x))
      return(y)
    })
    
    # bin color factor Z1
    my.calc.fac.Z1 = cut(my.calc$x,breaks=seq(0,100,by=10),labels=1:10,include.lowest=TRUE)
    names(my.calc.fac.Z1) = my.calc$x
    
    my.lengths = aggregate(data[,idx.Z1], by=list(fXY), length)
    
    ### combine all frequencies in one table
    my.calc = cbind(my.calc,fac.Z1=as.numeric(my.calc.fac.Z1))
    my.calc = cbind(my.calc,ncells=my.lengths$x)
  }
  this$my.calc = my.calc
  this$my.calc.fac.Z1 = my.calc.fac.Z1
  ########## DONE CALCULATE FREQUENCIES 
  
  ########## PLOT DOUBLE FREQUENCY BINS
  brackets.open = c("(","[")
  for (x in 1:nrow(tab)) {
    for (y in 1:ncol(tab)) {
      if ( tab[rownames(tab)[x],colnames(tab)[y]]>=mincells ) {
        brackets.idx.x = brackets.idx.y = 1
        if (x==1) brackets.idx.x = 2
        if (y==1) brackets.idx.y = 2
        
        fact = as.factor(paste(brackets.open[brackets.idx.x],rownames(tab)[x],',',as.numeric(rownames(tab)[x])+binsize,'] ',
                               brackets.open[brackets.idx.y],rownames(tab)[y],',',as.numeric(rownames(tab)[y])+binsize,']',sep=''))
        
        idx=which(as.character(fact)==as.character(my.calc$Group.1))
        
        tab[x,y] = my.calc[idx,'x']
      } else if ( tab[x,y]>0 & calc!="freq") {
        tab[x,y] = min.range.Z1
      } else {
        tab[x,y] = NA
      }
    }
  }
  
  this$tab3D = round(tab,3)
}