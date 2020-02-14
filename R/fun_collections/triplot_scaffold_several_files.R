#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019

### uncomment if not created yet
# fcs = new.env()

######## CREATE TRIPLOT BIN SCAFFOLD FROM SEVERAL SAMPLES #######
##### @param
# project.idx.list  list of project idxs according to file idx list
# file.idx.list     list of files to conclude in the scaffold
# feat.X            name of feature X
# feat.Y            name of feature Y
# binsize           bin size
# mincells          minimum amount of cells in a bin
# plot.range        plot range, first x axis, second y axis, default= c(2,12,2,12)
# sepa="."          separater for the feature names 
fcs$bintriploT_construct <- function(
  project.idx.list, 
  file.idx.list, 
  feat.X, 
  feat.Y, 
  binsize, 
  mincells, 
  plot.range=c(2,12,2,12), 
  sepa=".") 
{ 
    this=fcs 

    if (FALSE) {
      project.idx.list = fcs$project.idx.list
      file.idx.list = fcs$file.idx.list
      feat.X = "CD16"
      feat.Y = "CD14"
      binsize = 0.2
      mincells = 5
      plot.range = c(0,8,0,8)
    }

    ### remove all signs and write anything with capitals
    feat.X.clean = gsub("[^[:alnum:]]","",toupper(feat.X))
    feat.Y.clean = gsub("[^[:alnum:]]","",toupper(feat.Y))
    
    xmin.val = plot.range[1]
    xmax.val = plot.range[2]
    ymin.val = plot.range[3]
    ymax.val = plot.range[4]
    
    ########## CREATE BIN CONSTRUCT FIRST
    ### NO FREQUENCY CALCULATION HERE
    loop_i = 1
    ### START LOOP i FOR FILE INDEX i
    for (i in 1:length(file.idx.list)) {
        data = this$loaddata(project.idx.list[i],file.idx.list[i])

        features = colnames(data) 

        if (sepa == ".") {
          ### remove all signs and write anything with capitals
          features.clean = gsub("[^[:alnum:]]","",make.unique(unlist(lapply(features,function(x) {
              len = length(strsplit(x,"[.]")[[1]])
              y = toupper(strsplit(x,"[.]")[[1]][1])
              paste(y, collapse=".")
          }))))
        } else {
          ### remove all signs and write anything with capitals
          features.clean = gsub("[^[:alnum:]]","",make.unique(unlist(lapply(features,function(x) {
          len = length(strsplit(x,"[_]")[[1]])
          y = toupper(strsplit(x,"[_]")[[1]][2])
          paste(y, collapse=".")
        }))))
        }
        

        idx.X = which(features.clean == feat.X.clean)
        idx.Y = which(features.clean == feat.Y.clean)

        if ( is.na(idx.X) | is.na(idx.Y) ) stop(sprintf("Either feat.X (%s) or feat.Y (%s) not in file index %s",feat.X, feat.Y, i))

        fX = cut(data[,idx.X],breaks=seq(xmin.val,xmax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
        fY = cut(data[,idx.Y],breaks=seq(ymin.val,ymax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
        if ( loop_i == 1 ) {
            ### if first for loop run
            ### construct bin table with number of cells per bin
            global.tab = table(fX,fY)
            colnames(global.tab)=seq(ymin.val,ymax.val-binsize,by=binsize)
            rownames(global.tab)=seq(xmin.val,xmax.val-binsize,by=binsize)

            printf("loop #%s: bins(global.tab>=mincells)=%s",loop_i,length(which(global.tab>=mincells)))
        } else {
            ### construct bin table with number of cells per bin
            tab = table(fX,fY)
            colnames(tab)=seq(ymin.val,ymax.val-binsize,by=binsize)
            rownames(tab)=seq(xmin.val,xmax.val-binsize,by=binsize)
            
            updated.bins = 0
            for (x in rownames(tab)) {
              for (y in colnames(tab)) {
                # if there is a bin where in current file there are more cells than global
                if ( (tab[x,y]>global.tab[x,y])){#} & (tab[x,y]>=mincells) & (global.tab[x,y]<mincells) ) {
                  global.tab[x,y] = tab[x,y]
                  updated.bins = updated.bins + 1
                  #printf("#%s: x=%s y=%s #=%s",updated.bins,x,y,tab[x,y])
                }
              }
            }
            printf("loop #%s: bins(tab>=mincells)=%s, updated.bins=%s, global.bins=%s",loop_i,length(which(tab>=mincells)),updated.bins,length(which(global.tab>=mincells)))
        }
        loop_i = loop_i + 1
    }
    fXY=as.factor(paste(fX,fY))

    ### start plot frame
    par(oma=c(0.5,1,2,1),mar=c(3,3,4,2))
    plot(1,type='n',frame.plot=FALSE,xlim=c(xmin.val,xmax.val+10*binsize),axes=FALSE
        ,ylim=c(ymin.val-2.5*binsize,ymax.val+5*binsize),xlab=NA,ylab=NA,cex.lab=1,cex.axis=1,mgp=c(1.7,0.4,0))
    box(lwd=0.8,col="darkgrey")

    ### draw an axis on the bottom
    ### draw an axis on the left
    asinh.scale = c(asinh(-1000),asinh(-100),asinh(-10),asinh(0),asinh(10),asinh(100),asinh(1000),asinh(10000),asinh(100000),asinh(1000000),asinh(10000000),asinh(100000000))
    axis(side=1, at=asinh.scale,labels=c("1e-3","1e-2","1e-1","0","1e1","1e2","1e3","1e4","1e5","1e6","1e7","1e8")
        ,las=1,cex.axis=1,mgp=c(1.7,0.4,0),col="darkgrey")
    axis(side=2, at=asinh.scale,labels=c("1e-3","1e-2","1e-1","0","1e1","1e2","1e3","1e4","1e5","1e6","1e7","1e8")
        ,las=3,cex.axis=1,mgp=c(1.7,0.4,0),col="darkgrey")

    ### grid
    xgrid.steps=seq((xmin.val),(xmax.val),by=2)
    ygrid.steps=seq((ymin.val),(ymax.val),by=2)
    abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
    
    ##### plot bin construct in grey first
    for (x in rownames(global.tab)) {
        for (y in colnames(global.tab)) {
            if ( global.tab[x,y]>=mincells ) {
                rect(x,y,as.numeric(x)+binsize,as.numeric(y)+binsize,col="lightgrey",border=NA)
            }
        }
    }

    this$global.tab = global.tab

    ########## ADD DATABASE AND FILE NAMES
    title(main=sprintf("DB: %s; bin construct of %s files:",sub(".sqlite3","",fcs$db.name),length(fcs$file.idx.list))
        ,line=3,cex.main=0.7,adj=0)
    title(main=sprintf("%s",paste(fcs$shortenFilename(fcs$file.list[fcs$file.idx.list]),collapse=" / "))
        ,line=2,cex.main=0.7,adj=0)
}
