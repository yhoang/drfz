#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019

### uncomment if not created yet
# fcs = new.env()


######## CREATE TRIPLOT BINS WITH FREQUENCY OF DOUBLE POSITIVES ####
##### @param
# file.idx.list     file index
# feat.X            name of feature X
# feat.Y            name of feature Y
# feat.Z1           name of feature Z1
# feat.Z2           name of feature Z2
# col              color palette for frequency of Z1+/Z2+
# binsize           bin size
# mincells          minimum amount of cells in a bin
# plot.range        plot range, first x axis, second y axis, default= c(2,12,2,12)  
# sepa="."          separater for the feature names 
fcs$bintriploT_freq_doublepos <- function(
  project.idx, 
  file.idx, 
  feat.X, 
  feat.Y, 
  feat.Z1, 
  feat.Z2, 
  col, 
  binsize, 
  mincells, 
  maxfreq=100, 
  plot.range=c(2,12,2,12), sepa=".") 
{ 
    this=fcs 

    if (FALSE) {
      project.idx = fcs$project.idx
      file.idx = fcs$file.idx
      feat.X = fcs$feat.X
      feat.Y = fcs$feat.Y
      feat.Z1 = fcs$feat.Z1
      feat.Z2 = fcs$feat.Z2
      binsize = 0.2
      mincells = 10
      col = color
      maxfreq = maxfreq
      plot.range = c(1,12,1,12)
      sepa="."
    }

    ### remove all signs and write anything with capitals
    feat.X.clean = gsub("[^[:alnum:]]","",toupper(feat.X))
    feat.Y.clean = gsub("[^[:alnum:]]","",toupper(feat.Y))
    feat.Z1.clean = gsub("[^[:alnum:]]","",toupper(feat.Z1))
    feat.Z2.clean = gsub("[^[:alnum:]]","",toupper(feat.Z2))

    ### axes range
    xmin.val = plot.range[1]
    xmax.val = plot.range[2]
    ymin.val = plot.range[3]
    ymax.val = plot.range[4]

    ############ LOAD DATA AND CUTOFFS
    data = this$loaddata(project.idx,file.idx)
    cutoffs = as.numeric(this$current.staintable[which(this$current.staintable$file_ID==file.idx),5])
    
    #cutoffs[idx.Z1]=10
    #cutoffs[idx.Z2]=7.5
    #cutoffs[idx.Z2]=8.663

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
    idx.Z1 = which(features.clean == feat.Z1.clean)
    idx.Z2 = which(features.clean == feat.Z2.clean)
  
    printf("idx X=%s Y=%s Z1=%s Z2=%s",idx.X,idx.Y,idx.Z1,idx.Z2)
    printf("cutoffs X=%s Y=%s Z1=%s Z2=%s",cutoffs[idx.X],cutoffs[idx.Y],cutoffs[idx.Z1],cutoffs[idx.Z2])
    if ( is.na(idx.X) | is.na(idx.Y) | is.na(idx.Z1) | is.na(idx.Z2) ) stop(sprintf("Either feat.X (%s) or feat.Y (%s) not in file index %s",feat.X, feat.Y, i))

    feat.Z1 = paste(sprintf("%s(%s)",feat.Z1,cutoffs[idx.Z1]))
    feat.Z2 = paste(sprintf("%s(%s)",feat.Z2,cutoffs[idx.Z2]))
    
    ### calculate quadrant percentages
    # cut all cells which are not producing cells
    ncells.total = nrow(data)
    
    tdata.q1 = data[which( data[,idx.X]<cutoffs[idx.X] &  data[,idx.Y]<cutoffs[idx.Y] ),]
    tdata.q2 = data[which( data[,idx.X]>=cutoffs[idx.X] &  data[,idx.Y]<cutoffs[idx.Y] ),]
    tdata.q3 = data[which( data[,idx.X]>=cutoffs[idx.X] &  data[,idx.Y]>=cutoffs[idx.Y] ),]
    tdata.q4 = data[which( data[,idx.X]<cutoffs[idx.X] &  data[,idx.Y]>=cutoffs[idx.Y] ),]
    
    this$tdata.q1 = tdata.q1
    
    q1.prodcells.num = nrow(tdata.q1[which(tdata.q1[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q1[,idx.Z2]>=cutoffs[idx.Z2]),])
    q2.prodcells.num = nrow(tdata.q2[which(tdata.q2[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q2[,idx.Z2]>=cutoffs[idx.Z2]),])
    q3.prodcells.num = nrow(tdata.q3[which(tdata.q3[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q3[,idx.Z2]>=cutoffs[idx.Z2]),])
    q4.prodcells.num = nrow(tdata.q4[which(tdata.q4[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q4[,idx.Z2]>=cutoffs[idx.Z2]),])
    q1.prodcells = round( (100 * q1.prodcells.num / nrow(tdata.q1) ), 2)
    q2.prodcells = round( (100 * q2.prodcells.num / nrow(tdata.q2) ), 2)
    q3.prodcells = round( (100 * q3.prodcells.num / nrow(tdata.q3) ), 2)
    q4.prodcells = round( (100 * q4.prodcells.num / nrow(tdata.q4) ), 2)
    
    q1.prodcells.total = round( 100 * nrow(tdata.q1[which(tdata.q1[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q1[,idx.Z2]>=cutoffs[idx.Z2]),]) / ncells.total,2)
    q2.prodcells.total = round( 100 * nrow(tdata.q2[which(tdata.q2[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q2[,idx.Z2]>=cutoffs[idx.Z2]),]) / ncells.total,2)
    q3.prodcells.total = round( 100 * nrow(tdata.q3[which(tdata.q3[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q3[,idx.Z2]>=cutoffs[idx.Z2]),]) / ncells.total,2)
    q4.prodcells.total = round( 100 * nrow(tdata.q4[which(tdata.q4[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q4[,idx.Z2]>=cutoffs[idx.Z2]),]) / ncells.total,2)
    
    printf("cells quadrant:                %s %s %s %s",nrow(tdata.q1),nrow(tdata.q2),nrow(tdata.q3),nrow(tdata.q4))
    printf("cells double prod in quadrant: %s %s %s %s",q1.prodcells.num,q2.prodcells.num,q3.prodcells.num,q4.prodcells.num)
    printf("%% cells double prod in quadrant: %s %s %s %s",q1.prodcells,q2.prodcells,q3.prodcells,q4.prodcells)
    printf("cells double prod in quadrant to total: %s %s %s %s",q1.prodcells.total,q2.prodcells.total,q3.prodcells.total,q4.prodcells.total)

    ### construct bin table with number of cells per bin
    fX = cut(data[,idx.X],breaks=seq(xmin.val,xmax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
    fY = cut(data[,idx.Y],breaks=seq(ymin.val,ymax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
    fXY=as.factor(paste(fX,fY))

    this$fX=fX
    this$fY=fY

    tab = table(fX,fY)
    colnames(tab)=seq(ymin.val,ymax.val-binsize,by=binsize)
    rownames(tab)=seq(xmin.val,xmax.val-binsize,by=binsize)
    tab3D = tab3D_z2 = tab

    printf("length(bins with tab>=mincells)=%s",length(which(tab>=mincells)))
    
    ########## CALCULATE FREQUENCIES
    ### frequency of feature Z1
    my.calc = aggregate(data[,idx.Z1], 
        by=list(fXY),
        function(x) {
            y = round( 100 * length(which(x >= cutoffs[idx.Z1])) / length(x))
            return(y)
    })
    # bin color factor Z1
    my.calc.fac.Z1 = cut(my.calc$x,breaks=seq(0,100,by=10),labels=1:10,include.lowest=TRUE)
    names(my.calc.fac.Z1) = my.calc$x

    ### frequency of feature Z2        
    freq.Z2 = aggregate(data[,idx.Z2],
        by=list(fXY),
        function(x) {
            y = round( 100 * length(which(x >= cutoffs[idx.Z2])) / length(x))
            return(y)
    })
    # bin color factor Z2
    my.calc.fac.Z2 = cut(freq.Z2$x,breaks=seq(0,100,by=10),labels=1:10,include.lowest=TRUE)
    names(my.calc.fac.Z2) = freq.Z2$x
    
    ### get only data table where Z1 and Z2 is produced
    data.double = data[ which( (data[,idx.Z1]>=cutoffs[idx.Z1]) & (data[,idx.Z2]>=cutoffs[idx.Z2]) ),]
    ### construct bin table with number of cells per bin
    fX.double = cut(data.double[,idx.X],breaks=seq(xmin.val,xmax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
    fY.double = cut(data.double[,idx.Y],breaks=seq(ymin.val,ymax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
    fXY.double = as.factor(paste(fX.double,fY.double))

    length.double = aggregate(data.double[,1],by=list(fXY.double),length)
    rownames(length.double)=length.double$Group.1

    length.all = aggregate(data[,idx.Z1], by=list(fXY), length)
    rownames(length.all)=length.all$Group.1

    freq.double = merge(length.all,length.double,by="row.names",all.x=TRUE)
    freq.double = freq.double[,-c(2,4)]
    freq.double = cbind( freq.double, round(freq.double[,3]/freq.double[,2] * 100))
    # bin color factor double producer Z1+Z2
    my.calc.fac.double = cut(freq.double[,4],breaks=seq(0,maxfreq,by=maxfreq/10),labels=1:10,include.lowest=TRUE)
    names(my.calc.fac.double) = freq.double[,4]


    ### combine all frequencies in one table
    my.calc = cbind(my.calc,fac.Z1=as.numeric(my.calc.fac.Z1))
    my.calc = cbind(my.calc,freq.Z2=freq.Z2$x)
    my.calc = cbind(my.calc,fac.Z2=as.numeric(my.calc.fac.Z2))
    my.calc = cbind(my.calc,freq.double=freq.double[,4])
    my.calc = cbind(my.calc,fac.double=as.numeric(my.calc.fac.double))
    my.calc = cbind(my.calc,ncells=length.all$x)
    my.calc = cbind(my.calc,ncells.double=freq.double[,3])

    maxfreq.real = max(this$my.calc[which(!is.na(this$my.calc$freq.double) & this$my.calc$ncells>9),6])
    printf("MAXIMUM FREQUENCY = %s",maxfreq.real)
    if ( maxfreq.real > maxfreq ) print("    !!!!! MAXFREQ SHOULD BE HIGHER !!!!!")

    this$my.calc = my.calc
    this$my.calc.fac.Z1 = my.calc.fac.Z1
    this$my.calc.fac.Z2 = my.calc.fac.Z2
    this$my.calc.fac.double = my.calc.fac.double

    ########## DONE CALCULATE FREQUENCIES 

    cols = fcs$colors(col)
    cols.heat = cols[,as.numeric(my.calc.fac.double)]
    this$cols.heat = cols.heat

    brackets.open = c("(","[")
    ########## PLOT DOUBLE FREQUENCY BINS
    for (x in rownames(tab)) {
        for (y in colnames(tab)) {
            if ( tab[x,y]>=mincells ) {
                
                #fact = as.factor(paste('(',x,',',as.numeric(x)+binsize,'] ','(',y,',',as.numeric(y)+binsize,']',sep=''))
                brackets.idx.x = brackets.idx.y = 1
                if (x==0) brackets.idx.x = 2
                if (y==0) brackets.idx.y = 2
                      
                fact = as.factor(paste(brackets.open[brackets.idx.x],x,',',as.numeric(x)+binsize,'] ',
                            brackets.open[brackets.idx.y],y,',',as.numeric(y)+binsize,']',sep=''))

                idx = which(as.character(fact)==as.character(my.calc$Group.1))
                if (length(cols.heat[,idx])!=0) {
                  if ( !is.na(cols.heat[1,idx]) ) {  
                    
                      #printf("fact=%s idx=%s col=%s",fact,idx,paste(cols.heat[,idx],collapse=" "))
                      rect(x,y,as.numeric(x)+binsize,as.numeric(y)+binsize
                          #,col=cols[my.calc[idx,'fac.double']]
                          ,col=eval(parse(text=paste0("rgb(",paste0(cols.heat[,idx],collapse=","),",maxColorValue=255)")))
                          ,border=NA)
                      
                      tab3D[x,y] = my.calc[idx,'x']
                      tab3D_z2[x,y] = my.calc[idx,'freq.Z2']
                  }
                }
            } else {
                tab3D[x,y] = tab3D_z2[x,y] = NA
            }
        }
    }
    this$tab3D.Z1 = tab3D
    this$tab3D.Z2 = tab3D_z2
    ########### DONE PLOT BINS

    this$addProdline(cutoffs[c(idx.X,idx.Y)])

    ########## PLOT LEGEND
    start.legend = par()$usr[2] - 11.5*binsize       
    for (step in 1:10) {
        ###### from bottom to top
        ### legend Z1
        rect(xleft = start.legend+step*binsize
            ,ybottom = par()$usr[3]+0.10*(par()$usr[3]+par()$usr[4])-4*binsize
            ,xright = start.legend+(step+1)*binsize
            ,ytop = par()$usr[3]+0.10*(par()$usr[3]+par()$usr[4])-2*binsize
            #,col=cols[step]
            ,col=eval(parse(text=paste0("rgb(",paste0(cols[,step],collapse=","),",maxColorValue=255)")))
            ,border=NA
        )
    }

    ### x axis label
    text(x = 0.5*(par()$usr[1]+par()$usr[2])
        ,y = par()$usr[3] - 1
        ,label = sprintf("%s(%s)",feat.X,cutoffs[idx.X])
        ,xpd = TRUE
    )
    ### y axis label
    text(x = par()$usr[1] - 0.8
        ,y = 0.5*(par()$usr[3]+par()$usr[4])
        ,label = sprintf("%s(%s)",feat.Y,cutoffs[idx.Y])
        ,xpd = TRUE
        ,srt=90
    )
    ######## from bottom to top
    ### legend title Z1
    text(x = start.legend + 0.8*binsize
        ,y = par()$usr[3]+0.075*(par()$usr[3]+par()$usr[4])-3*binsize
        ,label=feat.Z1,cex=0.7,pos=2
    )
    ### legend title Z2
    text(x = start.legend + 0.8*binsize
        ,y = par()$usr[3]+0.075*(par()$usr[3]+par()$usr[4])-1*binsize
        ,label=feat.Z2,cex=0.7,pos=2
    )
    ########## DONE PLOT LEGEND

    calc="freq"
    ########## ADD FILE NAME
    if (calc=="freq"){
        title(main=sprintf("File: %s(MFrange=%s,MFreal=%s,mincells=%s)",fcs$shortenFilename(fcs$file.list[file.idx]),maxfreq,maxfreq.real,mincells)
            ,line=1,cex.main=0.7,adj=0)
    } else {
        title(main=sprintf("File: %s",fcs$shortenFilename(fcs$file.list[file.idx]))
            ,line=1,cex.main=0.7,adj=0)
    }
}
