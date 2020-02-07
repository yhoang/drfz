#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2018


####################################################################
######## LOAD FUNCTIONS WITHOUT ACTUALLY CALLING THEM ##############
##### @param
# project.idx   project index
# file.idx      file index
##### returns data table
fcs$loaddata <- function(
  project.idx,
  file.idx)  
{
    this=fcs

    ### list all tables in database
    this$table.list=dbListTables(this$conn)

    ### get every index which has metadata
    metadata.idx=grep("markerIdentity|colnameIndex|fileIdentity|fileIndex|UserHistory|Classification|equipmentInfo|fileComments|SPILL",this$table.list)

    this$metadata.list = vector()
    df.num = 0
    for (i in metadata.idx){
        df.num = df.num + 1 
        this$metadata.list[df.num] = this$table.list[i]
    }

    ### now only the table names are listed which have our data (fluorescence intensity for each feature)
    this$project.list = this$table.list[-metadata.idx]
    ### if there are several projects in the database, choose one
    this$current.project = this$project.list[project.idx]

    this$current.filetable = this$getDFtable(paste0(this$current.project,"_fileIdentity"))
    this$file.list = this$current.filetable[,2]

    this$current.file = this$file.list[file.idx]
    this$current.staintable = this$getDFtable(paste0(this$current.project,"_markerIdentity"))
    this$current.vars = this$getVariables(index=file.idx)

    data = this$getFile(this$current.project,file.idx)

    data
}

###### helpful print function
printf <- function(...) invisible(print(sprintf(...)))

######## GET DATA TABLE ############################################
##### @param
# table         project name
# fileidx       file index
# stain         staining vector
fcs$getFile <- function (
  table,
  fileidx,
  stain=NA) 
{
    this=fcs
    
    # get table from database
    data=dbGetQuery(this$conn,paste0("SELECT * FROM ",table," WHERE file_ID == '",fileidx,"'"))
    printf("do getFile: %s cells from table='%s' with fileidx=%s",nrow(data),table,fileidx)

    # cut column "file_ID" and ignore columns with NAs
    data$file_ID=NULL
    col.NA=NULL
    for (i in 1:ncol(data)) {
        if (any(is.na(data[,i]))) col.NA=c(col.NA,i)
    }
    if (!is.null(col.NA)) data=data[,-col.NA]

    data = asinh(data) 

    if ( !is.na(stain) ) { 
            colnames(data) = stain 
    } else { 
        colnames(data) = this$current.vars 
    }

    #print("Table header")
    #print(head(data,2))
    print("Feature names:")
    print(colnames(data))

    cat("Done getFile\n\n")

    this$data = data

    data
}
 
######## GET FEATURE SHORT NAMES ###################################
##### @param
# staintable    manual stain table (optional)
# index         file index (optional)    
fcs$getVariables <- function (staintable=NA,index=NA) {
    this=fcs

    if (is.na(staintable)) staintable=this$current.staintable
    
    if (is.na(index)) {
        if (exists("current.file",envir=fcs)) {
            file=tclvalue(tkget(this$tkchoosefile))
        } else {
            file=this$file.list[1]
        }
        index=this$current.filetable[which(this$current.filetable[,2]==file),1]
    }
    printf("do getVariables from table=%s: file idx=%s",this$current.project,index)

    vars=staintable[which(staintable[,1]==index),4]
    
    #printf("vars=%s",paste(vars,collapse=" "))
    print(vars)
    cat("Done getVariables\n\n")

    vars
}

######## GET DATA FRAME (META DATA) ################################
##### @param
# table         project name
fcs$getDFtable <- function (table) {
    this=fcs

    table.df=dbGetQuery(this$conn,paste("SELECT * FROM ", table))

    table.df
}

######## CLEAN FILE NAME ###########################################
##### @param
# name          file name
fcs$shortenFilename <- function(name) {
    ### shortens file name
    name = sub(".csv$","",name)
    name = sub(".fcs$","",name)
    name = sub("^temp\\d\\d_","",name)
    name = sub("export_","",name)
    name = sub("Stain_*","S\\1",name)
    name = sub("stain_","S\\1",name)
    name = sub("Stain*","S\\1",name)
    name = sub("stain*","S\\1",name)
    name = sub("Surface*","S\\1",name)
    name = sub("surface*","S\\1",name)
    name = sub("Stimulated_","Stim_",name)
    name = sub("stimulated_","stim_",name)
    name = sub("^Gruppe","Grp",name)
    name = sub("Gruppe","grp",name)
    name = sub("gruppe","grp",name)
    name = sub("group","grp",name)
    name = sub("Day","d",name)
    name = sub("day","d",name)
    name = sub("Antibodies","Ab",name)
    name = sub("Antibody","Ab",name)
    name = sub("Untreated","Untr",name)
    name = sub("untreated","untr",name)
    name = sub("PMA_Iono","PMAIono",name)
    
    name = sub("(\\d{4})_(\\d{2})_(\\d{2})","\\1\\2\\3",name,perl=TRUE)

    ### if file name starts with 0-2 letters and follows with 6-8 digits
    # (short or long version of date) YY/MM/DD or YYYY/MM/DD)
    # then code only with last 3 digits
    name = sub("^\\w\\w\\d{3,5}(\\d)(\\d)(\\d)_","\\1\\2\\3_",name,perl=TRUE)

    name
}

######## CHOOSE COLORS #############################################
# color         name of color options
fcs$colors <- function(color) {
    ### light to dark

    if (color == "green") {
        # green
        cols = c("#F4FBEF","#E0F3D3","#C3E1AD","#ACD88D","#93CC6B","#7DC24C","#71B441","#63AB30","#56972B","#4F8726")
        cols = col2rgb(cols)
    } else if (color == "orange") {
        # orange
        cols = c("#FFF1DC","#FEE2B9","#FED69A","#FDC168","#FCB140","#FBA018","#F29202","#E78C04","#DB8400","#C97A03")
        cols = col2rgb(cols)
    } else if (color == "blue") {
        # blue
        cols = c("#F6F9FE","#ECF2FD","#D5E2FC","#CBDCFB","#BDD2F9","#B3CCF9","#A9C5F9","#9EBCF5","#88AFF2","#709DF0")
        cols = col2rgb(cols)
    } else if (color == "lila") {
        # lila
        cols = c("#4f2a4f","#60335f","#713b70","#824482","#934c93","#a554a4","#b062af","#b973b9","#c283c2","#cb94ca")
        cols = cols[10:1]
        cols = col2rgb(cols)
    } else if (color == "turquoise") {
        # turquoise
        cols = c("#2a444e","#33535f","#3b6170","#437081","#4c7f93","#548da4","#629ab0","#72a5b9","#82b0c2","#93bbca")
        cols = cols[10:1]
        cols = col2rgb(cols)
    } else if (color == "green2") {
        # green2
        cols = c("#284a3a","#305b47","#396c54","#417d61","#498e6e","#51a07b","#5dad88","#6eb795","#7ec0a1","#8fc8ae")
        cols = cols[10:1]
        cols = col2rgb(cols)
    } else if (color == "red") {
        # red
        cols = c("#4d2929","#5e3232","#6f3a3b","#804343","#924b4b","#a35353","#af6061","#b87171","#c18182","#ca9292")
        cols = cols[10:1]
        cols = col2rgb(cols)
    } else if (color == "blue2") {
        # blue2
        cols = c("#2a2f4e","#33395f","#3b4270","#434c81","#4c5593","#545fa4","#626cb0","#727cb9","#828bc2","#939bca")
        cols = cols[10:1]
        cols = col2rgb(cols)
    } else if (color == "pink") {
        # pink
        cols = c("#4d2935","#5e3241","#6f3a4c","#804357","#924b63","#a3536e","#af607b","#b87189","#c18197","#ca92a5")
        cols = cols[10:1]
        cols = col2rgb(cols)
    } else if (color == "red2") {
        # red2
        cols = c("#FFBEB8","#FE948A","#FE7367","#FE4E3E","#FE1B06","#DF1401","#D11100","#BC1101","#BC1101","#930D01")
        cols = col2rgb(cols)
    } else if (color == "rainbow") {
        cols = c("#001AC0","#00A7D0","#00F0EC","#00F884","#81FF42","#D1FF00","#FFFF00","#FFDB00FF","#FF6D00FF","#FF0000FF")
    } else if (color == "heat") {
        cols = heat.colors(11)
        cols = cols[10:1]
    } else if (color == "terrain") {
        cols = terrain.colors(11)
        cols = cols[10:1]
    }

    cols
}



########################### DATA MANIPULATION

######## ADD CUTOFF LINES ##########################################
##### @param
# cutoffs       cutoff vector (x,y)
fcs$addProdline <- function (cutoffs=c(0,0)) {
    this=fcs
    # add production cutoff line if provided
    if ( cutoffs[1] > 0 ) {
        abline(v=cutoffs[1],col="darkgrey")
    }
    if ( cutoffs[2] > 0 ) {
        abline(h=cutoffs[2],col="darkgrey")
    }
}



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
        ,las=1,cex.axis=0.8,mgp=c(1.7,0.4,0),col="darkgrey")
    axis(side=2, at=asinh.scale,labels=c("1e-3","1e-2","1e-1","0","1e1","1e2","1e3","1e4","1e5","1e6","1e7","1e8")
        ,las=3,cex.axis=0.8,mgp=c(1.7,0.4,0),col="darkgrey")

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

######## CREATE TRIPLOT BINS WITH TWO FREQUENCIES/MFI IN ONE FILE ##
##### @param
# file.idx.list     file index
# feat.X            name of feature X
# feat.Y            name of feature Y
# feat.Z1           name of feature Z1
# feat.Z2           name of feature Z2
# calc              calculation method
# col1              monochrome color palette for frequency of Z1
# col2              monochrome color palette for frequency of Z2
# binsize           bin size
# mincells          minimum amount of cells in a bin
# plot.range        plot range, first x axis, second y axis, default= c(2,12,2,12)
fcs$bintriploT_double_onefile <- function(
  project.idx, 
  file.idx, 
  feat.X, 
  feat.Y, 
  feat.Z1, 
  feat.Z2, 
  calc, 
  col1, 
  col2, 
  binsize, 
  mincells, 
  plot.range=c(2,12,2,12)) 
{ 
    this=fcs 

    if (FALSE) {
        project.idx=1
        file.idx = 2
        feat.X = "PD-1"
        feat.Y = "IFNg"
        feat.Z1 = "IL21"
        feat.Z2 = "Bcl6"
        calc = "MFI"
        binsize = 0.2
        mincells = 10
        col1 = this$col2
        col2 = this$col1
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
    idx.Z2 = which(features.clean == feat.Z2.clean)

    if ( is.na(idx.X) | is.na(idx.Y) | is.na(idx.Z1) | is.na(idx.Z2) ) stop(sprintf("Either feat.X (%s) or feat.Y (%s) not in file index %s",feat.X, feat.Y, i))

    feat.Z1 = paste(sprintf("%s(%s)",feat.Z1,cutoffs[idx.Z1]))
    feat.Z2 = paste(sprintf("%s(%s)",feat.Z2,cutoffs[idx.Z2]))
    
    ### construct bin table with number of cells per bin
    fX = cut(data[,idx.X],breaks=seq(xmin.val,xmax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
    fY = cut(data[,idx.Y],breaks=seq(ymin.val,ymax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
    fXY=as.factor(paste(fX,fY))

    tab = table(fX,fY)
    colnames(tab)=seq(ymin.val,ymax.val-binsize,by=binsize)
    rownames(tab)=seq(xmin.val,xmax.val-binsize,by=binsize)

    tab3D.Z1 = tab3D.Z2 = tab

    printf("length(bins with tab>=mincells)=%s",length(which(tab>=mincells)))

    if ( calc == "MFI" ) {
        my.lengths = aggregate(data[,idx.Z1],by=list(fXY),length)
        idx.len = which(my.lengths$x>=mincells)

        ### start with MFI calculation of Z1
        my.calc = aggregate(data[,idx.Z1],by=list(fXY),mean)

        min.range.Z1 = floor(min(my.calc[idx.len,'x'])*10)/10
        max.range.Z1 = ceiling(max(my.calc[idx.len,'x'])*10)/10

        # get steps for Z1
        step=round(diff(range(max.range.Z1,min.range.Z1))/10,2) 
        steps.Z1=seq(min.range.Z1,max.range.Z1,by=step)
        
        # bin color factor Z2
        my.calc.fac.Z1=cut(my.calc$x,breaks=steps.Z1,labels=1:10,include.lowest=TRUE)
        names(my.calc.fac.Z1) = my.calc$x

        ### MFI of feature Z2
        MFI.Z2 = aggregate(data[,idx.Z2],by=list(fXY),mean)

        min.range.Z2 = floor(min(MFI.Z2[idx.len,'x'])*10)/10
        max.range.Z2 = ceiling(max(MFI.Z2[idx.len,'x'])*10)/10

        # get steps for Z2
        step=round(diff(range(max.range.Z2,min.range.Z2))/10,2) 
        steps.Z2=seq(min.range.Z2,max.range.Z2,by=step)
        
        # bin color factor Z2
        my.calc.fac.Z2=cut(MFI.Z2$x,breaks=steps.Z2,labels=1:10,include.lowest=TRUE)
        names(my.calc.fac.Z2) = MFI.Z2$x

        ### combine all frequencies in one table
        my.calc = cbind(my.calc,fac.Z1=as.numeric(my.calc.fac.Z1))
        my.calc = cbind(my.calc,calc.Z2=MFI.Z2$x)
        my.calc = cbind(my.calc,fac.Z2=as.numeric(my.calc.fac.Z2))
        my.calc = cbind(my.calc,ncells=my.lengths)

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

        ### frequency of feature Z2
        freq.Z2 = aggregate(data[,idx.Z2],by=list(fXY),function(x) {
            y= round( 100 * length(which(x >= cutoffs[idx.Z2])) / length(x))
            return(y)
        })
        # bin color factor Z2
        my.calc.fac.Z2 = cut(freq.Z2$x,breaks=seq(0,100,by=10),labels=1:10,include.lowest=TRUE)
        names(my.calc.fac.Z2) = freq.Z2$x


        my.lengths = aggregate(data[,idx.Z1], by=list(fXY), length)

        ### combine all frequencies in one table
        my.calc = cbind(my.calc,fac.Z1=as.numeric(my.calc.fac.Z1))
        my.calc = cbind(my.calc,calc.Z2=freq.Z2$x)
        my.calc = cbind(my.calc,fac.Z2=as.numeric(my.calc.fac.Z2))
        my.calc = cbind(my.calc,ncells=my.lengths$x)
    }
    this$my.calc = my.calc
    this$my.calc.fac.Z1 = my.calc.fac.Z1
    this$my.calc.fac.Z2 = my.calc.fac.Z2
    ########## DONE CALCULATE FREQUENCIES 

    ########## INITIATE COLORS 
    ### orange: light to dark
    cols.heat.Z1 = fcs$colors(col1)

    ### green: light to dark
    cols.heat.Z2 = fcs$colors(col2)

    cols.heat.Z1Z2 = round((cols.heat.Z2+cols.heat.Z1)/2)

    ### get color index vector
    col.Z1 = cols.heat.Z1[,my.calc.fac.Z1]
    col.Z2 = cols.heat.Z2[,my.calc.fac.Z2]
    ### mix colors
    col.mix = round((col.Z1+col.Z2)/2)
    ########## DONE INITITATE COLORS

    ########## PLOT DOUBLE FREQUENCY BINS
    for (x in rownames(tab)) {
        for (y in colnames(tab)) {
            if ( tab[x,y]>=mincells ) {
                fact=as.factor(paste('(',x,',',as.numeric(x)+binsize,'] ','(',y,',',as.numeric(y)+binsize,']',sep=''))
                idx=which(as.character(fact)==as.character(my.calc$Group.1))
                rect(x,y,as.numeric(x)+binsize,as.numeric(y)+binsize,col=eval(parse(text=paste0("rgb(",paste0(col.mix[,idx],collapse=","),",maxColorValue=255)"))),border=NA)
                
                tab3D.Z1[x,y] = my.calc[idx,'x']
                tab3D.Z2[x,y] = my.calc[idx,'calc.Z2']
            } else {
                tab3D.Z1[x,y] = tab3D.Z2[x,y] = NA
            }
        }
    }
    this$tab3D.Z1 = round(tab3D.Z1,3)
    this$tab3D.Z2 = round(tab3D.Z2,3)
    
    brackets.open = c("(","[")
    if ( calc == "freq" ) {
        ### MARK DOUBLE POSITIVES
        # MAKES SENSE ONLY IN FREQ
        for (x in rownames(tab)) {
            for (y in colnames(tab)) {
                if ( tab[x,y]>=mincells ) {
                    #fact=as.factor(paste('(',x,',',as.numeric(x)+binsize,'] ','(',y,',',as.numeric(y)+binsize,']',sep=''))
                    brackets.idx.x = brackets.idx.y = 1
                    if (x==0) brackets.idx.x = 2
                    if (y==0) brackets.idx.y = 2
                          
                    fact = as.factor(paste(brackets.open[brackets.idx.x],x,',',as.numeric(x)+binsize,'] ',
                                brackets.open[brackets.idx.y],y,',',as.numeric(y)+binsize,']',sep=''))

                    idx=which(as.character(fact)==as.character(my.calc$Group.1))
                    ### if bin has more than 50% Z1 and 50% Z2, then make a border
                    if (my.calc$fac.Z1[idx]>5 & my.calc$fac.Z2[idx]>5) {
                        rect(x,y,as.numeric(x)+binsize,as.numeric(y)+binsize
                            ,col=eval(parse(text=paste0("rgb(",paste0(col.mix[,idx],collapse=","),",maxColorValue=255)")))
                            ,border="white"
                        )
                    }
                }
            }
        }
    }
    ########### DONE PLOT BINS

    this$addProdline(cutoffs[c(idx.X,idx.Y)])

    ########## PLOT LEGEND
    start.legend = par()$usr[2] - 12*binsize       
    for (step in 1:10) {
        ###### from bottom to top
        ### legend Z1
        rect(xleft = start.legend+step*binsize
            ,ybottom = par()$usr[3]+0.10*(par()$usr[3]+par()$usr[4])+0*binsize
            ,xright = start.legend+(step+1)*binsize
            ,ytop = par()$usr[3]+0.10*(par()$usr[3]+par()$usr[4])+1*binsize
            ,col=eval(parse(text=paste0("rgb(",paste0(cols.heat.Z1[,step],collapse=","),",maxColorValue=255)"))),border=NA,pos=2
        )
        ### legend Z1/Z2 Mix
        rect(xleft = start.legend+step*binsize
            ,ybottom = par()$usr[3]+0.10*(par()$usr[3]+par()$usr[4])+2*binsize
            ,xright = start.legend+(step+1)*binsize
            ,ytop = par()$usr[3]+0.10*(par()$usr[3]+par()$usr[4])+3*binsize
            ,col=eval(parse(text=paste0("rgb(",paste0(cols.heat.Z1Z2[,step],collapse=","),",maxColorValue=255)"))),border=NA,pos=2
        )
        ### legend Z2
        rect(xleft = start.legend+step*binsize
            ,ybottom = par()$usr[3]+0.10*(par()$usr[3]+par()$usr[4])+4*binsize
            ,xright = start.legend+(step+1)*binsize
            ,ytop = par()$usr[3]+0.10*(par()$usr[3]+par()$usr[4])+5*binsize
            ,col=eval(parse(text=paste0("rgb(",paste0(cols.heat.Z2[,step],collapse=","),",maxColorValue=255)"))),border=NA,pos=2
        )
    }

    ### x axis label
    text(x = 0.4*(par()$usr[1]+par()$usr[2])
        ,y = par()$usr[3] - 1.5
        ,label = sprintf("%s(%s)",feat.X,cutoffs[idx.X])
        ,xpd = TRUE
    )
    ### y axis label
    text(x = par()$usr[1] - 1.3
        ,y = 0.3*(par()$usr[3]+par()$usr[4])
        ,label = sprintf("%s(%s)",feat.Y,cutoffs[idx.Y])
        ,xpd = TRUE
        ,srt=90
    )
    ######## from bottom to top
    ### legend title Z1
    text(x = start.legend + 0.5*binsize
        ,y = par()$usr[3]+0.105*(par()$usr[3]+par()$usr[4])
        ,label=feat.Z1,cex=0.7,pos=2
    )
    ### legend title Z2
    text(x = start.legend + 0.5*binsize
        ,y = par()$usr[3]+0.105*(par()$usr[3]+par()$usr[4])+4*binsize
        ,label=feat.Z2,cex=0.7,pos=2
    )
    ########## DONE PLOT LEGEND

    ########## ADD FILE NAME
    title(main=sprintf("File: %s",fcs$shortenFilename(fcs$file.list[file.idx]))
        ,line=1,cex.main=0.7,adj=0)
}


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


######## CREATE TRIPLOT BINS WITH FREQUENCY OF DOUBLE Z AND DECIDE TO USE NEG OR POS POPULATION ####
##### @param
# file.idx.list     file index
# feat.X            name of feature X
# feat.Y            name of feature Y
# feat.Z1           name of feature Z1
# feat.Z2           name of feature Z2
# popul             population of interest: 1 for positive, 0 for negative, default=c(1,1)
# col               color palette for frequency of Z1+/Z2+
# binsize           bin size
# mincells          minimum amount of cells in a bin
# plot.range        plot range, first x axis, second y axis, default= c(2,12,2,12)  
# sepa="."          separater for the feature names 
fcs$bintriploT_freq_doubleZ <- function(
  project.idx, 
  file.idx, 
  feat.X, 
  feat.Y, 
  feat.Z1, 
  feat.Z2,
  popul=c(1,1),
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
    popul=c(1,0)
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
  
  ### calculate quadrant percentages
  # cut all cells which are not producing cells
  ncells.total = nrow(data)
  
  tdata.q1 = data[which( data[,idx.X]<cutoffs[idx.X] &  data[,idx.Y]<cutoffs[idx.Y] ),]
  tdata.q2 = data[which( data[,idx.X]>=cutoffs[idx.X] &  data[,idx.Y]<cutoffs[idx.Y] ),]
  tdata.q3 = data[which( data[,idx.X]>=cutoffs[idx.X] &  data[,idx.Y]>=cutoffs[idx.Y] ),]
  tdata.q4 = data[which( data[,idx.X]<cutoffs[idx.X] &  data[,idx.Y]>=cutoffs[idx.Y] ),]
  
  this$tdata.q1 = tdata.q1
  
  
  ######## NEW: SWITCH SIGNS FOR POPULATION OF INTEREST
  if (popul[1]==1) {
    sign.Z1 = c(">=","<")
    feat.label.Z1 = paste0(feat.Z1,"+")
  } else {
    sign.Z1 = c("<=",">")
    feat.label.Z1 = paste0(feat.Z1,"-")
  }
  
  if (popul[2]==1) {
    sign.Z2 = c(">=","<")
    feat.label.Z2 = paste0(feat.Z2,"+")
  } else {
    sign.Z2 = c("<=",">")
    feat.label.Z2 = paste0(feat.Z2,"-")
  }
  feat.label.Z1 = paste(sprintf("%s(%s)",feat.label.Z1,cutoffs[idx.Z1]))
  feat.label.Z2 = paste(sprintf("%s(%s)",feat.label.Z2,cutoffs[idx.Z2]))
  
  # q1.prodcells.num
  eval(parse(text=paste("q1.prodcells.num = nrow(tdata.q1[which(tdata.q1[,idx.Z1]",sign.Z1[1],"cutoffs[idx.Z1] ",
                        "& tdata.q1[,idx.Z2]",sign.Z2[1],"cutoffs[idx.Z2]),])",sep="")))
  # q2.prodcells.num
  eval(parse(text=paste("q2.prodcells.num = nrow(tdata.q2[which(tdata.q2[,idx.Z1]",sign.Z1[1],"cutoffs[idx.Z1] ",
                        "& tdata.q2[,idx.Z2]",sign.Z2[1],"cutoffs[idx.Z2]),])",sep="")))
  # q3.prodcells.num
  eval(parse(text=paste("q3.prodcells.num = nrow(tdata.q1[which(tdata.q3[,idx.Z1]",sign.Z1[1],"cutoffs[idx.Z1] ",
                        "& tdata.q3[,idx.Z2]",sign.Z2[1],"cutoffs[idx.Z2]),])",sep="")))
  # q4.prodcells.num
  eval(parse(text=paste("q4.prodcells.num = nrow(tdata.q1[which(tdata.q4[,idx.Z1]",sign.Z1[1],"cutoffs[idx.Z1] ",
                        "& tdata.q4[,idx.Z2]",sign.Z2[1],"cutoffs[idx.Z2]),])",sep="")))
  # q1.prodcells.num = nrow(tdata.q1[which(tdata.q1[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q1[,idx.Z2]<=cutoffs[idx.Z2]),])
  # q2.prodcells.num = nrow(tdata.q2[which(tdata.q2[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q2[,idx.Z2]<=cutoffs[idx.Z2]),])
  # q3.prodcells.num = nrow(tdata.q3[which(tdata.q3[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q3[,idx.Z2]<=cutoffs[idx.Z2]),])
  # q4.prodcells.num = nrow(tdata.q4[which(tdata.q4[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q4[,idx.Z2]<=cutoffs[idx.Z2]),])
  q1.prodcells = round( (100 * q1.prodcells.num / nrow(tdata.q1) ), 2)
  q2.prodcells = round( (100 * q2.prodcells.num / nrow(tdata.q2) ), 2)
  q3.prodcells = round( (100 * q3.prodcells.num / nrow(tdata.q3) ), 2)
  q4.prodcells = round( (100 * q4.prodcells.num / nrow(tdata.q4) ), 2)
  
  #q1.prodcells.total
  eval(parse(text=paste("q1.prodcells.total = round( 100 * nrow(tdata.q1[which(tdata.q1[,idx.Z1]",sign.Z1[1],"cutoffs[idx.Z1] ",
                        "& tdata.q1[,idx.Z2]",sign.Z1[1],"cutoffs[idx.Z2]),]) / ncells.total,2)",sep="")))
  #q2.prodcells.total
  eval(parse(text=paste("q2.prodcells.total = round( 100 * nrow(tdata.q2[which(tdata.q2[,idx.Z1]",sign.Z1[1],"cutoffs[idx.Z1] ",
                        "& tdata.q2[,idx.Z2]",sign.Z1[1],"cutoffs[idx.Z2]),]) / ncells.total,2)",sep="")))
  #q3.prodcells.total
  eval(parse(text=paste("q3.prodcells.total = round( 100 * nrow(tdata.q3[which(tdata.q3[,idx.Z1]",sign.Z1[1],"cutoffs[idx.Z1] ",
                        "& tdata.q3[,idx.Z2]",sign.Z1[1],"cutoffs[idx.Z2]),]) / ncells.total,2)",sep="")))
  #q4.prodcells.total
  eval(parse(text=paste("q4.prodcells.total = round( 100 * nrow(tdata.q4[which(tdata.q4[,idx.Z1]",sign.Z1[1],"cutoffs[idx.Z1] ",
                        "& tdata.q4[,idx.Z2]",sign.Z1[1],"cutoffs[idx.Z2]),]) / ncells.total,2)",sep="")))
  # q1.prodcells.total = round( 100 * nrow(tdata.q1[which(tdata.q1[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q1[,idx.Z2]>=cutoffs[idx.Z2]),]) / ncells.total,2)
  # q2.prodcells.total = round( 100 * nrow(tdata.q2[which(tdata.q2[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q2[,idx.Z2]>=cutoffs[idx.Z2]),]) / ncells.total,2)
  # q3.prodcells.total = round( 100 * nrow(tdata.q3[which(tdata.q3[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q3[,idx.Z2]>=cutoffs[idx.Z2]),]) / ncells.total,2)
  # q4.prodcells.total = round( 100 * nrow(tdata.q4[which(tdata.q4[,idx.Z1]>=cutoffs[idx.Z1] & tdata.q4[,idx.Z2]>=cutoffs[idx.Z2]),]) / ncells.total,2)
  
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
  eval(parse(text=paste("my.calc = aggregate(data[,idx.Z1], 
                      by=list(fXY),
                      function(x) {
                        y = round( 100 * length(which(x", sign.Z1[1]," cutoffs[idx.Z1])) / length(x))
                        return(y)
                      })",sep="")))
  # bin color factor Z1
  my.calc.fac.Z1 = cut(my.calc$x,breaks=seq(0,100,by=10),labels=1:10,include.lowest=TRUE)
  names(my.calc.fac.Z1) = my.calc$x
  
  ### frequency of feature Z2        
  eval(parse(text=paste("freq.Z2 = aggregate(data[,idx.Z2],
                      by=list(fXY),
                      function(x) {
                        y = round( 100 * length(which(x", sign.Z2[1]," cutoffs[idx.Z2])) / length(x))
                        return(y)
                      })",sep="")))
  # bin color factor Z2
  my.calc.fac.Z2 = cut(freq.Z2$x,breaks=seq(0,100,by=10),labels=1:10,include.lowest=TRUE)
  names(my.calc.fac.Z2) = freq.Z2$x
  
  ### get only data table where Z1 and Z2 is produced
  eval(parse(text=paste("data.double = data[ which( (data[,idx.Z1]", sign.Z1[1]," cutoffs[idx.Z1]) ",
                        "& (data[,idx.Z2]", sign.Z2[1]," cutoffs[idx.Z2]) ),]",sep="")))
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
       ,y = par()$usr[3] - 1.0
       ,label = sprintf("%s(%s)",feat.X,cutoffs[idx.X])
       ,xpd = TRUE
       ,cex = 1.2
  )
  ### y axis label
  text(x = par()$usr[1] - 1.0
       ,y = 0.5*(par()$usr[3]+par()$usr[4])
       ,label = sprintf("%s(%s)",feat.Y,cutoffs[idx.Y])
       ,xpd = TRUE
       ,srt=90
       ,cex = 1.2
  )
  ######## from bottom to top
  ### legend title Z1
  text(x = start.legend + 0.8*binsize
       ,y = par()$usr[3]+0.08*(par()$usr[3]+par()$usr[4])-3*binsize
       ,label=feat.label.Z1,cex=0.75,pos=2
  )
  ### legend title Z2
  text(x = start.legend + 0.8*binsize
       ,y = par()$usr[3]+0.085*(par()$usr[3]+par()$usr[4])-1*binsize
       ,label=feat.label.Z2,cex=0.75,pos=2
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

######## CREATE TRIPLOT TABLE this$tab3D w/o PLOTTING ##############
##### @param
# file.idx.list     file index
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
  project.idx, 
  file.idx, 
  feat.X, 
  feat.Y, 
  feat.Z1, 
  calc, 
  col1, 
  col2, 
  binsize, 
  mincells, 
  plot.range=c(2,12,2,12)) 
{ 
    this=fcs 

    if (FALSE) {
        project.idx=1
        file.idx = 2
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
    data = this$loaddata(project.idx,file.idx)
    cutoffs = as.numeric(this$current.staintable[which(this$current.staintable$file_ID==file.idx),5])

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
    for (x in rownames(tab)) {
        for (y in colnames(tab)) {
            if ( tab[x,y]>=mincells ) {
                brackets.idx.x = brackets.idx.y = 1
                if (x==0) brackets.idx.x = 2
                if (y==0) brackets.idx.y = 2
                      
                fact = as.factor(paste(brackets.open[brackets.idx.x],x,',',as.numeric(x)+binsize,'] ',
                            brackets.open[brackets.idx.y],y,',',as.numeric(y)+binsize,']',sep=''))

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


####################################################################




if (FALSE) {
  ### light to dark
  cols.green = c("#F4FBEF","#E0F3D3","#C3E1AD","#ACD88D","#93CC6B","#7DC24C","#71B441","#63AB30","#56972B","#4F8726")
  cols.green = col2rgb(cols.green)
  
  cols.orange = c("#FFF1DC","#FEE2B9","#FED69A","#FDC168","#FCB140","#FBA018","#F29202","#E78C04","#DB8400","#C97A03")
  cols.orange = col2rgb(cols.orange)
  
  cols.greenorange = round((cols.green+cols.orange)/2)
  
  cols.blue = c("#F6F9FE","#ECF2FD","#D5E2FC","#CBDCFB","#BDD2F9","#B3CCF9","#A9C5F9","#9EBCF5","#88AFF2","#709DF0")
  cols.blue = col2rgb(cols.blue)
  
  col2rgb(heat.colors(10))
  
  col.Z1 = cols.green[,my.calc.fac.Z1]
  col.Z2 = cols.orange[,my.calc.fac.Z2]
  # mix colors
  col.mix = round((col.Z1+col.Z2)/2)
  
  cols.lila = c("#4f2a4f","#60335f","#713b70","#824482","#934c93","#a554a4","#b062af","#b973b9","#c283c2","#cb94ca")
  cols.lila = cols.lila[10:1]
  cols.lila = col2rgb(cols.lila)
  
  cols.turquoise = c("#2a444e","#33535f","#3b6170","#437081","#4c7f93","#548da4","#629ab0","#72a5b9","#82b0c2","#93bbca")
  cols.turquoise = cols.turquoise[10:1]
  cols.turquoise = col2rgb(cols.turquoise)
  
  
  cols.green2 = c("#284a3a","#305b47","#396c54","#417d61","#498e6e","#51a07b","#5dad88","#6eb795","#7ec0a1","#8fc8ae")
  cols.green2 = cols.green2[10:1]
  cols.green2 = col2rgb(cols.green2)
  
  cols.red = c("#dc292a","#5e3232","#6f3a3b","#804343","#924b4b","#a35353","#af6061","#b87171","#c18182","#ca9292")
  cols.red = cols.red[10:1]
  cols.red = col2rgb(cols.red)
  
  cols.blue2 = c("#2a2f4e","#33395f","#3b4270","#434c81","#4c5593","#545fa4","#626cb0","#727cb9","#828bc2","#939bca")
  cols.blue2 = cols.blue2[10:1]
  cols.blue2 = col2rgb(cols.blue2)
  
  cols.pink = c("#4d2935","#5e3241","#6f3a4c","#804357","#924b63","#a3536e","#af607b","#b87189","#c18197","#ca92a5")
  cols.pink = cols.pink[10:1]
  cols.pink = col2rgb(cols.pink)
  
}




