######## CREATE TRIPLOT TABLE this$tab3D w/o PLOTTING ##############
##### @param
# file.idx.list     file index
# feat.X            name of feature X
# feat.Y            name of feature Y
# feat.Z1           name of feature Z1
# calc              calculation method
# cofactor          asinh cofactor asinh(x/cof), at Good it was 0.2
# binsize           bin size, default = 0.2
# mincells          minimum amount of cells in a bin, default = 10
# plot.range        plot range, first x axis, second y axis, default= c(0.1,10,0.1,10)
fcs$bintriploT_table <- function(
  project.idx, 
  file.idx, 
  feat.X, 
  feat.Y, 
  feat.Z1, 
  calc,
  cofactor=0.2,
  binsize=0.2, 
  mincells=10, 
  plot.range=c(0.1,10,0.1,10)) 
{ 
  this=fcs 
  
  if (FALSE) {
    project.idx=1
    file.idx = 2
    feat.X = "PD-1"
    feat.Y = "IFNg"
    feat.Z1 = "IL21"
    calc = "MFI"
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
  data = this$loaddata(project.idx,file.idx,cof=cofactor)
  ### takes cutoff saved in database
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
  
  # this$tab = tab
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
  
  ########## INSERT INFORMATION INTO TRIPOT TABLE 'tab'
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
      # } else if ( tab[x,y]>0 & calc!="freq") {
      #   tab[x,y] = min.range.Z1
      } else {
        tab[x,y] = NA
      }
    }
  }
  
  # this$tab3D = round(tab,3)
  round(tab,3)
}



####################################################################
######## LOAD FUNCTIONS WITHOUT ACTUALLY CALLING THEM ##############
##### @param
# project.idx   project index
# file.idx      file index
# cof           asinh cofactor
##### returns data table
fcs$loaddata <- function(
  project.idx,
  file.idx,
  cof)  
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
  
  data = this$getFile(this$current.project,file.idx,cofactor=cof)
  
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


######## GET DATA TABLE ############################################
##### @param
# table         project name
# fileidx       file index
# stain         staining vector
# cofactor      cofactor to transform asinh(x/cofactor)
fcs$getFile <- function (
  table,
  fileidx,
  stain=NA,
  cofactor = 1) 
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
  
  data = asinh(data/cofactor) 
  
  if ( !is.na(stain) ) { 
    colnames(data) = stain 
  } else { 
    colnames(data) = this$current.vars 
  }
  
  print("Feature names:")
  print(colnames(data))
  
  cat("Done getFile\n\n")
  
  this$data = data
  
  data
}




printf <- function(...) invisible(print(sprintf(...)))
