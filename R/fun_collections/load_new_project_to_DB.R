#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019

##################################################################################################
######## LOAD NEW PROJECT TO DATABASE RB_20191002_Good2018.sqlite3                     ###########
######## WITH DIFFICULT DATA (variable col count or labeling)                          ###########
######## including four tables                                                         ###########
######## data table named [project.name]                                               ###########
######## [project.name]_fileIdentity table                                             ###########
######## [project.name]_markerIdentity table                                           ###########
######## [project.name]_equipmentInfo table                                            ###########
##################################################################################################

### create new environment
rm(list = ls())

### load libraries
# tools           :   Library for md5sum().
# dplyr           :   Fast combining rows/cols
# flowCore        :   read FCS files
# RSQLite         :   interact with database
# rChoiceDialogs  :   select files with GUI
libraries = c("flowCore","RSQLite","rChoiceDialogs","dplyr","tools")
lapply(libraries,require, character.only = TRUE)


########### START: INPUT - change if necessary
db.path=file.path("","data","databases")
db.name="RB_20191002_Good2018.sqlite3"

fcs.folder = file.path("","data","Cytometry","Good_nat_med_2018_BCP_ALL","fcs_for_PRI")
# project.name = "Basal"
project.name = "BCR"
# project.name = "IL7"
# project.name = "Pervanadate"
# project.name = "TSLP"

selectedColnames = as.vector(read.table(file.path(fcs.folder,"columns_sorted.txt"))[,1])
### standardized column names
colnames.stand = gsub("[-]", "_",toupper(selectedColnames))
########### END: INPUT 



### helpful print function
printf <- function(...) invisible(print(sprintf(...)))

DB = dbConnect(SQLite(), dbname = file.path(db.path,db.name))
print(dbListTables(DB))
# dbRemoveTable(DB,"Basal")


### SELECT FILES FOR THIS PROJECT
# Direct selection of files. Only fcs files are shown in opening window. 
files = jchoose.files(default = fcs.folder , filters = "fcs")

########### START: INPUT DATA TABLE "[project.name]"
### read intensities in one file "read.all"
cells.total = 0
for (i in 1:length(files)) {
  read = read.FCS(files[i])
  keys = keyword(read)
  cells.total = cells.total+as.numeric(keys$`$TOT`)
  printf("reading file #%s: %s with %s entries",i,basename(files)[i],keys$`$TOT`)
  
  colnames.bio = col.id = metal.vec = vector()
  for (j in 1:100) {
    name.in.fcs = eval(parse(text=paste0("keys$`$P",j,"S`")))
    metal = eval(parse(text=paste0("keys$`$P",j,"N`")))
    if (length(name.in.fcs)==0) break;
    
    ### START: INPUT - change if necessary
    if (name.in.fcs=="HLADR") {name.in.fcs="HLA_DR"}#;printf("changed! %s",name.in.fcs)}
    if (name.in.fcs=="FITC_myeloid") {name.in.fcs="CD33_CD16"}#;printf("changed! %s",name.in.fcs)}
    if (name.in.fcs=="CD33-CD16") {name.in.fcs="CD33_CD16"}#;printf("changed! %s",name.in.fcs)}
    if (name.in.fcs=="Gd157" || name.in.fcs=="(Gd157)Di") {name.in.fcs="CD132"}#;printf("changed! %s",name.in.fcs)}
    if (name.in.fcs=="Kappa_lambda") {name.in.fcs="kappa_lambda"}
    if (name.in.fcs=="pAkt") {name.in.fcs="pAKT"}
    if (name.in.fcs=="Pax5") {name.in.fcs="PAX5"}
    if (name.in.fcs=="pCreb") {name.in.fcs="pCREB"}
    if (name.in.fcs=="pErk") {name.in.fcs="pERK"}
    if (name.in.fcs=="pIkaros") {name.in.fcs="pIKAROS"}
    if (name.in.fcs=="tIkaros") {name.in.fcs="tIKAROS"}
    if (name.in.fcs=="TSLPR") {name.in.fcs="TSLPr"}
    name.in.fcs = gsub("[-]", "_",name.in.fcs)
    tmp.names = toupper(name.in.fcs)
    ### END: INPUT - change if necessary
    
    
    if (tmp.names %in% colnames.stand) {
      # save column index in right order
      col.id= c(col.id,j)
      # save column names in right order
      colnames.bio= c(colnames.bio,name.in.fcs);
      # save metal names in right order
      metal.vec = c(metal.vec,metal)
      
    }
  }
  if ( length(col.id)<39 ) {printf("Something wrong with file #%s",i);break;}
  
  ### sort order
  new.order = order(colnames.bio)
  col.id = col.id[new.order]
  colnames.bio = colnames.bio[new.order]
  metal.vec = metal.vec[new.order]
  
  ### read exprs data
  read.in = as.data.frame(exprs(read))[col.id]
  colnames(read.in) = colnames.bio
  
  ### Generate empty first column and fill it with  file ID.
  file.id = as.data.frame(rep(i, nrow(read.in)))
  names(file.id) ="file_ID"
  
  # if first file
  if (i==1) {
    read.all = bind_cols(file.id, read.in)
    colnames.file = metal.vec
  } else {
    read.in = bind_cols(file.id, read.in)
    read.all = bind_rows(read.all,read.in)
    colnames.file = c(colnames.file,metal.vec)
  }
}
dim(read.all)
#write.table(as.data.frame(colnames.bio),file="/data/Cytometry/Good_nat_med_2018_BCP_ALL/fcs_for_PRI/colums_sorted.txt",sep="\t",
#quote=FALSE,row.names=FALSE,col.names=FALSE)

### WRITE IN DATABASE IF PROJECT DOES NOT EXIST YET
if( !dbExistsTable(conn = DB, project.name)) {
  # Generates standardized colnames accordingly to rows of default table with file 
  # colnames and shortnames.
  colnames.DB = NULL
  for (l in 1:(ncol(read.in)-1)) {
    colnames.DB = c(colnames.DB, paste("col", sprintf("%02i", l), sep = ""))
  }
  names(colnames.DB) = "database_colnames"
  # Changes colnames of current file. (Only happens once -> first time the table
  # is created within database. During the following runs it only gets extended.)
  colnames(read.all) = c("file_ID", colnames.DB)
  dbWriteTable(conn = DB, project.name, read.all)#,overwrite=TRUE)
  
  print(paste("File (index ", length(files), ") ", i, " of ", length(files), " has been added to table ", 
        project.name, " of the database.", sep = ""))
} else {
  print("PROJECT ALREADY EXISTS! END.")
}

### INDEXING IN DATABASE FOR FASTER READING
# After complete insertion the column 'file_ID' is indexed. This quickens the access for later use but it takes a while here.
sql.query = paste0("create index ", project.name, "_idxfile on ", project.name, " (file_ID)")
dbGetQuery(DB, sql.query)

### CHECK TABLE
dbGetQuery(DB, paste0("select * from ",project.name," limit 2"))
########### END: INPUT DATA TABLE "[project.name]"


######### START: INPUT TABLE "[project.name]_equipmentInfo"
######### Include keywords of fcs file.
keywords = unlist(keys)
# Change names of keywords.
names(keywords) = gsub("^\\$","", names(keywords))
# Compiles dataframe with current file_ID and keywords. Convert factors to values within table.
for ( i in 1:length(files)) {
  keyTmp = data.frame(file_ID = rep(i, length(keywords)), 
                      Keywords = names(keywords), 
                      Values = unname(keywords))
  keyTmp[] = lapply(keyTmp, as.character)
  if (i==1) {
    keyTable = keyTmp
  } else {
    keyTable = bind_rows(keyTable,keyTmp)
  }
}
head(keyTable,2)
### WRITE TABLE IN DATABASE
dbWriteTable(conn = DB, paste0(project.name, "_equipmentInfo"), keyTable)
### CHECK TABLE
dbGetQuery(DB, paste0("select * from ",project.name,"_equipmentInfo limit 2"))
######### END: INPUT TABLE "[project.name]_equipmentInfo"



######### START: INPUT TABLE "[project.name]_markerIdentity"
# col01 = file_ID
# col02 = database_colnames
# col03 = file_colname
# col04 = shortname
###

# Generates standardized colnames accordingly to rows of default table with file 
# colnames and shortnames.
colnames.DB = NULL
for (l in 1:(ncol(read.in)-1)) {
  colnames.DB = c(colnames.DB, paste("col", sprintf("%02i", l), sep = ""))
}
names(colnames.DB) = "database_colnames"

for ( i in 1:length(files)) {
  colNames = bind_cols(
    data.frame(file_ID=rep(i,length(colnames.DB))),
    data.frame(database_colnames=colnames.DB),
    data.frame(shortname=colnames.bio)
    )
  if (i==1) tmpNames = colNames
  else tmpNames = bind_rows(tmpNames,colNames)
}

markerTable = bind_cols(tmpNames,data.frame(file_colname=colnames.file))
markerTable = markerTable[,c(1,2,4,3)]
head(markerTable,2)

### WRITE TABLE IN DATABASE
dbWriteTable(conn = DB, paste0(project.name, "_markerIdentity"), markerTable)
### CHECK TABLE
dbGetQuery(DB, paste0("select * from ",project.name,"_markerIdentity limit 2"))
######### END: INPUT TABLE "[project.name]_markerIdentity"

  


######### START: INPUT TABLE "[project.name]_fileIdentity"
# col01 = file_ID
# col02 = filename
# col03 = Date_of_Measurement
# col04 = Staining_Information
# col05 = md5_hash
###

# Compiles fileID accordingly to loaded files.
file_ID = 1:length(files)
# Basename extracts the filename of a path. It goes through all paths within an argument.
filename = basename(files)
# col03
Date = rep(NA, length(files))
# StainID given by user. Not needed for further calculations just for user to see which fcs file
# has which stain. Provided by user during file upload (shiny app).
stain_ID = rep("none", length(files))
# Calculates 32-byte MD5 hashes for all loaded files.
md5_hash = as.vector(md5sum(files))
# Combines four variables in one data frame.
fileTable = data.frame(file_ID, filename, 
                       Date_of_Measurement=Date,
                       Staining_Information=stain_ID, md5_hash)

head(fileTable,2)
### WRITE TABLE IN DATABASE
dbWriteTable(DB, paste0(project.name,"_fileIdentity"), fileTable)
### CHECK TABLE
dbGetQuery(DB, paste0("select * from ",project.name,"_fileIdentity limit 2"))

### Include date of computation in fileIdentity table. NA is replaced by date of 
### current fcs file.
# sql.query= paste("update ", project.name, "_fileIdentity set Date_of_Measurement = '", 
#                  keyTable[grep("^DATE", keyTable[, "Keywords"]), "Values"][1], 
#                  "' where file_ID = ", files, sep = "")
# dbGetQuery(DB, sql.query)
######### END: INPUT TABLE "[project.name]_fileIdentity"
  

##### FINAL CHECK
total.projects = dbListTables(DB)
project.idx = grep(project.name,total.projects)
print(total.projects[project.idx])
colnames(dbGetQuery(DB, paste0("select * from ",project.name," limit 2")))

###### CLOSE DB
dbDisconnect(DB)

