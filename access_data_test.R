#!/usr/bin/R

rm(list = ls())
fcs = new.env()

library(RSQLite)

fcs$db.path=file.path("","data","databases")
fcs$db.name="EM_20200210_EME002_WTCD4.sqlite3"

##################
####### PARAMETERS
### calculation method
calc = "MFI"
### bin size
binsize = 0.2
### minimum cell count in a bin
mincells = 10
cofactor = 0.2

source("create_triploT_matrx.R")

### CONNECT TO DATABASE
fcs$conn = dbConnect(SQLite(), dbname = file.path(fcs$db.path,fcs$db.name))

############################################################
############# LOAD ALL NECESSARY TABLES ####################
fcs$loaddata()
########## after running fcs$loaddata() you have ###########
# fcs$table.list                list all tables in database
# fcs$project.list              list all projects in database
# fcs$current.project           project chosen with index project.idx
# fcs$file.list                 list all files in one project
# fcs$data                      this is your data
# fcs$current.vars              column names of your data fcs$data
# fcs$current.staintable        lists all your column names in project
# fcs$current.staintable[,5]    lists your cutoffs in database according to your table and file index
############################################################

fcs$project.list

project.idx <- which(fcs$project.list=="Staining2")

head(fcs$loaddata(project.idx, file.idx = 1, cof = 0.2))

fcs$file.list
file.idx <- which(fcs$file.list=="CD4_Sap_STII.fcs")
head(fcs$loaddata(project.idx, file.idx = file.idx, cof = 0.2))

## look for markers of the file
fcs$current.vars

## look for cutoffs saved in database
fcs$current.staintable[,5]

data <- fcs$loaddata(project.idx, file.idx = file.idx, cof = 0.2)

fcs$bintriploT_table(
                     project.idx, 
                     file.idx, 
                     feat.X="CXCR5", 
                     feat.Y="CXCR3", 
                     feat.Z1="IL21", 
                     calc="MFI",
                     cofactor=0.2,
                     binsize=0.2, 
                     mincells=10, 
                     plot.range=c(0.1,10,0.1,10)
                     )


