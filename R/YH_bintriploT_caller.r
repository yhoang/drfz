#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2018

### load library
library(RSQLite)
### delete everything first
# and then create new environment
#dbDisconnect(fcs$conn)
rm(list = ls())
fcs = new.env()


############################################################
############################################################
#################### YOUR INPUT NEEDED #####################

##################
###### DATABASE
### set path where database is located
fcs$db.path=file.path("","scratch","db")
fcs$db.path=file.path("","data","databases")
### set database name
fcs$db.name="SG_20171220_ThSubsets-NZBxW.sqlite3"
#fcs$db.name="SG_20170501_NZBxWCD4Subsets.sqlite3"
#fcs$db.name="SG_20170529_NZBxWIFNg+PD1+Signaling.sqlite3"
#fcs$db.name="SG_20170915_NZBxWTCellSubsets.sqlite3"
fcs$db.name="SG_20170915_NZBxW-TSubsets.sqlite3"
fcs$db.name="SG_20170627_NZBxWTfh-Th1.sqlite3"
#fcs$db.name="SG_20170704_NZBxW-Tfh.sqlite3"

##################
###### DATA SELECTION
### choose your project index and file index
# you can look up projects with
# fcs$project.list
# you can look up files with
# fcs$file.list
# see part 'LOAD ALL NECESSARY TABLES

##################
####### PARAMETERS
### calculation method
calc = "freq"
### bin size
binsize = 0.2
### minimum cell count in a bin
mincells = 10


############################################################
##################### LOAD FUNCTIONS #######################
source("YH_bintriploT_functions.r")
############################################################

############################################################
### CONNECT TO DATABASE
fcs$conn = dbConnect(SQLite(), dbname = file.path(fcs$db.path,fcs$db.name))
############################################################


############################################################
############# LOAD ALL NECESSARY TABLES ####################
#fcs$loaddata()
########## after running fcs$loaddata() you have ###########
# fcs$project.list              list all projects in database
# fcs$current.project           project chosen with index project.idx
# fcs$file.list                 list all files in one project
# fcs$data                      this is your data
# fcs$current.vars              column names of your data fcs$data
# fcs$current.staintable        lists all your column names in project
# fcs$current.staintable[,5]    lists your cutoffs in database according to your table and file index
############################################################

### call this if(FALSE) part manually first
if(FALSE) {
#if(TRUE) {
source("YH_bintriploT_functions.r")

fcs$project.idx.list = c(1,1,1,1)
fcs$file.idx.list = c(2,6,10,13)
fcs$bintriploT_construct(
    project.idx.list = fcs$project.idx.list
    ,file.idx.list = fcs$file.idx.list
    ,feat.X = "PD-1"
    ,feat.Y = "IFNg"
    ,binsize=0.2
    ,mincells=20
)
fcs$project.idx = 1
fcs$file.idx = 16
fcs$feat.Z1 = "CXCR5"
fcs$feat.Z2 = "Tbet"
fcs$col2="rainbow"

source("YH_bintriploT_functions.r")

fcs$bintriploT_freq_doublepos(
    project.idx = fcs$project.idx
    ,file.idx = fcs$file.idx
    ,feat.X = "PD-1"
    ,feat.Y = "IFNg"
    ,feat.Z1 = fcs$feat.Z1
    ,feat.Z2 = fcs$feat.Z2
    ,binsize = 0.2
    ,mincells = 20
    ,col = fcs$col2
)
dev.copy(png,sprintf("%s_%s_%s_double.png",
    fcs$shortenFilename(fcs$file.list[fcs$file.idx]),fcs$feat.Z1,fcs$feat.Z2))
dev.off()





### colors
# orange
# green
# blue
# lila
# turquoise
# green2
# red
# red2
# blue2
# pink
source("YH_bintriploT_functions.r")
fcs$col1="blue"
fcs$col2="red2"
fcs$project.idx = 1
fcs$file.idx = 16
fcs$feat.Z1 = "CXCR5"
fcs$feat.Z2 = "Tbet"
fcs$bintriploT_double_onefile(
    project.idx = fcs$project.idx
    ,file.idx = fcs$file.idx
    ,feat.X = "PD-1"
    ,feat.Y = "IFNg"
    ,feat.Z1 = fcs$feat.Z1
    ,feat.Z2 = fcs$feat.Z2
    ,calc = "freq"
    ,binsize = 0.2
    ,mincells = 20
    ,col1 = fcs$col2
    ,col2 = fcs$col1
)
dev.copy(png,sprintf("%s_%s_%s_%s%s.png",fcs$shortenFilename(fcs$file.list[fcs$file.idx]),fcs$feat.Z1,fcs$feat.Z2,fcs$col1,fcs$col2))
dev.off()




#####################################
fcs$db.name="SG_20171220_ThSubsets-NZBxW.sqlite3"
fcs$project.idx.list = c(3)
fcs$file.idx.list = c(1)
fcs$bintriploT_construct(
    project.idx.list = fcs$project.idx.list
    ,file.idx.list = fcs$file.idx.list
    ,feat.X = "PD-1"
    ,feat.Y = "IFNg"
    ,binsize = 0.2
    ,mincells = 20
)

fcs$project.idx = 3
fcs$file.idx = 1
fcs$feat.Z1 = "CXCR5"
fcs$feat.Z2 = "CXCR3"
fcs$col2="rainbow"

source("YH_bintriploT_functions.r")

fcs$bintriploT_freq_doublepos(
    project.idx = fcs$project.idx
    ,file.idx = fcs$file.idx
    ,feat.X = "PD-1"
    ,feat.Y = "IFNg"
    ,feat.Z1 = fcs$feat.Z1
    ,feat.Z2 = fcs$feat.Z2
    ,binsize = 0.2
    ,mincells = 20
    ,col = fcs$col2
)

dev.copy(png,sprintf("%s_%s_%s_double.png",
    fcs$shortenFilename(fcs$file.list[fcs$file.idx]),fcs$feat.Z1,fcs$feat.Z2))
dev.off()



library(plotly)
x=seq(0,12,by=0.2)
y=x
plot_ly() %>% add_surface(x=x,y=y,z=this$tab)


}





#######################################
############## PRI manuscript Fig 7 E
# 627_4hP+I_011_CD44_IL21_CXCR5_double80_green.pdf
# 627_4hP+I_011_CD44_IL21_Bcl6_double80_green.pdf
# 627_4hP+I_011_CD44_ICOS_IL21_double80_green.pdf
# 627_4hP+I_011_CD44_ICOS_IL21_double50_green.pdf
# 627_4hP+I_011_CD44_CXCR5_IL21_double50_green.pdf
# 627_4hP+I_011_CD44_Bcl6_IL21_double50_green.pdf
# 627_4hP+I_011_CD44_Bcl6_CXCR5_double80_green.pdf
rm(list = ls())
fcs = new.env()
source("YH_bintriploT_functions.r")

####### PARAMETERS
### calculation method
calc = "freq"
### bin size
binsize = 0.2
### minimum cell count in a bin
mincells = 10

fcs$db.path=file.path("","data","databases")
fcs$db.name="SG_20170627_NZBxWTfh-Th1.sqlite3"

############################################################
### CONNECT TO DATABASE
fcs$conn = dbConnect(SQLite(), dbname = file.path(fcs$db.path,fcs$db.name))



fcs$project.idx.list = 1
fcs$file.idx.list = 2
fcs$bintriploT_construct(
    project.idx.list = fcs$project.idx.list
    ,file.idx.list = fcs$file.idx.list
    ,feat.X = "PD-1"
    ,feat.Y = "IFNg"
    ,binsize = 0.2
    ,mincells = 10
)

fcs$project.idx = 1
fcs$file.idx = 2
fcs$feat.Z1 = "cxcr5"
fcs$feat.Z2 = "bcl6"
source("YH_bintriploT_functions.r")

color="green"
maxfreq = 80
fcs$bintriploT_freq_doublepos(
    project.idx = fcs$project.idx
    ,file.idx = fcs$file.idx
    ,feat.X = "PD-1"
    ,feat.Y = "IFNg"
    ,feat.Z1 = fcs$feat.Z1
    ,feat.Z2 = fcs$feat.Z2
    ,binsize = 0.2
    ,mincells = 10
    ,col = color
    ,maxfreq = maxfreq
)

dev.copy(pdf,sprintf("%s_%s_%s_double%s_%s.pdf",
    fcs$shortenFilename(fcs$file.list[fcs$file.idx]),fcs$feat.Z1,fcs$feat.Z2,maxfreq,color))
dev.off()


##################################################
source("YH_bintriploT_functions.r")
fcs$col1="blue"
fcs$col2="red2"
fcs$project.idx = 1
fcs$file.idx = 2
fcs$feat.Z1 = "IL21"
fcs$feat.Z2 = "Bcl6"
fcs$calc = "MFI"
fcs$bintriploT_double_onefile(
    project.idx = fcs$project.idx
    ,file.idx = fcs$file.idx
    ,feat.X = "PD-1"
    ,feat.Y = "IFNg"
    ,feat.Z1 = fcs$feat.Z1
    ,feat.Z2 = fcs$feat.Z2
    ,calc = fcs$calc
    ,binsize = 0.2
    ,mincells = 10
    ,col1 = fcs$col2
    ,col2 = fcs$col1
)

x=seq(2,12,by=0.2)
y=seq(2,12,by=0.2)
plot_ly() %>% add_surface(x=x,y=y,z=t(this$tab3D.Z1),reversescale=T)
plot_ly() %>% add_surface(x=x,y=y,z=t(this$tab3D.Z2),reversescale=T)


##################################################
source("YH_bintriploT_functions.r")
fcs$project.idx = 1
fcs$file.idx = 2
fcs$feat.Z1 = "Bcl6"
fcs$feat.Z1 = "IL21"
fcs$calc = "freq"
fcs$bintriploT_table(
    project.idx = fcs$project.idx
    ,file.idx = fcs$file.idx
    ,feat.X = "PD-1"
    ,feat.Y = "IFNg"
    ,feat.Z1 = fcs$feat.Z1
    ,calc = fcs$calc
    ,binsize = 0.2
    ,mincells = 10
)

x=seq(2,12,by=0.2)
y=seq(2,12,by=0.2)
scene = list(xaxis=list(title= "PD-1")
    ,yaxis=list(title="IFN-g")
    ,zaxis=list(title="freq",range=c(0,77))
    #,camera=list(
    #    up=list(x=0,y=0,z=1)
    #    ,center=list(x=0,y=0,z=0)
    #    ,eye=list(x=0.5,y=1.2,z=1.5))
    ,cameraposition=list(c(-0.1, 0.4, -1, 0.33), c(0.0, 0, 0.0), 2)
)


# heat
plot_ly(showscale = TRUE) %>% add_surface(x=x,y=y,z=t(fcs$tab3D), cmin=0, cmax=82,colors="YlOrRd",reversescale=F) %>% layout(scene=scene)

#spectral
plot_ly(showscale = TRUE) %>% add_surface(x=x,y=y,z=t(fcs$tab3D), cmin=0, cmax=82,colors="Spectral",reversescale=T) %>% layout(scene=scene)

#ygn
plot_ly(showscale = TRUE) %>% add_surface(x=x,y=y,z=t(fcs$tab3D), cmin=0, cmax=82,colors="YlGnBu",reversescale=F) %>% layout(scene=scene)

#yg
plot_ly(showscale = TRUE) %>% add_surface(x=x,y=y,z=t(fcs$tab3D), cmin=0, cmax=82,colors="YlGn",reversescale=F) %>% layout(scene=scene)

#green
plot_ly(showscale = TRUE) %>% add_surface(x=x,y=y,z=t(fcs$tab3D), cmin=0, cmax=82,colors="Greens",reversescale=F) %>% layout(scene=scene)

plotly_IMAGE(p,format="pdf",outfile=sprintf("%s_greens.pdf",fcs$feat.Z1))

plot_ly(showscale = TRUE) %>% add_surface(x=x,y=y,z=t(fcs$tab3D), cmin=0, cmax=82,colorscale=list(c(0,"rgb(255, 0, 0)"), list(1, "rgb(0, 255, 0)"))) %>% layout(scene=scene)




############## PRI manuscript Fig 7 G
fcs$db.name="SG_20170627_NZBxWTfh-Th1.sqlite3"
#fcs$db.name="SG_20170704_NZBxW-Tfh.sqlite3"
fcs$conn = dbConnect(SQLite(), dbname = file.path(fcs$db.path,fcs$db.name))

source("YH_bintriploT_functions.r")
fcs$project.idx = 1
fcs$file.idx = 2
fcs$feat.Z1 = "Bcl6"
#fcs$feat.Z1 = "IL21"
fcs$calc = "freq"
fcs$bintriploT_table(
    project.idx = fcs$project.idx
    ,file.idx = fcs$file.idx
    ,feat.X = "PD-1"
    ,feat.Y = "IFNg"
    ,feat.Z1 = fcs$feat.Z1
    ,calc = fcs$calc
    ,binsize = 0.2
    ,mincells = 10
)
write.csv(fcs$tab3D,file=sprintf("%s_3d.csv",fcs$feat.Z1),row.names=rownames(fcs$tab3D))

x=seq(2,12,by=0.2)
y=seq(2,12,by=0.2)
scene = list(xaxis=list(title= "PD-1")
    ,yaxis=list(title="IFN-g")
    ,zaxis=list(title="IL-21",range=c(0,77))
    #,camera=list(
    #    up=list(x=0,y=0,z=1)
    #    ,center=list(x=0,y=0,z=0)
    #    ,eye=list(x=0.5,y=1.2,z=1.5))
    ,cameraposition=list(c(-0.1, 0.4, -1, 0.33), c(0.0, 0, 0.0), 2)
)
lighting = list(
    ambient = 0.6
    ,diffuse = 0.3
    ,roughness = 0.9
    ,specular = 0.3
    ,fresnel = 0.2
)


#green
plot_ly(showscale = TRUE) %>% add_surface(x=x,y=y,z=t(fcs$tab3D), cmin=0, cmax=82,colors="Greens",reversescale=F, lighting=lighting) %>% layout(scene=scene)

plot_ly(showscale = TRUE) %>% add_surface(x=x,y=y,z=t(fcs$tab3D), cmin=0, cmax=82,colors="Greens",reversescale=F) %>% layout(scene=scene)
##########################


############# FOR PRI-Manuscript Fig5
#X: TNFa
#Y: IFNg
#Z: IL2 + IL21.
#Die Datenbank nennt sich SGJP_20180524_NZBxWIL21+Cytokines. Gibt es ein Skript, mit dem ich die erstellen kann, 
#oder ist es besser, wenn Du die Abbildungen machst?

source("YH_bintriploT_functions.r")
fcs$db.name = "SGJP_20180524_NZBxWIL21+Cytokines.sqlite3"
### CONNECT TO DATABASE
fcs$conn = dbConnect(SQLite(), dbname = file.path(fcs$db.path,fcs$db.name))


fcs$feat.X = "TNFA"
fcs$feat.Y = "Ifng"
fcs$feat.Z1 = "il2"
fcs$feat.Z2 = "il21"

fcs$project.idx.list = 1
fcs$file.idx.list = 2
fcs$bintriploT_construct(
  project.idx.list = fcs$project.idx.list
  ,file.idx.list = fcs$file.idx.list
  ,feat.X = fcs$feat.X
  ,feat.Y = fcs$feat.Y
  ,binsize = 0.2
  ,mincells = 5
)

fcs$project.idx = 1
fcs$file.idx = 2
source("YH_bintriploT_functions.r")

color="green"
maxfreq = 100
fcs$bintriploT_freq_doublepos(
  project.idx = fcs$project.idx
  ,file.idx = fcs$file.idx
  ,feat.X = fcs$feat.X
  ,feat.Y = fcs$feat.Y
  ,feat.Z1 = fcs$feat.Z1
  ,feat.Z2 = fcs$feat.Z2
  ,binsize = 0.2
  ,mincells = 5
  ,col = color
  ,maxfreq = maxfreq
)

dev.copy(pdf,sprintf("%s_%s_%s_double%s_%s.pdf",
                     fcs$shortenFilename(fcs$file.list[fcs$file.idx]),fcs$feat.Z1,fcs$feat.Z2,maxfreq,color))
dev.off()
#########################################################








### for Andreas Radbruch
{
fcs$db.name="SG_20170627_NZBxWTfh-Th1.sqlite3"
source("YH_bintriploT_functions.r")
fcs$project.idx = 1
fcs$file.idx = 2
fcs$feat.Z1 = "Bcl6"
fcs$feat.Z1 = "CXCR5"
fcs$calc = "freq"
fcs$bintriploT_table(
    project.idx = fcs$project.idx
    ,file.idx = fcs$file.idx
    ,feat.X = "PD-1"
    ,feat.Y = "IL21"
    ,feat.Z1 = fcs$feat.Z1
    ,calc = fcs$calc
    ,binsize = 0.2
    ,mincells = 10
)

x=seq(2,12,by=0.2)
y=seq(2,12,by=0.2)
scene = list(xaxis=list(title= "PD-1")
    ,yaxis=list(title="IL21")
    ,zaxis=list(title="Freq(CXCR5)",range=c(0,70))
    #,camera=list(
    #    up=list(x=0,y=0,z=1)
    #    ,center=list(x=0,y=0,z=0)
    #    ,eye=list(x=0.5,y=1.2,z=1.5))
    ,cameraposition=list(c(-0.1, 0.4, -1, 0.33), c(0.0, 0, 0.0), 2)
)

#green
plot_ly(showscale = TRUE) %>% add_surface(x=x,y=y,z=t(fcs$tab3D), cmin=0, cmax=82,colors="Greens",reversescale=F) %>% layout(scene=scene)



fcs$db.name="SG_20170627_NZBxWTfh-Th1.sqlite3"
#fcs$db.name="SG_20170704_NZBxW-Tfh.sqlite3"
fcs$conn = dbConnect(SQLite(), dbname = file.path(fcs$db.path,fcs$db.name))
fcs$project.idx.list = 1
fcs$file.idx.list = 2
fcs$bintriploT_construct(
    project.idx.list = fcs$project.idx.list
    ,file.idx.list = fcs$file.idx.list
    ,feat.X = "PD-1"
    ,feat.Y = "IL21"
    ,binsize = 0.2
    ,mincells = 10
)

fcs$project.idx = 1
fcs$file.idx = 2
fcs$feat.Z1 = "CXCR5"
fcs$feat.Z2 = "Bcl6"
fcs$col2="red"
source("YH_bintriploT_functions.r")

color="green"
maxfreq = 70
fcs$bintriploT_freq_doublepos(
    project.idx = fcs$project.idx
    ,file.idx = fcs$file.idx
    ,feat.X = "PD-1"
    ,feat.Y = "IL21"
    ,feat.Z1 = fcs$feat.Z1
    ,feat.Z2 = fcs$feat.Z2
    ,binsize = 0.2
    ,mincells = 10
    ,col = color
    ,maxfreq = maxfreq
)
dev.copy(pdf,sprintf("%s_%s_%s_double%s_%s.pdf",
    fcs$shortenFilename(fcs$file.list[fcs$file.idx]),fcs$feat.Z1,fcs$feat.Z2,maxfreq,color))
dev.off()



}








































