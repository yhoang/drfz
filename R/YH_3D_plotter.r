#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019

# create new environment
rm(list = ls())
fcs = new.env()
param = new.env()
fcs$parent.env=ls()
fcs$version="v0.39b"


############# workspace
### "office": DRFZ Office
### "Ria"   : Ria Laptop
### "YH"    : YH Office



workspace="YH"
#workspace="rev"
# workspace="office"
#workspace="Ria"
#workspace="Praktika"



### load public libraries - DO NOT TOUCH!
if (workspace=="office") {
  library(RSQLite,quietly=TRUE,lib.loc = file.path("Y:","AG_Baumgrass","AG-PRI","R","R-3.5.2","library"))
  library(tcltk2,quietly=TRUE,lib.loc = file.path("Y:","AG_Baumgrass","AG-PRI","R","R-3.5.2","library"))
  library(R.devices,quietly=TRUE,lib.loc = file.path("Y:","AG_Baumgrass","AG-PRI","R","R-3.5.2","library"))
} else {
  library(RSQLite,quietly=TRUE)#,lib.loc = file.path("","home","yhoang","R","library"))
  library(tcltk2,quietly=TRUE)#,lib.loc = file.path("","home","yhoang","R","library"))
  library(R.devices,quietly=TRUE)#,lib.loc = file.path("","home","yhoang","R","library"))
#used only to detect number of cores on machine
#library(parallel,quietly=TRUE)
#package used for parallele computing
#library(snow,quietly=TRUE)
}


if (workspace=="office") {
	setwd(file.path("Y:","AG_Baumgrass","AG-PRI","PRIanalyzer"))
	fcs$db.path=file.path("Y:","AG_Baumgrass","AG-PRI","DB")
	fcs$db.name=""
	fcs$working=FALSE
	#fcs$working=TRUE

	### load private library
	source(file.path("Y:","AG_Baumgrass","AG-PRI","PRIanalyzer","tcl","libbwidget.r"))
	### set template location
	fcs$png.file=file.path("Y:","AG_Baumgrass","AG-PRI","PRIanalyzer","tcl","tmp.png")
	fcs$template.file=file.path("Y:","AG_Baumgrass","AG-PRI","PRIanalyzer","tcl","template.png")
} else if (workspace=="Ria") {
	setwd("C:/Users/ag-baumgras/Documents/R")
	fcs$db.path=file.path("data")
	fcs$db.name=""
	fcs$working=FALSE
	### load private library
	source(file.path("tcl","libbwidget.r"))
	### set template location
	fcs$png.file=file.path("tcl","tmp.png")
	fcs$template.file=file.path("tcl","template.png")
} else if (workspace=="YH") {
	setwd(file.path("","scratch","drfz","R"))
	fcs$db.path=file.path("","data","databases")
	fcs$db.name="SG_20161030_HumanCytokines_CD4mem.sqlite3"
	fcs$db.name="SG_20170517_NZBxW-TconPD1IFNgReceptors.sqlite3"
	fcs$db.name="SG_20170821_NZBxW-SubsetDistinction.sqlite3"
	#fcs$db.name="SG_20170821_NZBxW-TCellSubstets.sqlite3"
	fcs$db.name="SG_20170529_NZBxWIFNg+PD1+Signaling_old.sqlite3"
	fcs$db.name="SG_20170915_NZBxW-TSubsets.sqlite3"
	fcs$db.name="SG_20171220_ThSubsets-NZBxW.sqlite3"
	fcs$db.name="CK_20170123_PBMC_Melanoma.sqlite3"
	#fcs$db.name="AM_20180227_panel-redesgin_MS-project.sqlite3"
	fcs$db.name="CP_20171030_humanmelanoma.sqlite3"
	fcs$db.name="AM_20180227_panel-redesgin_MS-project.sqlite3"
	fcs$db.name="AM_20180111_MS-exp-4-CIS10-CIS11-CIS12-HC05-HC07-HC08_CD4+.sqlite3"
	fcs$db.name="CP_20171102_weidata3.sqlite3"
	fcs$db.name="YH_20170123_PBMC_Melanoma.sqlite3"
	fcs$db.name="YH_20180423_Myeloid_CD14CD16.sqlite3"
	#fcs$db.name="SGJP_20180524_NZBxWIL21+Cytokines.sqlite3"
	#fcs$db.name="SG_20170704_NZBxW-Tfh.sqlite3"
	#fcs$db.name="SG_20161030_HumanCytokines_CD4mem.sqlite3"
	#fcs$db.name="SG_20171220_ThSubsets-NZBxW.sqlite3"
	#fcs$db.name="SG_20170704_NZBxW-Tfh.sqlite3"
	#fcs$db.name="SG_20170501_NZBxWCD4Subsets.sqlite3"
	#fcs$db.name="SG_20170529_NZBxWIFNg+PD1+Signaling.sqlite3"
	fcs$db.name="AR_20180111_Chevrier-selection.sqlite3"
	#fcs$db.name="SG_20170627_NZBxWTfh-Th1.sqlite3"
	#fcs$db.name="AR_20171109_test_NZBxWTCellSubsets.sqlite3"
	#fcs$db.name="SG_20170915_NZBxWTCellSubsets.sqlite3"
	#fcs$db.name="AM_20170914_MS-01.sqlite3"
	#fcs$db.name="AM_20170928_MSIICIS4CIS5HC2.sqlite3"
	#fcs$db.name="AM_20170524_Cytokines_CD4_Yen.sqlite3"
	fcs$db.name="SG_20180612_NZBxWCytokines.sqlite3"
	fcs$db.name="AM_20180323_MS-exp-6-CIS17-CIS18-CIS20-CIS21-CIS22-CIS23-HC13-HC14-HC15-HC16-BF1-BF2_CD4+.sqlite3"
	fcs$db.name="SG_20181026_NZBxW_TCSubsets.sqlite3"
	fcs$db.name="TL_20160721_NZBxW-Kinetics.sqlite3"
	fcs$db.name="AR_20170817_Spitzer_bloodsamples_day3.sqlite3"
	#fcs$db.name="AM_20180312_MS-exp-5-CIS13-CIS14-CIS15-CIS16-HC09-HC10-HC11-HC12_CD4+.sqlite3"
	#fcs$db.name="AR_20170817_Spitzer_spleensamples.sqlite3"
	#fcs$db.name="AR_20170907_Spitzer_lymphnodesamples.sqlite3"
	fcs$db.name="SG_20170915_NZBxW-TSubsets.sqlite3"
	fcs$db.name="SG_20171220_ThSubsets-NZBxW.sqlite3"
	# fcs$db.name="SG_20180130_NZBxW-ThSubsets.sqlite3"
	# fcs$db.name="SG_20170704_NZBxW-Tfh.sqlite3"
	fcs$db.name="SG_20181012_NZBxWTfh-like.sqlite3"
	fcs$db.name="SG_20181026_NZBxW_TCSubsets.sqlite3"
	fcs$db.name="SG_20170627_NZBxWTfh-Th1.sqlite3"
	fcs$db.name="RB_20191002_Good2018.sqlite3"
	fcs$db.name="SG_20180413_NZBxW-IL21_CD44.sqlite3"
	fcs$db.name="RB_20191002_Good2018.sqlite3"
	# fcs$db.name="YH_20190524_Tordesillas2016.sqlite3"
	
	fcs$table.dir = "/scratch/drfz/Spitzer2017/glmnet_bio_keep/"
	fcs$working=TRUE
	#fcs$working=FALSE

	### load private library
	source(file.path("tcl","libbwidget.r"))
	### set template location
	fcs$png.file=file.path("tcl","tmp.png")
	fcs$template.file=file.path("tcl","template.png")
} else if (workspace=="rev") {
	setwd(file.path("","scratch","drfz","R"))
	fcs$db.path=file.path("","scratch","db")
	fcs$db.name="SG_20161030_HumanCytokines_CD4mem.sqlite3"
	#fcs$db.name="CP_20171030_humanmelanoma.sqlite3"
	fcs$working=TRUEdensity=TRUE

	### load private library
	source(file.path("tcl","libbwidget.r"))
	### set template location
	#fcs$png.file=file.path("tcl","tmp.png")
	#fcs$template.file=file.path("tcl","template.png")
}  else if (workspace=="Praktika") {
	setwd(file.path("Y:","AG_Baumgrass","Praktika2017_Grp1"))
	fcs$db.path=file.path("Y:","AG_Baumgrass","Praktika2017_Grp1","DB")
	fcs$db.name=""
	fcs$working=FALSE
	
	### load private library
	source(file.path("Y:","AG_Baumgrass","Praktika2017_Grp1","tcl","libbwidget.r"))
	### set template location
	fcs$png.file=file.path("Y:","AG_Baumgrass","Praktika2017_Grp1","tcl","tmp.png")
	fcs$template.file=file.path("Y:","AG_Baumgrass","Praktika2017_Grp1","tcl","template.png")
}

### global variables - DO NOT TOUCH!
options(warn=-1,scipen=999)
fcs$dataframes.name = vector()
fcs$dataframes1 = fcs$dataframes2 = as.data.frame(matrix(NA,1,1))
fcs$temptable.name = vector()
fcs$df.num = 2
fcs$plot.num = 0
fcs$plot.windows = vector()
fcs$plotter.di.graph = 0
fcs$plotter.di.num = 0
fcs$plotter.tri.num = 0
fcs$distepy = 0
fcs$origin.ncells = 0
fcs$ncell.sel = 0
fcs$ncell.perc = 0
fcs$cbVal = list()
fcs$cbcutoffperc=list()
fcs$saved.cutoffs = list()
fcs$saved.checks = list()
fcs$images = list()
fcs$plot.attr = list()

fcs$rbasinh = tclVar("1")

fcs$asinh = list(
  scale = c(asinh(-1000),asinh(-100),asinh(-10),asinh(0),asinh(10),asinh(100),asinh(1000),asinh(10000),asinh(100000),asinh(1000000),asinh(10000000),asinh(100000000)),
  label = c("1e-3","1e-2","1e-1","0","1e1","1e2","1e3","1e4","1e5","1e6","1e7","1e8"),
  step = 2,
  scipen = -1,   # when to display numbers as 1e5 rather than 10000
  range = c("0","12","0","12"), # current x and y range
  binsize = 0.2,
  binsizes = as.character(c(0.1,0.2,0.4,0.5,0.6)),
  mincount = 10
)

fcs$biex = list(
  scale = c(0,1000,2000,3000,4000,5000,6000),
  label = c(0,1000,2000,3000,4000,5000,6000),
  step = 500,
  scipen = 999,   # when to display numbers as 1e5 rather than 10000
  range = c("0","4000","0","4000"), # current x and y range
  binsize = 50,
  binsizes = as.character(c(25,50,100,200,400)),
  mincount = 5
)

fcs$coords.info = vector()
fcs$plotWindowExists = FALSE
fcs$densityWindowExists = FALSE
fcs$temp.data = FALSE
fcs$ttprojectsopen = FALSE
fcs$OverviewGate = FALSE
fcs$diploTGate = FALSE
fcs$current.xaxis = c(0,14)
fcs$current.project = ""
fcs$current.trans = "asinh"
fcs$current.projectnum = fcs$selected.projectnum = 1
fcs$current.filenum = fcs$selected.filenum = 1
fcs$current.tab = "0"
fcs$current.var1.diploT = ""
fcs$current.diploT.x = 2
fcs$current.cofactor = 1
fcs$temp.num = 0
fcs$min.counts = as.character(c(2,5,10,15,20,25,50))
fcs$min.counts.diploT = as.character(c(10,20,50,100,200,300,500))
fcs$binSizes = as.character(c(0.1,0.2,0.4,0.5,0.6))
fcs$binSizes.biex = as.character(c(25,50,100,200))
fcs$preproc.types = c("Choose option","Normalize selected files","Calculating cutoffs")
fcs$trim.num = 0.0001
fcs$plot.percentage = TRUE

### old
#fcs$col.rainbow = c("#021BC0","#021BC0","#0092FFFF","#00FFFFFF","#0DE378","#4dff4d","#B6FF00FF","#FFFF00","#FFDB00FF","#FF6D00FF","#FF0000FF","#FF0000FF")
### newIFNg.XMG12
fcs$col.rainbow = c("#001AC0","#001AC0","#00A7D0","#00F0EC","#00F884","#81FF42","#D1FF00","#FFFF00","#FFDB00FF","#FF6D00FF","#FF0000FF","#FF0000FF")
#fcs$col.rainbow.pale = c("#606DC0","#606DC0","#69BDD2","#78F1EF","#7CF8BE","#C1FFA1","#E9FF80","#FFFF80","#FFEF80","#FFB780","#FF8080")
fcs$col.rainbow.pale = c("#868EC0","#868EC0","#93C5D2","#A8F1EF","#AEF8D6","#DEFFCD","#F2FFB3","#FFFFB3","#FFFACD","#FFD4B3","#FFB3B3")
#																				mittelblau,hellblau,türkisgrün,mittelgrün,hellgrün

fcs$col.blackwhite = c("gray90","gray90","gray80","gray70","gray60","gray50","gray40","gray30","gray20","gray10","gray0","gray0")
fcs$col.blackred = c("gray95","gray95","gray80","gray65","gray50","gray35","darkred","darkred")

#fcs$ncores = detectCores() - 1

### initiate parameters
fcs$vx1=tclVar("")
fcs$vx2=tclVar("")
fcs$vy1=tclVar("")
fcs$vy2=tclVar("")

fcs$legend.space = 0
fcs$legend.space.di = 0
if (!grepl("linux",sessionInfo()$R.version$os)) {
  fcs$legend.space = 0.2
  fcs$legend.space.di = 0.03
}
### global variables - END

### if changes in database needs to bmatrix::e done
#dbGetQuery(DB, "ALTER TABLE Surface_Marker_SG_fileIndex RENAME TO Surface_Marker_SG_fileIdentity")
#dbGetQuery(DB, "ALTER TABLE Surface_Marker_SG_colnameIndex RENAME TO Surface_Marker_SG_markerIdentity")
#dbGetQuery(DB, "ALTER TABLE Cytokines_SG_fileIndex RENAME TO Cytokines_SG_fileIdentity")
#dbGetQuery(DB, "ALTER TABLE Cytokines_SG_colnameIndex RENAME TO Cytokines_SG_markerIdentity")

#rm(list = ls())
paramfile="myPRIparam.rda"
if (file.exists(paramfile)) {
	load(paramfile)
	print(sprintf("diploT :: var1=%s var2=%s",param$currVarDi1,param$currVarDi2))
	print(sprintf("triploT :: var1=%s var2=%s var3=%s",param$currVarTri1,param$currVarTri2,param$currVarTri3))
	fcs$vnrow=tclVar(param$nrow)
	fcs$vncol=tclVar(param$ncol)

	fcs$vminvalX=tclVar(param$minvalX)
	fcs$vmaxvalX=tclVar(param$maxvalX)
	fcs$vminvalY=tclVar(param$minvalY)
	fcs$vmaxvalY=tclVar(param$maxvalY)
	fcs$vminMSI=tclVar(param$minMSI)
	fcs$vmaxMSI=tclVar(param$maxMSI)
  
	fcs$rbasinh = param$cofactor
	fcs$rbcalc = tclVar(param$rbcalc)
	fcs$rbtrans = tclVar(param$rbtrans)
	fcs$cbtfeatA = tclVar(param$cbtfeatA)
	fcs$cbtfeatB = tclVar(param$cbtfeatB)
	fcs$cbtfeatC = tclVar(param$cbtfeatC)
	fcs$cbtaddFilename = tclVar(param$cbtaddFilename)
	fcs$cbtaddDate = tclVar(param$cbtaddDate)
	fcs$cbtdisplayA = tclVar(param$cbtdisplayA)
	fcs$cbtshowGrid = tclVar(param$cbtshowGrid)
	fcs$cbtshowPercentage = tclVar(param$cbtshowPercentage)
	fcs$cbttrimming = tclVar(param$cbttrimming)
	fcs$cbtmanRect = tclVar(param$cbtmanRect)
	fcs$cbtdynRange = tclVar(param$cbtdynRange)
	fcs$cbtshowLegend = tclVar(param$cbtshowLegend)
	fcs$cbtshowMinBins = tclVar(param$cbtshowMinBins)
	fcs$cbtgateData = tclVar(param$cbtgateData)
	fcs$cbtautoRect = tclVar(param$cbtautoRect)
	
} else {
	param$currVarDi1 = 1
	param$currVarDi2 = 2
	param$binSizesDiPos = 2
	param$minCountDiPos = 2

	param$currVarTri1 = 1
	param$currVarTri2 = 2
	param$currVarTri3 = 3
	param$binSizesTriPos = 2
	param$minCountTriPos = 3

	fcs$vnrow=tclVar("2")
	fcs$vncol=tclVar("4")

	fcs$vminvalX=tclVar("0")
	fcs$vmaxvalX=tclVar("14")
	fcs$vminvalY=tclVar("0")
	fcs$vmaxvalY=tclVar("14")
	fcs$vminMSI=tclVar("0")
	fcs$vmaxMSI=tclVar("0")
	
	fcs$rbasinh = "1"
	fcs$rbcalc = "MSI"
	fcs$rbtrans = "asinh"
	fcs$cbtfeatA = "0"
	fcs$cbtfeatB = "0"
	fcs$cbtfeatC = "0"
	fcs$cbtaddFilename = "1"
	fcs$cbtaddDate = "0"
	fcs$cbtdisplayA = "0"
	fcs$cbtshowGrid = "0"
	fcs$cbttrimming = "0"
	fcs$cbtmanRect = "0"
	fcs$cbtdynRange = "1"
	fcs$cbtshowLegend = "1"
	fcs$cbtshowPercentage = "1"
	fcs$cbtshowMinBins = "0"
	fcs$cbtgateData = "0"
	fcs$cbtautoRect = "0"
	
}


### helpful print function
printf <- function(...) invisible(print(sprintf(...)))

fcs$GUImain <- function () {
	this=fcs
	this$total.projects=dbListTables(this$conn)
	
	this$dataframes.name = vector()
	this$dataframes.name[1] = "Rectangle_Data_Raw"
	this$dataframes.name[2] = "Rectangle_Data_Transformed"

	# sort metadata from tables first
	# index for metadata
	idx=grep("markerIdentity|colnameIndex|fileIdentity|fileIndex|UserHistory|Classification|equipmentInfo|fileComments|SPILL",this$total.projects)

	for (i in idx){
		this$df.num = this$df.num + 1 

		nam <- paste("dataframes",this$df.num, sep="")
		df= this$getDFtable(this$total.projects[i])
		assign(nam, df, envir = fcs) 
		this$dataframes.name[this$df.num] = this$total.projects[i]
	}

	this$total.projects = this$total.projects[-idx]
	this$current.project = this$total.projects[1]

	### get indices for filenames and table
	idx.table=grep(paste0("^",this$total.projects[1],"_"),this$dataframes.name)
	idx.file=grep("fileIdentity|fileIndex",this$dataframes.name)
	if( grepl("Identity",this$dataframes.name[idx.file[1]]) ) this$fileid_name = "_fileIdentity"
	else this$fileid_name = "_fileIndex"
	idx.stain=grep("markerIdentity|colnameIndex",this$dataframes.name)[1]
	if( grepl("Identity",this$dataframes.name[idx.stain]) ) this$markerid_name = "_markerIdentity"
	else this$markerid_name = "_colnameIndex"
	

	### get filenames for table
	this$current.filetable = this$getDFtable(paste0(this$current.project,this$fileid_name))
	this$current.filenames = this$current.filetable[,2]

	### find longest filename to display in combobox
	this$max.nchar = 40
	for (i in this$current.filenames) {
		this$max.nchar = max(nchar(i),this$max.nchar)
	}
	this$max.nchar = this$max.nchar + 12

	this$current.staintable = this$getDFtable(paste0(this$current.project,this$markerid_name))
	this$current.filenum = this$current.filetable[which(this$current.filetable[,2]==this$current.filenames[1]),1]
	#this$recent.fileindex = this$current.filenum
	this$selected.filenum = this$current.filenum
	this$selected.vars = this$getVariables()
	this$current.vars = this$selected.vars

	if ( ncol(this$current.staintable)==6 ) {
		### if cutoffs and checks were already saved
		# get all cutoffs and checks
		for ( file.id in unique(this$current.filetable[,1]) ) {
			this$new.cutoffs[[file.id]] = this$current.staintable[which(this$current.staintable[,1]==file.id),5]
			this$new.checks[[file.id]] = this$current.staintable[which(this$current.staintable[,1]==file.id),6]
		}
	} else if ( ncol(this$current.staintable)==5 ) {
		### if only cutoffs were saved
		# get all cutoffs
		for ( file.id in unique(this$current.filetable[,1]) ) {
			this$new.cutoffs[[file.id]] = this$current.staintable[which(this$current.staintable[,1]==file.id),5]
			this$new.checks[[file.id]] = rep("0",nrow(this$current.staintable[which(this$current.staintable[,1]==file.id),]))
		}
	} else {
		### no cutoffs listed in this project
		# initiate cutoff and check list
		for ( file.id in unique(this$current.filetable[,1]) ) {
			this$new.cutoffs[[file.id]] = rep("0",nrow(this$current.staintable[which(this$current.staintable[,1]==file.id),]))
			this$new.checks[[file.id]] = rep("0",nrow(this$current.staintable[which(this$current.staintable[,1]==file.id),]))
		}
	}
	this$saved.cutoffs = this$new.cutoffs
	this$saved.checks = this$new.checks
	
	### get current cutoffs for our current file
	this$selected.cutoffs = this$new.cutoffs[[this$current.filenum]]
	this$current.cutoffs = this$selected.cutoffs
	this$current.checks = this$new.checks[[this$current.filenum]]
	this$recent.checks = this$current.checks

	### start GUI
	this$tt <- tktoplevel()

	### set window title
	tkwm.title(this$tt,"PRI-ana")

	### create topMenu
	topMenu <- tkmenu(this$tt)
		tkconfigure(this$tt, menu=topMenu)
		if (this$working) {
			tkadd(topMenu, "command",label="Refresh",command=this$exitandstart)
			tkadd(topMenu, "command",label="Quit",command=this$exit)
			#tkadd(topMenu, "command", label="Debugger",command=this$debug)
			tkadd(topMenu, "command", label=sprintf("Open project (%s)",length(this$total.projects)), command = this$openProject)
			tkadd(topMenu, "command", label="Open new database", command = this$openDB)
			fileMenu2 <- tkmenu(topMenu)
				tkadd(fileMenu2,"command",label="Cutoffs calculation",command=function() {this$GUIpreproc(mode=1)})
				tkadd(fileMenu2,"command",label="Normalization",command=function() {this$GUIpreproc(mode=2)})
				tkadd(fileMenu2, "command", label="Show legend colors",command=function() {this$plotLegend()})
			tkadd(topMenu, "cascade", label="Preprocession",menu=fileMenu2)
		} else {
			fileMenu <- tkmenu(topMenu)
				tkadd(fileMenu, "command", label="Open project", command = this$openProject)
				tkadd(fileMenu, "command", label="Open database", command = this$openDB)
				tkadd(fileMenu, "command", label="Quit",command=this$exit)
				tkadd(fileMenu, "command", label="Show legend colors",command=function() {this$plotLegend()})
			tkadd(topMenu, "cascade", label=sprintf("File (%s)",length(this$total.projects)),menu=fileMenu)

			fileMenu2 <- tkmenu(topMenu)
				tkadd(fileMenu2,"command",label="Cutoffs calculation",command=function() {this$GUIpreproc(mode=1)})
				tkadd(fileMenu2,"command",label="Normalization",command=function() {this$GUIpreproc(mode=2)})
			tkadd(topMenu, "cascade", label="Preprocession",menu=fileMenu2)
		}
		tkadd(topMenu, "command", label=sprintf("%s - DB: %s",this$version,sub(".sqlite3","",this$db.name)))


	len.vars= length(this$selected.vars)
	this$max.nchar.vars = max(nchar(this$selected.vars))
	
	this$width.panel = 315 + this$max.nchar.vars
  
	# start frame windows ----------------------------------------------------
	ttf_main = tkframe(this$tt)
	
	# top frame with file names and tabs --------------------------------------------------------
  ttf_maintop = tkframe(ttf_main)	
	
	### select/deselect buttons
	ttfButtons=tkframe(ttf_maintop)#,width=this$width.panel)
	### set label frame select markers of interest
	ttfselButtons=tkframe(ttfButtons)
	btselAll=tkbutton(ttfselButtons,text="Select all vars",command=function() {
	  for ( i in 1:len.vars) {
	    tclvalue(this$cbVal[[i]])="1"
	  }
	})
	btdeselAll=tkbutton(ttfselButtons,text="Deselect all vars",command=function() {
	  for ( i in 1:len ) {
	    tclvalue(this$cbVal[[i]])="0"
	  }
	})
	tkgrid(btselAll,btdeselAll,padx=6,pady=2,sticky="w")
	tkgrid(ttfselButtons)
	###
	
	### cutoff buttons
	ttfloadsave=tkframe(ttfButtons)
	btuseSavedCutoffs=tkbutton(ttfloadsave,text="Load cutoffs",command=function() {this$refreshCutoffs(saved=TRUE)})
	btsaveCutoffs2all=tkbutton(ttfloadsave,text="Save cutoffs to all",command=function() {this$saveCutoffsToAll()})
	tkgrid(btuseSavedCutoffs,btsaveCutoffs2all,padx=5,sticky="we")
	tk2tip(btuseSavedCutoffs,"Load from database. Only if available.")
	tk2tip(btsaveCutoffs2all,"Save cutoffs to all files in project.")
	tkgrid(ttfloadsave,pady=2,sticky="we")
	btresetCutoffs=tkbutton(ttfButtons,text="Reset cutoffs",command=function() {this$refreshCutoffs(reset=TRUE)})
	tkgrid(btresetCutoffs,pady=2,padx=5,sticky="we")
	# tkgrid(ttfloadsave,padx=10)
	
	# select files --------------------------------------------------------
	ttf_filenames_asinh = tkframe(ttf_maintop)
	ttf_filenames = tkframe(ttf_filenames_asinh)
	ttfchoosefile = tkframe(ttf_filenames,height=20)
	if (this$working) {
	  this$tkchoosefile = ttkcombobox(ttfchoosefile,values=this$current.filenames,width=this$max.nchar-4)
	  tkset(this$tkchoosefile,this$current.filenames[1])
	  btrefreshPanel=tkbutton(ttfchoosefile,text="GUI",command=function() {this$refreshPanel()})
	  tkgrid(btrefreshPanel,tklabel(ttfchoosefile,text="File:"),this$tkchoosefile,padx=5,sticky="w")
	} else {
	  this$tkchoosefile = ttkcombobox(ttfchoosefile,values=this$current.filenames,width=this$max.nchar)
	  tkset(this$tkchoosefile,this$current.filenames[1])
	  btrefreshPanel=tkbutton(ttfchoosefile,text="GUI",command=function() {this$refreshPanel()})
	  tkgrid(btrefreshPanel,tklabel(ttfchoosefile,text="File:"),this$tkchoosefile,padx=5,sticky="w")
	}
	tkbind(this$tkchoosefile,"<<ComboboxSelected>>",function(){this$refreshComboboxVars(tclvalue(tkget(this$tkchoosefile)))})
	tkgrid(ttfchoosefile,pady=1)
	
	ttlfcomment=ttklabelframe(ttf_filenames,text="Title text")
	this$title <- tk2text(ttlfcomment,width=this$max.nchar-15, height=2)
	displayfile = this$shortenFilename(this$current.filenames[1],title=TRUE)
	tkinsert(this$title,"1.0",displayfile)
	btaddtitle= tkbutton(ttlfcomment,text="Add title",command=function() {this$addTitle(mode="tri")})
	btadddate= tkbutton(ttlfcomment,text="Add date",command=function() {this$addTitle(mode="date")})
	tkgrid(this$title,btaddtitle,btadddate,padx=6,pady=2)
	tkgrid(ttlfcomment,pady=2)
	
	# asinh radiobuttons --------------------------------------------------------
	ttf_radioasinh = tkframe(ttf_filenames_asinh)
	tkgrid(tklabel(ttf_radioasinh,text="asinh cof"),sticky="e")
	tk2tip(ttf_radioasinh,"arcsinh() cofactor")
	tfradios = tkframe(ttf_radioasinh)
	rb1=tkradiobutton(tfradios,variable=this$rbasinh,value="0.1")
	rb2=tkradiobutton(tfradios,variable=this$rbasinh,value="0.2")
	rb3=tkradiobutton(tfradios,variable=this$rbasinh,value="1")
	rb4=tkradiobutton(tfradios,variable=this$rbasinh,value="5")
	rb5=tkradiobutton(tfradios,variable=this$rbasinh,value="150")
	tkgrid(rb1,tklabel(tfradios,text="1/10"))
	tkgrid(rb2,tklabel(tfradios,text="1/5"))
	tkgrid(rb3,tklabel(tfradios,text="1"),sticky="we")
	tkgrid(rb4,tklabel(tfradios,text="5"),sticky="we")
	tkgrid(rb5,tklabel(tfradios,text="150"),sticky="we")
	tk2tip(rb1,"asinh(x*10)")
	tk2tip(rb2,"asinh(x*5)")
	tk2tip(rb3,"asinh(x)")
	tk2tip(rb4,"asinh(x/5)")
	tk2tip(rb5,"asinh(x/150)")
	tkgrid(tfradios)
	
	tkgrid(ttf_filenames, ttf_radioasinh,sticky="we",padx=2)
	tkgrid(ttfButtons,ttf_filenames_asinh,padx=10,sticky="w")
	tkgrid(ttf_maintop,sticky="w")
	###

  # bottom frame ------------------------------------------------------------------------------
	ttf_mainbottom = tkframe(ttf_main)
	this$ttf_markers = tkframe(ttf_mainbottom)
	### scrollable frame
	this$sw <- tkwidget(this$ttf_markers,"ScrolledWindow",relief="sunken",borderwidth=2)
	scrollFrame <- tkwidget(this$sw,"ScrollableFrame",height=630,width=this$width.panel)
	tcl(this$sw,"setwidget",scrollFrame)
	subfID <- tclvalue(tcl(scrollFrame,"getframe"))
	#subfID <- tkgetframe(scrollFrame)$ID

	# label project name
	this$project.lab <- tcl("label",paste(subfID,".lab",sep=""),text=sprintf(">>%s",this$shortenFilename(this$current.project)))
	tkgrid(this$project.lab,sticky="w",padx=5)

	i = j = 1
	while (i < (len.vars*4) ) {
	  	### initiate
		this$vcutoffs[[j]] = tclVar(this$selected.cutoffs[j])
	  	this$cbVal[[j]]=tclVar(this$current.checks[j])
		this$cbcutoffperc[[j]] = tclVar("0")
		
		cutoffcb=tcl("checkbutton",paste(subfID,".",i,sep=""),variable=this$cbVal[[j]],text=this$selected.vars[j])
	  	cutofflabel=tcl("label",paste(subfID,".",i+1,sep=""),text="cutoff: ")
	  	cutoffentry=tcl("entry",paste(subfID,".",i+2,sep=""),width=6,textvariable=this$vcutoffs[[j]])
		cutoffperccb=tcl("checkbutton",paste(subfID,".",i+3,sep=""),variable=this$cbcutoffperc[[j]],text="in %")
	  	tkgrid(cutoffcb,cutofflabel,cutoffentry,cutoffperccb,sticky="w",pady=1)
	  	#tkbind(cutoffcb,"<FocusIn>",function() tcl(scrollFrame,"see",cutoffcb) )
	  	i = i + 4
	  	j = j + 1
	}
	tkgrid(this$sw)
	###
	
	# notebook panels --------------------------------------------------------
	ttf_tabs <- tkframe(ttf_mainbottom)
	this$panelTabs=ttknotebook(ttf_tabs)
		tab1=ttkframe(this$panelTabs)
		tab2=ttkframe(this$panelTabs)
		tab3=ttkframe(this$panelTabs)
		tab4=ttkframe(this$panelTabs)

		tkadd(this$panelTabs,tab1,text="n-diploTs")
		tkadd(this$panelTabs,tab2,text="n-triploTs")
		tkadd(this$panelTabs,tab3,text="table info")
		tkadd(this$panelTabs,tab4,text="log")
		tkbind(this$panelTabs,"<<NotebookTabChanged>>",this$clickTab)

		# tab: diploT
		this$GUIdiploT(tab1)
		# tab: triploT
		this$GUItriploT(tab2)
		# tab: table info
		this$GUIdataframe(tab3)
		# tab: log
		this$GUIlog(tab4)
	tkgrid(this$panelTabs,sticky="w",padx=5)
	###

	tkgrid(this$ttf_markers,ttf_tabs)
	tkgrid(ttf_mainbottom,pady=2,sticky="ns")
	tkgrid(ttf_main,pady=2)

	this$refreshComboboxStart(1)
	
	this$openProjectYet = TRUE
	
	### if there are more than one project, chose a project first
	if(length(this$total.projects)>1) {
	  this$openProjectYet = FALSE
	  printf("Chose your project first.")
	  this$openProject()
	  tkraise(this$ttprojects)
	  
	  ### if not clicked on a project yet
	  # but clicked on main window already
	  # bind mouse click with function
	  if (!this$openProjectYet) {
	    # tkbindtext(this$tt, "<Button-1>", function(...) {
	    #   tkmessageBox(title = "Chose project", message = "Chose your project first.", icon = "info", type = "ok")
	    #   })
	    tkbind(this$tt, "<Button-1>", function(...) {
	      tkmessageBox(title = "Chose project", message = "Chose your project first.", icon = "info", type = "ok");
	      tkraise(this$ttprojects)
	    })
	    tkraise(this$ttprojects)
	  } 
	} else {
	  tkraise(this$tt)
	}
}

fcs$GUIdiploT <- function(tab) {
	this=fcs

	if (this$working) printf("w: GUIdiploT")
	
	var.length = 38
	if (this$max.nchar.vars>var.length) var.length=this$max.nchar.vars

	# set frame and content of combobox with file, mincount, binSize
	ttftopframe = tkframe(tab)
	ttfcombo=tkframe(ttftopframe)
	# combobox: binSize
	this$binSizedi=ttkcombobox(ttfcombo,values=this$binSizes,width=var.length)
	tkset(this$binSizedi,this$binSizes[param$binSizesDiPos])
	# combobox: mincount cells
	this$minCountDi=ttkcombobox(ttfcombo,values=this$min.counts.diploT,width=var.length)
	tkset(this$minCountDi,this$min.counts.diploT[param$minCountDiPos])
	# combobox: param
	# if FSC/SSC are listed
	#tmp = max(grep("SC",this$selected.vars)) 
	#tmp2 = min(grep(c("Time"),this$selected.vars))
	#if (tmp2>tmp & !is.infinite(tmp2) & tmp2<(length(this$selected.vars)-3)) tmp=tmp2
	#if ( ((tmp+3)>length(this$selected.vars)) | is.infinite(tmp)) tmp = 0
	this$cbvar1di=ttkcombobox(ttfcombo,values=this$selected.vars,width=var.length)
	tkset(this$cbvar1di,this$selected.vars[param$currVarDi1])
	this$cbvar2di=ttkcombobox(ttfcombo,values=this$selected.vars,width=var.length)
	tkset(this$cbvar2di,this$selected.vars[param$currVarDi2])
	# call comboboxleft
	tkgrid(tklabel(ttfcombo,text="binSize:"),this$binSizedi,sticky="e",padx=5)
	tkgrid(tklabel(ttfcombo,text="minCount:"),this$minCountDi,sticky="e",padx=5)
	tkgrid(tklabel(ttfcombo,text="Feature A:"),this$cbvar1di,sticky="e",padx=5)
	tkgrid(tklabel(ttfcombo,text="Feature B:"),this$cbvar2di,sticky="e",padx=5)
	# tkgrid(ttfcombo,padx=5,pady=1,sticky="w")

	# plot button
	ttfbut=tkframe(ttftopframe)
	btplot= tkbutton(ttfbut,text="Plot diploT",command=function() {this$dodiploT(mode="single")})
	tkgrid(btplot,sticky="snwe",padx=10,pady=1)
	#tkgrid(btplot,sticky="snwe",padx=10,pady=1)
	tkconfigure(btplot,height=2,width=25)
	bthist=tkbutton(ttfbut,text="Plot histograms",command=function(){this$plotHistograms(plotter="di",pdf=FALSE)})
	tkgrid(bthist,sticky="we",padx=10,pady=1)
	tkgrid(ttfcombo,ttfbut,padx=5,pady=1,sticky="w")
	tkgrid(ttftopframe)
	
	##### set rect frame ttfvals with additional options min/max FI values  
	ttfvals=tkframe(tab)
	# set label frame min/max FI in rect frame "ttfvals"
	ttlfFI=ttklabelframe(ttfvals,text="Plot options")
	ttfrange = tkframe(ttlfFI)
	# # X axis
	# minvalXdi=tkentry(ttfrange,width=6,textvariable=this$vminvalX)
	# maxvalXdi=tkentry(ttfrange,width=6,textvariable=this$vmaxvalX)
	# tkgrid(tklabel(ttfrange,text="Range x-axis (min/max):"),minvalXdi,maxvalXdi,padx=7,pady=1,sticky="w")
	# dynamic range
	ttfdyn = tkframe(ttfrange)
	ttfdynL = tkframe(ttfdyn)
	tkgrid(tklabel(ttfdynL,text="Range x-axis (min/max):",pady=1,padx=15))
	cbtdynrange=tkcheckbutton(ttfdynL,variable=this$cbtdynRange,text="Dynamic MSI range")
	tkgrid(cbtdynrange,sticky="nws",pady=10,padx=15)
	tk2tip(cbtdynrange,"Uncheck for manual input")
	# second param
	ttfdynR = tkframe(ttfdyn)
	minvalXdi=tkentry(ttfdynR,width=6,textvariable=this$vminvalX)
	maxvalXdi=tkentry(ttfdynR,width=6,textvariable=this$vmaxvalX)
	tkgrid(minvalXdi,maxvalXdi,padx=5,pady=1,sticky="w")
	minMSI=tkentry(ttfdynR,width=4,textvariable=this$vminMSI)
	maxMSI=tkentry(ttfdynR,width=4,textvariable=this$vmaxMSI)
	tkgrid(tklabel(ttfdynR,text="min(MSI):"),minMSI,padx=5,pady=1,sticky="e")
	tkgrid(tklabel(ttfdynR,text="max(MSI):"),maxMSI,padx=5,pady=1,sticky="e")
	tkgrid(ttfdynL,ttfdynR)
	tkgrid(ttfdyn)
	tkgrid(ttfrange)
	# add checkbuttons for plot attributss
	ttfplot = tkframe(ttlfFI)
	ttfgridnew = tkframe(ttfplot)
	cbtgrid=tkcheckbutton(ttfgridnew,variable=this$cbtshowGrid,text="Grid")
	tkgrid(cbtgrid,sticky="w")
	tk2tip(cbtgrid,"Check to add grid.")
	cbtaddfilename=tkcheckbutton(ttfgridnew,variable=this$cbtaddFilename,text="File name")
	tkgrid(cbtaddfilename,sticky="w")
	tk2tip(cbtaddfilename,"Check to add file name at every diploT.")
	cbtaddDate=tkcheckbutton(ttfgridnew,variable=this$cbtaddDate,text="Date")
	tkgrid(cbtaddDate,sticky="w")
	cbtdisplayA=tkcheckbutton(ttfgridnew,variable=this$cbtdisplayA,text="Hide A-")
	tk2tip(cbtdisplayA,"Check if no distraction of (Feature A-) is desired.")
	tkgrid(cbtdisplayA,sticky="w")
	ttfsave = tkframe(ttfplot)
	btnewplot = tkbutton(ttfsave,text="New diploT window",command=function() {this$newdiploT(26)})
	tkgrid(btnewplot,padx=5)
	btsavePDF = tkbutton(ttfsave,text="Save active window",command=function() {this$saveWindow()})
	tkgrid(btsavePDF,padx=5,pady=1)
	tk2tip(btsavePDF,"Save as pdf")
	tkgrid(ttfgridnew,ttfsave,padx=10)
	tkgrid(ttfplot)
	# option panel
	ttlfinfo=ttklabelframe(ttfvals,text="Calculation options")
	# add radio buttons for transformation type
	ttfradios=tkframe(ttlfinfo)
	ttfradiosright1=tkframe(ttfradios)
	rb1=tkradiobutton(ttfradiosright1,variable=this$rbtrans,value="asinh")
	rb2=tkradiobutton(ttfradiosright1,variable=this$rbtrans,value="none")
	#rb2=tkradiobutton(ttfradiosright1,variable=this$rbtrans,value="biex")
	tkgrid(rb1,tklabel(ttfradiosright1,text="asinh"),rb2,tklabel(ttfradiosright1,text="none"))
	tkgrid(tklabel(ttfradios,text="Transformation:"),ttfradiosright1,sticky="w",pady=1)
	tkbind(rb1, "<Button-1>", function() { this$changeFI.range(1) })
	#tkbind(rb2, "<Button-1>", function() { this$changeFI.range(2) })

	# add radio buttons for calculation type
	ttfradiosright2=tkframe(ttfradios)
	rb1=tkradiobutton(ttfradiosright2,variable=this$rbcalc,value="MSI")
	rbl1=tklabel(ttfradiosright2,text="MSI")
	rb2=tkradiobutton(ttfradiosright2,variable=this$rbcalc,value="freq")
	rbl2=tklabel(ttfradiosright2,text="frequency")
	rb3=tkradiobutton(ttfradiosright2,variable=this$rbcalc,value="density")
	rbl3=tklabel(ttfradiosright2,text="density")
	rb4=tkradiobutton(ttfradiosright2,variable=this$rbcalc,value="SD")
	rbl4=tklabel(ttfradiosright2,text="SD")
	rb5=tkradiobutton(ttfradiosright2,variable=this$rbcalc,value="SEM")
	rbl5=tklabel(ttfradiosright2,text="SEM")
	rb6=tkradiobutton(ttfradiosright2,variable=this$rbcalc,value="RSEM")
	rbl6=tklabel(ttfradiosright2,text="RSEM")
	tkgrid(rb1,rbl1,rb2,rbl2,sticky="w")
	tk2tip(rb1,"Mean fluorescence intensity: If cutoff is set, it only considers producing cells of feature Y")
	tk2tip(rbl1,"Mean fluorescence intensity: If cutoff is set, it only considers producing cells of feature Y")
	tk2tip(rb2,"Frequency of producing cells")
	tk2tip(rbl2,"Frequency of producing cells")
	tkgrid(rb3,rbl3,rb4,rbl4,sticky="w")
	tk2tip(rb3,"Cell count")
	tk2tip(rbl3,"Cell count")
	tk2tip(rb4,"Standard deviation")
	tk2tip(rbl4,"Standard deviation")
	tkgrid(rb5,rbl5,rb6,rbl6,sticky="w")
	tk2tip(rb5,"Standard error of mean")
	tk2tip(rbl5,"Standard error of mean")
	tk2tip(rb6,"Relative standard error of mean")
	tk2tip(rbl6,"Relative standard error of mean")
	tkgrid(tklabel(ttfradios,text="Statistical method:"),ttfradiosright2,sticky="nw",pady=1)
	#tk2tip(ttfradiosright2,"If cutoff is set, it only considers producing cells of feature 3.")
	tkgrid(ttfradios,padx=10,pady=1)
	ttfbuttons=tkframe(ttlfinfo)
	btremovdub=tkbutton(ttfbuttons,text="Remove Doublets",command=function(){this$preprocData()})
	ttfcheckbuttons=tkframe(ttfbuttons)
	cbttrimming=tkcheckbutton(ttfcheckbuttons,variable=this$cbttrimming,text="Trim first")
	tk2tip(cbttrimming,sprintf("Trims %s of each column.",this$trim.num))
	cbtshowMinBins=tkcheckbutton(ttfcheckbuttons,variable=this$cbtshowMinBins,text="Bins <minCount")
	tk2tip(cbtshowMinBins,"Show bins in pale color.")
	tkgrid(cbttrimming,sticky="w",pady=1)
	tkgrid(cbtshowMinBins,sticky="w")
	tkgrid(btremovdub,ttfcheckbuttons,padx=5,pady=1,sticky="we")
	tkgrid(ttfbuttons,sticky="w")
	tkgrid(ttlfFI,ttlfinfo,padx=5,pady=1,sticky="wens")
	tkgrid(ttfvals,pady=1)
	#####

	# plot button "overview"
	ttfplotten = tkframe(tab)
	btplot1= tkbutton(ttfplotten,text="Plot diploT-Overview with fixed Feature A",command=function() {this$dodiploT(mode="overview")})
	tk2tip(btplot1,"max. 35 diploTs")
	btplot2= tkbutton(ttfplotten,text="Plot diGraph-Overview with fixed Feature A",command=function() {this$dodiGraph(mode="overview")})
	tk2tip(btplot2,"max. 10 diGraphs")
	tkgrid(btplot1,btplot2,sticky="snwe",padx=15,pady=5)
	tkconfigure(btplot1,height=2)
	btplot3= tkbutton(ttfplotten,text="Plot diploT-Overview total",command=function() {this$dodiploTtotal()})
	btplot4= tkbutton(ttfplotten,text="Plot diGraph-Overview0 with fixed Feature A",command=function() {this$dodiGraph(mode="overview",set.zero=TRUE)})
	tk2tip(btplot4,"max. 10 diGraphs")
	tkgrid(btplot3,btplot4,sticky="snwe",padx=15,pady=5)
	tkconfigure(btplot3,height=2)
	tkgrid(ttfplotten)

	# frame: Data Info
	lfinfo=ttklabelframe(tab,text="Data info")
	this$ncell.gui.di=tklabel(lfinfo,text=as.character(this$origin.ncells))
	this$ncell.sel.gui.di=tklabel(lfinfo,text=as.character(this$ncell.sel))
	this$ncell.perc.gui.di=tklabel(lfinfo,text=as.character(this$ncell.perc))
	tkgrid(tklabel(lfinfo,text="# of total cells:"),this$ncell.gui.di,sticky="w",padx=5)
	tkgrid(tklabel(lfinfo,text="# of selected cells:"),this$ncell.sel.gui.di,tklabel(lfinfo,text="("),this$ncell.perc.gui.di,tklabel(lfinfo,text="%)"),sticky="w",padx=5)
	# call frame Data Info
	tkgrid(lfinfo,sticky="s",columnspan=2,padx=5)
	tkgrid.configure(lfinfo,sticky="swe",pady=5,padx=5)  
}

fcs$GUItriploT <- function(tab){
	this=fcs

	if (this$working) printf("w: GUItriploT")

	# set frame and content of combobox with file, mincount, binSize
	ttftopframe = tkframe(tab)
	ttfcombo=tkframe(ttftopframe)
	# content comboboxleft
	
	var.length = 39
	if (this$max.nchar.vars>var.length) var.length=this$max.nchar.vars
	this$minCountTri=ttkcombobox(ttfcombo,values=this$min.counts,width=var.length)
	tkset(this$minCountTri,this$min.counts[param$minCountTriPos])
	this$binSize=ttkcombobox(ttfcombo,values=this$binSizes,width=var.length)
	tkset(this$binSize,this$binSizes[param$binSizesTriPos])

	# if FSC/SSC are listed
	#tmp = max(grep("SC",this$selected.vars)) 
	#tmp2 = min(grep(c("Time"),this$selected.vars))
	#if (tmp2>tmp & !is.infinite(tmp2) & tmp2<(length(this$selected.vars)-3)) tmp=tmp2
	#if ( ((tmp+3)>length(this$selected.vars)) | is.infinite(tmp)) tmp = 0

	this$cbvar1=ttkcombobox(ttfcombo,values=this$selected.vars,width=var.length)
	tkset(this$cbvar1,this$selected.vars[param$currVarTri1])
	this$cbvar2=ttkcombobox(ttfcombo,values=this$selected.vars,width=var.length)
	tkset(this$cbvar2,this$selected.vars[param$currVarTri2])
	this$cbvar3=ttkcombobox(ttfcombo,values=this$selected.vars,width=var.length)
	tkset(this$cbvar3,this$selected.vars[param$currVarTri3])
	cbtfeatA=tkcheckbutton(ttfcombo,variable=this$cbtfeatA)
	cbtfeatB=tkcheckbutton(ttfcombo,variable=this$cbtfeatB)
	cbtfeatC=tkcheckbutton(ttfcombo,variable=this$cbtfeatC)
	# call comboboxleft
	tkgrid(tklabel(ttfcombo,text="binSize:"),this$binSize,sticky="e",padx=5)
	tkgrid(tklabel(ttfcombo,text="minCount:"),this$minCountTri,tklabel(ttfcombo,text="Fix"),sticky="e",padx=5)
	tkgrid(tklabel(ttfcombo,text="Feature A:"),this$cbvar1,cbtfeatA,sticky="e",padx=5)
	tk2tip(cbtfeatA,"Check to fix this feature for triploT-Overview")
	tkgrid(tklabel(ttfcombo,text="Feature B:"),this$cbvar2,cbtfeatB,sticky="e",padx=5)
	tk2tip(cbtfeatB,"Check to fix this feature for triploT-Overview")
	tkgrid(tklabel(ttfcombo,text="Feature C:"),this$cbvar3,cbtfeatC,sticky="e",padx=5)
	tk2tip(cbtfeatC,"Check to fix this feature for triploT-Overview")
	# tkgrid(ttfcombo,padx=25,pady=1,sticky="w")
	
	#### rect buttons "Plot"
	ttfbut=tkframe(ttftopframe)
	btplot= tkbutton(ttfbut,text="Plot triploT",command=this$dotriploT)
	btplot2= tkbutton(ttfbut,text="Plot for all files",command=this$dotriploTfiles)
	yhbtplot2= tkbutton(ttfbut,text="Read Table",command=this$dotriploTtable)
	bthist=tkbutton(ttfbut,text="Plot histograms",command=function(){this$plotHistograms(plotter="tri",pdf=FALSE)})
	if (this$working) {
	  tkgrid(btplot,sticky="snwe",padx=2,pady=1)
	  tkgrid(btplot2,sticky="snwe",padx=2,pady=1)
	  tkgrid(yhbtplot2,sticky="snwe",padx=2,pady=1)
	  tkgrid(bthist,sticky="snwe",padx=2,pady=1)
	} else {
	  tkgrid(btplot,btplot2,bthist,sticky="snwe",padx=5,pady=1)
	}
	tkconfigure(btplot,height=2,width=25)
	#btdens=tkbutton(ttfbut,text="Plot density",command=function(){this$plotDensities(plotter="tri",pdf=FALSE)})
	#tkgrid(btdens,sticky="we",pady=1)
	tkgrid(ttfcombo,ttfbut,padx=5,sticky="w")
	tkgrid(ttftopframe)


	##### set rect frame ttfvals with additional options min/max FI values  
	ttfvals=tkframe(tab)
	# set label frame min/max FI in rect frame "ttfvals"
	ttlfFI=ttklabelframe(ttfvals,text="Plot options")
	# checkbox axes range label -------------------------------------------
	ttfRange = tkframe(ttlfFI)
	ttfRangeLabel = tkframe(ttfRange)
	tkgrid(tklabel(ttfRangeLabel,text="Range x-axis (min/max):"),padx=10,pady=1)
	tkgrid(tklabel(ttfRangeLabel,text="Range y-axis (min/max):"),padx=10,pady=1)
	cbtdynrange=tkcheckbutton(ttfRangeLabel,variable=this$cbtdynRange,text="Dynamic MSI range")
	tkgrid(cbtdynrange,sticky="ews",padx=10,pady=10)
	tk2tip(cbtdynrange,"Uncheck for manual input")
	# checkbox axes range value -------------------------------------------
	ttfRangeValue = tkframe(ttfRange)
	# X axis
	this$minvalX=tkentry(ttfRangeValue,width=6,textvariable=this$vminvalX)
	this$maxvalX=tkentry(ttfRangeValue,width=6,textvariable=this$vmaxvalX)
	tkgrid(this$minvalX,this$maxvalX,padx=5,pady=1,sticky="e")
	# Y axis
	this$minvalY=tkentry(ttfRangeValue,width=6,textvariable=this$vminvalY)
	this$maxvalY=tkentry(ttfRangeValue,width=6,textvariable=this$vmaxvalY)
	tkgrid(this$minvalY,this$maxvalY,padx=5,pady=1,sticky="e")
	# Z range
	this$minMSI=tkentry(ttfRangeValue,width=4,textvariable=this$vminMSI)
	this$maxMSI=tkentry(ttfRangeValue,width=4,textvariable=this$vmaxMSI)
	tkgrid(tklabel(ttfRangeValue,text="min(MSI):"),this$minMSI,padx=5,pady=1,sticky="e")
	tkgrid(tklabel(ttfRangeValue,text="max(MSI):"),this$maxMSI,padx=5,pady=1,sticky="e")
	tkgrid(ttfRangeLabel,ttfRangeValue, sticky="w")
	tkgrid(ttfRange)
# checkbox show grid --------------------------------------------------
	ttfcheckPlot = tkframe(ttlfFI)
	cbtshowgrid=tkcheckbutton(ttfcheckPlot,variable=this$cbtshowGrid,text="Grid")
	tk2tip(cbtshowgrid,"Check to plot grid lines")
	cbtpercentage = tkcheckbutton(ttfcheckPlot, variable=this$cbtshowPercentage,text="Percentages")
	tk2tip(cbtpercentage,"Check to plot percentages")
	cbtshowlegend=tkcheckbutton(ttfcheckPlot,variable=this$cbtshowLegend,text="Legend")
	tk2tip(cbtshowlegend,"Uncheck to hide legend on plot")
	cbtaddDate=tkcheckbutton(ttfcheckPlot,variable=this$cbtaddDate,text="Date")
	tkgrid(cbtshowgrid,cbtpercentage,cbtshowlegend,cbtaddDate,padx=5,pady=8,sticky="we")
	tkgrid(ttfcheckPlot)
	# option panel

# radio buttons for transformation type -----------------------------------
	ttlfinfo=ttklabelframe(ttfvals,text="Calculation options")
	ttfradios=tkframe(ttlfinfo)
	ttfradiosright1=tkframe(ttfradios)
	rb1=tkradiobutton(ttfradiosright1,variable=this$rbtrans,value="asinh")
	rb2=tkradiobutton(ttfradiosright1,variable=this$rbtrans,value="none")
	tkgrid(rb1,tklabel(ttfradiosright1,text="asinh"),
		rb2,tklabel(ttfradiosright1,text="none"),sticky="we")
	tkgrid(tklabel(ttfradios,text="Transformation:"),ttfradiosright1,sticky="w")
	tkbind(rb1, "<Button-1>", function() { this$changeFI.range(1) })
	#tkbind(rb2, "<Button-1>", function() { this$changeFI.range(2) })

	# add radio buttons for calculation type
	ttfradiosright2=tkframe(ttfradios)
	rb1=tkradiobutton(ttfradiosright2,variable=this$rbcalc,value="MSI")
	rbl1=tklabel(ttfradiosright2,text="MSI")
	rb2=tkradiobutton(ttfradiosright2,variable=this$rbcalc,value="MSI(+)")
	rbl2=tklabel(ttfradiosright2,text="MSI(+)")
	rb3=tkradiobutton(ttfradiosright2,variable=this$rbcalc,value="density")
	rbl3=tklabel(ttfradiosright2,text="density")
	rb4=tkradiobutton(ttfradiosright2,variable=this$rbcalc,value="freq")
	rbl4=tklabel(ttfradiosright2,text="frequency")
	rb5=tkradiobutton(ttfradiosright2,variable=this$rbcalc,value="SD")
	rbl5=tklabel(ttfradiosright2,text="SD")
	rb6=tkradiobutton(ttfradiosright2,variable=this$rbcalc,value="SEM")
	rbl6=tklabel(ttfradiosright2,text="SEM")
	rb7=tkradiobutton(ttfradiosright2,variable=this$rbcalc,value="RSEM")
	rbl7=tklabel(ttfradiosright2,text="RSEM")
	tkgrid(rb1,rbl1,rb2,rbl2,sticky="w")
	tk2tip(rb1,"Mean signal intensity")
	tk2tip(rbl1,"Mean signal intensity")
	tk2tip(rb2,"Mean signal intensity of Feature C producing cells")
	tk2tip(rbl2,"Mean signal intensity of Feature C producing cells")
	tkgrid(rb3,rbl3,rb4,rbl4,sticky="w")
	tk2tip(rb3,"Cell count")
	tk2tip(rbl3,"Cell count")
	tk2tip(rb4,"Frequency of Feature C producing cells")
	tk2tip(rbl4,"Frequency of Feature C producing cells")
	tkgrid(rb5,rbl5,rb6,rbl6,sticky="w")
	tk2tip(rb5,"Standard deviation")
	tk2tip(rbl5,"Standard deviation")
	tk2tip(rb6,"Standard error of mean")
	tk2tip(rbl6,"Standard error of mean")
	tkgrid(rb7,rbl7,sticky="w")
	tk2tip(rb7,"Relative standard error of mean")
	tk2tip(rbl7,"Relative standard error of mean")
	tkgrid(tklabel(ttfradios,text="Statistical method:"),ttfradiosright2,sticky="nw")
	tkgrid(ttfradios,padx=10)
	ttfbuttons=tkframe(ttlfinfo)
	btremovdub=tkbutton(ttfbuttons,text="Remove Doublets",command=function(){this$preprocData()})
	tk2tip(btremovdub,"Only possible with FSH/SSH channels")
	cbttrimming=tkcheckbutton(ttfbuttons,variable=this$cbttrimming,text="Trim first")
	tk2tip(cbttrimming,sprintf("Trims %s of each column.",this$trim.num))
	cbtshowMinBins=tkcheckbutton(ttfbuttons,variable=this$cbtshowMinBins,text="Bins <minCount")
	tk2tip(cbtshowMinBins,"Show bins in pale color.")
	tkgrid(btremovdub,cbttrimming,cbtshowMinBins,padx=5,pady=1,sticky="we")
	tkgrid(ttfbuttons,sticky="w")
	tkgrid(ttlfFI,ttlfinfo,padx=5,pady=1,sticky="nsew")
	tkgrid(ttfvals)
	#####

	##### set rect frame ttfgraph with Graphics/Comment and Set Rectangle   
	ttfgraph=tkframe(tab)
	##### set rect frame ttfgraphcomment with Graphics and Comment
	ttfgraphcomment=tkframe(ttfgraph)
	# frame: Graphics
	btOverviewXY = tkbutton(ttfgraphcomment,text="Plot triploT-Overview with fixed Features A and B", command = function() {this$dotriploTOverviewXY()})
	btOverviewX = tkbutton(ttfgraphcomment,text="Plot triploT-Overview with fixed Feature A", command = function() {this$dotriploTOverviewX()})
	btOverview = tkbutton(ttfgraphcomment,text="Plot triploT-Overview total", command = function() {this$dotriploTOverview()})
	
	tkgrid(btOverviewXY,sticky="snwe",pady=1,padx=5)
	tkconfigure(btOverviewXY,height=2)
	tk2tip(btOverviewX,"Select your Feature A and Y above and check all features to the left which you would like to plot as heat bins.")
	tkgrid(btOverviewX,sticky="snwe",pady=1,padx=5)
	tkconfigure(btOverviewX,height=2)	
	tk2tip(btOverviewX,"Select your Feature A above and check all features to the left which you would like to plot.")
	tkgrid(btOverview,sticky="snwe",pady=1,padx=5)
	tkconfigure(btOverview,height=2)
	tk2tip(btOverview,"Check your features to the left. If more than 16 markers are selected, several PDFs will be created.")
	ttlfgraph=ttklabelframe(ttfgraphcomment,text="Graphics")
	ttfrowcol=tkframe(ttlfgraph)
	this$ncol=tkentry(ttfrowcol,width=2,textvariable=this$vncol)
	this$nrow=tkentry(ttfrowcol,width=2,textvariable=this$vnrow)
	# call buttons and frame Graphics
	tkgrid(tklabel(ttfrowcol,text="Rows:"),this$nrow,padx=10,sticky="e")
	tkgrid(tklabel(ttfrowcol,text="Columns:"),this$ncol,padx=10,sticky="e")
	tkgrid(tkbutton(ttfrowcol,text="New plot window",command=this$newPlotWindow))
	ttfgridsave = tkframe(ttlfgraph)
	ttfsave = tkframe(ttfgridsave)
	btsavePNG = tkbutton(ttfsave,text="Save last plot as PNG",command=function(){this$savePNG(single=TRUE)})
	btsavePDF = tkbutton(ttfsave,text="Save active window",command=function(){this$saveWindow(type="triploT")})
	tkgrid(btsavePNG,pady=1,sticky="we")
	tkgrid(btsavePDF,pady=1,sticky="we")
	tk2tip(btsavePDF,"Save as pdf")
	tk2tip(btsavePNG,"Save as jpg")
	tkgrid(ttfsave,padx=3)
	tkgrid(ttfrowcol,ttfgridsave,padx=5,pady=1)
	# set active window buttons
	#tkgrid(tkbutton(ttlfgraph,text="Set prev window as active",command=function(){dev.set(dev.prev())}),
	#		tkbutton(ttlfgraph,text="Set next window as active",command=function(){dev.set(dev.next())}),padx=5,pady=1)
	tkgrid(ttlfgraph,padx=5,pady=1)


	### set label frame "Set Rectangle"
	ttlfRect=ttklabelframe(ttfgraph,text="Set rectangle")
	ttfdef=tkframe(ttlfRect)
	butdef=tkbutton(ttfdef,text="Click rectangle",command=this$doRect)
	# checkbox manual rectangle
	cbtmanrect=tkcheckbutton(ttfdef,variable=this$cbtmanRect,text="Manual rectangle")
	tkgrid(butdef, cbtmanrect,padx=10,pady=1)
	tk2tip(butdef,"Define rectangle by mouse clicks")
	tk2tip(cbtmanrect,"Define rectangle with coordinates below")
	tkgrid(ttfdef)
	#
	ttfcoords=tkframe(ttlfRect)
	# choose coords from x axis
	ttfcoords.x=tkframe(ttfcoords)
	x1=tkentry(ttfcoords.x,width=4,textvariable=this$vx1)
	x2=tkentry(ttfcoords.x,width=4,textvariable=this$vx2)
	tkgrid(tklabel(ttfcoords.x,text="x1:"),x1,padx=6,sticky="e",pady=1)
	tkgrid(tklabel(ttfcoords.x,text="x2:"),x2,padx=6,sticky="e",pady=1)
	# choose coords from y ayis
	ttfcoords.y=tkframe(ttfcoords)
	y1=tkentry(ttfcoords.y,width=4,textvariable=this$vy1)
	y2=tkentry(ttfcoords.y,width=4,textvariable=this$vy2)
	tkgrid(tklabel(ttfcoords.y,text="y1:"),y1,padx=8,sticky="e",pady=1)
	tkgrid(tklabel(ttfcoords.y,text="y2:"),y2,padx=8,sticky="e",pady=1)
	tkgrid(ttfcoords.x,ttfcoords.y)
	tkgrid(ttfcoords)
	#
	ttfcheckboxes=tkframe(ttlfRect)
	# checkbox tmp.data
	cbttmpdata=tkcheckbutton(ttfcheckboxes,variable=this$cbtgateData,text="Gate data")
	#ä#tkgrid(,padx=5,pady=1,sticky="n")
	tk2tip(cbttmpdata,"Save gated data temporarily")
	# call button "Rect" and "Redo Rect"
	ttfaddInfo=tkframe(ttfcheckboxes)
	# checkbox auto rectangle
	cbtautorect=tkcheckbutton(ttfaddInfo,variable=this$cbtautoRect,text="Auto add info")
	tk2tip(cbtautorect,"Add rect and cell info automatically")
	tkgrid(cbtautorect,padx=5,sticky="n")
	tkgrid(tkbutton(ttfaddInfo,text="Add rect info",command=function(){this$addRectInfo(setcex=1.0)}),pady=1)
	tkgrid(tkbutton(ttfaddInfo,text="Add cell info",command=function(){this$addCellInfo(setcex=1.0)}))
	tkgrid(cbttmpdata,ttfaddInfo,padx=5)
	tkgrid(ttfcheckboxes,padx=5,pady=1,sticky="nsew")
	##
	btplotrect=tkbutton(ttlfRect,text="Plot rectangle data only",command=this$dotriploTRectData)
	tkgrid(btplotrect,sticky="we",pady=5,padx=5)
	tkconfigure(btplotrect,height=2)
	tkgrid(ttfgraphcomment,ttlfRect,padx=5,sticky="we")
	## 
	tkgrid(ttfgraph,padx=5)
	##### set rect frame ttfgraph with Graphics and rectangle end
	
	# frame: Data Info
	lfinfo=ttklabelframe(tab,text="Data info")
	this$ncell.gui=tklabel(lfinfo,text=as.character(this$origin.ncells))
	this$ncell.sel.gui=tklabel(lfinfo,text=as.character(this$ncell.sel))
	this$ncell.perc.gui=tklabel(lfinfo,text=as.character(this$ncell.perc))
	tkgrid(tklabel(lfinfo,text="# of total cells:"),this$ncell.gui,sticky="w",padx=5)
	tkgrid(tklabel(lfinfo,text="# of selected cells:"),this$ncell.sel.gui,tklabel(lfinfo,text="("),this$ncell.perc.gui,tklabel(lfinfo,text="%)"),sticky="w",padx=5)
	# call frame Data Info
	tkgrid(lfinfo,columnspan=2,padx=5,sticky="we")
	tkgrid.configure(lfinfo,sticky="ew")  
}

fcs$GUIdataframe <- function(tab) {
	this=fcs
	
	# create paneled windows
	# resize the window and move the panel separator with the mouse
	this$pw2 <- tk2panedwindow(tab, orient = "horizontal")
	 
	##### add tree widget
	# create scrolling function
	xScr <- tkscrollbar(this$pw2, command=function(...)
		tkxview(this$treeDF, ...), orient = "horizontal")
	yScr <- tkscrollbar(this$pw2, command=function(...)
		tkyview(this$treeDF, ...), orient = "vertical")
		
	# create tree grid
	this$treeDF <- tkwidget(this$pw2, "Tree",
		xscrollcommand=function(...) tkset(xScr, ...),
		yscrollcommand=function(...) tkset(yScr, ...),
		width = 15, height =20,
		opencmd=function (...) { this$mytkopentree(this$treeDF,...) },
		closecmd=function (...) { this$mytkclosetree(this$treeDF,...)}
	)
	
	# call treewidget
	tkgrid(this$treeDF,yScr)
	tkgrid.configure(this$treeDF,stick="nsew")
	tkgrid.configure(yScr,stick="nsw")
	tkgrid(xScr)
	
	# insert tree content
	for (i in 1:length(this$parent.env)) {
		children=ls(eval(parse(text=paste(this$parent.env[i]))))

		for (j in 1:length(children)) {
			if ( (class(eval(parse(text=paste("this$",children[j],sep=""))))[1] == "data.frame") & (grepl("dataframes",children[j])) ) {
				# if class of children is a data.frame 
				# AND
				# if data.frame named dataframes
				df.num=as.numeric(strsplit(children[j],"[^0-9]+")[[1]][2])
				tkinsert(this$treeDF,"end","root",children[j], text=this$dataframes.name[df.num],image=imageMatrix)
			}
		}
	}
	
	# add tree widget to panel
	tkadd(this$pw2, this$treeDF)

	# bind mouse click with function
	tkbindtext(this$treeDF, "<Button-1>", function(...) { this$clickinsertDF(this$treeDF,...)})
	##### tree widget stop
	
	######## add table widget on the right
	this$treetable=ttktreeview(tab,columns=5,show="headings")
	
	# create scroll bar
	yscr <- tkscrollbar(this$treetable, command=function(...) tkyview(this$treetable,...))
	tkconfigure(this$treetable,yscrollcommand=function(...) tkset(yscr,...))
	tkpack(yscr, side="right", fill="both")
	
	xscr <- tkscrollbar(this$treetable, command=function(...) tkxview(this$treetable, ...), orient = "horizontal")
	tkconfigure(this$treetable,xscrollcommand=function(...) tkset(xscr,...))
	tkpack(xscr, side="bottom", fill="both")
	
	tkadd(this$pw2, this$treetable)
	##### tree table stop
	
	# show text paneled window
	tkpack(this$pw2, fill="both", expand=TRUE)

	# add buttons
	lfbut=ttkframe(tab)
	btn1=ttkbutton(lfbut,text="Save Rect_Data_Raw",command=function() {this$saveDF(1)},width=25)
	tk2tip(btn1,"Save as csv table.")
	btn2=ttkbutton(lfbut,text="Save Rect_Data_Transformed",command=function() {this$saveDF(2)},width=25)
	tk2tip(btn2,"Save as csv table.")
	tkpack(btn1,btn2,padx=5,pady=1,side="left")
	tkpack(lfbut,side='top',fill='y',expand=FALSE)
}

fcs$GUIlog <- function(tab) {
	this=fcs

	# set frame and content of combobox with file, mincount, binSize
	ttfmain=tkframe(tab)
	
	text <- tk2text(ttfmain,width=88, height=35)
	log=readLines("PRIanalyzer.log")
	for (i in log) {
		tkinsert(text,"end",i)
		tkinsert(text,"end","\n")
	}
	tkgrid(text)
	tkconfigure(text,state="disabled")
	tkgrid(ttfmain)
}

fcs$GUIpreproc <- function(mode) {
	this=fcs

	if (exists("ttpreproc",where=eval(parse(text="fcs")))) {tkdestroy(this$ttpreproc)}

	this$ttpreproc=tktoplevel()
	
	### mode 1 = automatic gating
	### mode 2 = normalization
	if (mode == 1) {
		mode.title="Calculating cutoffs"
		mode.size="1200x800"
	} else if (mode == 2) {
		mode.title="Normalization"
		mode.size="1100x600"
	}


	# set window title
	tkwm.title(this$ttpreproc,sprintf("Data preprocession - %s",mode.title))
	#set window size
	tkwm.geometry(this$ttpreproc,mode.size)

	# create topMenu
	topMenu <- tkmenu(this$ttpreproc)
		tkconfigure(this$ttpreproc, menu=topMenu)
		fileMenu <- tkmenu(topMenu, tearoff=FALSE)
		tkadd(topMenu, "command", label="Close",command=function(){this$closePreproc(mode)})

	### insert preprocession mode
	this$GUIpreprocMode(mode)
}

fcs$GUIpreprocMode <- function(mode) {
	this = fcs
	printf("w: GUI PreProcTab :: mode %s",mode)
	
	### mode 1 = automatic gating
	### mode 2 = normalization
	if (mode == 1) {
		# create info panel
		this$ttlfautogate=tkframe(this$ttpreproc)

		### combobox choose file
		ttfchoosefile=tkframe(this$ttlfautogate)
		this$current.filePrep=ttkcombobox(ttfchoosefile,values=this$current.filenames,width=this$max.nchar)
		tkset(this$current.filePrep,this$current.filenames[as.numeric(this$current.filenum)])
		tkbind(this$current.filePrep,"<<ComboboxSelected>>",function() {this$refreshComboboxVars(tclvalue(tkget(this$current.filePrep)))})
		tkgrid(tklabel(ttfchoosefile,text="File:  "),this$current.filePrep,sticky="w")
		tkgrid(ttfchoosefile,pady=10)

		### add checkbuttons for data preparation
		ttfcheckboxes=ttklabelframe(this$ttlfautogate,text="Calculation options")
		this$cbtremovedub=tclVar("0")
		cbtremovedub=tkcheckbutton(ttfcheckboxes,variable=this$cbtremovedub,text="Remove doublets")
		tk2tip(cbtremovedub,"Check for removing doublets.")
		this$cbttrimming=tclVar("0")
		cbttrimming=tkcheckbutton(ttfcheckboxes,variable=this$cbttrimming,text="Trimming")
		tk2tip(cbttrimming,"Check for trimming 0.02% of each marker.")
		tkgrid(cbttrimming,cbtremovedub,sticky="w",pady=1,padx=10)
		tkgrid(ttfcheckboxes,pady=10)

		### add plot/start buttons
		ttfbuttons=tkframe(this$ttlfautogate)
		bthist=tkbutton(ttfbuttons,width=25,height=3,text="Plot histograms",command=function(){this$plotHistograms(plotter="preproc",pdf=FALSE)})
		btautogate=tkbutton(ttfbuttons,width=25,height=3,text="Calculate cutoffs",command=function(){this$autoSetCutoff()})
		btinteractgate=tkbutton(ttfbuttons,width=25,height=3,text="Set cutoffs interactively",command=function(){this$interactGate()})
		tkgrid(bthist,btautogate,btinteractgate,padx=20,pady=10)
		### add buttons to select all markers or none
		len=length(this$selected.vars)
		#for ( i in 1:len ) {
		#	this$cbValprep[[i]]=tclVar("0")
		#}
		btselAll=tkbutton(ttfbuttons,text="Select all",command=function() {
			len=length(this$selected.vars)
			for ( i in 1:len ) {
				tclvalue(this$cbVal[[i]]) = "1"
			}
		})
		btdeselAll=tkbutton(ttfbuttons,text="Deselect all",command=function() {
			len=length(this$selected.vars)
			for ( i in 1:len ) {
				tclvalue(this$cbVal[[i]]) = "0"
			}
		})
		tkgrid(btselAll,btdeselAll)
		tkgrid(ttfbuttons)

		### add markers
		this$ttlfselMarker=ttklabelframe(this$ttlfautogate,text="Select marker:")
		ttfcol1=tkframe(this$ttlfselMarker)
		ttfcol2=tkframe(this$ttlfselMarker)
		ttfcol3=tkframe(this$ttlfselMarker)
		ttfcol4=tkframe(this$ttlfselMarker)
		residue = 0
		if (len%%4>0) residue=1
		quarter=floor(1*len/4)+ residue
		for ( i in 1:quarter ) {
			cb=tkcheckbutton(ttfcol1,variable=this$cbVal[[i]],text=this$selected.vars[i])
			cutoffentry=tkentry(ttfcol1,width=6,textvariable=this$vcutoffs[[i]])
			tkgrid(cb,tklabel(ttfcol1,text="cutoff: "),cutoffentry,padx=2,pady=1,sticky="w")
		}
		for ( j in (i+1):(i+quarter) ) {
			cb=tkcheckbutton(ttfcol2,variable=this$cbVal[[j]],text=this$selected.vars[j])
			cutoffentry=tkentry(ttfcol2,width=6,textvariable=this$vcutoffs[[j]])
			tkgrid(cb,tklabel(ttfcol2,text="cutoff: "),cutoffentry,padx=2,pady=1,sticky="w")
		}
		for ( k in (j+1):(j+quarter) ) {
			cb=tkcheckbutton(ttfcol3,variable=this$cbVal[[k]],text=this$selected.vars[k])
			cutoffentry=tkentry(ttfcol3,width=6,textvariable=this$vcutoffs[[k]])
			tkgrid(cb,tklabel(ttfcol3,text="cutoff: "),cutoffentry,padx=2,pady=1,sticky="w")
		}

		for ( l in (k+1):len ) {
			cb=tkcheckbutton(ttfcol4,variable=this$cbVal[[l]],text=this$selected.vars[l])
			cutoffentry=tkentry(ttfcol4,width=6,textvariable=this$vcutoffs[[l]])
			tkgrid(cb,tklabel(ttfcol4,text="cutoff: "),cutoffentry,padx=2,pady=1,sticky="w")
		}
		tkgrid(ttfcol1,ttfcol2,ttfcol3,ttfcol4,padx=20,pady=1,sticky="we")
		tkgrid(this$ttlfselMarker,padx=10)
		
		tkgrid(this$ttlfautogate,pady=1)
	} else if (mode == 2) {
		# create info panel
		this$ttlfnormalize=tkframe(this$ttpreproc)

		### add buttons select/deselect all
		ttfbuttons=tkframe(this$ttlfnormalize)
		# add start normalization button
		btnormalize=tkbutton(ttfbuttons,height=3,text="Start normalization",command=function(){this$normalizeFiles()})
		tkgrid(btnormalize,pady=10,sticky="snew")

		len=length(this$current.filenames)
		for ( i in 1:len ) {
			this$cbValFile[[i]]=tclVar("0")
		}
		btselAll=tkbutton(ttfbuttons,text="Select all",command=function() {
			for ( i in 1:len ) {
				tclvalue(this$cbValprepFile[[i]]) = "1"
			}
		})
		btdeselAll=tkbutton(ttfbuttons,text="Deselect all",command=function() {
				for ( i in 1:len ) {
						tclvalue(this$cbValprepFile[[i]]) = "0"
				}
		})
		tkgrid(btselAll,btdeselAll,padx=20,pady=1)
		tkgrid(ttfbuttons)

		### add files to select
		ttfselFile=ttklabelframe(this$ttlfnormalize,text="Select files:")
		ttfcol1=tkframe(ttfselFile)
		ttfcol2=tkframe(ttfselFile)
		ttfcol3=tkframe(ttfselFile)
		for ( i in 1:(floor(len/3)+(len%%3)) ) {
			cb=tkcheckbutton(ttfcol1,variable=this$cbValprepFile[[i]],text=this$current.filenames[i])
			tkgrid(cb,padx=2,pady=1,sticky="w")
		}
		for ( j in (i+1):(i+floor(len/3)) ) {
			cb=tkcheckbutton(ttfcol2,variable=this$cbValprepFile[[j]],text=this$current.filenames[j])
			tkgrid(cb,padx=2,pady=1,sticky="w")
		}
		for ( k in (j+1):len ) {
			cb=tkcheckbutton(ttfcol3,variable=this$cbValprepFile[[k]],text=this$current.filenames[k])
			tkgrid(cb,padx=2,pady=1,sticky="w")
		}
		tkgrid(ttfcol1,ttfcol2,ttfcol3,padx=20,pady=1,sticky="wen")
		tkgrid(ttfselFile,pady=10,padx=5)


		tkgrid(this$ttlfnormalize)
	}   
}

fcs$refreshPanel <- function() {
	this=fcs

	len=length(this$selected.cutoffs)

	#### refresh checkbox markers version 3
	scrollFrame <- tkwidget(this$sw,"ScrollableFrame",height=630,width=this$width.panel)
	tcl(this$sw,"setwidget",scrollFrame)
	subfID <- tclvalue(tcl(scrollFrame,"getframe"))
	# label project name
	lab <- tcl("label",paste(subfID,".lab",sep=""),text=sprintf("%s",this$shortenFilename(this$current.project)))
	tkgrid(lab,sticky="w",padx=5)

	i = j = 1
	while (i < (len*4) ) {
		cutoffcb=tcl("checkbutton",paste(subfID,".",i,sep=""),variable=this$cbVal[[j]],text=this$selected.vars[j])
		cutofflabel=tcl("label",paste(subfID,".",i+1,sep=""),text="cutoff: ")
		cutoffentry=tcl("entry",paste(subfID,".",i+2,sep=""),width=6,textvariable=this$vcutoffs[[j]])
		cutoffperccb=tcl("checkbutton",paste(subfID,".",i+3,sep=""),variable=this$cbcutoffperc[[j]],text="in %")
		tkgrid(cutoffcb,cutofflabel,cutoffentry,cutoffperccb,sticky="w",pady=1)
		i = i + 4
		j = j + 1
	}
}

fcs$openProject <- function() {
	this=fcs

	printf("w: do openProject")

	this$ttprojects=tktoplevel()
	this$ttprojectsopen = TRUE
	
	# set window title
	tkwm.title(this$ttprojects,"Project List")
	#set window size
	tkwm.geometry(this$ttprojects,"400x600")

	##### add table widget
	listProjects=ttktreeview(this$ttprojects,columns=1,show="headings")
	
	# create scroll bar
	yscr <- tkscrollbar(listProjects, command=function(...) tkyview(listProjects,...))
	tkconfigure(listProjects,yscrollcommand=function(...) tkset(yscr,...))
	tkpack(yscr, side="right", fill="both")
	xscr <- tkscrollbar(listProjects, command=function(...) tkxview(listProjects, ...), orient = "horizontal")
	tkconfigure(listProjects,xscrollcommand=function(...) tkset(xscr,...))
	tkpack(xscr, side="bottom", fill="both")

	shade = c("none","gray")
	idx.file=grep("fileIdentity|fileIndex",this$dataframes.name)
	tcl(listProjects,"column",1,width=this$max.nchar)
	for (i in 1:length(this$total.projects)){
		size=nrow(this$getDFtable(this$dataframes.name[idx.file[i]]))
		tkinsert(listProjects,"","end",values=as.tclObj(paste(this$total.projects[i],"(",size,")",sep="")),tag=shade[i%%2+1])
	}
	#tcl(this$treetable,"configure",background="gray95")
	tktag.configure(this$treetable,"gray",background="gray95")
	tkpack(listProjects,fill="both",expand=TRUE)
	tk2tip(listProjects,"Double Click to Open")

	# get row index
	#this$projectChildren=tcl(listProjects,"children","")
	#print(sapply(as.character(this$projectChildren), function(i) tclvalue(tkindex(listProjects,i))))

	# double click on project    
	tkbind(listProjects, "<Double-Button-1>", function(W,x,y) {
		sel= as.character(tcl(W, "identify","row",x,y))
		this$selected.project = strsplit(tclvalue(tcl(W, "item", sel, "-values")),"\\(")[[1]][1]
		project.num = which(this$total.projects==this$selected.project)
		#if ( project.num != this$current.projectnum) {
		#	print(paste("Project selected:",this$selected.project))
			this$refreshComboboxStart(project.num)
		#}
		tkwm.withdraw(tkwinfo("toplevel",W))
	})
}

fcs$openDB <- function() {
	this=fcs

	printf("w: do openDB")
	file <- tclvalue(tkgetOpenFile(initialdir=fcs$db.path,defaultextension="sqlite3"))
	this$exit()
	this$connectDb(file)
	this$db.name=file
	this$GUImain()
}

fcs$refreshComboboxStart <- function(num) {
	# refresh File Combobox
	this = fcs

	printf("w: do refreshComboboxStart with project #%s",num)


	this$selected.projectnum = num
	this$selected.project = this$total.projects[num]

	### ask to save cutoffs
	new = unlist(this$new.cutoffs)
	new = new[!is.na(new)]
	saved = unlist(this$saved.cutoffs)
	saved = saved[!is.na(saved)]
	both <- new == saved

	new2 = unlist(this$new.checks)
	new2 = new2[!is.na(new2)]
	saved2 = unlist(this$saved.checks)
	saved2 = saved2[!is.na(saved2)]
	both2 <- new2 == saved2

	if ( !all(both) & !all(both2)) {
		answer=tkmessageBox(title = "Cutoffs..", message = "Save cutoffs first?", icon = "info", type = "yesno")
		if (tclvalue(answer) == "yes") this$saveCutoffsToDB()
	}

	# get filenames for table
	this$current.filetable = this$getDFtable(paste0(this$selected.project,this$fileid_name))
	this$current.filenames = this$current.filetable[,2]
	if (this$working) {
		names=this$current.filenames
		if (length(names)>5) names = names[1:5]
		print("w: First files in project:")
		print(names)
	}
	#this$current.filenames2 = this$current.filenames
	this$current.staintable = this$getDFtable(paste0(this$selected.project,this$markerid_name))
	this$selected.filenum = this$current.filetable[which(this$current.filetable[,2]==this$current.filenames[1]),1]

	### load cutoffs
	if ( ncol(this$current.staintable)==6 ) {
		### if cutoffs were already saved
		# get total cutoffs
		for ( file.id in unique(this$current.filetable[,1]) ) {
			this$new.cutoffs[[file.id]] = this$current.staintable[which(this$current.staintable[,1]==file.id),5]
			this$new.checks[[file.id]] = this$current.staintable[which(this$current.staintable[,1]==file.id),6]
		}
		this$saved.cutoffs = this$new.cutoffs
		this$saved.checks = this$new.checks
	} else if ( ncol(this$current.staintable)==5 ) {
		# get total cutoffs
		for ( file.id in unique(this$current.filetable[,1]) ) {
			this$new.cutoffs[[file.id]] = this$current.staintable[which(this$current.staintable[,1]==file.id),5]
			this$new.checks[[file.id]] = rep("0",nrow(this$current.staintable[which(this$current.staintable[,1]==file.id),]))
		}
		this$saved.cutoffs = this$new.cutoffs
		this$saved.checks = list()
	} else {
		print("Initiate new cutoff list")
		### no cutoffs listed in this project
		# initiate total cutoff list
		for ( file.id in unique(this$current.filetable[,1]) ) {
			this$new.cutoffs[[file.id]] = rep("0",nrow(this$current.staintable[which(this$current.staintable[,1]==file.id),]))
			this$new.checks[[file.id]] = rep("0",nrow(this$current.staintable[which(this$current.staintable[,1]==file.id),]))
		}
		this$saved.cutoffs = list()
		this$saved.checks = list()
	}


	tkconfigure(this$tkchoosefile,values=this$current.filenames)
	tkset(this$tkchoosefile,this$current.filenames[1])


	### refresh filenames in title
	displayfile = this$shortenFilename(this$current.filenames[1],title=TRUE)
	tkdelete(this$title,"1.0","end")
	tkinsert(this$title,"1.0",displayfile)

	### set variables from project in tab di tri and tetraplot
	this$selected.vars = this$getVariables(table=this$selected.projectnum,index=this$selected.filenum)

	### if FSC/SSC are listed
	# set vars after the FSC/SSCs
	tmp = max(grep("SC",this$selected.vars)) 
	tmp2 = min(grep(c("Time"),this$selected.vars))
	if (tmp2>tmp & !is.infinite(tmp2) & tmp2<(length(this$selected.vars)-3)) tmp=tmp2
	if ( ((tmp+3)>length(this$selected.vars)) | is.infinite(tmp)) tmp = 0

	tkconfigure(this$cbvar1,values=this$selected.vars)
	if (length(param$currVarTri1) == 0) param$currVarTri1 = 1 + tmp
	else if (param$currVarTri1<tmp) param$currVarTri1 = 1 + tmp
	tkset(this$cbvar1,this$selected.vars[param$currVarTri1])

	tkconfigure(this$cbvar2,values=this$selected.vars)
	if (length(param$currVarTri2) == 0) param$currVarTri2 = 2 + tmp
	else if (param$currVarTri2<tmp) param$currVarTri2 = 2 + tmp
	tkset(this$cbvar2,this$selected.vars[param$currVarTri2])

	tkconfigure(this$cbvar3,values=this$selected.vars)
	if (length(param$currVarTri3) == 0) param$currVarTri3 = 3 + tmp
	else if(param$currVarTri3<tmp) param$currVarTri3 = 3 + tmp
	tkset(this$cbvar3,this$selected.vars[param$currVarTri3])

	tkconfigure(this$cbvar1di,values=this$selected.vars)
	if (length(param$currVarDi1) == 0) param$currVarDi1 = 1 + tmp
	else if (param$currVarDi1<tmp) param$currVarDi1 = 1 + tmp
	tkset(this$cbvar1di,this$selected.vars[param$currVarDi1])

	tkconfigure(this$cbvar2di,values=this$selected.vars)
	if (length(param$currVarDi2) == 0) param$currVarDi2 = 2 + tmp
	else if (param$currVarDi2<tmp) param$currVarDi2 = 2 + tmp
	tkset(this$cbvar2di,this$selected.vars[param$currVarDi2])

	# renew markers and marker position
	this$refreshMarker(openProject=TRUE)
	
	### project chosen, now no function if you click on main window
	this$openProjectYet = TRUE
	tkbind(this$tt, "<Button-1>", function(...) {})
}

fcs$refreshMarker <- function(tab=FALSE,openProject=FALSE) {
	# renew checkboxes "markers" in tab "Overview"
	this = fcs

	len=length(this$selected.vars)
	this$max.nchar.vars = max(nchar(this$selected.vars))

	if (tab) {
		### tab 0 = diploT
		### tab 1 = triploT

		### tab numbering starts with "0"
		tab.old=this$current.tab
		tab.new=tclvalue(tkindex(this$panelTabs,"current"))

		printf("w: do refreshMarker: tab=%s tab.new=%s tab.old=%s",tab,tab.new,tab.old)

		### set x-/y-/z-features the same within the other tabs
		if (tab.old == "0" & exists("cbvar1",env=fcs) ) {
			### tab 0 = diploT

			### set combobox
			var1 = this$checkMarker(tclvalue(tkget(this$cbvar1di)))
			tkset(this$cbvar1,var1)

			var2 = this$checkMarker(tclvalue(tkget(this$cbvar2di)))
			tkset(this$cbvar2,var2)

			### if z variable is the same as in x and y, then go one var up
			# in other words:
			# if cbvar3 == cbvar1/2, then set cbvar3 randomly but different
			idx=which(tclvalue(tkget(this$cbvar1))==this$selected.vars)
			idy=which(tclvalue(tkget(this$cbvar2))==this$selected.vars)
			idz=which(tclvalue(tkget(this$cbvar3))==this$selected.vars)

			if ( length(which(idz==c(idx,idy))) > 0 | length(idz) == 0 ) {
				### if var index z is identical to var index x or y
				# take another random var index for z
				t = 1:len
				t = t[-c(idx,idy)]
				random = sample(t,1)
				tkset(this$cbvar3,this$selected.vars[random])
				printf("w: changed z-feature from %s to %s",this$selected.vars[idz],this$selected.vars[random])
			}

			### initiate
			this$coords.info = vector()
		} else if (tab.old == "1") {
			### tab 1 = triploT

			### set combobox
			var1 = this$checkMarker(tclvalue(tkget(this$cbvar1)))
			tkset(this$cbvar1di,var1)

			var2 = this$checkMarker(tclvalue(tkget(this$cbvar2)))
			tkset(this$cbvar2di,var2)
		} 

	} else {

		### initiate
		this$idx=list()

		printf("w: do refreshMarker: len(vars)=%s fileidx=%s tab=%s openProject=%s",len,this$selected.filenum,tab,openProject)

		### get var names before the "." and make them capitals
		vars.old = this$shortenMarkername(this$current.vars)
		vars.new = this$shortenMarkername(this$selected.vars)

		### go through all new selected.vars
		# if any of selected.vars is in current.vars
		# save check and cutoff
		# i : new position of var
		# idx[[i]][1] : new position of var
		# idx[[i]][2] : var checked or not ("0"/"1")
		# idx[[i]][3] : cutoff value of var
		for ( i in 1:len ) {
			this$idx[[i]]=tmp=which(vars.new[i]==vars.old)
			if ( length(tmp) > 0 ) {
				this$idx[[i]][2]=tclvalue(this$cbVal[[tmp[1]]])
				this$idx[[i]][3]=tclvalue(this$vcutoffs[[tmp[1]]])
			}
		}

		# save cutoffs from old file
		if (!openProject) this$refreshCutoffs()
		this$refreshCutoffs(saved=TRUE)
		len=length(this$selected.vars)

		#### refresh checkbox markers version 3
		scrollFrame <- tkwidget(this$sw,"ScrollableFrame",height=630,width=this$width.panel)
		tcl(this$sw,"setwidget",scrollFrame)
		subfID <- tclvalue(tcl(scrollFrame,"getframe"))
		# label project name
		lab <- tcl("label",paste(subfID,".lab",sep=""),text=sprintf("%s",this$shortenFilename(this$selected.project)))
		tkgrid(lab,sticky="w",padx=5)

		i = j = 1
		while (i < (len*4) ) {
			#printf("j=%s vars=%s cutoffs=%s subfID=%s",j,this$selected.vars[j],this$selected.cutoffs[j],paste(subfID,".",i+3,sep=""))
			#tclvalue(this$cbVal[[j]]) = "0"
			#tclvalue(this$vcutoffs[[j]]) = this$selected.cutoffs[j]
			#tclvalue(this$cbcutoffperc[[j]]) = "0"
				
			cutoffcb=tcl("checkbutton",paste(subfID,".",i,sep=""),variable=this$cbVal[[j]],text=this$selected.vars[j])
			cutofflabel=tcl("label",paste(subfID,".",i+1,sep=""),text="cutoff: ")
			cutoffentry=tcl("entry",paste(subfID,".",i+2,sep=""),width=6,textvariable=this$vcutoffs[[j]])
			cutoffperccb=tcl("checkbutton",paste(subfID,".",i+3,sep=""),variable=this$cbcutoffperc[[j]],text="in %")
			tkgrid(cutoffcb,cutofflabel,cutoffentry,cutoffperccb,sticky="w",pady=1)
			i = i + 4
			j = j + 1
		}
		#tkbind(cutoffcb,"<FocusIn>",function() tcl(scrollFrame,"see",cutoffcb) )
		#tkgrid(this$sw)
		###

		### then go through idx list and set check/cutoff 
		for ( i in 1:len) {
			### if one of the var where in the old file
			# and cutoff was not zero
			if ( length(this$idx[[i]])>0 & as.numeric(this$idx[[i]][3]) != 0 ) {
				position = as.numeric(this$idx[[i]][1])
				#if (position != i) {
				tclvalue(this$cbVal[[i]])=this$idx[[i]][2]
				tclvalue(this$vcutoffs[[i]])=this$idx[[i]][3]
				if (this$working & tclvalue(this$vcutoffs[[i]])!="0" ) {
					printf("w: changed var=%s from position %s to %s with value=%s, checked=%s",
						this$selected.vars[i],position,i,tclvalue(this$vcutoffs[[position]]),tclvalue(this$cbVal[[position]])
					)
				}
				#}
				
			}
		}

		#tkgrid(this$sw)
	}
	printf("w: Done refreshMarker")
}

fcs$refreshCutoffs <- function (fileindex=NA,current=FALSE,reset=FALSE,saved=FALSE) {
		this=fcs

		if (is.na(fileindex)) {
			#if (exists("tkchoosefile",envir=fcs)) {
			#	file=tclvalue(tkget(this$tkchoosefile))
			#} else {    
			#	file=this$current.filenames[1]
			#}
			#fileindex=this$current.filetable[which(this$current.filetable[,2]==file),1]
			fileindex=this$selected.filenum
		}
		len = length(this$selected.cutoffs)
		len = length(this$selected.vars)

		printf("w: do refreshCutoffs: table=%s fileidx=%s current=%s reset=%s load_saved=%s len=%s",this$selected.project,fileindex,current,reset,saved,len)

		if (reset) {
			## rest cutoffs in GUI and in array new.cutoffs
			this$selected.cutoffs = rep(0,len)
			this$new.cutoffs[[this$selected.filenum]] = this$selected.cutoffs
			this$selected.checks = rep("0",len)
			this$new.checks[[this$selected.filenum]] = this$selected.checks
		} else if (current) {
			### save current cutoffs in array new.cutoffs
			for (i in 1:len) {
				this$selected.cutoffs[i] = this$checkDigits(cutoff_id=i)
				this$selected.checks[i] = tclvalue(this$cbVal[[i]])
			}
		  
		  
			this$new.cutoffs[[this$current.filenum]] = this$selected.cutoffs
			this$new.checks[[this$current.filenum]] = this$selected.checks
		} else if (saved) {
			if (length(this$saved.cutoffs)>0) {
			  ### load saved cutoffs from database if available
			  if (any(this$saved.cutoffs[[fileindex]]!="0")) {
  			  ### AND
  			  ### if saved cutoffs arent only zeros "0"
  				this$selected.cutoffs = this$saved.cutoffs[[fileindex]]
  				this$selected.cutoffs[is.na(this$selected.cutoffs)] = 0
  				this$selected.checks = this$saved.checks[[fileindex]]
  				this$selected.checks[is.na(this$selected.checks)] = "0"
			  }
			} else {
			  printf("No cutoffs in database saved.")
			  if (any(this$new.cutoffs[[fileindex]]!="0")) {
			    printf("Loading cutoffs from current session.")
  				this$selected.cutoffs = this$new.cutoffs[[fileindex]]
  				this$current.checks = this$new.checks[[fileindex]]
			  }
			}
		} else {
			len = length(this$selected.vars)
			### first save cutoffs in recent file 
			for ( i in 1:len ) {
				this$current.cutoffs[i] = this$checkDigits(cutoff_id=i)
				this$current.checks[i] = tclvalue(this$cbVal[[i]])
			}
			this$new.cutoffs[[this$current.filenum]] = this$current.cutoffs
			this$new.checks[[this$current.filenum]] = this$current.checks

			### then load cutoffs from current chosen file
			if (any(this$new.cutoffs[[fileindex]]!="0")) {
  			this$selected.cutoffs = this$new.cutoffs[[fileindex]]
  			this$selected.checks = this$new.checks[[fileindex]]
			}
		}

		### if present cutoff length is longer than recent 
		# initiate additional positions 
		len = length(this$selected.vars)
		if (len > length(this$vcutoffs)) {
			for ( i in len:length(this$vcutoffs)) {
			  this$cbVal[[i]] = tclVar("0")
				this$vcutoffs[[i]] = tclVar(this$selected.cutoffs[i])
				this$cbcutoffperc[[i]] = tclVar("0")
			}
		}

		### display cutoffs
		for ( i in 1:len ) {
			tclvalue(this$vcutoffs[[i]]) = as.character(this$selected.cutoffs[i])
			tclvalue(this$cbVal[[i]]) = this$selected.checks[i]
		}

		printf("w: Done refreshCutoffs: %s cutoffs:%s",len,paste(this$selected.cutoffs,collapse=" "))
		this$current.cutoffs = this$selected.cutoffs
		this$current.checks = this$selected.checks
}

fcs$refreshComboboxVars <- function(file) {
	# refresh var1 var2 var3 combobox in tabs
	this = fcs

	displayfile = this$shortenFilename(file,title=TRUE)

	### refresh filenames in title
	tkdelete(this$title,"1.0","end")
	tkinsert(this$title,"1.0",displayfile)

	if ( grepl("^temp\\d\\d", file) ) {
		fileindex = 1
	} else {
		fileindex = this$current.filetable[which(this$current.filetable[,2]==file),1]
		#this$selected.project = this$total.projects[this$selected.projectnum]
		#this$selected.filenum = fileindex
	}



	# if stain index changed set stain variables and refresh
	if (fileindex != this$current.filenum) {

		printf("w: do refreshComboboxVars: fileidx=%s file=%s",fileindex,file)
	  
	  this$selected.filenum = fileindex
		#this$current.vars = this$selected.vars
		this$selected.vars = this$getVariables(index=fileindex)
		#this$current.cutoffs = this$selected.cutoffs

		#this$recent.fileindex = this$selected.filenum
		#this$selected.filenum = fileindex
  
  	# if FSC/SSC are listed
  	tmp = max(grep("SC",this$selected.vars)) 
  	tmp2 = min(grep(c("Time"),this$selected.vars))
  	if (tmp2>tmp & !is.infinite(tmp2) & tmp2<(length(this$selected.vars)-3)) tmp=tmp2
  	if ( ((tmp+3)>length(this$selected.vars)) | is.infinite(tmp)) tmp = 0

  
  	######## triploT
  	#### compare with old vars: if match then let it put
  	var1 = this$checkMarker(tclvalue(tkget(this$cbvar1)))
  	## if var1 is not available in new file
  	if ( length(var1)==0 ) {
  		var1 = this$selected.vars[tmp+1]
  		printf("w: var1-tri changed to %s",var1)
  	}
  	tkconfigure(this$cbvar1,values=this$selected.vars)
  	tkset(this$cbvar1,var1)
  	var2 = this$checkMarker(tclvalue(tkget(this$cbvar2)))
  	## if var2 is not available in new file
  	if ( length(var2)==0 ) {
  		var2 = this$selected.vars[tmp+2]
  		printf("w: var2-tri changed to %s",var2)
  	}
  	tkconfigure(this$cbvar2,values=this$selected.vars)
  	tkset(this$cbvar2,var2)
  	var3 = this$checkMarker(tclvalue(tkget(this$cbvar3)))
  	## if var3 is not available in new file
  	if ( length(var3)==0 ) {
  		var3 = this$selected.vars[tmp+3]
  		printf("w: var3-tri changed to %s",var3)
  	}
  	tkconfigure(this$cbvar3,values=this$selected.vars)
  	tkset(this$cbvar3,var3)
  
  
  	######## diploT
  	#### compare with old vars: if match then let it put
  	var1 = this$checkMarker(tclvalue(tkget(this$cbvar1di)))
  	## if var1 is not available in new file
  	if ( length(var1)==0 ) {
  		var1 = this$selected.vars[tmp+1]
  		printf("w: var1-di changed to %s",var1)
  	}
  	tkconfigure(this$cbvar1di,values=this$selected.vars)
  	tkset(this$cbvar1di,var1)
  	## if var2 is not available in new file
  	if ( length(var2)==0 ) {
  		var2 = this$selected.vars[tmp+2]
  		printf("w: var2-di changed to %s",var2)
  	}
  	var2 = this$checkMarker(tclvalue(tkget(this$cbvar2di)))
  	tkconfigure(this$cbvar2di,values=this$selected.vars)
  	tkset(this$cbvar2di,var2)
  			
  
  	this$refreshMarker()
  
  	# if change file, undo "gate data"
  	# so that temporary data is not set
  	tclvalue(this$cbtgateData) = "0"
  
  	printf("w: Done refreshComboboxVars")
	}
}

fcs$closePreproc <- function(mode) {
	this=fcs

	### ask for apply procession
	# refresh files and marker according to procession window
	if (mode==1) this$refreshComboboxVars(tclvalue(tkget(this$current.filePrep)))
	# close window
	tkdestroy(this$ttpreproc)
}

fcs$plotDensities <- function(plotter,pdf) {this=fcs

	this$refreshPlotters()

	if (!pdf) tkconfigure(this$tt, cursor = "watch")

	printf("w: do plotDensities: plotter=%s, pdf=%s", plotter, pdf)
	len=length(this$selected.vars)
	dev.label=paste("density",plotter,sep=".")

	if ( plotter =="tri") {
		# mode plotter == triploT
		dev.label=paste("plotter","tri",this$plotter.tri.num,sep=".")
		
		if ( length(which(dev.label==names(devList()))) == 0 ) {            
			ncol=as.numeric(tclvalue(this$vncol))
			nrow=as.numeric(tclvalue(this$vnrow))
			this$plotter.tri.num = this$plotter.tri.num + 1
			dev.label=paste("plotter","tri",this$plotter.tri.num,sep=".")
			devNew(type="x11",title="Densities for triploTs",width=ncol*4,height=nrow*4,label=dev.label)
			# mar in points, mai in inches
			# oma adds title lines
			# order: bottom, left, top, and right
			par(mfrow=c(nrow,ncol),oma=c(1,1,3,1),mar=c(3,3,5,2))
			this$plot.windows = c(this$plot.windows,dev.label)
		} else {
			devSet(devList()[which(dev.label==names(devList()))])
		}

		file=tclvalue(tkget(this$tkchoosefile))
		displayfile = this$shortenFilename(file)
		var1 = this$checkMarker(tclvalue(tkget(this$cbvar1)))
		var2 = this$checkMarker(tclvalue(tkget(this$cbvar2)))
		tkset(this$cbvar2,var2)
		vars = c(var1,var2)

		### if features are not in sample
		if ( length(var1)==0 ){
			tkmessageBox(title = "An error has occured!",
				message = "Check your features.")
			stop(sprintf("Feature A=%s is not existent in that file.",tclvalue(tkget(this$cbvar1))))
		} else {
			tkset(this$cbvar1,var1)
		}
		if ( length(var2)==0 ){
			tkmessageBox(title = "An error has occured!",
				message = "Check your features.")
			stop(sprintf("Feature B=%s is not existent in that file.",tclvalue(tkget(this$cbvar2))))
		} else {
			tkset(this$cbvar2,var2)
		}


		var1.idx = which(this$selected.vars==var1)
		var2.idx = which(this$selected.vars==var2)

		cutoff1 = this$checkDigits(cutoff_id=var1.idx)
		cutoff2 = this$checkDigits(cutoff_id=var2.idx)
		cutoffs = c(cutoff1,cutoff2)
		colvec = c(var1.idx,var2.idx)
		len.colvec = length(colvec)

		cbperc1 = tclvalue(this$cbcutoffperc[[var1.idx]])
		cbperc2 = tclvalue(this$cbcutoffperc[[var2.idx]])
		cbperc = c(cbperc1,cbperc2)
	} else if (plotter == "triover") {
		# mode triploT overview
		file=tclvalue(tkget(this$tkchoosefile))

		if ( this$OverviewGate ) {
			displayfile = this$shortenFilename(this$plot.attr[[1]]$file.name)}
		else {
			displayfile = this$shortenFilename(file)
		}

		# column and cutoff number vector
		vars = vector()
		colvec = vector()
		cbperc = vector()
		cutoffs = vector()
		for (i in 1:len) {
			if ( tclvalue(this$cbVal[[i]]) == "1") {
				colvec = c(colvec,i);
				cutoffs = c(cutoffs,this$checkDigits(cutoff_id=i))
				cbperc = c(cbperc,tclvalue(this$cbcutoffperc[[i]]))
				vars = c(vars,this$selected.vars[i])
			}
		}
		len.colvec = length(colvec)

		if ( len.colvec < 1 ) {
			tkmessageBox(title = "An error has occured!",
					message = "Please select at least one markers.", icon = "error", type = "ok")
			stop("Missing values")
		}

		if (pdf) {
			print(paste0("Creating page 2/",len.colvec+2,".."))

		  			if (len.colvec > 8) par(mfrow=c(len.colvec-1,len.colvec-2),oma=c(2,1,3,1),mar=c(3.5,4,5,2.5))
			else par(mfrow=c(8,6),oma=c(2,1,3,1),mar=c(3.5,4,5,2.5))
		} else if ( length(which(dev.label==names(devList()))) == 0 ) {
			devNew(type="x11",title="Densities for triploTs overview",width=3*4,height=3*3.8,label=dev.label)
			par(mfrow=c(3,3),oma=c(1,1,3,1),mar=c(3,3,5,1))
			#devSet(grep(dev.label,names(devList()))+1)
			this$plot.windows = c(this$plot.windows,dev.label)
		} else {
			#devSet(grep(dev.label,names(devList()))+1)
			devSet(devList()[which(dev.label==names(devList()))])
		}
	}

	if (pdf & this$OverviewGate) {
		table=this$temptable.name[this$temp.num]
		file.idx = this$selected.filenum
	} else if ( grepl("^temp\\d\\d", file) ) {
		table = file
		file.idx = 1
	} else {
		table = this$current.project
		file.idx = this$current.filetable[which(this$current.filetable[,2]==file),1]
	}

	if (pdf & this$OverviewGate) table=this$temptable.name[this$temp.num]

	timeSTART = Sys.time()
	if ( is.null(this$data) | this$current.project != table | this$current.filenum != file.idx |
	     this$current.cofactor != as.numeric(tclvalue(this$rbasinh)) ) {
		if ( this$working ) print("Time loading data:")
		this$getFile(table,file.idx)
		if ( this$working ) print(Sys.time()-timeSTART)
	}

	checkTRANS = tclvalue(this$rbtrans)
	checkGATED = tclvalue(this$cbtgateData)
	checkGRID = tclvalue(this$cbtshowGrid)
	checkTRIMMING = tclvalue(this$cbttrimming)
	checkCALC = "density"

	xminval = as.numeric(tkget(this$minvalX))
	xmaxval = as.numeric(tkget(this$maxvalX))
	yminval = as.numeric(tkget(this$minvalY))
	ymaxval = as.numeric(tkget(this$maxvalY))

	data = this$data[colvec]

	if ( checkTRANS == "asinh" ) {
		scale = this$asinh$scale
		label = this$asinh$label
		grid.step = this$asinh$step
	} else {
		scale = this$biex$scale
		label = this$biex$label
		grid.step = this$biex$step
	} 

	### check cutoffs
	for ( i in 1:len.colvec ) {
		if ( cbperc[i] == "1" ) {
			tdata=data[,i]
			cutoffs[i] = this$calcCUTOFF(tdata,cutoffs[i],vars[i],colvec[i])
		}
	}

	binSize=as.numeric(tkget(this$binSize))
	mincount=as.numeric(tkget(this$minCountTri))
	#cutoffs[3]=0
	t = 0

	for ( v1 in 1:len.colvec ) {
		for ( v2 in (v1+1):len.colvec ) {

			## skip loop if one of the axes are the same
			if ( v1 == v2 | v2>len.colvec) next
			

			#### big start
			set.cex=1.1
			set.cex.axes=1.0
			set.mgp=c(1.9, 0.7, 0)
			if ( cutoffs[1] > 0 ) title.axis = sprintf("%s (%s)",v1,cutoffs[1])
			else title.axis = v1
			if ( cutoffs[2] > 0 ) title.axis = c(title.axis,sprintf("%s (%s)",v2,cutoffs[2]))
			else title.axis = c(title.axis,v2)

			plot(1,type='n',frame.plot=FALSE,xlim=c(xminval,xmaxval+10*binSize),axes=FALSE,
				ylim=c(yminval-2.5*binSize,ymaxval+5*binSize),xlab=title.axis[1],ylab=title.axis[2],cex.lab=set.cex,cex.axis=0.5*set.cex.axes,mgp=set.mgp)
			box(lwd=0.8,col="darkgrey")

			### draw axis on the bottom and on the left
			axis(side=1, at=scale,labels=labels,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
			axis(side=2, at=scale,labels=labels,las=3,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")

			### add grid
			if(checkGRID=="1") {
			  xgrid.steps=seq(0,(xmaxval),by=grid.step)
			  ygrid.steps=seq(0,(ymaxval),by=grid.step)
			  abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
			}

			#select columns to plot
			tdata=as.matrix(data[,c(v1,v2,v2)])

			### calc quadrants in total
			ncells = nrow(tdata)
			# q1 Quadrant unten links
			# q2 Quadrant unten rechts
			# q3 Quadrant oben rechts
			# q4 Quadrant oben links
			if ( cutoffs[v1] > 0 & cutoffs[v2] > 0 ) {
				tdata.q1 = tdata[which( tdata[,1]<cutoffs[v1] &  tdata[,2]<cutoffs[v2] ),3]
				tdata.q2 = tdata[which( tdata[,1]>=cutoffs[v1] &  tdata[,2]<cutoffs[v2] ),3]
				tdata.q3 = tdata[which( tdata[,1]>=cutoffs[v1] &  tdata[,2]>=cutoffs[v2] ),3]
				tdata.q4 = tdata[which( tdata[,1]<cutoffs[v1] &  tdata[,2]>=cutoffs[v2] ),3]

				this$q1.total = abs(100 * length( tdata.q1 ) / ncells)
				this$q2.total = abs(100 * length( tdata.q2 ) / ncells)
				this$q3.total = abs(100 * length( tdata.q3 ) / ncells)
				this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
			}			

			# get only the cells which are greater than 0
			tdata.zero = tdata[ tdata[,1]>=0 & tdata[,2]>=0, ]
			ncells.zero = nrow(tdata.zero)

			### calc quadrants with only positive values
			# q1 Quadrant unten links
			# q2 Quadrant unten rechts
			# q3 Quadrant oben rechts
			# q4 Quadrant oben links
			if ( cutoffs[v1] > 0 & cutoffs[v2] > 0 ) {
				q1.zero = abs(100 * length( which( tdata.zero[,1]<cutoffs[v1] &  tdata.zero[,2]<cutoffs[v2] ) ) / ncells.zero)      
				q2.zero = abs(100 * length( which( tdata.zero[,1]>=cutoffs[v1] &  tdata.zero[,2]<cutoffs[v2] ) ) / ncells.zero)                    
				q3.zero = abs(100 * length( which( tdata.zero[,1]>=cutoffs[v1] &  tdata.zero[,2]>=cutoffs[v2] ) ) / ncells.zero)                    
				q4.zero = abs(100 - q1.zero - q2.zero - q3.zero)
			}

			if (checkGATED != "1") this$origin.ncells=ncells; this$coords=list()

			this$bintriplot(tdata,c(cutoffs[c(v1,v2)],0),binSize=binSize,mincells=mincount,density=TRUE)

			if (checkGATED != "1") {
				this$ncell.sel=this$origin.ncells
				this$ncell.perc=round(this$ncell.sel/this$origin.ncells*100,2)
				tkconfigure(this$ncell.gui,text=as.character(this$origin.ncells))
				tkconfigure(this$ncell.sel.gui,text=as.character(this$ncell.sel))
				tkconfigure(this$ncell.perc.gui,text=as.character(this$ncell.perc))
				tkconfigure(this$ncell.gui.di,text=as.character(this$origin.ncells))
				tkconfigure(this$ncell.sel.gui.di,text=as.character(this$ncell.sel))
				tkconfigure(this$ncell.perc.gui.di,text=as.character(this$ncell.perc))
			}

			### add title for single plot
			#if (checkGATED == "1") firstLine=sprintf("%s(%0.1f%%): density",displayfile,this$ncell.perc)
			if ( this$OverviewGate ) firstLine=sprintf("%s*(%0.1f%%): density/cof=%s",displayfile,this$ncell.perc,this$current.cofactor)
			else firstLine=sprintf("%s: density/cof=%s",displayfile,this$current.cofactor)
			title(main=firstLine,line=3.2,cex.main=0.9,adj=0)

			#secondLine=sprintf("%s-%s(%0.1f%%); min/max=%s/%s; binSize=%s,#bins=%s",ncells,ncells.zero,(ncells.zero/ncells*100),mincount,this$maxcells,binSize,this$bincount)
			secondLine=sprintf("cells(min/max)=%s/%s",mincount,this$maxcells)
			title(main=secondLine,line=2.4,cex.main=0.9,adj=0)
											
			#thirdLine=""
			#if ( cutoffs[v1] > 0 & cutoffs[v2] > 0 ) {
			#	thirdLine=sprintf("Q1=%0.1f/%0.1f Q2=%0.1f/%0.1f Q3=%0.1f/%0.1f Q4=%0.1f/%0.1f",q1.zero,this$q1.total,q2.zero,this$q2.total,q3.zero,this$q3.total,q4.zero,this$q4.total)
			#}
			thirdLine=sprintf("%s-%s(%0.1f%%); binSize=%s,#bins=%s",
			     ncells,ncells.zero,(ncells.zero/ncells*100),binSize,this$bincount)

							
			fifthLine=""
			if (checkGATED == "1" | grepl("^temp\\d+",file)) {
				if ( length(this$coords.info)>2 ) {
					fourthLine=sprintf("%s; %s",this$coords.info[1],this$coords.info[2])
					fifthLine=sprintf("%s; %s",this$coords.info[3],this$coords.info[4])
				} else {
					fourthLine=sprintf("%s",paste(this$coords.info,collapse="; "))
				}
				title(main=fourthLine,line=0.9,cex.main=0.6,adj=0)
				title(main=fifthLine,line=0.2,cex.main=0.6,adj=0)
			}

			tkconfigure(this$tt, cursor = "left_ptr")

			t = t + 1
		}
		#}
		if ( plotter=="triover" & (t==1) ) break
	}
	#devSet(devList()[which(dev.label==names(devList()))])
	
	if (pdf) {
		#add title for page in 3D-Overview tab
		#title=as.character(tclvalue(tkget(this$title,"1.0","end-1c")))
		#mtext(title, outer = TRUE, cex = 1.5,line=1.3,pos=2)
		#if (length(this$coords.info>0)) mtext(sprintf("(%s)",paste(this$coords.info,collapse="; ")),outer=TRUE,cex = 0.8)
		#dev.off()
		print("w: Done plotting densities in PDF.")
	}

	############# start plot history 
	if (plotter == "tri") {
		### only in plotter = "tri"
		
		this$plot.num = this$plot.num + 1
		### push plot attributes one down
		for (i in this$plot.num:1) {
			this$plot.attr[i+1] = this$plot.attr[i]
		}
		
		tdata=as.matrix(data[,c(1,2,2)])
		### pngs for history
		png(file=this$png.file,width=140,height=145,bg="transparent")
		# mar in points, mai in inches
		# oma adds title lines
		# order: bottom, left, top, and right
		par(mar=c(1.1,1.1,1.3,0.3))
		if ( cutoffs[1] > 0 & cutoffs[2] > 0 ) title.axis = c(sprintf("%s (%s)",vars[1],cutoffs[1]),sprintf("%s (%s)",vars[2],cutoffs[2]))
		else title.axis = c(vars[1],vars[2])
		plot(1,type='n',frame.plot=FALSE,xlim=c(xminval,xmaxval+10*binSize),axes=FALSE,
				ylim=c(yminval-2.5*binSize,ymaxval+5*binSize),xlab=title.axis[1],ylab=title.axis[2],cex.lab=set.cex*0.9,cex.axis=0.5*set.cex.axes,mgp=c(0.1,0,0))
		box(lwd=0.8,col="darkgrey")

		### add grid
		if(checkGRID=="1") {
		  xgrid.steps=seq(0,(xmaxval),by=grid.step)
		  abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
		}

		this$bintriplot(tdata,c(cutoffs[c(v1,v2)],0),set.cex=0.8,set.cex.axes=0.2,set.mgp=c(0.2,0,0),binSize=binSize,png=TRUE,density=TRUE)
		title(main=sprintf("%s: %s",displayfile,checkCALC),line=0.5,cex.main=0.8,adj=0)
		dev.off()
		###

		### plot history
		this$images[[this$plot.num]]=tclVar()
		tkimage.create("photo",this$images[[this$plot.num]],file=this$png.file)

		if (FALSE) {
		### save plot attributes
		this$plot.attr[[1]]=list(
			table=table,
			origin.table=this$current.project,
			file.name=displayfile,
			file.idx=file.idx,
			stain.idx=this$current.filenum,
			origin.ncells=this$origin.ncells,
			ncells.zero=ncells.zero,
			sel.ncells=this$ncell.sel,
			binSize=binSize,
			bincount=this$bincount,
			mincount=mincount,
			maxcount=this$maxcells,
			vars=c(vars[1],vars[2],vars[2]),
			FI=as.character(c(xminval,xmaxval,yminval,ymaxval)),
			MSI=c(tclvalue(this$vminMSI),tclvalue(this$vmaxMSI)),
			#coords=c(this$coords$x[1],this$coords$x[2],this$coords$y[1],this$coords$y[2]),
			coords=this$coords.info,
			dynrange=tclvalue(this$cbtdynRange),
			calc=checkCALC,
			grid=checkGRID,
			transformation=checkTRANS,
			xcutoff=cutoffs[1],
			ycutoff=cutoffs[2],
			zcutoff=NULL,
			CellInfo=NULL,
			RectInfo=NULL
		)

		if (tclvalue(this$cbtautoRect)=="1") {
				this$plot.attr[[1]]$CellInfo=TRUE
				this$plot.attr[[1]]$RectInfo=TRUE
		}
		}

		### save FI ranges for this transformation type
		this$changeFI.range(mode=3)

	}
	############# stop plot history

	tkconfigure(this$tt, cursor = "left_ptr")
}

fcs$plotHistograms <- function(plotter,pdf) {
	this=fcs
	printf("w: do plotHistograms: plotter=%s, pdf=%s",plotter,pdf)
	
	if ( plotter=="preproc" ) {
	  file=tclvalue(tkget(this$current.filePrep))
	} else {
	  file=tclvalue(tkget(this$tkchoosefile))
	}
	len=length(this$selected.vars)

	dev.label=paste("histogram",plotter,sep=".")

	### get column vector
	colvec=vector()
	cutoffs=vector()
	cbperc=vector()
	vars=vector()
	for (i in 1:len) {
		if ( tclvalue(this$cbVal[[i]]) == "1") {
			colvec=c(colvec,i)
			cutoffs = c(cutoffs,this$checkDigits(cutoff_id=i))
			cbperc = c(cbperc,tclvalue(this$cbcutoffperc[[i]]))
			vars = c(vars,this$selected.vars[i])
		}
	}
	len.colvec = length(colvec)

	if (len.colvec == 0) {

	}

	#if (plotter=="preproc") {

	
	#} else {
	if ( plotter!="preproc" ) {
		this$refreshPlotters()
		this$refreshCutoffs(current=TRUE)

		if (plotter=="di") {
			var1 = this$checkMarker(tclvalue(tkget(this$cbvar1di)))
			var1.idx = which(this$selected.vars==var1)
			tkset(this$cbvar1di,var1)
			var2 = this$checkMarker(tclvalue(tkget(this$cbvar2di)))
			var2.idx = which(this$selected.vars==var2)
			tkset(this$cbvar2di,var2)

			### if features are not in sample
			if ( length(var1)==0 | length(var2)==0 ){
				tkmessageBox(title = "An error has occured!",
					message = "Check your features.")
				stop("One of the features are not existent.")
			}

			### remove double features
			double.idx = which(colvec==var1.idx | colvec==var2.idx)
			colvec = colvec[-double.idx]
			cutoffs = cutoffs[-double.idx]

			colvec = c(var1.idx,var2.idx,colvec)
			cutoffs = c(
				this$checkDigits(cutoff_id=var1.idx),
				this$checkDigits(cutoff_id=var2.idx),
				cutoffs
			)
		} else if (plotter=="tri") {
			var1 = this$checkMarker(tclvalue(tkget(this$cbvar1)))
			var1.idx = which(this$selected.vars==var1)
			tkset(this$cbvar1,var1)
			var2 = this$checkMarker(tclvalue(tkget(this$cbvar2)))
			var2.idx = which(this$selected.vars==var2)
			tkset(this$cbvar2,var2)
			var3 = this$checkMarker(tclvalue(tkget(this$cbvar3)))
			var3.idx = which(this$selected.vars==var3)
			tkset(this$cbvar3,var3)

			### if features are not in sample
			if ( length(var1)==0 | length(var2)==0 | length(var3)==0 ){
				tkmessageBox(title = "An error has occured!",
					message = "Check your features.")
				stop("One of the features are not existent.")
			}

			### remove double features
			double.idx = which(colvec==var1.idx | colvec==var2.idx | colvec==var3.idx)
			colvec = colvec[-double.idx]
			cutoffs = cutoffs[-double.idx]

			colvec = c(
				var1.idx,
				var2.idx,
				var3.idx,
				colvec
			)
			cutoffs = c(
				this$checkDigits(cutoff_id=var1.idx),
				this$checkDigits(cutoff_id=var2.idx),
				this$checkDigits(cutoff_id=var3.idx),
				cutoffs
			)
		}

		#colvec = sort(unique(colvec))
		len.colvec=length(colvec)
	}


	if ( length(which(dev.label==names(devList()))) != 0 ) {
		devSet(devList()[which(dev.label==names(devList()))])
	} else if (!pdf) {
		if ( len.colvec>36 ) {
			devNew(type="x11",title="Histograms",width=7*2.8,height=ceiling(len.colvec/7)*2.0,label=dev.label)
			par(mfrow=c(ceiling(len.colvec/7),7),oma=c(1.5,1,1.5,1),mar=c(2,2.4,1.5,1),mgp=c(1.3,0.4,0))
		} else if ( len.colvec>25 ) {
			devNew(type="x11",title="Histograms",width=6*2.5,height=ceiling(len.colvec/6)*2.0,label=dev.label)
			par(mfrow=c(ceiling(len.colvec/6),6),oma=c(1.5,1,1.5,1),mar=c(2,2.4,1.5,1),mgp=c(1.3,0.4,0))
		} else if ( len.colvec>16 ) {
			devNew(type="x11",title="Histograms",width=5*3,height=ceiling(len.colvec/5)*2.3,label=dev.label)
			par(mfrow=c(ceiling(len.colvec/5),5),oma=c(1.5,1,1.5,1),mar=c(2,2.4,1.5,1),mgp=c(1.3,0.4,0))
		} else if ( len.colvec>9 ) {
			devNew(type="x11",title="Histograms",width=4*4,height=ceiling(len.colvec/4)*2.5,label=dev.label)
			par(mfrow=c(ceiling(len.colvec/4),4),oma=c(1.5,1,1.5,1),mar=c(2,2.4,1.5,1),mgp=c(1.3,0.4,0))
		} else {
			devNew(type="x11",title="Histograms",width=3*4,height=ceiling(len.colvec/3)*3,label=dev.label)
			par(mfrow=c(ceiling(len.colvec/3),3),oma=c(1.5,1,1.5,1),mar=c(2,2.4,1.5,1),mgp=c(1.3,0.4,0))
		}
		this$plot.windows = c(this$plot.windows,dev.label)
	}

	set.cex = 1.1
	if (pdf) {
		if (FALSE) {
			if (len.colvec/5>2) {
				#print("5")
				par(mfrow=c(7,ceiling(len.colvec/7)),oma=c(3,3,4,3),mar=c(5,5,7,5))
				set.cex = 2.5
			} else if (len.colvec/4>2) {
				#print("4")
				par(mfrow=c(6,ceiling(len.colvec/6)),oma=c(3,3,4,3),mar=c(5,5,7,5))
				set.cex = 2.5
			} else if (len.colvec/3>2) {
				#print("3")
				par(mfrow=c(5,ceiling(len.colvec/5)),oma=c(3,3,4,3),mar=c(5,5,7,5))
				set.cex = 2
			} else {
				par(mfrow=c(4,ceiling(len.colvec/4)),oma=c(3,3,4,3),mar=c(3,3,5,3))
				set.cex = 1.5
			}
		}
		print("hallo")
	} else if ( length(which(dev.label==names(devList()))) == 0 ) {
		devNew(type="x11",title=sprintf("Histograms for %sploTs",plotter),width=3*4,height=ceiling(len.colvec/3)*3.1,label=dev.label)
		if (len.colvec/5>2) par(mfrow=c(ceiling(len.colvec/5),5),oma=c(2.5,2,2.5,2),mar=c(3,2,3,2))
		else if (len.colvec/4>2) par(mfrow=c(ceiling(len.colvec/4),4),oma=c(2.5,2,2.5,2),mar=c(3,2,3,2))
		else  par(mfrow=c(ceiling(len.colvec/3),3),oma=c(2.5,2,2.5,2),mar=c(3,2,3,2))
		this$plot.windows = c(this$plot.windows,dev.label)
	} else {
		devSet(devList()[which(dev.label==names(devList()))])
	}

	if (pdf & this$OverviewGate) {
		table=this$temptable.name[this$temp.num]
		file.idx = this$selected.filenum
	} else if ( grepl("^temp\\d\\d", file) ) {
		table = file
		file.idx = 1
	} else {
		table = this$selected.project
		file.idx = this$current.filetable[which(this$current.filetable[,2]==file),1]
	}

	timeSTART = Sys.time()

	if ( is.null(this$data) | this$current.project != table | this$current.filenum != file.idx |
	     this$current.cofactor != as.numeric(tclvalue(this$rbasinh)) ) {
		if ( this$working ) print("Time loading data:")
		this$getFile(table,file.idx)
		if ( this$working ) print(Sys.time()-timeSTART)
	}

	data = this$data
	for ( i in colvec ) {
		xminval=as.numeric(tkget(this$minvalX))
		xmaxval=as.numeric(tkget(this$maxvalX))

		tdata=data[,i]
		d=density(tdata)

		if ( min(d$x) > xminval ) xminval = min(d$x)
		if ( max(d$x) < xmaxval ) xmaxval = max(d$x)
		
		### check for inifinity
		if (is.infinite(xminval) | is.na(xminval)) xminval = 0
		if (is.infinite(xmaxval) | is.na(xminval)) xmaxval = 8
		
		if (max(d$y)>10) {
			# if density too narrow
			# plot frequency histogram instead of density
			hist(tdata,freq=F,breaks=50,main=paste(colnames(data)[i]))
		} else {
			hist(tdata,freq=F,breaks=50,main=paste(colnames(data)[i]),xlab="",
				#col = "palegreen2",lwd=0.5)
				col = "#66bd63",lwd=0.5)
		  #lines(density(tdata), lwd = 2, col = "brown")
		  lines(density(tdata), lwd = 2, col = "#b2182b")
			#plot(d,main=paste(colnames(data)[i]),xlim=c(xminval,xmaxval),xlab="",cex=set.cex)
			#polygon(d, col="palegreen2", border="palegreen2")
      
			#"#b2182b","#006837"
			#"#f46d43","#66bd63"
		}
		# vertical lines
		for (x in ceiling(xminval) : floor(xmaxval) ) {
			abline(v=x,col="darkgrey")
		}

		### cutoff line and percentage
		if ( this$selected.cutoffs[i] != 0 ) {
			abline(v=this$selected.cutoffs[i],lty=2,lwd=2)
			tperc = round( length(tdata[which(tdata>=this$selected.cutoffs[i])]) / length(tdata)*100,1)
			### display percentage of production
			text(this$selected.cutoffs[i], max(d$y)-0.05*max(d$y),label=sprintf("%0.1f%%",tperc),cex=set.cex,pos=4,xpd=TRUE)
			### display asinh-value of cutoff
			text(this$selected.cutoffs[i], min(d$y)+0.05*max(d$y),label=sprintf("[%s]",this$selected.cutoffs[i]),cex=set.cex,pos=4,xpd=TRUE)
		}
	}

	if (pdf) {
		#add title for page in 3D-Overview tab
		title=as.character(tclvalue(tkget(this$title,"1.0","end-1c")))
		mtext(title, outer = TRUE, cex = 1.5, line=1.0, pos=2)
		if (length(this$coords.info>0)) mtext(sprintf("(%s)",paste(this$coords.info,collapse=";")),outer=TRUE,cex = set.cex,xpd=TRUE)
		#dev.off()
		print("w: Done plotting histograms in PDF.")
	} else {
		devSet(devList()[which(dev.label==names(devList()))])
	}
}

fcs$plotLegend <- function () {
	this=fcs

	dev.label="legend"
	dev.cur=dev.cur()
	if ( length(which(dev.label==names(devList()))) != 0 ) {
			devOff(devList()[which(dev.label==names(devList()))])
	} 

	devNew(type="x11",title="Color legend",width=4,height=0.7,label=dev.label)
	par(oma=c(1.3,0.5,1.0,0.5),mar=c(0.5,0,0,0))
	plot(1,type="n",frame.plot=F,xlim=c(1,11),axes=FALSE,
			ylim=c(0,1),xlab="Legend colors",ylab="",cex.lab=1.0,cex.axis=1.0)
	axis(side=1, at=c(1.35,6,10.55),
			labels=c("low","medium","high"),las=1,cex.axis=1.0,tick=F,mgp=c(1.0,0.3,0))
	for ( i in 1:10 ) {
			rect(i,0,i+1,1,col=this$col.rainbow[i+1],border=F)
	}

	devSet(dev.cur) 
}

if (FALSE) {
fcs$plotLegend <- function () {
	this=fcs

	printf("w: do plotLegend")

	this$ttprojects=tktoplevel()
	this$ttprojectsopen = TRUE
	
	# set window title
	tkwm.title(this$ttprojects,"Project List")
	#set window size
	tkwm.geometry(this$ttprojects,"400x600")

	##### add table widget
	listProjects=ttktreeview(this$ttprojects,columns=1,show="headings")
	
	# create scroll bar
	yscr <- tkscrollbar(listProjects, command=function(...) tkyview(listProjects,...))
	tkconfigure(listProjects,yscrollcommand=function(...) tkset(yscr,...))
	tkpack(yscr, side="right", fill="both")
	xscr <- tkscrollbar(listProjects, command=function(...) tkxview(listProjects, ...), orient = "horizontal")
	tkconfigure(listProjects,xscrollcommand=function(...) tkset(xscr,...))
	tkpack(xscr, side="bottom", fill="both")

	shade = c("none","gray")
	idx.file=grep("fileIdentity|fileIndex",this$dataframes.name)
	tcl(listProjects,"column",1,width=this$max.nchar)
	for (i in 1:length(this$total.projects)){
		size=nrow(this$getDFtable(this$dataframes.name[idx.file[i]]))
		tkinsert(listProjects,"","end",values=as.tclObj(paste(this$total.projects[i],"(",size,")",sep="")),tag=shade[i%%2+1])
	}
	#tcl(this$treetable,"configure",background="gray95")
	tktag.configure(this$treetable,"gray",background="gray95")
	tkpack(listProjects,fill="both",expand=TRUE)
	tk2tip(listProjects,"Double Click to Open")

	# get row index
	#this$projectChildren=tcl(listProjects,"children","")
	#print(sapply(as.character(this$projectChildren), function(i) tclvalue(tkindex(listProjects,i))))

	# double click on project    
	tkbind(listProjects, "<Double-Button-1>", function(W,x,y) {
		sel= as.character(tcl(W, "identify","row",x,y))
		this$selected.project = strsplit(tclvalue(tcl(W, "item", sel, "-values")),"\\(")[[1]][1]
		project.num = which(this$total.projects==this$selected.project)
		if ( project.num != this$current.projectnum) {
			print(paste("Project selected:",this$selected.project))
			this$refreshComboboxStart(project.num)
		}
		tkwm.withdraw(tkwinfo("toplevel",W))
	})
}
}

fcs$refreshPlotters <- function() {
	### function to delete plotter names if window is closed
	this=fcs

	if ( FALSE ) {#this$working) {
		printf("w: current dev names: %s",paste(names(devList()),collapse=" "))
		printf("w: old plotter names: %s",paste(this$plot.windows,collapse=" "))
	}
	if ( length(this$plot.windows)>0 ) {
		dev.names=names(devList)
		dev.del.idx=vector()
		for ( i in 1:length(this$plot.windows) ) {
			if ( length(which(this$plot.windows[i]==names(devList()))) == 0 ) {
				dev.del.idx= c(dev.del.idx,i)
			}
		}
		if ( length(dev.del.idx>0) ) this$plot.windows=this$plot.windows[-dev.del.idx]
	}
}

fcs$dodiploTtotal <- function() {
	this=fcs
	printf("w: do dodiploTtotal")

	checkTRANS = tclvalue(this$rbtrans)
	checkGRID = tclvalue(this$cbtshowGrid)

	file=tclvalue(tkget(this$tkchoosefile))
	binSize=as.numeric(tkget(this$binSizedi))
	mincells=as.numeric(tkget(this$minCountDi))
	
	xminval=as.numeric(tkget(this$minvalX))
	xmaxval=as.numeric(tkget(this$maxvalX))

	len.vars = length(this$selected.vars)
	colvec = vector()
	for (i in 1:len.vars) {
		if ( tclvalue(this$cbVal[[i]]) == "1") colvec=c(colvec,i)
	}
	len.colvec = length(colvec)



	if ( len.colvec < 2 ) {
		tkconfigure(this$tt, cursor = "left_ptr")
		tkmessageBox(title = "An error has occured!",
			message = "Please select at least two marker.", icon = "error", type = "ok")
		stop("Please select at least two marker.")
	}

	set.cex.axes = 0.75
	set.mgp = c(1.1, 0.3, 0)


	# if path didnt save yet
	if (!exists("lastpath",where=eval(parse(text="fcs")))){
		this$lastpath = getwd()
	}
	filename = this$shortenFilename(file)
	#filepath = tkgetSaveFile(defaultextension=".pdf", initialfile=file.path(this$lastpath,filename))
	file = tkgetSaveFile(defaultextension=".pdf", initialdir= this$lastpath, initialfile=filename)

	### remove file name to get file path
	filepath = unlist(strsplit(tclvalue(file),"/"))
	filepath = filepath[-length(filepath)]
	filepath = paste(filepath,collapse="/")

	if (this$working) printf("w path=%s filename=%s file=%s",filepath,filename,file)


	#devNew(type="pdf")
	#pdf(file=filepath)
	#devEval("pdf",name=filepath, {
	toPDF(filename, path=filepath, {
	if (len.colvec < 5) {
		par(mfrow=c(round(len.colvec/2),2),oma=c(0.5,0.5,2,0),mar=c(3,5,1,2))
		plot.num.per.page = round(len.colvec/2)*2
	} else if (len.colvec < 12) {
		par(mfrow=c(4,2),oma=c(0.5,0.5,2,0),mar=c(3,5,1,2))
		plot.num.per.page = 4*2
	} else if (len.colvec < 20) {
		par(mfrow=c(3,2),oma=c(0.5,0.5,2,0),mar=c(3,5,1,2))
		plot.num.per.page = 3*2
	} else {
		par(mfrow=c(2,2),oma=c(0.5,0.5,2,0),mar=c(3,5,1,2))
		plot.num.per.page = 2*2
	}

	cex.size = 12/len.colvec
	if (len.colvec <12) cex.size = 1

	it = 0
	for ( i in colvec ) {
		plot(1,type='n',frame.plot=FALSE,xlim=c(xminval-0.5,xmaxval+2),axes=FALSE,
			ylim=c(0,len.colvec),xlab=sprintf("%s (%s,binSize=%s/mincells=%s,cutoff=%s)",this$selected.vars[i],checkTRANS,binSize,mincells,tclvalue(this$vcutoffs[[i]])),
			ylab="",cex.lab=cex.size*0.7,cex.axis=cex.size,mgp=set.mgp)

		### no box for now
		# box(lwd=0.8,col="darkgrey")

		if ( checkTRANS == "asinh") {
			scale = this$asinh$scale
			labels = this$asinh$label
			grid.step = this$asinh$step
		} else {
			scale = this$biex$scale
			labels = this$biex$label
			grid.step = this$biex$step
		} 

		### draw axis on the bottom
		axis(side=1, at=scale,labels=labels,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")

		### grid
		if(checkGRID=="1") {
			xgrid.steps=seq(0,(xmaxval),by=grid.step)
			abline(v=xgrid.steps,col="grey",lty=3)
		}

		this$distepy = 0
		this$axis.label = rep("",len.colvec+4)
		this$dodiploT(mode="overview",totalID=i,total.cex=0.7)
		
		it = it + 1
		if (it == plot.num.per.page) {
			this$addTitle()
			it = 0
		}
	}

	this$addTitle()	

	})
	#dev.off()

	if (this$working) printf("Done.")
	else tkmessageBox(title = "PDF saved in.. ",
			message = paste(file), icon = "info", type = "ok")

	this$lastpath = filepath

	printf("w: Done dodiploTtotal")
}

fcs$dodiploT <- function(mode,totalID=NA,total.cex=NA) {
	this=fcs
	printf("w: dodiploT: mode=%s totalID=%s",mode,totalID)

	checkTRANS = tclvalue(this$rbtrans)
	if (checkTRANS =="") checkTRANS = tclvalue(this$rbtrans) = "asinh"
	checkCALC = tclvalue(this$rbcalc)
	if (checkCALC == "") {
	  checkCALC = tclvalue(this$rbcalc) = "MSI"
	  this$rbcalc = tclVar("MSI")
	}
	checkGRID = tclvalue(this$cbtshowGrid)
	checkDYNRANGE = tclvalue(this$cbtdynRange)
	checkTRIMMING = tclvalue(this$cbttrimming)
	if ( checkTRIMMING == "1" ) {
	  this$preprocData(mode="trim")
	}

	if (checkDYNRANGE != "1") {
		min.MSI=this$checkDigits(num=tkget(this$minMSI))
		max.MSI=this$checkDigits(num=tkget(this$maxMSI))
		if (max.MSI <= min.MSI) {
		  tmp = min.MSI 
		  min.MSI = max.MSI
		  max.MSI = tmp
		}
	}

	if ( checkTRANS == "asinh" ) {
	  scale = this$asinh$scale
	  label = this$asinh$label
	  grid.step = this$asinh$step
	} else {
	  scale = this$biex$scale
	  label = this$biex$label
	  grid.step = this$biex$step
	} 

	file=tclvalue(tkget(this$tkchoosefile))
	binSize=as.numeric(tkget(this$binSizedi))

	xminval=as.numeric(tkget(this$minvalX))
	xmaxval=as.numeric(tkget(this$maxvalX))

	len.vars = length(this$selected.vars)
	if ( is.na(totalID) ) {
		var1=this$checkMarker(tclvalue(tkget(this$cbvar1di)))
		tkset(this$cbvar1di,var1)
		var1.num = which(this$selected.vars==var1)

		### if features are not in sample
		if ( length(var1)==0 ){
			tkmessageBox(title = "An error has occured!",
				message = "Check your features again.")
			stop("One of the features are not existent.")
		}
		# set plot window
		new.plot=TRUE
		this$refreshPlotters()
		dev.label=paste("plotter","di",this$plotter.di.num,sep=".")
		if ( length(which(dev.label==names(devList()))) != 0 ) {
				devSet(devList()[which(dev.label==names(devList()))])
				new.plot=FALSE
		}
	} else {
		var1=this$selected.vars[totalID]
		var1.num = totalID
		new.plot=FALSE
	}
	cutoff1=this$checkDigits(cutoff_id=var1.num)

	tkconfigure(this$tt, cursor = "watch")
	
	# column and cutoff number vector
	colvec=vector()
	cutoffs=vector()
	cbperc=vector()
	vars=vector()
	if (mode == "single") {
		var2=this$checkMarker(tclvalue(tkget(this$cbvar2di)))
		tkset(this$cbvar2di,var2)

		### if features are not in sample
		if ( length(var2)==0 ){
			tkmessageBox(title = "An error has occured!",
				message = "Check your features again.")
			stop("One of the features are not existent.")
		}

		var2.num = which(this$selected.vars==var2)
		cutoff2=this$checkDigits(cutoff_id=var2.num)

		colvec=c(var2.num,var1.num)
		cutoffs=c(cutoff2,cutoff1)

		cbperc1 = tclvalue(this$cbcutoffperc[[var1.num]])
		cbperc2 = tclvalue(this$cbcutoffperc[[var2.num]])
		cbperc = c(cbperc2,cbperc1)
		vars = c(var2,var1)
		len = 2
	} else if (mode =="overview") {
		# if plotting overview
		for (i in 1:len.vars) {
			if ( tclvalue(this$cbVal[[i]]) == "1") {
				colvec=c(colvec,i);
				cbperc=c(cbperc,tclvalue(this$cbcutoffperc[[i]]))
				cutoffs=c(cutoffs,this$checkDigits(cutoff_id=i))
				vars=c(vars,this$selected.vars[i])
			}
		}
		colvec=c(colvec,var1.num)
		cutoffs=c(cutoffs,cutoff1)
		vars=c(vars,var1)
		len=length(colvec)

		if ( len == 1 ) {
			tkconfigure(this$tt, cursor = "left_ptr")
			tkmessageBox(title = "An error has occured!",
				message = "Please select at least one marker.", icon = "error", type = "ok")
			stop("Please select at least one marker.")
		}
	}

	if (this$working) printf("w: dodiploT vars=%s",paste(vars,collapse=" "))
	
	### import data   
	table = this$selected.project
	if ( grepl("temp",file) ) {
		file.idx = 1
		table=file
	} else {
		file.idx=this$current.filetable[which(this$current.filetable[,2]==file),1]
		this$selected.filenum = file.idx
	}
					
	# do not trim and remove doublets if data is gated
	if (this$diploTGate) {
		# if data is gated, just recall
		data=this$getData(table,file.idx,columns=colvec)
	}  else if ( this$current.project==table & this$current.filenum==file.idx & !is.null(this$data) & 
	             this$current.cofactor == as.numeric(tclvalue(this$rbasinh))) {
		data = this$data
	} else {
		this$getFile(table,file.idx)
		data=this$data
	}

	# new plot window if
	len.default=26
	if (len>len.default) len.default=len
	if (len.default>35) len.default=35
	if ( is.na(totalID) & (
		length(devList())==0
		| new.plot
		| (mode=="overview" & (this$distepy+len.default)>len.default+3) 
		| this$distepy>(len.default+3)
		| var1 != this$current.var1.diploT
		| !all(c(xminval,xmaxval) == this$current.xaxis)) 
		) {
		this$newdiploT(len=len.default)
		this$current.xaxis=c(xminval,xmaxval)
		this$current.var1.diploT = var1
	}
	
	### add grid
	if (FALSE) {
		if(checkGRID=="1") {
		xgrid.steps=seq(0,(xmaxval),by=grid.step)
		abline(v=xgrid.steps,col="grey",lty=3)
		}
	}
	

	### check cutoffs
	for ( i in 1:(len-1) ) {
		if ( cbperc[i] == "1" ) cutoffs[i] = this$calcCUTOFF(data[,colvec[i]],cutoffs[i],vars[i],colvec[i])
	}
	
	### plot bins
	# tmp.plot for a new plot window
	tmp.plot = FALSE
	pos.y = 0
	for ( i in 1:(len-1) ) {
		if (this$distepy > len.default) this$newdiploT(len.default)
		if (pos.y > len.default) {
			this$newdiploT(len.default)
			tmp.plot = TRUE
		}
		
		if (tmp.plot) {
			pos.y = i-len.default-3
			printf("TRUE pos.y = %s distepy=%s len.default=%s ",pos.y,this$distepy,len.default)
			#this$distepy = 0
		} else {
			pos.y = this$distepy+i-1
		}

		if (mode=="single") {
			this$bindiplot(data[,c(colvec[i],var1.num)],cutoffs[c(i,len)],bin.pos=this$distepy)
			this$distepy = this$distepy + 1
			
			this$origin.ncells = nrow(data)
			if (checkCALC != "freq") {
			  this$ncell.sel = nrow(data[which(data[,1]> cutoffs[1]),])
			  this$ncell.perc=round(this$ncell.sel/this$origin.ncells*100,2)
			  tkconfigure(this$ncell.sel.gui.di,text=as.character(this$ncell.sel))
			  tkconfigure(this$ncell.perc.gui.di,text=as.character(this$ncell.perc))
			}
			tkconfigure(this$ncell.gui.di,text=as.character(this$origin.ncells))
			
		} else {
			#pos.y = this$distepy+i-1
			#printf("pos.y=%s",pos.y)
			if ( !is.na(totalID) ) this$bindiplot(data[,c(colvec[i],var1.num)],cutoffs[c(i,len)],bin.pos=pos.y,totalcex=total.cex)
			else this$bindiplot(data[,c(colvec[i],var1.num)],cutoffs[c(i,len)],bin.pos=pos.y)
		}
	}

	if (mode=="overview") this$distepy = this$distepy+i
	
	tkconfigure(this$tt, cursor = "left_ptr")

	printf("w: Done dodiploT distepy=%s",this$distepy)
}

fcs$newdiploT <- function(len=26) {
	this=fcs

	checkTRANS = tclvalue(this$rbtrans)
	checkGRID = tclvalue(this$cbtshowGrid)
	checkDATE = tclvalue(this$cbtaddDate)

	xminval=as.numeric(tkget(this$minvalX))
	xmaxval=as.numeric(tkget(this$maxvalX))
	binSize=as.numeric(tkget(this$binSizedi))
	mincells=as.numeric(tkget(this$minCountDi))
	len.vars=length(this$selected.vars)
	var1=this$checkMarker(tclvalue(tkget(this$cbvar1di)))
	var1.idx=which(this$selected.vars==var1)

	### if features are not in sample
	if ( length(var1)==0 ){
		tkmessageBox(title = "An error has occured!",
			message = "Check your features.")
		stop("One of the features are not existent.")
	}



	# start plot
	# mgp: A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. 
	# The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
	set.cex.axes = 0.75
	set.mgp = c(1.3, 0.2, 0)

	this$plotter.di.num = this$plotter.di.num + 1
	dev.label=paste("plotter","di",this$plotter.di.num,sep=".")
	this$plot.windows = c(this$plot.windows,dev.label)
	
	devNew(type="x11",title=sprintf("n-diploTs  x-range=[%.1f,%.1f] x-axis=%s",xminval,xmaxval,var1),width=8,height=len*0.28,label=dev.label)
	par(oma=c(0,0.5,2,0),mar=c(2.5,11,1,2))
	plot(1,type='n',frame.plot=FALSE,xlim=c(xminval-0.5,xmaxval+2),axes=FALSE,
			ylim=c(0.5,len),xlab=sprintf("%s (%s,binSize=%s/mincells=%s,cutoff=%s)",var1,checkTRANS,binSize,mincells,tclvalue(this$vcutoffs[[var1.idx]])),ylab="",cex.lab=1.0,cex.axis=0.9,mgp=set.mgp)
	
	### no box for now
	# box(lwd=0.8,col="darkgrey")

	if ( checkTRANS == "asinh") {
		scale = this$asinh$scale
		labels = this$asinh$label
		grid.step = this$asinh$step
	} else {
		scale = this$biex$scale
		labels = this$biex$label
		grid.step = this$biex$step
	} 

	### draw axis on the bottom and on the left
	axis(side=1, at=scale,labels=labels,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")

	### grid
	if(checkGRID=="1") {
		xgrid.steps=seq(0,(xmaxval),by=grid.step)
		abline(v=xgrid.steps,col="grey",lty=3)
	}

	### add date
	if ( checkDATE=="1" ) {
		date = gsub("-","",Sys.Date())
		title(main=date,outer=T,line=1,cex.main=1.0,adj=1)
	}

	this$distepy = 0
	this$axis.label = rep("",len+4)

	printf("w: Done newdiploT: distepy=%s var1=%s var1.idx=%s total plots=%s",this$distepy,var1,var1.idx,len)
}

fcs$bindiplot <- function(
  data,              # table with two columns, col01=Feat. B, col02=Feat. A!
  cutoffs,           # vector with two values, cutoffs[1] for Feat. B, cutoffs[2] for Feat. A
  bin.pos,           # y-position where to plot the diplot
  set.cex=1.1,
  set.cex.axes=1.0,
  totalcex=1.0
){ 
	this=fcs

	#if (this$working) printf("w: do bindiplot: plot #%s, cutoffs=%s set.cex=%s",y,paste(cutoffs,collapse=" "),set.cex)
	empty=FALSE

	data.lengths=nrow(data)
	data.old = data

	### get variables
	xmin.val=as.numeric(tkget(this$minvalX))
	xmax.val=as.numeric(tkget(this$maxvalX))
	min.MSI=this$checkDigits(num=tkget(this$minMSI))
	max.MSI=this$checkDigits(num=tkget(this$maxMSI))
	binSize=as.numeric(tkget(this$binSizedi))
	mincells=as.numeric(tkget(this$minCountDi))
	len.vars=length(this$selected.vars)
	set.mgp = c(1.1, 0.3, 0)

	### checkbutton options
	checkTRANS = tclvalue(this$rbtrans)
	checkCALC = tclvalue(this$rbcalc)
	checkADDFILE = tclvalue(this$cbtaddFilename)
	checkDYNRANGE = tclvalue(this$cbtdynRange)
	checkSHOWMINBIN = tclvalue(this$cbtshowMinBins)
	checkDISPLAYA = tclvalue(this$cbtdisplayA)

	if (checkTRANS=="asinh") scipen = this$asinh$scipen
	else scipen = this$biex$scipen

	options(scipen=scipen)

	# legend from blue to red
	cols = this$col.rainbow
	cols.pale = this$col.rainbow.pale
	
	### file name
	title.text=as.character(tclvalue(tkget(this$title,"1.0","end-1c")))
	displayfile=this$shortenFilename(title.text)
	
	### write metadata in csv table
	metadatafile = sprintf("%s_PRIdivis_metadata.csv",this$current.project)
	if(!file.exists(metadatafile)) {
	  header = c("date","sample", "cofactor","calc"
	             ,"feat.A", "feat.B"
	             ,"minVal.B","maxVal.B","absRange.B"
	             ,"cutoff.A","cutoff.B"
	             ,"neg.total","pos.total"
	             ,"neg.prodcells","pos.prodcells")
	  write.table(t(header), metadatafile, sep = "\t", row.names=F, col.names=F)
	}
	current.date = format(Sys.Date(),"%y%m%d")

	### if cutoff is set for Feat. B, cut data from cutoff first
	if (cutoffs[1]>0 & checkCALC!="freq") {
	  data=data[which(data[,1]> cutoffs[1]),]
	  printf(">>>>>>>>>>>>>>>>> DATA CUTTED to %s",nrow(data))
	}

	### if data is empty, only display 'empty data'
	if ( nrow(data) == 0 ) {
		if (this$distepy==0) text(x=xmax.val*1/2,y=bin.pos+0.2,labels="[no bins to display]"
			,cex=0.5*set.cex,adj=1,mgp=set.mgp,xpd=TRUE)
		else text(x=this$current.diploT.x-1,y=bin.pos+0.2,labels="[no bins to display]"
			,cex=0.5*set.cex,adj=1,mgp=set.mgp,xpd=TRUE)
		empty=TRUE
	} else {
		### set negative values of colnum 1 to zero
		if ( checkCALC != "density" ) data[which(data[,1]<0),1] = 0

		### fX divides the range of data[,2] into intervals and codes the values in data[,2] according to which interval they fall. 
		# The leftmost interval corresponds to level one, the next leftmost to level two and so on.
		fX=cut(data[,2],breaks=seq(xmin.val,xmax.val,by=binSize),include.lowest=TRUE,dig.lab=5)

		if ( checkCALC == "density") {
			### number of cells in bin
			my.calc=aggregate(data[,2],by=list(fX),length)
		} else if ( checkCALC == "SD" ) {
			my.calc=aggregate(data[,1],by=list(fX),sd)
			cols=this$col.blackwhite
		} else if ( checkCALC == "SEM" ) {
			my.calc=aggregate(data[,1],by=list(fX),function(x) {
				SEM = sd(x)/sqrt(length(x))
				### if normally distributed, 95,4 % of the cells should lie inside the interval mean +/- SEM
				# interval_min = mean(x) - SEM
				# interval_max = mean(x) + SEM 
				SEM
			})
			cols=this$col.blackwhite
		} else if ( checkCALC == "RSEM" ) {
			my.calc=aggregate(data[,1],by=list(fX),function(x) {
				RSEM = sd(x)/sqrt(length(x))
				RSEM/mean(x)*100
			})
			cols=this$col.blackred
		} else if ( checkCALC == "freq" ) {
		  my.calc = aggregate(data[,1],by=list(fX),function(x) {
		    y= round( 100 * length(which(x >= cutoffs[1])) / length(x))
		    return(y)
			})
		} else {#} if ( grepl("MSI",checkCALC) ) {
			my.calc=aggregate(data[,1],by=list(fX),mean)
		} 

		my.lengths=aggregate(data[,1],by=list(fX),length)
		my.calc=cbind(my.calc,ncells=my.lengths$x)

		### get color steps and color "factors" my.calc.fac
		idx=which(my.calc$ncells>=mincells)
		idx.lower=which(my.calc$ncells<mincells & my.calc$ncells>0)
		range=""
		col.minmax="black"
		if (checkCALC == "freq") {
			### bin color factor
			my.calc.fac=cut(my.calc$x,breaks=seq(0,100,by=10),labels=1:10,include.lowest=TRUE)
			levels(my.calc.fac)=c(0,levels(my.calc.fac),11,1)
		} else if (checkCALC == "RSEM") {
			### bin color factor
			my.calc.fac=cut(my.calc$x,breaks=seq(0,50,by=5),labels=1:10,include.lowest=TRUE)
			levels(my.calc.fac)=c(0,levels(my.calc.fac),11,12)
			this$my.calc.fac2 = my.calc.fac
			for ( i in 1:length(my.calc.fac) ){
				if ( !is.na(my.calc$x[i]) ) {
					if (my.calc$x[i]>=50 & my.calc$ncells[i]>=mincells) {
						my.calc.fac[i] = 11
					}
				}
			}
		} else {
			min.y=floor(min(my.calc[idx,'x'])*10)/10
			max.y=ceiling(max(my.calc[idx,'x'])*10)/10
      
			### IF NO BINS TO DISPLAY
			if ( is.infinite(min.y) ) {
			  printf("I AM EMPTY, min.y = %s",min.y)
				empty=TRUE
				if (this$distepy==0) text(x=xmax.val*1/2,y=bin.pos+0.2,labels="[no bins to display]",cex=0.5*set.cex,adj=1,mgp=set.mgp,xpd=TRUE)
				else text(x=this$current.diploT.x-1,y=bin.pos+0.2,labels="[no bins to display]",cex=0.5*set.cex,adj=1,mgp=set.mgp,xpd=TRUE)
			} else {
				if ( grepl("MSI",checkCALC) & diff(c(min.y,max.y)) <= 0.5 ) col.minmax="red"
				if ( checkCALC=="density" & max.y < 100 ) col.minmax="red"

				### get steps for legend and plot
				if ( checkCALC == "density" | checkDYNRANGE=="1" ) {
					min.range=min.y
					max.range=max.y
				} else {
					min.range=min.MSI
					max.range=max.MSI
					range=sprintf(", man range: %0.2f-%0.2f",min.range,max.range)
				}
			}

			if (!empty) {
				### get steps
				step=round(diff(range(max.range,min.range))/10,2) 
				steps=seq(min.range,max.range,length.out=11)

				### bin color factor
				# printf("steps=%s",paste(steps))
				my.calc.fac=cut(my.calc$x,breaks=steps,labels=2:11,include.lowest=TRUE)

				levels(my.calc.fac)=c(0,levels(my.calc.fac),11,12)
				### if x < min.range
				my.calc.fac[which(my.calc$x<steps[1] & my.calc$ncells>=mincells)]=0
				### if x > max.range
				my.calc.fac[which(my.calc$x>steps[11] & my.calc$ncells>=mincells)]=11
			}
		}
		
		### if there are bins to display
		if (!empty) {
			my.calc=cbind(my.calc,fac=as.numeric(my.calc.fac)+1)
			
			### get min, max and absolute bin range of feat. B in all bins with mincells
			range.B = range(my.calc$x[which(my.calc$ncells>=mincells)])
			minVal.B = round(range.B[1],3)
			maxVal.B  = round(range.B[2],3)
			absRange = round(diff(range.B),3)

			this$my.calc.fac = my.calc.fac
			this$my.calc = my.calc

			max.cells = 0
			min.cells = 10000
			maxX = 0
			minX = 10000
			
			############ bindiplot
			plotting = TRUE
			if (checkDISPLAYA=="1") plotting = FALSE
			for ( i in min(idx,idx.lower):max(idx,idx.lower) ) {
				sub1=strsplit(as.character(my.calc$Group.1[i]),",")[[1]][1]
				if (grepl("^\\(",sub1)) sub2=as.numeric(strsplit(sub1,"\\(")[[1]][2])
				else sub2=as.numeric(strsplit(sub1,"\\[")[[1]][2])
				
				### display only bins for Feat. A+
				if ( (sub2+binSize)>cutoffs[2]) {
				  plotting=TRUE
				  if (!exists("plot.xmin")) plot.xmin = sub2
				  if (i == idx[length(idx)]) plot.xmax=sub2
				}
				
				### bin configs
				# if i is in idx.lower then there were not enough cells in the bin
				# so plot the bin in pale
				bin.color = bin.border = bin.lwd = NA
				bin.border = NA
				if (any(i==idx.lower)) {
				  if (checkSHOWMINBIN =="1") {
				    if ( any(checkCALC==c("SEM","RSEM","SD")) ) {
				      bin.color = cols[my.calc$fac[i]]
				      bin.border = "white"
				      bin.lwd = 3
				    } else {
				      bin.color = cols.pale[my.calc$fac[i]]
				      bin.border = "white"
				      bin.lwd = 3
				    }
				  }
				} else {
				  bin.color = cols[my.calc$fac[i]]
				  bin.border = "darkgrey"
				  bin.lwd = 0.05
				}
				
				### plot bin
				if (plotting) rect(sub2,bin.pos,sub2+binSize,bin.pos+0.5,col=bin.color,border=bin.border,lwd=bin.lwd)
				
				if (my.calc$ncells[i]>max.cells) max.cells=my.calc$ncells[i]   
				if (my.calc$ncells[i]<min.cells) min.cells=my.calc$ncells[i]   
				
				if (my.calc$x[i]>maxX) maxX=my.calc$x[i]   
				if (my.calc$x[i]<minX) minX=my.calc$x[i]   
				
				if (i == idx[1]) plot.filename = sub2
			}
      ################
			this$current.diploT.x = sub2

			### set position where to text the percentages
			if (plot.xmin < xmin.val) {plot.xmin=xmin.val}
  
			### draw line sigments between pairs of points
			segments(cutoffs[2], bin.pos-0.12, x1 = cutoffs[2], y1 = bin.pos+0.55)
			### get percentage of cells left and right from cutoff x
			data.left = data.old[which(data.old[,2]<cutoffs[2]),]
			data.right = data.old[which(data.old[,2]>=cutoffs[2]),]
			cells.left = 100*(nrow(data.left))/data.lengths
			cells.right = 100-cells.left
			
			### and plot the percentages of cells left and right from cutoff x
			text(x=plot.xmin-binSize,y=bin.pos+0.63+this$legend.space.di,label=sprintf("%0.1f%%",cells.left),cex=0.48*set.cex,pos=4)
			text(x=plot.xmax+binSize+0.2,y=bin.pos+0.63+this$legend.space.di,label=sprintf("%0.1f%%",cells.right),cex=0.48*set.cex,pos=2)
			
			is.csv.content = FALSE
			if ( cutoffs[2] != 0 & checkCALC == "freq") {
			  ### plot cutoff line Feature A+/- and cell frequency
				
				### get percentage of cells which are productive on left and right side
				prod.left = 100 * (nrow(data.left[which(data.left[,1]>=cutoffs[1]),])/nrow(data.left))
				prod.right = 100 * (nrow(data.right[which(data.right[,1]>=cutoffs[1]),])/nrow(data.right))
				
				text(x=plot.xmin+binSize,y=bin.pos+0.63+this$legend.space.di,label=sprintf("%0.1f%%",prod.left),cex=0.48*set.cex,col="red",pos=2)
				text(x=plot.xmax+binSize-0.2,y=bin.pos+0.63+this$legend.space.di,label=sprintf("%0.1f%%",prod.right),cex=0.48*set.cex,col="red",pos=4)
				
				csv.content <- cbind(current.date,displayfile,this$current.cofactor,checkCALC
				                     ,colnames(data)[2],colnames(data)[1]
				                     ,minVal.B,maxVal.B,absRange
				                     ,cutoffs[1],cutoffs[2]
				                     ,round(cells.left,2),round(cells.right,2)
				                     ,round(prod.left,2 ),round(prod.right,2))
				write.table(csv.content, sprintf("%s_PRIdivis_metadata.csv",this$current.project), sep = "\t", col.names = F, row.names = F, append = T)
				printf("w: diploT meta1 written in %s.",sprintf("%s_PRIdivis_metadata.csv",this$current.project))
				
			}	else if ( cutoffs[2] != 0 & checkCALC!="density") {
				### display production percentage of Feature B+/- on each side of cutoff line
				if ( cutoffs[1] != 0 ) {
					prod.left = 100 * (nrow(data.left[which(data.left[,1]>=cutoffs[1]),])/data.lengths)
					prod.right = 100 * (nrow(data.right[which(data.right[,1]>=cutoffs[1]),])/data.lengths)

					text(x=plot.xmin+binSize,y=bin.pos+0.63+this$legend.space.di,label=sprintf("%0.1f%%",prod.left),cex=0.48*set.cex,col="chartreuse4",pos=2)
					text(x=plot.xmax+binSize-0.2,y=bin.pos+0.63+this$legend.space.di,label=sprintf("%0.1f%%",prod.right),cex=0.48*set.cex,col="chartreuse4",pos=4)
					
					csv.content <- cbind(current.date,displayfile,this$current.cofactor,checkCALC
					                     ,colnames(data)[2],colnames(data)[1]
					                     ,minVal.B,maxVal.B,absRange
					                     ,cutoffs[1],cutoffs[2]
					                     ,round(cells.left,2),round(cells.right,2)
					                     ,round(prod.left,2 ),round(prod.right,2))
					write.table(csv.content, sprintf("%s_PRIdivis_metadata.csv",this$current.project), sep = "\t", col.names = F, row.names = F, append = T)
					printf("w: diploT meta2 written in %s.",sprintf("%s_PRIdivis_metadata.csv",this$current.project))
				} 
			} else {
			  csv.content <- cbind(current.date,displayfile,this$current.cofactor,checkCALC
			                       ,colnames(data)[2],colnames(data)[1]
			                       ,minVal.B,maxVal.B,absRange)
			  write.table(csv.content, sprintf("%s_PRIdivis_metadata.csv",this$current.project), sep = "\t", col.names = F, row.names = F, append = T)
			  printf("w: diploT meta3 written in %s.",sprintf("%s_PRIdivis_metadata.csv",this$current.project))
			}

			### print min/max value
			if ( any(checkCALC ==c("density","freq","RSEM")) ) print.minmax=sprintf("min/max=%0.0f/%0.0f",minX,maxX)
			else print.minmax = sprintf("min/max=%0.2f/%0.2f%s",minX,maxX,range)
			# print to the right of plot area
			text(x=sub2+0.8,y=bin.pos-0.03,label=print.minmax,col=col.minmax,cex=0.48*set.cex*totalcex,pos=4,xpd=TRUE)
		}
	}

	### print calculation method and cutoff if available
	print.label = checkCALC
	if (checkCALC != "freq" & cutoffs[1] > 0) {
	  print.label=sprintf("%s(+), cutoff=%s",print.label,cutoffs[1])
	} else if (cutoffs[1] > 0) {
	  print.label = sprintf("%s, cutoff=%s",print.label,cutoffs[1])
	}
	
	if (empty) {
		if (this$distepy=="1") text(x=xmax.val*2/3,y=bin.pos+0.21,label=print.label,cex=0.6*set.cex*totalcex,pos=4,col="red",xpd=TRUE)
		else text(x=this$current.diploT.x,y=bin.pos+0.34,label=print.label,cex=0.6*set.cex*totalcex,pos=4,col="red",xpd=TRUE)
	} else {
		text(x=sub2+0.8,y=bin.pos+0.34,label=print.label,cex=0.6*set.cex*totalcex,pos=4)
	}

	### print second feature and file name on y axis horizontal
	if (checkADDFILE=="1") {
		this$axis.label[bin.pos+1]=sprintf("%s: ",displayfile)
	} 
	
	if ( cutoffs[1]>0 | checkCALC != "density" ) {
	# if ( all(checkCALC != c("density","SD","RSEM","SEM")) ) {
	  this$axis.label[bin.pos+1]=paste0(this$axis.label[bin.pos+1],colnames(data)[1])
	} 
	
	text(x=plot.filename,y=bin.pos+0.2,labels=this$axis.label[bin.pos+1],
			cex=0.6*set.cex,pos=2,mgp=set.mgp,xpd=TRUE)

}

fcs$dodiGraph <- function(mode, total.cex=NA, set.zero=FALSE) {
	this=fcs
	printf("w: dodiGraph: mode=%s set.zero=%s",mode,set.zero)

	### checkbutton options
	checkTRANS = tclvalue(this$rbtrans)
	if (checkTRANS == "") {
		checkTRANS = "asinh"
		tclvalue(this$rbtrans) = "asinh"
	}

	if ( checkTRANS == "asinh") {
		scale = this$asinh$scale
		labels = this$asinh$label
		grid.step = this$asinh$step
	} else {
		scale = this$biex$scale
		labels = this$biex$label
		grid.step = this$biex$step
	} 

	checkGRID = tclvalue(this$cbtshowGrid)
	checkGATED = tclvalue(this$cbtgateData)
	checkDATE = tclvalue(this$cbtaddDate)
	checkCALC = tclvalue(this$rbcalc)
	if (checkCALC == "") {
		checkCALC = "MSI"
		tclvalue(this$rbcalc) = "MSI"
	}
	checkADDFILE = tclvalue(this$cbtaddFilename)
	if (checkADDFILE == "") {
		checkADDFILE = "1"
		tclvalue(this$cbtaddFilename) = "1"
	}
	checkSHOWMINBIN = tclvalue(this$cbtshowMinBins)

	checkDYNRANGE = tclvalue(this$cbtdynRange)
	if (checkDYNRANGE == "") {
		checkDYNRANGE = "1"
		tclvalue(checkDYNRANGE) = "1"
	}
	if (checkDYNRANGE != "1") {
		min.MSI=this$checkDigits(num=tkget(this$minMSI))
		max.MSI=this$checkDigits(num=tkget(this$maxMSI))
		if (max.MSI <= min.MSI) {
		  tmp = min.MSI 
		  min.MSI = max.MSI
		  max.MSI = tmp
		}
	}

	file=tclvalue(tkget(this$tkchoosefile))
	xminval=as.numeric(tkget(this$minvalX))
	xmaxval=as.numeric(tkget(this$maxvalX))
	binSize=as.numeric(tkget(this$binSizedi))
	mincells=as.numeric(tkget(this$minCountDi))

	tkconfigure(this$tt, cursor = "watch")


	len.vars=length(this$selected.vars)
	var1=this$checkMarker(tclvalue(tkget(this$cbvar1di)))
	var1.idx=which(this$selected.vars==var1)
	### if features are not in sample
	if ( length(var1)==0 ){
		tkmessageBox(title = "An error has occured!",
			message = "Check your Feature A.")
		stop("One of the features are not existent.")
	}

	cutoff1=this$checkDigits(cutoff_id=var1.idx)

	
	# column and cutoff number vector
	colvec=vector()
	cutoffs=vector()
	cbperc=vector()
	vars=vector()
	if (mode == "single") {
		var2=this$checkMarker(tclvalue(tkget(this$cbvar2di)))
		### if features are not in sample
		if ( length(var2)==0 ){
			tkmessageBox(title = "An error has occured!",
				message = "Check your feature Y.")
			stop("One of the features are not existent.")
		}

		tkset(this$cbvar2di,var2)
		var2.num = which(this$selected.vars==var2)
		cutoff2=this$checkDigits(cutoff_id=var2.num)

		colvec=c(var2.num,var1.idx)
		cutoffs=c(cutoff2,cutoff1)

		cbperc1 = tclvalue(this$cbcutoffperc[[var1.idx]])
		cbperc2 = tclvalue(this$cbcutoffperc[[var2.num]])
		cbperc = c(cbperc2,cbperc1)
		vars = c(var2,var1)
		len = 2
	} else if (mode =="overview") {
		# if plotting overview
		for (i in 1:len.vars) {
			if ( tclvalue(this$cbVal[[i]]) == "1") {
				colvec=c(colvec,i);
				cbperc=c(cbperc,tclvalue(this$cbcutoffperc[[i]]))
				cutoffs=c(cutoffs,this$checkDigits(cutoff_id=i))
				vars=c(vars,this$selected.vars[i])
			}
		}

		### import data   
		table = this$selected.project
		if ( grepl("temp",file) ) {
			file.idx = 1
			table=file
		} else {
			file.idx=this$current.filetable[which(this$current.filetable[,2]==file),1]
			this$selected.filenum = file.idx
		}
						
		# do not trim and remove doublets if data is gated
		if (this$diploTGate) {
			# if data is gated, just recall
			data=this$getData(table,file.idx,columns=colvec)
		}  else if ( this$current.project==table & this$current.filenum==file.idx & !is.null(this$data)  
		             & this$current.cofactor == as.numeric(tclvalue(this$rbasinh))) {
			data = this$data
		} else {
			this$getFile(table,file.idx)
			data=this$data
		}
		ncells = nrow(data)


		### check cutoffs
		for ( i in 1:length(cbperc) ) {
			if ( cbperc[i] == "1" ) cutoffs[i] = this$calcCUTOFF(data[,colvec[i]],cutoffs[i],vars[i],colvec[i])
		}

		if (any(var1.idx==colvec)) {
			idx.col = which(colvec==var1.idx)
			colvec = colvec[-idx.col]
			cutoffs = cutoffs[-idx.col]
			vars = vars[-idx.col]
		}

		colvec=c(colvec,var1.idx)
		cutoffs=c(cutoffs,cutoff1)
		vars=c(vars,var1)
		len=length(colvec)
		colors = rainbow(len-1)

		if ( len == 1 ) {
			tkconfigure(this$tt, cursor = "left_ptr")
			tkmessageBox(title = "An error has occured!",
				message = "Please select at least one marker.", icon = "error", type = "ok")
			stop("Please select at least one marker.")
		}
	}

	if (this$working) printf("w: vars=%s",paste(vars,collapse=" "))
	

	if (checkGATED != "1") {
		this$origin.ncells=ncells
		this$ncell.sel=this$origin.ncells
		this$ncell.perc=round(this$ncell.sel/this$origin.ncells*100,2)
		tkconfigure(this$ncell.gui,text=as.character(this$origin.ncells))
		tkconfigure(this$ncell.sel.gui,text=as.character(this$ncell.sel))
		tkconfigure(this$ncell.perc.gui,text=as.character(this$ncell.perc))
		tkconfigure(this$ncell.gui.di,text=as.character(this$origin.ncells))
		tkconfigure(this$ncell.sel.gui.di,text=as.character(this$ncell.sel))
		tkconfigure(this$ncell.perc.gui.di,text=as.character(this$ncell.perc))
	} else {
		this$ncell.sel=ncells
		this$ncell.perc=round(this$ncell.sel/this$origin.ncells*100,2)
		tkconfigure(this$ncell.sel.gui,text=as.character(this$ncell.sel))
		tkconfigure(this$ncell.perc.gui,text=as.character(this$ncell.perc))
		tkconfigure(this$ncell.sel.gui.di,text=as.character(this$ncell.sel))
		tkconfigure(this$ncell.perc.gui.di,text=as.character(this$ncell.perc))
	}


	### create plot table
	my.plot = vector()


	### fX divides the range of data[,2] (x-axis) into intervals and codes the values in data[,2] according to which interval they fall. 
	# The leftmost interval corresponds to level one, the next leftmost to level two and so on.
	fX=cut(data[,colvec[len]],breaks=seq(xminval,xmaxval,by=binSize),include.lowest=TRUE,dig.lab=5)

	for ( i in 1:(len-1) ) {

		### if cutoff is set, cut data from cutoff 
		#if (cutoffs[i]>0 & checkCALC == "MSI+") {
		#	data = data[which(data[,i] > cutoffs[i]),]
		#}

		### if data is empty, only display 'empty data'
		if ( nrow(data) != 0 ) {
			
			### set negative values of colnum 1 to zero
			#if ( checkCALC != "density" ) data[which(data[,colvec[i]]<0),1] = 0

			if ( checkCALC == "density") {
				### number of cells in bin
				my.calc=aggregate(data[,colvec[i]],by=list(fX),length)
			} else if ( checkCALC == "SD" ) {
				my.calc=aggregate(data[,colvec[i]],by=list(fX),sd)
				cols=this$col.blackwhite
			} else if ( checkCALC == "SEM" ) {
				my.calc=aggregate(data[,colvec[i]],by=list(fX),function(x) {
					SEM = sd(x)/sqrt(length(x))
					### if normally distributed, 95,4 % of the cells should lie inside the interval mean +/- SEM
					# interval_min = mean(x) - SEM
					# interval_max = mean(x) + SEM 
					SEM
				})
				cols=this$col.blackwhite
			} else if ( checkCALC == "RSEM" ) {
				my.calc=aggregate(data[,colvec[i]],by=list(fX),function(x) {
					RSEM = sd(x)/sqrt(length(x))
					RSEM/mean(x)*100
				})
				cols=this$col.blackred
			} else if ( checkCALC == "freq" ) {
				my.calc = aggregate(data[,colvec[i]],by=list(fX),function(x) {
					y= round( 100 * length(which(x >= cutoffs[1])) / length(x))
					return(y)
					})
			} else {
				# MSI
				if (set.zero) {
					printf("%s::%s values set to zero.",vars[i],length(which(data[,colvec[i]]<0)))
					data[which(data[,colvec[i]]<0),colvec[i]] = 0
				}
				my.calc=aggregate(data[,colvec[i]],by=list(fX),mean)
			} 

		}

		if (i==1) my.plot = my.calc
		else my.plot = cbind(my.plot,my.calc[,2])
		colnames(my.plot)[i+1] = vars[i]
	}

	my.lengths=aggregate(data[,1],by=list(fX),length)
	my.plot=cbind(my.plot,ncells=my.lengths$x)

	this$my.plot = my.plot


	### parameters for plotting frame
	idx.ncells = which(my.plot$ncells>=mincells)
	this$idx.ncells = idx.ncells
	
	xmin = strsplit(as.character(my.plot$Group.1[idx.ncells[1]]),",")[[1]][1]
	if ( grepl("^\\(",xmin) ) xmin = as.numeric(strsplit(xmin,"\\(")[[1]][2])
	else xmin = as.numeric(strsplit(xmin,"\\[")[[1]][2])
	xmin = xmin - 2*binSize
	xmax = strsplit(as.character(my.plot$Group.1[idx.ncells[length(idx.ncells)]]),",")[[1]][1]
	xmax = as.numeric(strsplit(xmax,"\\(")[[1]][2])
	xmax = xmax + 2*binSize


	# start plot
	# mgp: A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. 
	# The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
	set.cex.axes = 0.75
	set.mgp = c(1.3, 0.4, 0)

	this$plotter.di.graph = this$plotter.di.graph + 1
	dev.label=paste("plotter","digraph",this$plotter.di.graph,sep=".")
	this$plot.windows = c(this$plot.windows,dev.label)
	
	devNew(type="x11",title=sprintf("n-diGraphs  x-range=[%.1f,%.1f] x-axis=%s",xmin,xmax,var1),width=9,height=8,label=dev.label)
	par(oma=c(0,1,2,0),mar=c(3,3,1,21))
	
	### get x ticks for line graphs
	x.tick = vector()
	for ( i in 1:dim(my.plot)[1] ) {
		sub1=strsplit(as.character(my.calc$Group.1[i]),",")[[1]][1]
		if (grepl("^\\(",sub1)) sub2=as.numeric(strsplit(sub1,"\\(")[[1]][2])
		else sub2=as.numeric(strsplit(sub1,"\\[")[[1]][2])
		x.tick = c(x.tick,sub2)	
	}

	### plot line graphs
	legend.width = legend.type = legend.col = legend.point = legend.var =  legend.label = legend.range = vector()
	for ( i in 2:len ) {
		if (diff(range(my.plot[idx.ncells,i]))<0.5) {
			width = 0.5
			type = 3
			point = 3 
		} else {
			width = 1.5
			type = 1
			point = 1
		}
		yrange = round(diff(range(my.plot[idx.ncells,i])),2)
		#lines(x.tick[idx.ncells],my.plot[idx.ncells,i], type="o",lty=type,col=colors[i-1],lwd=width,pch=point,cex=0.5)

		legend.width = c(legend.width,width)
		legend.type = c(legend.type,type)
		#legend.col = c(legend.col,colors[i-1])
		legend.point = c(legend.point,point)
		legend.range = c(legend.range,yrange)
		legend.var = c(legend.var,names(my.plot)[i])
		legend.label = c(legend.label,sprintf("%s (%s)",names(my.plot)[i],yrange))
	}
	legend.table = data.frame(
		var = legend.var
		,label = legend.label
		,range = legend.range
		#,col = legend.col
		,width = legend.width
		,type = legend.type
		,point = legend.point
		,stringsAsFactors=FALSE
	)

	plot.len = 10
	if (dim(legend.table)[1]<plot.len) plot.len = dim(legend.table)[1]


	legend.table = legend.table[order(legend.table$range,decreasing=TRUE),]
	legend.table = legend.table[1:plot.len,]
	colors = rainbow(plot.len+2)
	legend.table = cbind(legend.table,col=colors[1:(length(colors)-2)])
	
	y_range = range(my.plot[idx.ncells,c(legend.table$var)])

	### plot frame
	plot(1,type='n',frame.plot=FALSE,axes=FALSE
			,xlim=c(xmin-10*binSize,xmax)
			,ylim=c(y_range[1],y_range[2])
			,xlab=sprintf("%s (%s,binSize=%s/mincells=%s,cutoff=%s)",var1,checkTRANS,binSize,mincells,tclvalue(this$vcutoffs[[var1.idx]]))
			,ylab=checkCALC,cex.lab=1.0,cex.axis=0.7,mgp=set.mgp)

	### no box for now
	# box(lwd=0.8,col="darkgrey")

	### draw axis on the bottom and on the left
	axis(side=1, at=scale,labels=labels,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
	axis(side=2, at=scale,labels=labels,las=3,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")

	### add grid
	if (FALSE) {
		if(checkGRID=="1") {
		  xgrid.steps=seq(0,(xmaxval),by=grid.step)
		  abline(v=xgrid.steps,col="grey",lty=3)
		}
	}

	### add date
	if ( checkDATE=="1" ) {
		date = gsub("-","",Sys.Date())
		title(main=date,outer=T,line=1,cex.main=1.0,adj=1)
	}

	this$legend.table = legend.table
	this$x.tick = x.tick
	### ṕlot line graphs
	for ( i in 1:plot.len ) {
		lines(x = x.tick[idx.ncells]
			,y=my.plot[idx.ncells,legend.table$var[i]]
			,lty=legend.table$type[i]
			,lwd=legend.table$width[i]
			,pch=legend.table$point[i]
			,col=as.character(legend.table$col[i])
			,type="o"
			,cex=0.5
		)

		if ( (i%%2)==0 ) {
			text(x=x.tick[idx.ncells[1]]
				,y=my.plot[idx.ncells[1],legend.table$var[i]]
				,labels=legend.table$lab[i]
				,col=as.character(legend.table$col[i])
				,pos=2
				,xpd=TRUE
			)
		} else {
			text(x=x.tick[idx.ncells[length(idx.ncells)]]
				,y=my.plot[idx.ncells[length(idx.ncells)],legend.table$var[i]]
				,labels=legend.table$lab[i]
				,col=as.character(legend.table$col[i])
				,pos=4
				,xpd=TRUE
			)
		}

	}

	### create legend
	legend(
		x = par()$usr[2] + 3*0.1*(par()$usr[2]-par()$usr[1]),
		y = par()$usr[4],
		legend.table$lab,
		col=as.character(legend.table$col),
		lty=legend.table$type,
		lwd=legend.table$width,
		pch=legend.table$point,
		xpd=TRUE
	)


	tkconfigure(this$tt, cursor = "left_ptr")

	printf("w: Done dodiGraph")
}

fcs$dotriploT <- function() {
  # calculate and display plot
  this=fcs
  printf("w: do dotriploT")
  
  this$refreshPlotters()
  
  checkGATED = tclvalue(this$cbtgateData)
  checkCALC = tclvalue(this$rbcalc)
  checkTRANS = tclvalue(this$rbtrans)
  checkGRID = tclvalue(this$cbtshowGrid)
  checkDATE = tclvalue(this$cbtaddDate)
  checkTRIMMING = tclvalue(this$cbttrimming)
  if ( checkTRIMMING == "1" ) {
    this$preprocData(mode="trim")
  }
  checkDYNRANGE = tclvalue(this$cbtdynRange)
  if (checkDYNRANGE != "1") {
    min.MSI=this$checkDigits(num=tkget(this$minMSI))
    max.MSI=this$checkDigits(num=tkget(this$maxMSI))
    if (max.MSI <= min.MSI) {
      tmp = min.MSI 
      min.MSI = max.MSI
      max.MSI = tmp
    }
  }
  
  if ( checkTRANS == "asinh" ) {
    scale = this$asinh$scale
    label = this$asinh$label
    grid.step = this$asinh$step
  } else {
    scale = this$biex$scale
    label = this$biex$label
    grid.step = this$biex$step
  } 
  
  
  
  quadrants.col = "black"
  
  v1=this$checkMarker(tclvalue(tkget(this$cbvar1)))
  tkset(this$cbvar1,v1)
  v2=this$checkMarker(tclvalue(tkget(this$cbvar2)))
  tkset(this$cbvar2,v2)
  v3=this$checkMarker(tclvalue(tkget(this$cbvar3)))
  tkset(this$cbvar3,v3)
  vars = c(v1,v2,v3)
  
  ### if Feature A is not in sample
  if ( length(v1)==0 ){
    tkmessageBox(title = "An error has occured!",
                 message = "Check your Feature A.")
    stop("Feature A is not existent.")
  }
  ### if Feature A is not in sample
  if ( length(v2)==0 ){
    tkmessageBox(title = "An error has occured!",
                 message = "Check your feature B.")
    stop("Feature B is not existent.")
  }
  ### if Feature A is not in sample
  if ( length(v3)==0 ){
    tkmessageBox(title = "An error has occured!",
                 message = "Check your feature C.")
    stop("Feature C is not existent.")
  }
  
  ### if manual range for z-axis is checked but no input
  if ( tclvalue(this$cbtdynRange) == "0" & tclvalue(this$vmaxMSI) == "0" ) {
    tkmessageBox(title = "An error has occured!",
                 message = "You forgot to set maximum manual range for Feature C (It is still zero).",icon = "error", type = "ok")
    stop("Set maximum manual range for Feature C (It is still zero).")
  }
  
  cutoffz = this$checkDigits(cutoff_id=which(this$selected.vars==v3))
  
  ### if method is freq or MSI(+), cutoff(z) needs to be setted
  if ( cutoffz <= 0 & (checkCALC=="freq" | checkCALC=="MSI(+)") ) {
    tkmessageBox(title = "An error has occured in mode: MSI(+) or freq!",
                 message = "You forgot to set cutoff for Feature C.", icon = "error", type = "ok")
    stop("Missing production cutoff for Feature C.")
  }
  ### if axes ranges are not the same
  if ( !this$checkAxesRange() ) {
    tkmessageBox(title = "An error has occured!",
                 message = "Please set x and y axis with the same range.", icon = "error", type = "ok")
    stop("Set x and y axis with the same range.")
  }
  
  dev.label=paste("plotter","tri",this$plotter.tri.num,sep=".")
  if ( length(which(dev.label==names(devList()))) == 0 ) {
    this$plotter.tri.num = this$plotter.tri.num + 1
    dev.label=paste("plotter","tri",this$plotter.tri.num,sep=".")
    ncol=as.numeric(tclvalue(this$vncol))
    nrow=as.numeric(tclvalue(this$vnrow))
    devNew(type="x11",title="n-triploTs",width=ncol*3.4,height=nrow*3.7,label=dev.label)
    # mar in points, mai in inches
    # oma adds title lines
    # order: bottom, left, top, and right
    par(mfrow=c(nrow,ncol),oma=c(0.5,1,2,1),mar=c(3,3,4,2))
    
    this$plot.windows = c(this$plot.windows,dev.label)
  } else {
    devSet(devList()[which(dev.label==names(devList()))])
  }
  
  this$rect.lwd=1
  
  empty=FALSE
  file=tclvalue(tkget(this$tkchoosefile))
  file=unlist(strsplit(file,"\\; |\\;|\\. | ")[[1]])
  binSize = as.numeric(tkget(this$binSize))
  mincount = as.numeric(tkget(this$minCountTri))
  xminval = as.numeric(tkget(this$minvalX))
  xmaxval = as.numeric(tkget(this$maxvalX))
  yminval = as.numeric(tkget(this$minvalY))
  ymaxval = as.numeric(tkget(this$maxvalY))
  
  
  tkconfigure(this$tt, cursor = "watch")
  if ( (checkGATED == "1") & (this$temp.num > 0) ) {
  #   table = file
  #   file = this$plot.attr[[1]]$file.name
  #   file.idx = 1
  # } else if ( grepl("^temp",file) ) {
    table=file
    file.idx = 0
  } else {
    table = this$selected.project
    file.idx = this$current.filetable[which(this$current.filetable[,2]==file),1]
    this$selected.filenum = file.idx
    this$coords.info = vector()
  }
  
  displayfile = this$shortenFilename(file)
  
  file2=""
  #### new getFile if temporary file is not present
  timeSTART = Sys.time()
  if ( is.null(this$data) | this$current.project != table | this$current.filenum != file.idx | 
       this$current.trans != tclvalue(this$rbtrans) | this$current.cofactor != as.numeric(tclvalue(this$rbasinh))) {
    #if ( !exists("data",env=fcs) | this$current.project != table | this$current.filenum != file.idx | this$current.trans != tclvalue(this$rbtrans)) {
    if ( this$working ) print("Time loading data:")
    this$getFile(table,file.idx)
    if ( this$working ) print(Sys.time()-timeSTART)
  }
  ####
  
  tdata = this$data[vars]
  this$tdata = tdata
  
  cutoff_idx = which(this$selected.vars==v1)
  cutoff_idy = which(this$selected.vars==v2)
  cutoff_idz = which(this$selected.vars==v3)
  cutoffs = c(cutoff_idx,cutoff_idy,cutoff_idz)
  
  ### if percentage is checked
  # calculate cutoffs and set check button to zero
  for ( i in 1:length(cutoffs)) {
    if ( tclvalue(this$cbcutoffperc[[cutoffs[i]]]) == "1" ) {
      cutoffs[i] = this$calcCUTOFF(tdata[vars[i]],this$checkDigits(cutoff_id=cutoffs[i]),vars[i],cutoffs[i])
    } else {
      cutoffs[i] = this$checkDigits(cutoff_id=cutoffs[i])
    }
  }
  printf("w: do dotriploT :: file: %s",file)
  printf("w: do dotriploT :: table: %s",table)
  printf("w: do dotriploT :: file.idx=%s",file.idx)
  printf("w: do dotriploT :: cutoffs=%s",paste(cutoffs,collapse=" "))
  
  ### calculate cells which were not plotted 
  cells.overmaxFI = length(which( tdata[,1]>xmaxval | tdata[,2]>ymaxval ))
  cells.underminFI = length(which( tdata[,1]<xminval | tdata[,2]<yminval ))
  cells.overmaxFI.perc = round(100 * cells.overmaxFI / (dim(tdata)[1]-cells.underminFI))
  ### warn if more then 5% productive cells (q2+q3+q4) werent plotted
  if ( cells.overmaxFI.perc >= 5 & !this$working) {
    tkmessageBox(title = "Warning!",
                 message = sprintf("Your cells exceed %s%% of your plot max ranges. You might want to increase your max ranges.",cells.overmaxFI.perc), 
                 icon = "info", type = "ok")
  }
  
  timeSTART = Sys.time()
  if ( this$working ) print("Time loading plot:")
  
  # start plot
  # mgp: A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. 
  # The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
  set.cex=1.2
  set.cex.axes=1.0
  set.mgp=c(1.7, 0.4, 0)
  if ( cutoffs[1] > 0 ) title.axis = sprintf("%s (%s)",v1,cutoffs[1])
  else title.axis = v1
  if ( cutoffs[2] > 0 ) title.axis = c(title.axis,sprintf("%s (%s)",v2,cutoffs[2]))
  else title.axis = c(title.axis,v2)
  
  plot(1,type='n',frame.plot=FALSE,xlim=c(xminval,xmaxval+10*binSize),axes=FALSE,
       ylim=c(yminval-2.5*binSize,ymaxval+5*binSize),xlab=title.axis[1],ylab=title.axis[2],cex.lab=set.cex,cex.axis=0.5*set.cex.axes,mgp=set.mgp)
  box(lwd=0.8,col="darkgrey")
  
  ### draw axis on the bottom and on the left
  axis(side=1, at=scale,labels=label,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
  axis(side=2, at=scale,labels=label,las=3,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
  
  ### add grid
  if(checkGRID=="1") {
    xgrid.steps=seq(0,(xmaxval),by=grid.step)
    ygrid.steps=seq(0,(ymaxval),by=grid.step)
    abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
  }
  
  
  ### calc quadrants in total
  ncells = ncells.total = nrow(tdata)
  
  ### OLD tdata.zero where all rows with negative values where filtered
  #tdata.zero = tdata[apply(tdata[, -1], MARGIN = 1, function(row) {all(row > 0)}), ]
  ### NEW tdata.zero where ONLY rows where filtered if v1 or v2 are negative
  tdata.zero = tdata[which(tdata[,1]>=0 | tdata[,2]>=0),]
  ncells.zero=nrow(tdata.zero)
  ### cut all cells which are not producing cells
  if ( checkCALC == "MSI(+)" ) tdata.plus = tdata[which(tdata[,3]> cutoffs[3]),]
  
  # q1 Quadrant unten links
  # q2 Quadrant unten rechts
  # q3 Quadrant oben rechts
  # q4 Quadrant oben links
  if ( cutoffs[1] > 0 & cutoffs[2] > 0 ) {
    
    ### count cells in quadrant
    tdata.q1 = tdata[which( tdata[,1]<cutoffs[1] &  tdata[,2]<cutoffs[2] ),3]
    tdata.q2 = tdata[which( tdata[,1]>=cutoffs[1] &  tdata[,2]<cutoffs[2] ),3]
    tdata.q3 = tdata[which( tdata[,1]>=cutoffs[1] &  tdata[,2]>=cutoffs[2] ),3]
    tdata.q4 = tdata[which( tdata[,1]<cutoffs[1] &  tdata[,2]>=cutoffs[2] ),3]
    
    ### q[x].total [ink=black]
    ### percentage of cells in quadrant to total cells 
    ### or in MSI(+): percentage of cells in quadrant to total positive cells
    this$q1.total = abs(100 * length( tdata.q1 ) / ncells)
    this$q2.total = abs(100 * length( tdata.q2 ) / ncells)
    this$q3.total = abs(100 * length( tdata.q3 ) / ncells)
    this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
    
    if ( cutoffs[3] > 0 ) {
      ### number of cells which are producing cells in feature C
      ncells = nrow(tdata[which(tdata[,3]> cutoffs[3]),])
      
      ### q[x].prodcells [ink=red]
      ### percentage of cells which are positive for feature C in quadrant to total quadrant cells
      this$q1.prodcells = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[3])]) /length(tdata.q1)
      if (is.nan(this$q1.prodcells)) this$q1.prodcells = 0
      this$q2.prodcells = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[3])]) /length(tdata.q2)
      if (is.nan(this$q2.prodcells)) this$q2.prodcells = 0
      this$q3.prodcells = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[3])]) /length(tdata.q3)
      if (is.nan(this$q3.prodcells)) this$q3.prodcells = 0
      this$q4.prodcells = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[3])]) /length(tdata.q4)
      if (is.nan(this$q4.prodcells)) this$q4.prodcells = 0
      
      ### only do MSI plots on producing cells only
      if ( checkCALC == "MSI(+)" ) {
        ncells = nrow(tdata.plus)
        
        tdata.q1 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
        tdata.q2 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
        tdata.q3 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
        tdata.q4 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
        
        ### q[x].total [ink=blue]
        ### in MSI(+): percentage of cells in quadrant to total positive cells
        this$q1.total = abs(100 * length( which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
        this$q2.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
        this$q3.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ) ) / ncells)
        this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
      }
      
      ### q[x].prodcellsplus [ink=green]
      ### percentage of cells which are positive for feature C to total cells
      this$q1.prodcellsplus = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[3])]) / ncells.total
      if (is.nan(this$q1.prodcellsplus)) this$q1.prodcellsplus = 0
      this$q2.prodcellsplus = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[3])]) / ncells.total
      if (is.nan(this$q2.prodcellsplus)) this$q2.prodcellsplus = 0
      this$q3.prodcellsplus = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[3])]) / ncells.total
      if (is.nan(this$q3.prodcellsplus)) this$q3.prodcellsplus = 0
      this$q4.prodcellsplus = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[3])]) / ncells.total
      if (is.nan(this$q4.prodcellsplus)) this$q4.prodcellsplus = 0
      
      if (this$working) {
        print("w: ncells ::  q1  q2  q3  q4")
        printf("w: %s(total) :: %s %s %s %s",ncells.total,length(tdata.q1),length(tdata.q2),length(tdata.q3),length(tdata.q4))
        printf("w: %s(black/blue) :: %.1f %.1f %.1f %.1f",ncells,this$q1.total,this$q2.total,this$q3.total,this$q4.total)
        printf("w: prodcells(red) :: %.1f %.1f %.1f %.1f",this$q1.prodcells,this$q2.prodcells,this$q3.prodcells,this$q4.prodcells)
        printf("w: prodcellsplus(green) :: %.1f %.1f %.1f %.1f",this$q1.prodcellsplus,this$q2.prodcellsplus,this$q3.prodcellsplus,this$q4.prodcells)
      }
    }
  } 
  
  
  if (checkGATED != "1") this$origin.ncells=ncells; this$coords=list()
  
  if (checkCALC == "density") {
    this$bintriplot(data=tdata,cutoffs=cutoffs,density=TRUE, binSize=binSize,mincells=mincount)
  } else if (checkCALC == "MSI(+)") {
    this$bintriplot(data=tdata.plus,cutoffs=cutoffs,binSize=binSize,mincells=mincount, quadrants.color = "blue", data.origin=tdata)
  } else {
    this$bintriplot(data=tdata,cutoffs=cutoffs,binSize=binSize,mincells=mincount)
  }
  
  if (checkGATED != "1") {
    this$ncell.sel=this$origin.ncells
    this$ncell.perc=round(this$ncell.sel/this$origin.ncells*100,2)
    tkconfigure(this$ncell.gui,text=as.character(this$origin.ncells))
    tkconfigure(this$ncell.sel.gui,text=as.character(this$ncell.sel))
    tkconfigure(this$ncell.perc.gui,text=as.character(this$ncell.perc))
    tkconfigure(this$ncell.gui.di,text=as.character(this$origin.ncells))
    tkconfigure(this$ncell.sel.gui.di,text=as.character(this$ncell.sel))
    tkconfigure(this$ncell.perc.gui.di,text=as.character(this$ncell.perc))
  }
  
  ### add title for single plot
  if (checkGATED == "1") firstLine=sprintf("%s*(%0.1f%%): %s/cof=%s",displayfile,this$ncell.perc,checkCALC,this$current.cofactor)
  else firstLine=sprintf("%s: %s/cof=%s",displayfile,checkCALC,this$current.cofactor)
  title(main=firstLine,line=3.2,cex.main=0.9,adj=0)
  
  if (checkCALC == "freq" | grepl("MSI",checkCALC)) secondLine=sprintf("cells(min/max)=%s/%s; %s",mincount,this$maxcells,v3)
  else secondLine=sprintf("cells(min/max)=%s/%s",mincount,this$maxcells)
  title(main=secondLine,line=2.4,cex.main=0.9,adj=0)
  
  thirdLine=sprintf("%s-%s(%0.1f%%); binSize=%s,#bins=%s",ncells.total,ncells.zero,(ncells.zero/ncells.total*100),binSize,this$bincount)
  title(main=thirdLine,line=1.6,cex.main=0.7,adj=0)
  
  if (checkGATED == "1" | grepl("temp",file)) {
    if ( length(this$coords.info)>2 ) {
      fourthLine=sprintf("%s %s",this$coords.info[1],this$coords.info[2])
      fifthLine=sprintf("%s %s",this$coords.info[3],this$coords.info[4])
    } else {
      fourthLine=sprintf("%s",paste(this$coords.info,collapse=";"))
      fifthLine=""
    }
    title(main=fourthLine,line=0.9,cex.main=0.6,adj=0)
    title(main=fifthLine,line=0.2,cex.main=0.6,adj=0)
  }
  
  ############# start plot history
  this$plot.num = this$plot.num + 1
  ### push plot attributes one down
  # for (i in this$plot.num:1) {
  #   this$plot.attr[i+1] = this$plot.attr[i]
  # } 
  
  ### save FI ranges for this transformation type
  this$changeFI.range(mode=3)
  
  ############# stop plot history
  if ( this$working ) print(Sys.time()-timeSTART)
  
  tkconfigure(this$tt, cursor = "left_ptr")
  
  if ( checkDATE=="1" ) {
    date = gsub("-","",Sys.Date())
    title(main=date,outer=T,line=1,cex.main=1.3,adj=1)
  }
}

fcs$dotriploTfiles <- function(read=FALSE) {
  # calculate and display plot
  this=fcs
  printf("w: do dotriploT for all files")
  
  this$refreshPlotters()
  
  checkGATED = tclvalue(this$cbtgateData)
  checkCALC = tclvalue(this$rbcalc)
  checkTRANS = tclvalue(this$rbtrans)
  checkGRID = tclvalue(this$cbtshowGrid)
  checkDATE = tclvalue(this$cbtaddDate)
  checkTRIMMING = tclvalue(this$cbttrimming)
  if ( checkTRIMMING == "1" ) {
    this$preprocData(mode="trim")
  }
  checkDYNRANGE = tclvalue(this$cbtdynRange)
  if (checkDYNRANGE != "1") {
    min.MSI=this$checkDigits(num=tkget(this$minMSI))
    max.MSI=this$checkDigits(num=tkget(this$maxMSI))
    
    if (max.MSI <= min.MSI) {
      tmp = min.MSI 
      min.MSI = max.MSI
      max.MSI = tmp
    }
  }
  
  if ( checkTRANS == "asinh" ) {
    scale = this$asinh$scale
    label = this$asinh$label
    grid.step = this$asinh$step
  } else {
    scale = this$biex$scale
    label = this$biex$label
    grid.step = this$biex$step
  } 
  
  
  
  quadrants.col = "black"
  
  v1=this$checkMarker(tclvalue(tkget(this$cbvar1)))
  tkset(this$cbvar1,v1)
  v2=this$checkMarker(tclvalue(tkget(this$cbvar2)))
  tkset(this$cbvar2,v2)
  v3=this$checkMarker(tclvalue(tkget(this$cbvar3)))
  tkset(this$cbvar3,v3)
  vars = c(v1,v2,v3)
  
  ### if Feature A is not in sample
  if ( length(v1)==0 ){
    tkmessageBox(title = "An error has occured!",
                 message = "Check your Feature A.")
    stop("Feature A is not existent.")
  }
  ### if Feature A is not in sample
  if ( length(v2)==0 ){
    tkmessageBox(title = "An error has occured!",
                 message = "Check your feature B.")
    stop("Feature B is not existent.")
  }
  ### if Feature A is not in sample
  if ( length(v3)==0 ){
    tkmessageBox(title = "An error has occured!",
                 message = "Check your feature C.")
    stop("Feature C is not existent.")
  }
  
  ### if manual range for z-axis is checked but no input
  if ( tclvalue(this$cbtdynRange) == "0" & tclvalue(this$vmaxMSI) == "0" ) {
    tkmessageBox(title = "An error has occured!",
                 message = "You forgot to set maximum manual range for Feature C (It is still zero).",icon = "error", type = "ok")
    stop("Set maximum manual range for Feature C (It is still zero).")
  }
  
  cutoffz = this$checkDigits(cutoff_id=which(this$selected.vars==v3))
  
  ### if method is freq or MSI(+), cutoff(z) needs to be set
  if ( cutoffz <= 0 & (checkCALC=="freq" | checkCALC=="MSI(+)") ) {
    tkmessageBox(title = "An error has occured in mode: MSI(+) or freq!",
                 message = "You forgot to set cutoff for Feature C.", icon = "error", type = "ok")
    stop("Missing production cutoff for Feature C.")
  }
  ### if axes ranges are not the same
  if ( !this$checkAxesRange() ) {
    tkmessageBox(title = "An error has occured!",
                 message = "Please set x and y axis with the same range.", icon = "error", type = "ok")
    stop("Set x and y axis with the same range.")
  }
  
  this$rect.lwd=1
  
  empty=FALSE
  binSize = as.numeric(tkget(this$binSize))
  mincount = as.numeric(tkget(this$minCountTri))
  xminval = as.numeric(tkget(this$minvalX))
  xmaxval = as.numeric(tkget(this$maxvalX))
  yminval = as.numeric(tkget(this$minvalY))
  ymaxval = as.numeric(tkget(this$maxvalY))
  
  tkconfigure(this$tt, cursor = "watch")
  timeSTART = Sys.time()

# SET FILES SPITZER ---------------------------------------------------------------
  if (workspace == "YH") {
    max.nhorizplots = 4
    max.nvertiplots = 4
    
    file.idx.vec = c(1:length(this$current.filenames))
    if (grepl("blood",this$db.name)) {
      file.idx.vec = c(1,4,3,5,6,8,7,9,10,12,11,13,2)
    } else if (grepl("spleen",this$db.name)) {
      file.idx.vec = c(1,3:13,2)
    } else if (grepl("lymph",this$db.name)) {
      file.idx.vec = c(2:4,6:13)
    }
    
    pdf.file = sprintf("%s/CD4_TIN_BL_d3_triploTs_%s_%s_%s_%s_%s.pdf",this$table.dir,this$version,v1,v2,v3,checkCALC)
  } else {
    max.nhorizplots = as.numeric(tclvalue(this$vncol))
    max.nvertiplots = as.numeric(tclvalue(this$vnrow))
    file.idx.vec = 1:length(this$current.filenames)
    pdf.file = sprintf("%s/%s_%s_%s_%s_%s_triploTs_%s.pdf",getwd(),this$current.project,v1,v2,v3,checkCALC,this$version)
  }
  
  timeSTART = Sys.time()
  cat("\n\n>>>> Start triploTOverviewXY with total data files: \n\n")
  if (this$working) printf("w: %s - time started",timeSTART)
  
  toPDF(file=pdf.file, 
        path=this$saveinFolder, 
        title=sprintf("project %s - triploTs of %s(%s.%s.%s)",this$current.project,checkCALC,v1,v2,v3), 
        ### 
        width=3.21*max.nhorizplots,
        height=3.5*max.nvertiplots,
        pointsize=11,
        ###
        {
          label.cex = 1.2
          set.cex.axes = 1
          set.mgp = c(1.9,0.5,0)
          par(mfrow=c(max.nvertiplots,max.nhorizplots),oma=c(0.5,1,6,1),mar=c(3,4,5,1))
          
  for ( file.idx in file.idx.vec ) {
    
    table = this$selected.project
    displayfile = this$shortenFilename(this$current.filenames[file.idx])
    this$coords.info = vector()
    
    #### new getFile 
    this$getFile(table,file.idx)
    ####
    tdata = this$data[vars]
    this$tdata = tdata
    
    
    cutoff_idx = which(this$selected.vars==v1)
    cutoff_idy = which(this$selected.vars==v2)
    cutoff_idz = which(this$selected.vars==v3)
    cutoffs = c(cutoff_idx,cutoff_idy,cutoff_idz)
    
    ### if percentage is checked
    # calculate cutoffs and set check button to zero
    for ( i in 1:length(cutoffs)) {
      if ( tclvalue(this$cbcutoffperc[[cutoffs[i]]]) == "1" ) {
        cutoffs[i] = this$calcCUTOFF(tdata[vars[i]],this$checkDigits(cutoff_id=cutoffs[i]),vars[i],cutoffs[i])
      } else {
        cutoffs[i] = this$checkDigits(cutoff_id=cutoffs[i])
      }
    }
    printf("w: do dotriploT :: file: %s",displayfile)
    printf("w: do dotriploT :: table: %s",table)
    printf("w: do dotriploT :: file.idx=%s",file.idx)
    printf("w: do dotriploT :: cutoffs=%s",paste(cutoffs,collapse=" "))
    printf("w: do dotriploT :: vars=%s",paste(colnames(this$tdata),collapse=" "))
    
    ### calculate cells which were not plotted 
    cells.overmaxFI = length(which( tdata[,1]>xmaxval | tdata[,2]>ymaxval ))
    cells.underminFI = length(which( tdata[,1]<xminval | tdata[,2]<yminval ))
    cells.overmaxFI.perc = round(100 * cells.overmaxFI / (dim(tdata)[1]-cells.underminFI))
    ### warn if more then 5% productive cells (q2+q3+q4) werent plotted
    if ( cells.overmaxFI.perc >= 5 & !this$working) {
      tkmessageBox(title = "Warning!",
                   message = sprintf("Your cells exceed %s%% of your plot max ranges. You might want to increase your max ranges.",cells.overmaxFI.perc), 
                   icon = "info", type = "ok")
    }
    
    timeSTART = Sys.time()
    if ( this$working ) print("Time loading plot:")
    
    # start plot
    # mgp: A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. 
    # The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
    set.cex=1.2
    set.cex.axes=1.0
    set.mgp=c(1.7, 0.4, 0)
    if ( cutoffs[1] > 0 ) title.axis = sprintf("%s (%s)",v1,cutoffs[1])
    else title.axis = v1
    if ( cutoffs[2] > 0 ) title.axis = c(title.axis,sprintf("%s (%s)",v2,cutoffs[2]))
    else title.axis = c(title.axis,v2)
    
    plot(1,type='n',frame.plot=FALSE,xlim=c(xminval,xmaxval+10*binSize),axes=FALSE,
         ylim=c(yminval-2.5*binSize,ymaxval+5*binSize),xlab=title.axis[1],ylab=title.axis[2],cex.lab=set.cex,cex.axis=0.5*set.cex.axes,mgp=set.mgp)
    box(lwd=0.8,col="darkgrey")
    
    ### draw axis on the bottom and on the left
    axis(side=1, at=scale,labels=label,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
    axis(side=2, at=scale,labels=label,las=3,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
    
    ### add grid
    if(checkGRID=="1") {
      xgrid.steps=seq(0,(xmaxval),by=grid.step)
      ygrid.steps=seq(0,(ymaxval),by=grid.step)
      abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
    }
    
    
    ### calc quadrants in total
    ncells = ncells.total = nrow(tdata)
    
    ### OLD tdata.zero where all rows with negative values where filtered
    #tdata.zero = tdata[apply(tdata[, -1], MARGIN = 1, function(row) {all(row > 0)}), ]
    ### NEW tdata.zero where ONLY rows where filtered if v1 or v2 are negative
    tdata.zero = tdata[which(tdata[,1]>=0 | tdata[,2]>=0),]
    ncells.zero=nrow(tdata.zero)
    # q1 Quadrant unten links
    # q2 Quadrant unten rechts
    # q3 Quadrant oben rechts
    # q4 Quadrant oben links
    if ( cutoffs[1] > 0 & cutoffs[2] > 0 ) {
      
      ### count cells in quadrant
      tdata.q1 = tdata[which( tdata[,1]<cutoffs[1] &  tdata[,2]<cutoffs[2] ),3]
      tdata.q2 = tdata[which( tdata[,1]>=cutoffs[1] &  tdata[,2]<cutoffs[2] ),3]
      tdata.q3 = tdata[which( tdata[,1]>=cutoffs[1] &  tdata[,2]>=cutoffs[2] ),3]
      tdata.q4 = tdata[which( tdata[,1]<cutoffs[1] &  tdata[,2]>=cutoffs[2] ),3]
      
      ### q[x].total [ink=black]
      ### percentage of cells in quadrant to total cells 
      ### or in MSI(+): percentage of cells in quadrant to total positive cells
      this$q1.total = abs(100 * length( tdata.q1 ) / ncells)
      this$q2.total = abs(100 * length( tdata.q2 ) / ncells)
      this$q3.total = abs(100 * length( tdata.q3 ) / ncells)
      this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
      
      if ( cutoffs[3] > 0 ) {
        ### number of cells which are producing cells in feature Z
        ncells = nrow(tdata[which(tdata[,3]> cutoffs[3]),])
        
        ### q[x].prodcells [ink=red]
        ### percentage of cells which are positive for feature Z in quadrant to total quadrant cells
        this$q1.prodcells = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[3])]) /length(tdata.q1)
        if (is.nan(this$q1.prodcells)) this$q1.prodcells = 0
        this$q2.prodcells = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[3])]) /length(tdata.q2)
        if (is.nan(this$q2.prodcells)) this$q2.prodcells = 0
        this$q3.prodcells = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[3])]) /length(tdata.q3)
        if (is.nan(this$q3.prodcells)) this$q3.prodcells = 0
        this$q4.prodcells = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[3])]) /length(tdata.q4)
        if (is.nan(this$q4.prodcells)) this$q4.prodcells = 0
        
        ### only do MSI plots on producing cells only
        if ( checkCALC == "MSI(+)" ) {
          ### cut all cells which are not producing cells
          tdata.plus = tdata[which(tdata[,3]> cutoffs[3]),]
          ncells = nrow(tdata.plus)
          
          tdata.q1 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
          tdata.q2 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
          tdata.q3 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
          tdata.q4 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
          
          ### q[x].total [ink=blue]
          ### in MSI(+): percentage of cells in quadrant to total positive cells
          this$q1.total = abs(100 * length( which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
          this$q2.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
          this$q3.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ) ) / ncells)
          this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
        }
        
        ### q[x].prodcellsplus [ink=green]
        ### percentage of cells which are positive for feature Z to total cells
        this$q1.prodcellsplus = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[3])]) / ncells.total
        if (is.nan(this$q1.prodcellsplus)) this$q1.prodcellsplus = 0
        this$q2.prodcellsplus = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[3])]) / ncells.total
        if (is.nan(this$q2.prodcellsplus)) this$q2.prodcellsplus = 0
        this$q3.prodcellsplus = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[3])]) / ncells.total
        if (is.nan(this$q3.prodcellsplus)) this$q3.prodcellsplus = 0
        this$q4.prodcellsplus = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[3])]) / ncells.total
        if (is.nan(this$q4.prodcellsplus)) this$q4.prodcellsplus = 0
        
        if (this$working) {
          print("w: ncells ::  q1  q2  q3  q4")
          printf("w: %s(total) :: %s %s %s %s",ncells.total,length(tdata.q1),length(tdata.q2),length(tdata.q3),length(tdata.q4))
          printf("w: %s(black/blue) :: %.1f %.1f %.1f %.1f",ncells,this$q1.total,this$q2.total,this$q3.total,this$q4.total)
          printf("w: prodcells(red) :: %.1f %.1f %.1f %.1f",this$q1.prodcells,this$q2.prodcells,this$q3.prodcells,this$q4.prodcells)
          printf("w: prodcellsplus(green) :: %.1f %.1f %.1f %.1f",this$q1.prodcellsplus,this$q2.prodcellsplus,this$q3.prodcellsplus,this$q4.prodcells)
        }
      }
    } 
    
    
    if (checkGATED != "1") this$origin.ncells=ncells; this$coords=list()
    
    if (checkCALC == "density") {
      this$bintriplot(data=tdata,cutoffs=cutoffs,density=TRUE, binSize=binSize,mincells=mincount,file=displayfile)
    } else if (checkCALC == "MSI(+)") {
      this$bintriplot(data=tdata.plus,cutoffs=cutoffs,binSize=binSize,mincells=mincount, quadrants.color = "blue", data.origin=tdata,file=displayfile)
    } else {
      this$bintriplot(data=tdata,cutoffs=cutoffs,binSize=binSize,mincells=mincount,file=displayfile)
    }
    
    if (checkGATED != "1") {
      this$ncell.sel=this$origin.ncells
      this$ncell.perc=round(this$ncell.sel/this$origin.ncells*100,2)
      tkconfigure(this$ncell.gui,text=as.character(this$origin.ncells))
      tkconfigure(this$ncell.sel.gui,text=as.character(this$ncell.sel))
      tkconfigure(this$ncell.perc.gui,text=as.character(this$ncell.perc))
      tkconfigure(this$ncell.gui.di,text=as.character(this$origin.ncells))
      tkconfigure(this$ncell.sel.gui.di,text=as.character(this$ncell.sel))
      tkconfigure(this$ncell.perc.gui.di,text=as.character(this$ncell.perc))
    }
    
    ### add title for single plot
    if (checkGATED == "1") firstLine=sprintf("%s(%0.1f%%): %s/cof=%s",displayfile,this$ncell.perc,checkCALC,this$current.cofactor)
    else firstLine=sprintf("%s: %s/cof=%s",displayfile,checkCALC,this$current.cofactor)
    title(main=firstLine,line=3.2,cex.main=0.9,adj=0)
    
    if (checkCALC == "freq" | grepl("MSI",checkCALC)) secondLine=sprintf("cells(min/max)=%s/%s;%s",mincount,this$maxcells,v3)
    else secondLine=sprintf("cells(min/max)=%s/%s",mincount,this$maxcells)
    title(main=secondLine,line=2.4,cex.main=0.9,adj=0)
    
    thirdLine=sprintf("%s-%s(%0.1f%%); binSize=%s,#bins=%s",ncells.total,ncells.zero,(ncells.zero/ncells.total*100),binSize,this$bincount)
    title(main=thirdLine,line=1.6,cex.main=0.7,adj=0)
  }
        })
  
  ### save FI ranges for this transformation type
  this$changeFI.range(mode=3)
  
  ############# stop 
  if ( this$working ) print(Sys.time()-timeSTART)
  
  printf("Saved in %s.",pdf.file)
  
  tkconfigure(this$tt, cursor = "left_ptr")
  
  if ( checkDATE=="1" ) {
    date = gsub("-","",Sys.Date())
    title(main=date,outer=T,line=1,cex.main=1.3,adj=1)
  }
  
}

fcs$dotriploTtable <- function() {
  # calculate and display plot
  this=fcs
  printf("w: do dotriploT for all files first 10 triploTs")
  
  this$refreshPlotters()
  
  checkGATED = tclvalue(this$cbtgateData)
  checkCALC = tclvalue(this$rbcalc)
  checkTRANS = tclvalue(this$rbtrans)
  checkGRID = tclvalue(this$cbtshowGrid)
  checkDATE = tclvalue(this$cbtaddDate)
  checkTRIMMING = tclvalue(this$cbttrimming)
  if ( checkTRIMMING == "1" ) {
    this$preprocData(mode="trim")
  }
  checkDYNRANGE = tclvalue(this$cbtdynRange)
  if (checkDYNRANGE != "1") {
    min.MSI=this$checkDigits(num=tkget(this$minMSI))
    max.MSI=this$checkDigits(num=tkget(this$maxMSI))
    
    if (max.MSI <= min.MSI) {
      tmp = min.MSI 
      min.MSI = max.MSI
      max.MSI = tmp
    }
  }
  
  if ( checkTRANS == "asinh" ) {
    scale = this$asinh$scale
    label = this$asinh$label
    grid.step = this$asinh$step
  } else {
    scale = this$biex$scale
    label = this$biex$label
    grid.step = this$biex$step
  } 
  
  quadrants.col = "black"
  
  
  this$rect.lwd=1
  
  ### if manual range for z-axis is checked but no input
  if ( tclvalue(this$cbtdynRange) == "0" & tclvalue(this$vmaxMSI) == "0" ) {
    tkmessageBox(title = "An error has occured!",
                 message = "You forgot to set maximum manual range for Feature C (It is still zero).",icon = "error", type = "ok")
    stop("Set maximum manual range for Feature C (It is still zero).")
  }
  ### if axes ranges are not the same
  if ( !this$checkAxesRange() ) {
    tkmessageBox(title = "An error has occured!",
                 message = "Please set x and y axis with the same range.", icon = "error", type = "ok")
    stop("Set x and y axis with the same range.")
  }
  
  empty=FALSE
  binSize = as.numeric(tkget(this$binSize))
  mincount = as.numeric(tkget(this$minCountTri))
  xminval = as.numeric(tkget(this$minvalX))
  xmaxval = as.numeric(tkget(this$maxvalX))
  yminval = as.numeric(tkget(this$minvalY))
  ymaxval = as.numeric(tkget(this$maxvalY))
  
  tkconfigure(this$tt, cursor = "watch")
  timeSTART = Sys.time()
  

# ask for triplot values ---------------------------------------------------
  printf("w: do open triplot table")
  file <- tclvalue(tkgetOpenFile(initialdir=this$table.dir,defaultextension="csv"))
  triplot.table = read.table(file,fill=T,header=T,sep="\t")
  
  
  
  
  
  for ( i in 1:10 ) {
    v1 = this$checkMarker(as.character(triplot.table$X[i]))
    v2 = this$checkMarker(as.character(triplot.table$Y[i]))
    v3 = this$checkMarker(as.character(triplot.table$Z[i]))
    vars = c(v1,v2,v3)
    
    triplot.q = as.character(triplot.table$Q[i])
    
    # SET FILES SPITZER ---------------------------------------------------------------
    if (workspace == "YH") {
      max.nhorizplots = 4
      max.nvertiplots = 4
      
      file.idx.vec = c(1:length(this$current.filenames))
      if (grepl("blood",this$db.name)) {
        file.idx.vec = c(1,4,3,5,6,8,7,9,10,12,11,13,2)
        tissue = "BL"
      } else if (grepl("spleen",this$db.name)) {
        file.idx.vec = c(1,3:13,2)
        tissue = "SP"
      } else if (grepl("lymph",this$db.name)) {
        file.idx.vec = c(2:4,6:13)
        tissue = "LN"
      }
      
      pdf.file = sprintf("%s/CD4_TIN_%s_d3_triploTs_%s_%s_%s_%s_%s_%s.pdf",
                         this$table.dir,tissue,this$version,v1,v2,v3,triplot.q,checkCALC)
    } else {
      max.nhorizplots = as.numeric(tclvalue(this$vncol))
      max.nvertiplots = as.numeric(tclvalue(this$vnrow))
      file.idx.vec = 1:length(this$current.filenames)
      pdf.file = sprintf("%s/%s_%s_%s_%s_%s_triploTs_%s.pdf",getwd(),this$current.project,v1,v2,v3,checkCALC,this$version)
      tissue = ""
    }
    
    timeSTART = Sys.time()
    if (this$working) printf("w: %s - time started",timeSTART)
    
    toPDF(file=pdf.file, 
          path=this$saveinFolder, 
          title=sprintf("%s - triploTs of %s(%s.%s.%s) %s",tissue,checkCALC,v1,v2,v3,triplot.q), 
          ### 
          width=3.21*max.nhorizplots,
          height=3.5*max.nvertiplots,
          pointsize=11,
          ###
          {
            label.cex = 1.2
            set.cex.axes = 1
            set.mgp = c(1.9,0.5,0)
            par(mfrow=c(max.nvertiplots,max.nhorizplots),oma=c(0.5,1,6,1),mar=c(3,4,5,1))
            
            for ( file.idx in file.idx.vec ) {
              
              table = this$selected.project
              displayfile = this$shortenFilename(this$current.filenames[file.idx])
              this$coords.info = vector()
              
              #### new getFile 
              this$getFile(table,file.idx)
              ####
              tdata = this$data[vars]
              printf(">>>>>>%s",colnames(tdata))
              this$tdata = tdata
              
              
              cutoff_idx = which(this$selected.vars==v1)
              cutoff_idy = which(this$selected.vars==v2)
              cutoff_idz = which(this$selected.vars==v3)
              cutoffs = c(cutoff_idx,cutoff_idy,cutoff_idz)
              
              ### if percentage is checked
              # calculate cutoffs and set check button to zero
              for ( i in 1:length(cutoffs)) {
                if ( tclvalue(this$cbcutoffperc[[cutoffs[i]]]) == "1" ) {
                  cutoffs[i] = this$calcCUTOFF(tdata[vars[i]],this$checkDigits(cutoff_id=cutoffs[i]),vars[i],cutoffs[i])
                } else {
                  cutoffs[i] = this$checkDigits(cutoff_id=cutoffs[i])
                }
              }
              printf("w: do dotriploT :: file: %s",displayfile)
              printf("w: do dotriploT :: table: %s",table)
              printf("w: do dotriploT :: file.idx=%s",file.idx)
              printf("w: do dotriploT :: cutoffs=%s",paste(cutoffs,collapse=" "))
              printf("w: do dotriploT :: vars=%s",paste(colnames(this$tdata),collapse=" "))
              
              ### calculate cells which were not plotted 
              cells.overmaxFI = length(which( tdata[,1]>xmaxval | tdata[,2]>ymaxval ))
              cells.underminFI = length(which( tdata[,1]<xminval | tdata[,2]<yminval ))
              cells.overmaxFI.perc = round(100 * cells.overmaxFI / (dim(tdata)[1]-cells.underminFI))
              ### warn if more then 5% productive cells (q2+q3+q4) werent plotted
              if ( cells.overmaxFI.perc >= 5 & !this$working) {
                tkmessageBox(title = "Warning!",
                             message = sprintf("Your cells exceed %s%% of your plot max ranges. You might want to increase your max ranges.",cells.overmaxFI.perc), 
                             icon = "info", type = "ok")
              }
              
              timeSTART = Sys.time()
              if ( this$working ) print("Time loading plot:")
              
              # start plot
              # mgp: A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. 
              # The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
              set.cex=1.2
              set.cex.axes=1.0
              set.mgp=c(1.7, 0.4, 0)
              if ( cutoffs[1] > 0 ) title.axis = sprintf("%s (%s)",v1,cutoffs[1])
              else title.axis = v1
              if ( cutoffs[2] > 0 ) title.axis = c(title.axis,sprintf("%s (%s)",v2,cutoffs[2]))
              else title.axis = c(title.axis,v2)
              
              plot(1,type='n',frame.plot=FALSE,xlim=c(xminval,xmaxval+10*binSize),axes=FALSE,
                   ylim=c(yminval-2.5*binSize,ymaxval+5*binSize),xlab=title.axis[1],ylab=title.axis[2],cex.lab=set.cex,cex.axis=0.5*set.cex.axes,mgp=set.mgp)
              box(lwd=0.8,col="darkgrey")
              
              ### draw axis on the bottom and on the left
              axis(side=1, at=scale,labels=label,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
              axis(side=2, at=scale,labels=label,las=3,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
              
              ### add grid
              if(checkGRID=="1") {
                xgrid.steps=seq(0,(xmaxval),by=grid.step)
                ygrid.steps=seq(0,(ymaxval),by=grid.step)
                abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
              }
              
              
              ### calc quadrants in total
              ncells = ncells.total = nrow(tdata)
              
              ### OLD tdata.zero where all rows with negative values where filtered
              #tdata.zero = tdata[apply(tdata[, -1], MARGIN = 1, function(row) {all(row > 0)}), ]
              ### NEW tdata.zero where ONLY rows where filtered if v1 or v2 are negative
              tdata.zero = tdata[which(tdata[,1]>=0 | tdata[,2]>=0),]
              ncells.zero=nrow(tdata.zero)
              # q1 Quadrant unten links
              # q2 Quadrant unten rechts
              # q3 Quadrant oben rechts
              # q4 Quadrant oben links
              if ( cutoffs[1] > 0 & cutoffs[2] > 0 ) {
                
                ### count cells in quadrant
                tdata.q1 = tdata[which( tdata[,1]<cutoffs[1] &  tdata[,2]<cutoffs[2] ),3]
                tdata.q2 = tdata[which( tdata[,1]>=cutoffs[1] &  tdata[,2]<cutoffs[2] ),3]
                tdata.q3 = tdata[which( tdata[,1]>=cutoffs[1] &  tdata[,2]>=cutoffs[2] ),3]
                tdata.q4 = tdata[which( tdata[,1]<cutoffs[1] &  tdata[,2]>=cutoffs[2] ),3]
                
                ### q[x].total [ink=black]
                ### percentage of cells in quadrant to total cells 
                ### or in MSI(+): percentage of cells in quadrant to total positive cells
                this$q1.total = abs(100 * length( tdata.q1 ) / ncells)
                this$q2.total = abs(100 * length( tdata.q2 ) / ncells)
                this$q3.total = abs(100 * length( tdata.q3 ) / ncells)
                this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
                
                if ( cutoffs[3] > 0 ) {
                  ### number of cells which are producing cells in feature C
                  ncells = nrow(tdata[which(tdata[,3]> cutoffs[3]),])
                  
                  ### q[x].prodcells [ink=red]
                  ### percentage of cells which are positive for feature C in quadrant to total quadrant cells
                  this$q1.prodcells = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[3])]) /length(tdata.q1)
                  if (is.nan(this$q1.prodcells)) this$q1.prodcells = 0
                  this$q2.prodcells = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[3])]) /length(tdata.q2)
                  if (is.nan(this$q2.prodcells)) this$q2.prodcells = 0
                  this$q3.prodcells = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[3])]) /length(tdata.q3)
                  if (is.nan(this$q3.prodcells)) this$q3.prodcells = 0
                  this$q4.prodcells = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[3])]) /length(tdata.q4)
                  if (is.nan(this$q4.prodcells)) this$q4.prodcells = 0
                  
                  ### only do MSI plots on producing cells only
                  if ( checkCALC == "MSI(+)" ) {
                    ### cut all cells which are not producing cells
                    tdata.plus = tdata[which(tdata[,3]> cutoffs[3]),]
                    ncells = nrow(tdata.plus)
                    
                    tdata.q1 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
                    tdata.q2 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
                    tdata.q3 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
                    tdata.q4 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
                    
                    ### q[x].total [ink=blue]
                    ### in MSI(+): percentage of cells in quadrant to total positive cells
                    this$q1.total = abs(100 * length( which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
                    this$q2.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
                    this$q3.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ) ) / ncells)
                    this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
                  }
                  
                  ### q[x].prodcellsplus [ink=green]
                  ### percentage of cells which are positive for feature C to total cells
                  this$q1.prodcellsplus = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[3])]) / ncells.total
                  if (is.nan(this$q1.prodcellsplus)) this$q1.prodcellsplus = 0
                  this$q2.prodcellsplus = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[3])]) / ncells.total
                  if (is.nan(this$q2.prodcellsplus)) this$q2.prodcellsplus = 0
                  this$q3.prodcellsplus = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[3])]) / ncells.total
                  if (is.nan(this$q3.prodcellsplus)) this$q3.prodcellsplus = 0
                  this$q4.prodcellsplus = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[3])]) / ncells.total
                  if (is.nan(this$q4.prodcellsplus)) this$q4.prodcellsplus = 0
                  
                  if (this$working) {
                    print("w: ncells ::  q1  q2  q3  q4")
                    printf("w: %s(total) :: %s %s %s %s",ncells.total,length(tdata.q1),length(tdata.q2),length(tdata.q3),length(tdata.q4))
                    printf("w: %s(black/blue) :: %.1f %.1f %.1f %.1f",ncells,this$q1.total,this$q2.total,this$q3.total,this$q4.total)
                    printf("w: prodcells(red) :: %.1f %.1f %.1f %.1f",this$q1.prodcells,this$q2.prodcells,this$q3.prodcells,this$q4.prodcells)
                    printf("w: prodcellsplus(green) :: %.1f %.1f %.1f %.1f",this$q1.prodcellsplus,this$q2.prodcellsplus,this$q3.prodcellsplus,this$q4.prodcells)
                  }
                }
              } 
              
              
              if (checkGATED != "1") this$origin.ncells=ncells; this$coords=list()
              
              if (checkCALC == "density") {
                this$bintriplot(data=tdata,cutoffs=cutoffs,density=TRUE, binSize=binSize,mincells=mincount)
              } else if (checkCALC == "MSI(+)") {
                this$bintriplot(data=tdata.plus,cutoffs=cutoffs,binSize=binSize,mincells=mincount, quadrants.color = "blue", data.origin=tdata)
              } else {
                this$bintriplot(data=tdata,cutoffs=cutoffs,binSize=binSize,mincells=mincount)
              }
              
              if (checkGATED != "1") {
                this$ncell.sel=this$origin.ncells
                this$ncell.perc=round(this$ncell.sel/this$origin.ncells*100,2)
                tkconfigure(this$ncell.gui,text=as.character(this$origin.ncells))
                tkconfigure(this$ncell.sel.gui,text=as.character(this$ncell.sel))
                tkconfigure(this$ncell.perc.gui,text=as.character(this$ncell.perc))
                tkconfigure(this$ncell.gui.di,text=as.character(this$origin.ncells))
                tkconfigure(this$ncell.sel.gui.di,text=as.character(this$ncell.sel))
                tkconfigure(this$ncell.perc.gui.di,text=as.character(this$ncell.perc))
              }
              
              ### add title for single plot
              if (checkGATED == "1") firstLine=sprintf("%s(%0.1f%%): %s/cof=%s",displayfile,this$ncell.perc,checkCALC,this$current.cofactor)
              else firstLine=sprintf("%s: %s/cof=%s",displayfile,checkCALC,this$current.cofactor)
              title(main=firstLine,line=3.2,cex.main=0.9,adj=0)
              
              if (checkCALC == "freq" | grepl("MSI",checkCALC)) secondLine=sprintf("cells(min/max)=%s/%s;%s",mincount,this$maxcells,v3)
              else secondLine=sprintf("cells(min/max)=%s/%s",mincount,this$maxcells)
              title(main=secondLine,line=2.4,cex.main=0.9,adj=0)
              
              thirdLine=sprintf("%s-%s(%0.1f%%); binSize=%s,#bins=%s",ncells.total,ncells.zero,(ncells.zero/ncells.total*100),binSize,this$bincount)
              title(main=thirdLine,line=1.6,cex.main=0.7,adj=0)
            }
          })
    printf("Saved file in %s",pdf.file)
    
    ### save FI ranges for this transformation type
    this$changeFI.range(mode=3)
    
    ############# stop 
    if ( this$working ) print(Sys.time()-timeSTART)
    
    printf("Saved in %s.",pdf.file)
    
    tkconfigure(this$tt, cursor = "left_ptr")
    
    if ( checkDATE=="1" ) {
      date = gsub("-","",Sys.Date())
      title(main=date,outer=T,line=1,cex.main=1.3,adj=1)
    }
  }
}

fcs$dotriploTRectData <- function() {
	this=fcs
	printf("w: do dotriploTRectData")
	
	# get information: file, vars, cutoffs
	file=tclvalue(tkget(this$tkchoosefile))
	var1=this$checkMarker(tclvalue(tkget(this$cbvar1)))
	tkset(this$cbvar1,var1)
	var2=this$checkMarker(tclvalue(tkget(this$cbvar2)))
	tkset(this$cbvar2,var2)
	var3=this$checkMarker(tclvalue(tkget(this$cbvar3)))
	tkset(this$cbvar3,var3)
	vars=c(var1,var2,var3)

	### if features are not in sample
	if ( length(var1)==0 | length(var2)==0 | length(var3)==0 ){
		tkmessageBox(title = "An error has occured!",
			message = "Check your features.")
		stop("One of the features are not existent.")
	}

	cutoffx = this$checkDigits(cutoff_id=which(this$selected.vars==var1))
	cutoffy = this$checkDigits(cutoff_id=which(this$selected.vars==var2))
	cutoffz = this$checkDigits(cutoff_id=which(this$selected.vars==var3))
	cutoffs=c(cutoffx,cutoffy,cutoffz)
	
	### if axes ranges are not the same
	try ( if ( !this$checkAxesRange() ) {
			tkmessageBox(title = "An error has occured!",
					message = "Please set x and y axis with the same range.", icon = "error", type = "ok")
			stop("Set x and y axis with the same range.")
	})

	# set plot window
	dev.label=paste("plotter","tri",this$plotter.tri.num,sep=".")
	devSet(devList()[which(dev.label==names(devList()))])
	# set info label
	tkconfigure(this$ncell.gui,text=as.character(this$origin.ncells))
	tkconfigure(this$ncell.gui.di,text=as.character(this$origin.ncells))

	# get information: binSize, mincount, FIs
	binSize=as.numeric(tkget(this$binSize))
	mincount=as.numeric(tkget(this$minCountTri))
	xminval=as.numeric(tkget(this$minvalX))
	xmaxval=as.numeric(tkget(this$maxvalX))
	yminval=as.numeric(tkget(this$minvalY))
	ymaxval=as.numeric(tkget(this$maxvalY))
	
	checkGATED=tclvalue(this$cbtgateData)
	checkCALC = tclvalue(this$rbcalc)
	checkGRID = tclvalue(this$cbtshowGrid)
	checkTRANS=tclvalue(this$rbtrans)
	if ( checkTRANS == "asinh" ) {
	  scale = this$asinh$scale
	  label = this$asinh$label
	  grid.step = this$asinh$step
	} else {
	  scale = this$biex$scale
	  label = this$biex$label
	  grid.step = this$biex$step
	} 

	tkconfigure(this$tt, cursor = "watch")
	

	data = this$data
	if ( nrow(data) == 0 ) {
		tkmessageBox(title = "An error has occured!",
			message = "Please select an area with bins.", icon = "error", type = "ok")
		stop("No area with bins selected.")
	} else {
		# if temporary rect data should be saved
		if (checkGATED=="1") {
			this$temp.num = this$temp.num + 1

			file = unlist(strsplit(file,"temp[0-9]{2}_"))
			file = file[[length(file)]]
			if (this$working) printf("w: do dotriploTRectData: file: %s",file)
			
			### only allow underscores and dots for special characters
			temp.name = sprintf("temp%02i_%s",this$temp.num,file)
			temp.name = sub(".fcs$","",temp.name)
			temp.name = sub(".csv$","",temp.name)
			temp.name = gsub("[^[:alnum:]_]","",temp.name)

			if (this$working) printf("w: do dotriploTRectData: temp.name: %s",temp.name)
			
			this$temptable.name[this$temp.num]=temp.name
			col1=as.matrix(rep(1,nrow(data)))
			colnames(col1) = "file_ID"
			data=sinh(data)
			data=cbind(col1,data)
			data=transform(data,as.numeric("file_ID"))

			dbWriteTable(this$conn, temp.name, data, overwrite=TRUE)
			print(paste("writing temporary",nrow(data),"cells in database",temp.name,".."))
			

			this$current.filenames = c(temp.name,this$current.filenames)

			# set file name
			tkconfigure(this$tkchoosefile,values=this$current.filenames)
			tkset(this$tkchoosefile,temp.name)

			this$temp.data=TRUE

			table = temp.name
			file.idx = 0
			
		} else {
			table = this$selected.project
			file.idx = this$current.filetable[which(this$current.filetable[,2]==file),1]
		}
	  
	  
	  ### get rect data
	  if (this$temp.data) {
	    this$getInfo(temp.name,vars)
	  } else {
	    this$getInfo(file,vars)
	  }
	  
		### remove ending from filename
		displayfile = this$shortenFilename(file)

		timeSTART = Sys.time()
		if ( this$working ) print("Time loading plot:")

		# start plot
		# mgp: A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. 
		# The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
		set.cex=1.2
		set.cex.axes=1.0
		set.mgp=c(1.7, 0.4, 0)
	
		if ( cutoffs[1] > 0 ) title.axis = sprintf("%s (%s)",var1,cutoffs[1])
		else title.axis = var1
		if ( cutoffs[2] > 0 ) title.axis = c(title.axis,sprintf("%s (%s)",var2,cutoffs[2]))
		else title.axis = c(title.axis,var2)

		plot(1,type='n',frame.plot=FALSE,xlim=c(xminval,xmaxval+10*binSize),axes=FALSE,
				ylim=c(yminval-2.5*binSize,ymaxval+5*binSize),xlab=title.axis[1],ylab=title.axis[2],cex.lab=set.cex,cex.axis=0.5*set.cex.axes,mgp=set.mgp)
		box(lwd=0.8,col="darkgrey")

		### draw axis on the bottom and on the left
		axis(side=1, at=scale,labels=label,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
		axis(side=2, at=scale,labels=label,las=3,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")

		### add grid
		if(checkGRID=="1") {
		  xgrid.steps=seq(0,(xmaxval),by=grid.step)
		  ygrid.steps=seq(0,(ymaxval),by=grid.step)
		  abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
		}

		# this$getInfo(file,vars)

		# get rect data to plot
		if (checkGATED=="1") tdata = this$data[c(var1,var2,var3)]
		else tdata = this$rectdata[c(var1,var2,var3)]
		
		ncells = ncells.total = nrow(tdata)
		quadrants.col = "black"
		# q1 Quadrant unten links
		# q2 Quadrant unten rechts
		# q3 Quadrant oben rechts
		# q4 Quadrant oben links
		if ( cutoffs[1] > 0 & cutoffs[2] > 0 ) {

			### count cells in quadrant
			tdata.q1 = tdata[which( tdata[,1]<cutoffs[1] &  tdata[,2]<cutoffs[2] ),3]
			tdata.q2 = tdata[which( tdata[,1]>=cutoffs[1] &  tdata[,2]<cutoffs[2] ),3]
			tdata.q3 = tdata[which( tdata[,1]>=cutoffs[1] &  tdata[,2]>=cutoffs[2] ),3]
			tdata.q4 = tdata[which( tdata[,1]<cutoffs[1] &  tdata[,2]>=cutoffs[2] ),3]

			### q[x].total [ink=black]
			### percentage of cells in quadrant to total cells 
			### or in MSI(+): percentage of cells in quadrant to total positive cells
			this$q1.total = abs(100 * length( tdata.q1 ) / ncells)
			this$q2.total = abs(100 * length( tdata.q2 ) / ncells)
			this$q3.total = abs(100 * length( tdata.q3 ) / ncells)
			this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
		
			if (cutoffs[3] > 0 ) {
				### number of cells which are producing cells in feature C
				ncells = nrow(tdata[which(tdata[,3]> cutoffs[3]),])

				### q[x].prodcells [ink=red]
				### percentage of cells which are positive for feature C in quadrant to total quadrant cells
				this$q1.prodcells = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[3])]) /length(tdata.q1)
				if (is.nan(this$q1.prodcells)) this$q1.prodcells = 0
				this$q2.prodcells = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[3])]) /length(tdata.q2)
				if (is.nan(this$q2.prodcells)) this$q2.prodcells = 0
				this$q3.prodcells = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[3])]) /length(tdata.q3)
				if (is.nan(this$q3.prodcells)) this$q3.prodcells = 0
				this$q4.prodcells = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[3])]) /length(tdata.q4)
				if (is.nan(this$q4.prodcells)) this$q4.prodcells = 0

				### only do MSI plots on producing cells only
				if ( FALSE ) { #checkCALC == "MSI(+)" ) {
					#this$bintriplot(tdata,cutoffs,binSize=binSize,mincells=mincount,bg=TRUE)
					#quadrants.col = "blue" 

					### cut all cells which are not producing cells
					tdata = tdata[which(tdata[,3]> cutoffs[3]),]

					tdata.q1 = tdata[which( tdata[,1]<cutoffs[1] &  tdata[,2]<cutoffs[2] ),3]
					tdata.q2 = tdata[which( tdata[,1]>=cutoffs[1] &  tdata[,2]<cutoffs[2] ),3]
					tdata.q3 = tdata[which( tdata[,1]>=cutoffs[1] &  tdata[,2]>=cutoffs[2] ),3]
					tdata.q4 = tdata[which( tdata[,1]<cutoffs[1] &  tdata[,2]>=cutoffs[2] ),3]

					### q[x].total [ink=blue]
					### in MSI(+): percentage of cells in quadrant to total positive cells
					this$q1.total = abs(100 * length( which( tdata[,1]<cutoffs[1] &  tdata[,2]<cutoffs[2] ) ) / ncells)
					this$q2.total = abs(100 * length( which( tdata[,1]>=cutoffs[1] &  tdata[,2]<cutoffs[2] ) ) / ncells)
					this$q3.total = abs(100 * length( which( tdata[,1]>=cutoffs[1] &  tdata[,2]>=cutoffs[2] ) ) / ncells)
					this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
				}
				### only do MSI plots on producing cells only
				if ( checkCALC == "MSI(+)" ) {
				  ### cut all cells which are not producing cells
				  tdata.plus = tdata[which(tdata[,3]> cutoffs[3]),]
				  ncells = nrow(tdata.plus)
				  
				  tdata.q1 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
				  tdata.q2 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
				  tdata.q3 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
				  tdata.q4 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
				  
				  ### q[x].total [ink=blue]
				  ### in MSI(+): percentage of cells in quadrant to total positive cells
				  this$q1.total = abs(100 * length( which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
				  this$q2.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
				  this$q3.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ) ) / ncells)
				  this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
				}

				### q[x].prodcellsplus [ink=green]
				### percentage of cells which are positive for feature C to total cells
				this$q1.prodcellsplus = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[3])]) / ncells.total
				if (is.nan(this$q1.prodcellsplus)) this$q1.prodcellsplus = 0
				this$q2.prodcellsplus = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[3])]) / ncells.total
				if (is.nan(this$q2.prodcellsplus)) this$q2.prodcellsplus = 0
				this$q3.prodcellsplus = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[3])]) / ncells.total
				if (is.nan(this$q3.prodcellsplus)) this$q3.prodcellsplus = 0
				this$q4.prodcellsplus = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[3])]) / ncells.total
				if (is.nan(this$q4.prodcellsplus)) this$q4.prodcellsplus = 0

				if (this$working) {
					print("w: ncells ::  q1  q2  q3  q4")
					printf("w: %s(total) :: %s %s %s %s",ncells.total,length(tdata.q1),length(tdata.q2),length(tdata.q3),length(tdata.q4))
					printf("w: %s(black/blue) :: %.1f %.1f %.1f %.1f",ncells,this$q1.total,this$q2.total,this$q3.total,this$q4.total)
					printf("w: prodcells(red) :: %.1f %.1f %.1f %.1f",this$q1.prodcells,this$q2.prodcells,this$q3.prodcells,this$q4.prodcells)
					printf("w: prodcellsplus(green) :: %.1f %.1f %.1f %.1f",this$q1.prodcellsplus,this$q2.prodcellsplus,this$q3.prodcellsplus,this$q4.prodcells)
				}
			}
		} 


		# get only the cells which are greater than 0
		tdata.zero = tdata[ tdata[,1]>=0 & tdata[,2]>=0, ]
		ncells.zero = nrow(tdata.zero)

		### calc quadrants with only positive values
		# q1 Quadrant unten links
		# q2 Quadrant unten rechts
		# q3 Quadrant oben rechts
		# q4 Quadrant oben links
		if ( cutoffs[1] > 0 & cutoffs[2] > 0 ) {
			q1.zero = abs(100 * length( which( tdata.zero[,1]<cutoffs[1] &  tdata.zero[,2]<cutoffs[2] ) ) / ncells.zero)      
			q2.zero = abs(100 * length( which( tdata.zero[,1]>=cutoffs[1] &  tdata.zero[,2]<cutoffs[2] ) ) / ncells.zero)                    
			q3.zero = abs(100 * length( which( tdata.zero[,1]>=cutoffs[1] &  tdata.zero[,2]>=cutoffs[2] ) ) / ncells.zero)                    
			q4.zero = abs(100 - q1.zero - q2.zero - q3.zero)
		}

		if (checkGATED != "1") this$origin.ncells=ncells 
		tkconfigure(this$ncell.sel.gui,text=as.character(ncells))
		tkconfigure(this$ncell.sel.gui.di,text=as.character(ncells))

		# plot rect data
		if (checkCALC == "density") {
			this$bintriplot(data=tdata,cutoffs=cutoffs,density=TRUE, binSize=binSize,mincells=mincount)
		} else if (checkCALC == "MSI(+)") {
			this$bintriplot(data=tdata.plus,cutoffs=cutoffs,binSize=binSize,mincells=mincount, quadrants.color = "blue", data.origin=tdata)
		} else {
			this$bintriplot(data=tdata,cutoffs=cutoffs,binSize=binSize,mincells=mincount)
		}
		#this$bintriplot(tdata,cutoffs,binSize=binSize,mincells=mincount, quadrants.color = quadrants.col)

		# attr = this$plot.attr[[1]]
		this$ncell.perc=ncells/this$origin.ncells*100
		if ( checkGATED == "1") {
			### remove ending from filename
			firstLine=sprintf("%s*(%0.1f%%): %s/cof=%s",displayfile,this$ncell.perc,checkCALC,this$current.cofactor)
		} else {
			firstLine=sprintf("%s(%0.1f%%): %s/cof=%s",displayfile,this$ncell.perc,checkCALC,this$current.cofactor)
		}
		title(main=firstLine,line=3.2,cex.main=0.9,adj=0)

		#if (checkCALC == "freq" | checkCALC == "MSI(+)") secondLine=sprintf("cells(min/max)=%s/%s; cutoff=%s",mincount,this$maxcells, cutoffs[3])
		if (checkCALC == "freq" | grepl("MSI",checkCALC)) secondLine=sprintf("cells(min/max)=%s/%s; %s",mincount,this$maxcells,var3)
		else secondLine=sprintf("cells(min/max)=%s/%s",mincount,this$maxcells)
		title(main=secondLine,line=2.4,cex.main=0.9,adj=0)
						
		thirdLine=sprintf("%s-%s(%0.1f%%); binSize=%s,#bins=%s",ncells,ncells.zero,(ncells.zero/ncells*100),binSize,this$bincount)
		title(main=thirdLine,line=1.6,cex.main=0.7,adj=0)


		if (checkGATED == "1") {
			coords.string = paste(var1,"[",this$coords$x[1],",",this$coords$x[2],"]",var2,"[",this$coords$y[1],",",this$coords$y[2],"]",sep="")
			print(coords.string)
			this$coords.info = c(this$coords.info,coords.string)
		}

		if ( length(this$coords.info)>2 ) {
			fourthLine=sprintf("%s; %s",this$coords.info[1],this$coords.info[2])
			fifthLine=sprintf("%s; %s",this$coords.info[3],this$coords.info[4])
		} else {
			fourthLine=sprintf("%s",paste(this$coords.info,collapse=";"))
			fifthLine=""
		}
		title(main=fourthLine,line=0.9,cex.main=0.6,adj=0)
		title(main=fifthLine,line=0.2,cex.main=0.6,adj=0)


		this$plot.num = this$plot.num + 1
		# push plot attributes one up
		# for (i in this$plot.num:1) {
		# 		this$plot.attr[i+1] = this$plot.attr[i]
		# }

		if (FALSE) {
			# pngs for history
			png(file=this$png.file,width=140,height=145,bg="transparent")
			# mar in points, mai in inches
			# oma adds title lines
			# order: bottom, left, top, and right
			par(mar=c(1.1,1.1,1.3,0.3))
			if ( cutoffs[1] > 0 ) title.axis = sprintf("%s (%s)",var1,cutoffs[1])
			else title.axis = var1
			if ( cutoffs[2] > 0 ) title.axis = c(title.axis,sprintf("%s (%s)",var2,cutoffs[2]))
			else title.axis = c(title.axis,var2)

			plot(1,type='n',frame.plot=FALSE,xlim=c(xminval,xmaxval+10*binSize),axes=FALSE,
					ylim=c(yminval-2.5*binSize,ymaxval+5*binSize),xlab=title.axis[1],ylab=title.axis[2],cex.lab=set.cex*0.9,cex.axis=0.5*set.cex.axes,mgp=c(1.9,0.7,0))
			box(lwd=0.8,col="darkgrey")
			
			### add grid
			if(checkGRID=="1") {
			  xgrid.steps=seq(0,(xmaxval),by=grid.step)
			  ygrid.steps=seq(0,(ymaxval),by=grid.step)
			  abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
			}

			this$bintriplot(tdata,cutoffs,set.cex=0.8,set.cex.axes=0.2,set.mgp=c(1.9,0.7,0),binSize=binSize,png=TRUE, quadrants.color = quadrants.col)

			title(main=sprintf("%s: %s",displayfile,checkCALC),line=0.5,cex.main=0.8,adj=0)
			dev.off()

			this$images[[this$plot.num]]=tclVar()
			tkimage.create("photo",this$images[[this$plot.num]],file=this$png.file)
		}

		if (FALSE) {
		### save plot attributes
		this$plot.attr[[1]]=list(
			table=table,
			origin.table=this$current.project,
			file.name=displayfile,
			file.idx=file.idx,
			stain.idx=this$current.filenum,
			origin.ncells=this$origin.ncells,
			ncells.zero=ncells.zero,
			sel.ncells=ncells,
			binSize=binSize,
			bincount=this$bincount,
			mincount=mincount,
			maxcount=this$maxcells,
			vars=c(var1,var2,var3),
			FI=as.character(c(xminval,xmaxval,yminval,ymaxval)),
			MSI=c(as.double(tclvalue(this$vminMSI)),as.double(tclvalue(this$vmaxMSI))),
			coords=this$coords.info,
			dynrange=tclvalue(this$cbtdynRange),
			calc=checkCALC,
			grid=checkGRID,
			transformation=checkTRANS,
			xcutoff=cutoffs[1],
			ycutoff=cutoffs[2],
			zcutoff=cutoffs[3],
			CellInfo=NULL,
			RectInfo=NULL
		)
		if (tclvalue(this$cbtautoRect)=="1") {
			this$plot.attr[[1]]$CellInfo=TRUE
			this$plot.attr[[1]]$RectInfo=TRUE
		}
		}


		### save FI ranges for this transformation type
		# this$changeFI.range(mode=3)
	}
	tkconfigure(this$tt, cursor = "left_ptr")
}

fcs$dotriploTOverview <- function(table=NA) {
	this=fcs
	answer = tclVar("yes")
	freqans = tclVar("no")


	### if axes ranges are not the same
	if ( !this$checkAxesRange() ) {
		tkmessageBox(title = "An error has occured!",
			message = "Please set x and y axis with the same range.", icon = "error", type = "ok")
		stop("Set x and y axis with the same range.")
	}

	if (!this$working) answer=tkmessageBox(title = "Starting..", message = "Are you sure?", icon = "info", type = "yesno")

	if (tclvalue(answer) == "yes") {
		checkCALC = tclvalue(this$rbcalc)
		checkTRANS = tclvalue(this$rbtrans)
		checkGATED = tclvalue(this$cbtgateData)
		checkGRID = tclvalue(this$cbtshowGrid)
		checkTRIMMING = tclvalue(this$cbttrimming)

		if ( checkTRANS == "asinh" ) {
		  scale = this$asinh$scale
		  label = this$asinh$label
		  grid.step = this$asinh$step
		} else {
		  scale = this$biex$scale
		  label = this$biex$label
		  grid.step = this$biex$step
		} 

		quadrants.col = "black"

		len=length(this$selected.vars)

		# column and cutoff number vector
		colvec=vector()
		cutoffs=vector()
		for (i in 1:len) {
			if ( tclvalue(this$cbVal[[i]]) == "1") {
				colvec=c(colvec,i);
				cutoffs=c(cutoffs,this$checkDigits(cutoff_id=i))
			}
		}
		len.colvec = length(colvec)
		# title length
		titlelen = length(unlist(strsplit(tclvalue(tkget(this$title,"1.0","end-1c")),"")))
		titleans = tclVar("yes")

		if ( len.colvec < 3 ) {
			tkmessageBox(title = "An error has occured!",
				message = "Please select at least three markers.", icon = "error", type = "ok")
			stop("Please select at least three markers.")
		}
		if ( titlelen <= 1 & this$working==FALSE) {
			titleans = tkmessageBox(title = "Are you sure?",
				message = "Do not forget to set page title. Continue?", icon = "info", type = "yesno")
		} else if ( checkCALC == "freq" & !is.integer(which(cutoffs==0)) & this$working==FALSE ) {
			freqans = tkmessageBox(title = "Are you sure?",
				message = "Do not forget to set cutoffs. Continue?", icon = "info", type = "yesno")
		} else if (checkCALC=="density") {
		  tkmessageBox(title = "For your information.",
		               message = "It is not meaningful to plot method \"density\".", icon = "error", type = "ok")
		  stop("Chose another statistical method.")
		}

		if ( (tclvalue(titleans) == "yes") | tclvalue(freqans) == "yes" ) {
			if (this$working) {
				this$saveinFolder=getwd()
			} else {
				# if path didnt save yet
				if (!exists("lastpath",where=eval(parse(text="fcs")))){
					this$saveinFolder = tk_choose.dir(default = getwd(), caption = "Select directory to save the PDFs")
				} else {
					this$saveinFolder = tk_choose.dir(default = this$lastpath, caption = "Select directory to save the PDFs")
				}
			}

			if ( !is.na(this$saveinFolder) ) {
				tkconfigure(this$tt, cursor = "watch")

				this$lastpath = this$saveinFolder

				max.nhorizplots = 6
				max.nvertiplots = 8

				res.nhorizplots = (len.colvec-2)%%max.nhorizplots
				if (res.nhorizplots != 0) res.nhorizplots = max.nhorizplots - res.nhorizplots
				max.nplots = len.colvec*(len.colvec-1)*(len.colvec-2)
				max.npages = ceiling( (max.nplots+res.nhorizplots*(len.colvec-1)*len.colvec)/(max.nhorizplots*max.nvertiplots) ) # + 1

				### Create Progress Bar
				this$pb <- tkProgressBar(title="triploT-Overview", label=paste("Creating page 1/",max.npages,sep=""),
					min = 0, max = max.npages, initial = 1, width = 300)
				# set progress bar
				setTkProgressBar(this$pb,1,label=paste("Creating page 1/",max.npages,sep=""))
				
				timeSTART = Sys.time()
				if (!this$OverviewGate) cat("\n\n>>>> Start triploT-Overview with total data: \n\n")
				else cat("\n\n>>>> Start triploT overview with gated data: \n\n")
				if (this$working) print(sprintf("w: %s - time started",timeSTART))

				tmp.folder = file.path(this$saveinFolder,"tmp")
				if ( grepl("linux",sessionInfo()$R.version$os) ) {
					### on linux
					system(paste("mkdir -p -v ",tmp.folder,sep=""))
				} else {
					### on Windows
					dir.create(tmp.folder)
					printf(">>on Windows<< Creating folder: %s",tmp.folder)
					#system(paste("mkdir ",tmp.folder,sep=""))
					# or shell()?
				}

				file=tclvalue(tkget(this$tkchoosefile))
				if ( this$OverviewGate ) {
				displayfile = this$shortenFilename(this$plot.attr[[1]]$file.name,title=TRUE)
				} else {
						displayfile = this$shortenFilename(file)
				}

				table = this$total.projects[this$selected.projectnum]
				if ( grepl("temp",file) ) {
					file.idx = 1
					table=file
				} else {
					file.idx=this$current.filetable[which(this$current.filetable[,2]==file),1]
					this$selected.filenum = file.idx
				}

				if (this$working) printf("w: table=%s file.idx=%s",table,file.idx)

				binSize=as.numeric(tkget(this$binSize))
				mincount=as.numeric(tkget(this$minCountTri))
				xminval=as.numeric(tkget(this$minvalX))
				xmaxval=as.numeric(tkget(this$maxvalX))
				yminval=as.numeric(tkget(this$minvalY))
				ymaxval=as.numeric(tkget(this$maxvalY))

				#### Do not trim and remove doublets if data is gated
				if (this$OverviewGate) {
					# if data is gated, just recall
					data = this$getData(table,file.idx,columns=colvec)
				} else if ( checkTRIMMING == "1" ) {
					this$preprocData(mode="trim")
					data = this$data[,colvec]
				} else {
					this$getFile(table,file.idx)
					data = this$data[,colvec]
				}

				len = dim(data)[2]

				printf("Columns selected (%s): %s",dim(data)[2],paste(colvec,collapse=" "))
				printf("Features selected: %s",paste(colnames(data),collapse=" "))

				title=as.character(tclvalue(tkget(this$title,"1.0","end-1c")))

				### if less than 16 markers were selected
				# do all triplots in one pdf file
				# otherwise make one pdf file for one page
	            if (len.colvec<=16){

					toPDF(file=sprintf("%s_triploTOverview_%s%s_%s.pdf",displayfile,len.colvec,checkCALC,this$version), 
						path=this$saveinFolder, 
						title=sprintf("triploTOverview of %s",displayfile), 
						width=3.21*max.nhorizplots,
						height=3.5*max.nvertiplots,
						pointsize=11,
						{
							# plot histograms
							#this$plotHistograms(plotter="triover",pdf=TRUE)

							# set progress bar
							#setTkProgressBar(this$pb,2,label=paste("Creating page 2/",len.colvec+2,sep=""))
							# plot density plot
							#this$plotDensities(plotter="triover",pdf=TRUE)

							##### start triploT overview
							### set label and axes font sizes

							label.cex = 1.2
							set.cex.axes = 1
							set.mgp = c(1.9,0.5,0)
							par(mfrow=c(max.nvertiplots,max.nhorizplots),oma=c(0.5,1,6,1),mar=c(3,4,5,1))

							
							plot.idx = 0
							page.idx = 0
							current.page = 0
							for ( v1 in 1:len ) {
								for ( v2 in 1:len) {

									if ( v1 == v2 ) next

									for ( v3 in 1:len) {

										## skip loop if one of the axes are the same
										if ( v1 == v3 | v2 == v3 ) next
						
										#select columns to plot
										tdata=as.matrix(data[,c(v1,v2,v3)])
										this$cnames = colnames(tdata)

										if ( cutoffs[v1] > 0 ) title.axis = sprintf("%s (%s)",colnames(data)[v1],cutoffs[v1])
										else title.axis = colnames(data)[v1]
										if ( cutoffs[v2] > 0 ) title.axis = c(title.axis,sprintf("%s (%s)",colnames(data)[v2],cutoffs[v2]))
										else title.axis = c(title.axis,colnames(data)[v2])
										# start plot
										# mgp: A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. 
										# The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
										plot(1,type='n',frame.plot=FALSE,xlim=c(xminval,xmaxval+10*binSize),axes=FALSE
											,ylim=c(yminval-2.5*binSize,ymaxval+5*binSize),xlab=title.axis[1]
											,ylab=title.axis[2],cex.lab=label.cex,cex.axis=set.cex.axes,mgp=set.mgp)
										box(lwd=0.5,col="darkgrey")

										plot.idx = plot.idx + 1

										### draw axis on the bottom and on the left
										axis(side=1, at=scale,labels=label,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
										axis(side=2, at=scale,labels=label,las=3,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")

										### add grid
										if(checkGRID=="1") {
										  xgrid.steps=seq(0,(xmaxval),by=grid.step)
										  ygrid.steps=seq(0,(ymaxval),by=grid.step)
										  abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
										}


										ncells = ncells.total = nrow(tdata)
										# q1 Quadrant unten links
										# q2 Quadrant unten rechts
										# q3 Quadrant oben rechts
										# q4 Quadrant oben links
										if ( cutoffs[v1] > 0 & cutoffs[v2] > 0 ) {

											### count cells in quadrant
											tdata.q1 = tdata[which( tdata[,1]<cutoffs[v1] &  tdata[,2]<cutoffs[v2] ),3]
											tdata.q2 = tdata[which( tdata[,1]>=cutoffs[v1] &  tdata[,2]<cutoffs[v2] ),3]
											tdata.q3 = tdata[which( tdata[,1]>=cutoffs[v1] &  tdata[,2]>=cutoffs[v2] ),3]
											tdata.q4 = tdata[which( tdata[,1]<cutoffs[v1] &  tdata[,2]>=cutoffs[v2] ),3]

											### q[x].total [ink=black]
											### percentage of cells in quadrant to total cells 
											### or in MSI(+): percentage of cells in quadrant to total positive cells
											this$q1.total = abs(100 * length( tdata.q1 ) / ncells)
											this$q2.total = abs(100 * length( tdata.q2 ) / ncells)
											this$q3.total = abs(100 * length( tdata.q3 ) / ncells)
											this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
										
											if (cutoffs[v3] > 0 ) {
												### number of cells which are producing cells in feature C
												ncells = nrow(tdata[which(tdata[,3]> cutoffs[v3]),])

												### q[x].prodcells [ink=red]
												### percentage of cells which are positive for feature C in quadrant to total quadrant cells
												this$q1.prodcells = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[v3])]) /length(tdata.q1)
												if (is.nan(this$q1.prodcells)) this$q1.prodcells = 0
												this$q2.prodcells = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[v3])]) /length(tdata.q2)
												if (is.nan(this$q2.prodcells)) this$q2.prodcells = 0
												this$q3.prodcells = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[v3])]) /length(tdata.q3)
												if (is.nan(this$q3.prodcells)) this$q3.prodcells = 0
												this$q4.prodcells = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[v3])]) /length(tdata.q4)
												if (is.nan(this$q4.prodcells)) this$q4.prodcells = 0

												### only do MSI plots on producing cells only
												if ( checkCALC == "MSI(+)" ) {
												  ### cut all cells which are not producing cells
												  tdata.plus = tdata[which(tdata[,3]> cutoffs[3]),]
												  ncells = nrow(tdata.plus)
												  
												  tdata.q1 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
												  tdata.q2 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
												  tdata.q3 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
												  tdata.q4 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
												  
												  ### q[x].total [ink=blue]
												  ### in MSI(+): percentage of cells in quadrant to total positive cells
												  this$q1.total = abs(100 * length( which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
												  this$q2.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
												  this$q3.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ) ) / ncells)
												  this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
												}
												
												### q[x].prodcellsplus [ink=green]
												### percentage of cells which are positive for feature C to total cells
												this$q1.prodcellsplus = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[v3])]) / ncells.total
												if (is.nan(this$q1.prodcellsplus)) this$q1.prodcellsplus = 0
												this$q2.prodcellsplus = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[v3])]) / ncells.total
												if (is.nan(this$q2.prodcellsplus)) this$q2.prodcellsplus = 0
												this$q3.prodcellsplus = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[v3])]) / ncells.total
												if (is.nan(this$q3.prodcellsplus)) this$q3.prodcellsplus = 0
												this$q4.prodcellsplus = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[v3])]) / ncells.total
												if (is.nan(this$q4.prodcellsplus)) this$q4.prodcellsplus = 0
											}
										} 


										tdata.zero = tdata[apply(tdata[, -1], MARGIN = 1, function(row) {all(row > 0)}), ]
										ncells.zero=nrow(tdata.zero)

										### calc quadrants with only positive values
										# q1 Quadrant unten links
										# q2 Quadrant unten rechts
										# q3 Quadrant oben rechts
										# q4 Quadrant oben links
										if ( cutoffs[v1] > 0 & cutoffs[v2] > 0 ) {
											q1.zero = abs(100 * length( which( tdata.zero[,1]<cutoffs[v1] &  tdata.zero[,2]<cutoffs[v2] ) ) / ncells.zero)      
											q2.zero = abs(100 * length( which( tdata.zero[,1]>=cutoffs[v1] &  tdata.zero[,2]<cutoffs[v2] ) ) / ncells.zero)                    
											q3.zero = abs(100 * length( which( tdata.zero[,1]>=cutoffs[v1] &  tdata.zero[,2]>=cutoffs[v2] ) ) / ncells.zero)                    
											q4.zero = abs(100 - q1.zero - q2.zero - q3.zero)
										}

										#this$bintriplot(tdata,cutoffs[c(v1,v2,v3)],set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp
										#	,binSize=binSize,mincells=mincount, quadrants.color = quadrants.col)
										if (checkCALC == "density") {
											this$bintriplot(data=tdata,cutoffs=cutoffs[c(v1,v2,v3)],density=TRUE,binSize=binSize,mincells=mincount,overview=TRUE
												,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
										} else if (checkCALC == "MSI(+)") {
											this$bintriplot(data=tdata.plus,cutoffs=cutoffs[c(v1,v2,v3)],binSize=binSize,mincells=mincount,,overview=TRUE
												,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp,quadrants.color = "blue", data.origin=tdata)
										} else {
											this$bintriplot(data=tdata,cutoffs=cutoffs[c(v1,v2,v3)],binSize=binSize,mincells=mincount,overview=TRUE
												,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
										}


										if (this$OverviewGate) {this$ncell.perc = this$ncell.sel/this$origin.ncells*100}

										# add title for single plot
										if ( checkCALC == "freq" | grepl("MSI",checkCALC) ) {
											if ( this$OverviewGate ) firstLine=sprintf("%s(%0.1f%%): %s,%s/cof=%s",displayfile,this$ncell.perc,checkCALC,v3,this$current.cofactor)
											else firstLine=sprintf("%s: %s,cutoff=%s/cof=%s",displayfile,checkCALC,cutoffs[v3],this$current.cofactor)
										} else {
											if ( this$OverviewGate ) firstLine=sprintf("%s(%0.1f%%): %s/cof=%s",displayfile,this$ncell.perc,checkCALC,this$current.cofactor)
											else firstLine=sprintf("%s: %s/cof=%s",displayfile,checkCALC,this$current.cofactor)
										}
										title(main=firstLine,line=3.2,cex.main=1.0*label.cex,adj=0)

										secondLine=sprintf("%s-%s(%0.1f%%); min/max=%s/%s; binSize=%s,#bins=%s",ncells,ncells.zero,(ncells.zero/ncells*100),mincount,this$maxcells,binSize,this$bincount)
										title(main=secondLine,line=2.4,cex.main=0.7*label.cex,adj=0)

										thirdLine=""
										if ( cutoffs[v1] > 0 & cutoffs[v2] > 0 ) {
											thirdLine=sprintf("Q1=%0.1f/%0.1f Q2=%0.1f/%0.1f Q3=%0.1f/%0.1f Q4=%0.1f/%0.1f",q1.zero,this$q1.total,q2.zero,this$q2.total,q3.zero,this$q3.total,q4.zero,this$q4.total)
										}
										title(main=thirdLine,line=1.6,cex.main=0.7*label.cex,adj=0)

										if (checkGATED == "1" | grepl("temp",file)) {
											if ( length(this$coords.info)>2 ) {
												fourthLine=sprintf("%s; %s",this$coords.info[1],this$coords.info[2])
												fifthLine=sprintf("%s; %s",this$coords.info[3],this$coords.info[4])
											} else {
												fourthLine=sprintf("%s",paste(this$coords.info,collapse="; "))
												fifthLine=""
											}
											title(main=fourthLine,line=0.9,cex.main=0.6*label.cex,adj=0)
											title(main=fifthLine,line=0.1,cex.main=0.6*label.cex,adj=0)
										}

										### update progress bar
										if ( (plot.idx != 0) & ((plot.idx-1) %% (max.nhorizplots*max.nvertiplots) == 0) ) {
											page.idx = page.idx + 1
											info <- paste("Creating page ",page.idx,"/",max.npages,sep="")
											if (this$working) printf("%s with plot# %s",info,plot.idx)
											setTkProgressBar(this$pb, page.idx, label=info)
										}

										#add title for page in 3D-Overview tab
										if (current.page != page.idx) {
											mtext(title, outer = TRUE, cex = 1.5,line=1.0,pos=2,xpd=TRUE)
											if (length(this$coords.info>0)) mtext(sprintf("(%s)",paste(this$coords.info,collapse="; ")),outer=TRUE,cex = 0.8)
											#printf("Print on Page %s",page.idx)
											current.page = page.idx
										}
									}

									### plot empty plots
									if ( (plot.idx != 0) & (res.nhorizplots > 0) ) {
										for ( i in 1:res.nhorizplots ) {
											plot.new(); 
											plot.idx = plot.idx+1
										} 
									}

								}
							}

							### plot histograms
							this$plotHistograms(plotter="triover",pdf=TRUE)

							tkmessageBox(title = "Output of overview",
								message = paste ("PDF created in folder ",this$saveinFolder, sep=""))

						}
					)
				} else {
					toPDF(file=sprintf("%s_triploTOverview_%shist_%s.pdf",displayfile,len.colvec,this$version), 
						path=this$saveinFolder, 
						title=sprintf("Histograms of %s",displayfile), 

						width=3.21*max.nhorizplots,
						height=3.5*max.nvertiplots,
						pointsize=11,
						{
							# plot histograms
							this$plotHistograms(plotter="triover",pdf=TRUE)
						  
						  # set progress bar
							#setTkProgressBar(this$pb,2,label=paste("Creating page 2/",max.npages,sep=""))

							# plot density plot
							#this$plotDensities(plotter="triover",pdf=TRUE)
						}
					)
					
					for ( v1 in 1:len ) {
						printf("New PDF: #%s",v1)
						plot.idx = 0
						page.idx = 0

						### update progress bar
						#info <- paste("Creating page ",page.idx+1,"/",max.npages,sep="")
						#print(info)
						
						#setTkProgressBar(this$pb, v1+1, label=info)
						toPDF(file=sprintf("%s_triploTOverview_%s_%s%s_%s.pdf",displayfile,colnames(data)[v1],len.colvec-1,checkCALC,this$version), 
						#toPDF(file=sprintf("%s_triploTOverview_%s%s_%s_%s.pdf",displayfile,len.colvec,checkCALC,colnames(data)[v1],this$version), 
							path=this$saveinFolder, 
							title=sprintf("triploTOverview of %s(%s)",displayfile,colnames(data)[v1]), 

							width=3.21*max.nhorizplots,
							height=3.5*max.nvertiplots,
							pointsize=11,
							{

								### new
								label.cex = 1.2
								set.cex.axes = 1
								set.mgp = c(1.9,0.5,0)
								par(mfrow=c(max.nvertiplots,max.nhorizplots),oma=c(0.5,1,5,1),mar=c(3,4,5,1))
								###
								
								for ( v2 in 1:len) {
									if ( v1 == v2 ) next
									for ( v3 in 1:len) {

										## skip loop if one of the axes are the same
										if ( v1 == v3 | v2 ==v3 ) next
										plot.idx = plot.idx + 1
						
										### update progress bar
										#if (plot.idx %% (max.nhorizplots*max.nvertiplots) == 0 ) {
										#	page.idx = page.idx + 1
										#	info <- paste("Creating page ",page.idx,"/",max.npages,sep="")
										#	if (this$working) printf("%s with plot# %s",info,plot.idx)
										#	setTkProgressBar(this$pb, page.idx, label=info)
										#}

										#else {
										#select columns to plot
										tdata=as.matrix(data[,c(v1,v2,v3)])
										this$cnames = c(colnames(data)[v1],colnames(data)[v2],colnames(data)[v3])

										if ( cutoffs[v1] > 0 ) title.axis = sprintf("%s (%s)",colnames(data)[v1],cutoffs[v1])
										else title.axis = colnames(data)[v1]
										if ( cutoffs[v2] > 0 ) title.axis = c(title.axis,sprintf("%s (%s)",colnames(data)[v2],cutoffs[v2]))
										else title.axis = c(title.axis,colnames(data)[v2])

										# start plot
										# mgp: A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. 
										# The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
										plot(1,type='n',frame.plot=FALSE,xlim=c(xminval,xmaxval+10*binSize),axes=FALSE
											,ylim=c(yminval-2.5*binSize,ymaxval+5*binSize),xlab=title.axis[1]
											,ylab=title.axis[2],cex.lab=label.cex,cex.axis=set.cex.axes,mgp=set.mgp)
										box(lwd=0.5,col="darkgrey")
										
										#add title for page in 3D-Overview tab
										#mtext(title, outer = TRUE, cex = 1.5,line=1.3,pos=2,xpd=TRUE)
										#if (length(this$coords.info>0)) mtext(sprintf("(%s)",paste(this$coords.info,collapse="; ")),outer=TRUE,cex = 0.8)
										
										### draw axis on the bottom and on the left
										axis(side=1, at=scale,labels=label,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
										axis(side=2, at=scale,labels=label,las=3,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")

										### add grid
										if(checkGRID=="1") {
										  xgrid.steps=seq(0,(xmaxval),by=grid.step)
										  ygrid.steps=seq(0,(ymaxval),by=grid.step)
										  abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
										}


										### calc quadrants in total
										ncells = ncells.total = nrow(tdata)
										# q1 Quadrant unten links
										# q2 Quadrant unten rechts
										# q3 Quadrant oben rechts
										# q4 Quadrant oben links
										if ( cutoffs[v1] > 0 & cutoffs[v2] > 0 ) {

											### count cells in quadrant
											tdata.q1 = tdata[which( tdata[,1]<cutoffs[v1] &  tdata[,2]<cutoffs[v2] ),3]
											tdata.q2 = tdata[which( tdata[,1]>=cutoffs[v1] &  tdata[,2]<cutoffs[v2] ),3]
											tdata.q3 = tdata[which( tdata[,1]>=cutoffs[v1] &  tdata[,2]>=cutoffs[v2] ),3]
											tdata.q4 = tdata[which( tdata[,1]<cutoffs[v1] &  tdata[,2]>=cutoffs[v2] ),3]

											### q[x].total [ink=black]
											### percentage of cells in quadrant to total cells 
											### or in MSI(+): percentage of cells in quadrant to total positive cells
											this$q1.total = abs(100 * length( tdata.q1 ) / ncells)
											this$q2.total = abs(100 * length( tdata.q2 ) / ncells)
											this$q3.total = abs(100 * length( tdata.q3 ) / ncells)
											this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
										
											if (cutoffs[v3] > 0 ) {
												### number of cells which are producing cells in feature C
												ncells = nrow(tdata[which(tdata[,3]> cutoffs[v3]),])

												### q[x].prodcells [ink=red]
												### percentage of cells which are positive for feature C in quadrant to total quadrant cells
												this$q1.prodcells = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[v3])]) /length(tdata.q1)
												if (is.nan(this$q1.prodcells)) this$q1.prodcells = 0
												this$q2.prodcells = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[v3])]) /length(tdata.q2)
												if (is.nan(this$q2.prodcells)) this$q2.prodcells = 0
												this$q3.prodcells = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[v3])]) /length(tdata.q3)
												if (is.nan(this$q3.prodcells)) this$q3.prodcells = 0
												this$q4.prodcells = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[v3])]) /length(tdata.q4)
												if (is.nan(this$q4.prodcells)) this$q4.prodcells = 0

												### only do MSI plots on producing cells only
												if ( checkCALC == "MSI(+)" ) {
												  ### cut all cells which are not producing cells
												  tdata.plus = tdata[which(tdata[,3]> cutoffs[3]),]
												  ncells = nrow(tdata.plus)
												  
												  tdata.q1 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
												  tdata.q2 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
												  tdata.q3 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
												  tdata.q4 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
												  
												  ### q[x].total [ink=blue]
												  ### in MSI(+): percentage of cells in quadrant to total positive cells
												  this$q1.total = abs(100 * length( which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
												  this$q2.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
												  this$q3.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ) ) / ncells)
												  this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
												}

												### q[x].prodcellsplus [ink=green]
												### percentage of cells which are positive for feature C to total cells
												this$q1.prodcellsplus = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[v3])]) / ncells.total
												if (is.nan(this$q1.prodcellsplus)) this$q1.prodcellsplus = 0
												this$q2.prodcellsplus = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[v3])]) / ncells.total
												if (is.nan(this$q2.prodcellsplus)) this$q2.prodcellsplus = 0
												this$q3.prodcellsplus = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[v3])]) / ncells.total
												if (is.nan(this$q3.prodcellsplus)) this$q3.prodcellsplus = 0
												this$q4.prodcellsplus = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[v3])]) / ncells.total
												if (is.nan(this$q4.prodcellsplus)) this$q4.prodcellsplus = 0
											}
										} 

										tdata.zero = tdata[apply(tdata[, -1], MARGIN = 1, function(row) {all(row > 0)}), ]
										ncells.zero=nrow(tdata.zero)

										### calc quadrants with only positive values
										# q1 Quadrant unten links
										# q2 Quadrant unten rechts
										# q3 Quadrant oben rechts
										# q4 Quadrant oben links
										if ( cutoffs[v1] > 0 & cutoffs[v2] > 0 ) {
											q1.zero = abs(100 * length( which( tdata.zero[,1]<cutoffs[v1] &  tdata.zero[,2]<cutoffs[v2] ) ) / ncells.zero)      
											q2.zero = abs(100 * length( which( tdata.zero[,1]>=cutoffs[v1] &  tdata.zero[,2]<cutoffs[v2] ) ) / ncells.zero)                    
											q3.zero = abs(100 * length( which( tdata.zero[,1]>=cutoffs[v1] &  tdata.zero[,2]>=cutoffs[v2] ) ) / ncells.zero)                    
											q4.zero = abs(100 - q1.zero - q2.zero - q3.zero)
										}

										#this$bintriplot(tdata,cutoffs[c(v1,v2,v3)],set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp
										#	,binSize=binSize,mincells=mincount, quadrants.color = quadrants.col)
										if (checkCALC == "density") {
											this$bintriplot(data=tdata,cutoffs=cutoffs[c(v1,v2,v3)],density=TRUE,binSize=binSize,mincells=mincount,overview=TRUE
												,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
										} else if (checkCALC == "MSI(+)") {
											this$bintriplot(data=tdata.plus,cutoffs=cutoffs[c(v1,v2,v3)],binSize=binSize,mincells=mincount,,overview=TRUE
												,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp,quadrants.color = "blue", data.origin=tdata)
										} else {
											this$bintriplot(data=tdata,cutoffs=cutoffs[c(v1,v2,v3)],binSize=binSize,mincells=mincount,overview=TRUE
												,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
										}


										if (this$OverviewGate) {this$ncell.perc = this$ncell.sel/this$origin.ncells*100}

										# add title for single plot
										if ( checkCALC == "freq" | grepl("MSI",checkCALC) ) {
											if ( this$OverviewGate ) firstLine=sprintf("%s(%0.1f%%): %s,%s/cof=%s",
											                                           displayfile,this$ncell.perc,checkCALC,v3,this$current.cofactor)
											else firstLine=sprintf("%s: %s,cutoff=%s/cof=%s",
											                       displayfile,checkCALC,cutoffs[v3],this$current.cofactor)
										} else {
											if ( this$OverviewGate ) firstLine=sprintf("%s(%0.1f%%): %s/cof=%s",
											                                           displayfile,this$ncell.perc,checkCALC,this$current.cofactor)
											else firstLine=sprintf("%s: %s/cof=%s",
											                       displayfile,checkCALC,this$current.cofactor)
										}
										title(main=firstLine,line=3.2,cex.main=1.0*label.cex,adj=0)

										secondLine=sprintf("%s-%s(%0.1f%%); min/max=%s/%s; binSize=%s,#bins=%s",ncells,ncells.zero,(ncells.zero/ncells*100),mincount,this$maxcells,binSize,this$bincount)
										title(main=secondLine,line=2.4,cex.main=0.7*label.cex,adj=0)

										thirdLine=""
										if ( cutoffs[v1] > 0 & cutoffs[v2] > 0 ) {
											thirdLine=sprintf("Q1=%0.1f/%0.1f Q2=%0.1f/%0.1f Q3=%0.1f/%0.1f Q4=%0.1f/%0.1f"
												,q1.zero,this$q1.total,q2.zero,this$q2.total,q3.zero,this$q3.total,q4.zero,this$q4.total)
										}
										title(main=thirdLine,line=1.6,cex.main=0.7*label.cex,adj=0)

										if (checkGATED == "1" | grepl("temp",file)) {
											if ( length(this$coords.info)>2 ) {
												fourthLine=sprintf("%s; %s",this$coords.info[1],this$coords.info[2])
												fifthLine=sprintf("%s; %s",this$coords.info[3],this$coords.info[4])
											} else {
												fourthLine=sprintf("%s",paste(this$coords.info,collapse="; "))
												fifthLine=""
											}
											title(main=fourthLine,line=0.9,cex.main=0.6*label.cex,adj=0)
											title(main=fifthLine,line=0.1,cex.main=0.6*label.cex,adj=0)
										}
									}
									### update progress bar
									if ( (plot.idx != 0) & ((plot.idx-1) %% (max.nhorizplots*max.nvertiplots) == 0) ) {
										page.idx = page.idx + 1
										info <- paste("Creating page ",page.idx,"/",max.npages,sep="")
										if (this$working) printf("%s with plot# %s",info,plot.idx)
										setTkProgressBar(this$pb, page.idx, label=info)
									}

									#add title for page in 3D-Overview tab
									if (current.page != page.idx) {
										mtext(title, outer = TRUE, cex = 1.5,line=1.0,pos=2,xpd=TRUE)
										if (length(this$coords.info>0)) mtext(sprintf("(%s)",paste(this$coords.info,collapse="; ")),outer=TRUE,cex = 0.8)
										#printf("Print on Page %s",page.idx)
										current.page = page.idx
									}
									
									### plot emptry plots
									printf("last plot #:%s with X:%s Y:%s; plotting %s more empty plots",plot.idx,v1,v2,res.nhorizplots)
									if (res.nhorizplots > 0) {
										for ( i in 1:res.nhorizplots ) {
											plot(1,type='n',frame.plot=FALSE,axes=FALSE,ylab="")
											plot.idx = plot.idx+1
										}
									}
								}
							}
						)
					}

					tkmessageBox(title = "Output of overview",
						message = paste ("Plots created in folder ",this$saveinFolder, sep=""))
				}
				
				tkconfigure(this$tt, cursor = "left_ptr")

				cat("\n")
				# close progress bar
				close(this$pb)
				
				#if (checkPNG == "1") printf("PNGs are saved in %s.",file.path(this$saveinFolder,"png"))

				print(Sys.time()-timeSTART)

				cat("\n>>>> Done.\n")
			}
		}
	}	
	this$OverviewGate=FALSE
}


fcs$dotriploTOverview_ALL <- function(table=NA) {
  this=fcs
  answer = tclVar("yes")
  freqans = tclVar("no")
  
  
  ### if axes ranges are not the same
  if ( !this$checkAxesRange() ) {
    tkmessageBox(title = "An error has occured!",
                 message = "Please set x and y axis with the same range.", icon = "error", type = "ok")
    stop("Set x and y axis with the same range.")
  }
  
  if (!this$working) answer = tkmessageBox(title = "Starting..", 
                                           message = "Are you sure? Did you check the setting with fixed features?", 
                                           icon = "info", type = "yesno")
  
  if (tclvalue(answer) == "no") stop("OK.");
  
  # if (tclvalue(answer) == "yes") {
    checkCALC = tclvalue(this$rbcalc)
    checkTRANS = tclvalue(this$rbtrans)
    checkGATED = tclvalue(this$cbtgateData)
    checkGRID = tclvalue(this$cbtshowGrid)
    checkTRIMMING = tclvalue(this$cbttrimming)
    checkFeatA = tclvalue(this$cbtfeatA)
    checkFeatB = tclvalue(this$cbtfeatB)
    checkFeatC = tclvalue(this$cbtfeatC)
    
    if (checkCALC=="density") {
      tkmessageBox(title = "For your information.",
                   message = "It is not meaningful to plot method \"density\".", icon = "error", type = "ok")
      stop("Chose another statistical method.")
    }
    
    if ( checkTRANS == "asinh" ) {
      scale = this$asinh$scale
      label = this$asinh$label
      grid.step = this$asinh$step
    } else {
      scale = this$biex$scale
      label = this$biex$label
      grid.step = this$biex$step
    } 
    quadrants.col = "black"
    
    FEATURES.fix = c(checkFeatA,checkFeatB,checkFeatC)
    FEATURES.idx = FEATURES.cutoff = c(0,0,0)
    
    fixed.features = 0
    for ( i in 1:length(FEATURES.fix)) {
      if (FEATURES.fix[i]=="1") {
        tmp.var = this$checkMarker(tclvalue(tkget(eval(parse(text=paste0("this$cbvar",i))))))
        if ( length(tmp.var)==0 ) {
          tkmessageBox(title = "An error has occured!", message = sprintf("Check your Feature pos %s.",i))
          stop("This feature is not existent in this file.")
        } else {
          FEATURES.fix[i] = tmp.var
          FEATURES.idx[i] = which(this$selected.vars==tmp.var)
          FEATURES.cutoff[i] = this$checkDigits(cutoff_id=as.numeric(checkFEATURES[[i]][2]))
          fixed.features = fixed.features +1
          if (this$working) printf("Fixed Feature #%s",i)
        }
      }
    }
    # stop("Funzt!")
    
    ### get column and cutoff number vector
    len=length(this$selected.vars)
    colvec=vector()
    cutoffs=vector()
    for (i in 1:len) {
      if ( tclvalue(this$cbVal[[i]]) == "1") {
        colvec=c(colvec,i);
        cutoffs=c(cutoffs,this$checkDigits(cutoff_id=i))
      }
    }
    ### remove doubled checked features
    remove.idx = which(colvec==FEATURES.cutoff)
    colvec = colvec[-remove.idx]
    cutoffs = cutoffs[-remove.idx]
    
    len.colvec = length(colvec)
    
    ### title length
    titlelen = length(unlist(strsplit(tclvalue(tkget(this$title,"1.0","end-1c")),"")))
    titleans = tclVar("yes")
    
    if ( titlelen <= 1 & this$working==FALSE) {
      titleans = tkmessageBox(title = "Are you sure?",
                              message = "Do not forget to set page title. Continue?", icon = "info", type = "yesno")
    } else if ( checkCALC == "freq" & !is.integer(which(cutoffs==0)) & this$working==FALSE ) {
      freqans = tkmessageBox(title = "Are you sure?",
                             message = "Do not forget to set cutoffs. Continue?", icon = "info", type = "yesno")
    }
    
    if ( (tclvalue(titleans) == "no") | tclvalue(freqans) == "no" ) stop("OK.")
    
    if (this$working) {
      this$saveinFolder=getwd()
    } else {
      # if path havent been saved yet
      if (!exists("lastpath",where=eval(parse(text="fcs")))){
        this$saveinFolder = tk_choose.dir(default = getwd(), caption = "Select directory to save the PDFs")
      } else {
        this$saveinFolder = tk_choose.dir(default = this$lastpath, caption = "Select directory to save the PDFs")
      }
    }
    
    if ( !is.na(this$saveinFolder) ) {
      tkconfigure(this$tt, cursor = "watch")
      
      this$lastpath = this$saveinFolder
      
      max.nhorizplots = 6
      max.nvertiplots = 8
      
      res.nhorizplots = (len.colvec-2)%%max.nhorizplots
      if (res.nhorizplots != 0) res.nhorizplots = max.nhorizplots - res.nhorizplots
      max.nplots = len.colvec*(len.colvec-1)*(len.colvec-2)
      max.npages = ceiling( (max.nplots+res.nhorizplots*(len.colvec-1)*len.colvec)/(max.nhorizplots*max.nvertiplots) ) # + 1
      
      ### Create Progress Bar
      this$pb <- tkProgressBar(title="triploT-Overview", label=paste("Creating page 1/",max.npages,sep=""),
                               min = 0, max = max.npages, initial = 1, width = 300)
      # set progress bar
      setTkProgressBar(this$pb,1,label=paste("Creating page 1/",max.npages,sep=""))
      
      timeSTART = Sys.time()
      if (!this$OverviewGate) cat("\n\n>>>> Start triploT-Overview with total data: \n\n")
      else cat("\n\n>>>> Start triploT overview with gated data: \n\n")
      if (this$working) print(sprintf("w: %s - time started",timeSTART))
      
      tmp.folder = file.path(this$saveinFolder,"tmp")
      if ( grepl("linux",sessionInfo()$R.version$os) ) {
        ### on linux
        system(paste("mkdir -p -v ",tmp.folder,sep=""))
      } else {
        ### on Windows
        dir.create(tmp.folder)
        printf(">>on Windows<< Creating folder: %s",tmp.folder)
        #system(paste("mkdir ",tmp.folder,sep=""))
        # or shell()?
      }
      
      file=tclvalue(tkget(this$tkchoosefile))
      if ( this$OverviewGate ) {
        displayfile = this$shortenFilename(this$plot.attr[[1]]$file.name,title=TRUE)
      } else {
        displayfile = this$shortenFilename(file)
      }
      
      table = this$total.projects[this$selected.projectnum]
      if ( grepl("temp",file) ) {
        file.idx = 1
        table=file
      } else {
        file.idx=this$current.filetable[which(this$current.filetable[,2]==file),1]
        this$selected.filenum = file.idx
      }
      
      if (this$working) printf("w: table=%s file.idx=%s",table,file.idx)
      
      binSize=as.numeric(tkget(this$binSize))
      mincount=as.numeric(tkget(this$minCountTri))
      xminval=as.numeric(tkget(this$minvalX))
      xmaxval=as.numeric(tkget(this$maxvalX))
      yminval=as.numeric(tkget(this$minvalY))
      ymaxval=as.numeric(tkget(this$maxvalY))
      
      #### Do not trim and remove doublets if data is gated
      if (this$OverviewGate) {
        # if data is gated, just recall
        data = this$getData(table,file.idx,columns=colvec)
      } else if ( checkTRIMMING == "1" ) {
        this$preprocData(mode="trim")
        data = this$data[,colvec]
      } else {
        this$getFile(table,file.idx)
        data = this$data[,colvec]
      }
      
      len = dim(data)[2]
      
      printf("Columns selected (%s): %s",dim(data)[2],paste(colvec,collapse=" "))
      printf("Features selected: %s",paste(colnames(data),collapse=" "))
      
      title=as.character(tclvalue(tkget(this$title,"1.0","end-1c")))
      
      ### if less than 16 markers were selected
      # do all triplots in one pdf file
      # otherwise make one pdf file for one page
      if (len.colvec<=16){
        
        toPDF(file=sprintf("%s_triploTOverview_%s%s_%s.pdf",displayfile,len.colvec,checkCALC,this$version), 
              path=this$saveinFolder, 
              title=sprintf("triploTOverview of %s",displayfile), 
              width=3.21*max.nhorizplots,
              height=3.5*max.nvertiplots,
              pointsize=11,
              {
                # plot histograms
                #this$plotHistograms(plotter="triover",pdf=TRUE)
                
                # set progress bar
                #setTkProgressBar(this$pb,2,label=paste("Creating page 2/",len.colvec+2,sep=""))
                # plot density plot
                #this$plotDensities(plotter="triover",pdf=TRUE)
                
                ##### start triploT overview
                ### set label and axes font sizes
                
                label.cex = 1.2
                set.cex.axes = 1
                set.mgp = c(1.9,0.5,0)
                par(mfrow=c(max.nvertiplots,max.nhorizplots),oma=c(0.5,1,6,1),mar=c(3,4,5,1))
                
                
                plot.idx = 0
                page.idx = 0
                current.page = 0
                for ( v1 in 1:len ) {
                  for ( v2 in 1:len) {
                    
                    if ( v1 == v2 ) next
                    
                    for ( v3 in 1:len) {
                      
                      ## skip loop if one of the axes are the same
                      if ( v1 == v3 | v2 == v3 ) next
                      
                      #select columns to plot
                      tdata=as.matrix(data[,c(v1,v2,v3)])
                      this$cnames = colnames(tdata)
                      
                      if ( cutoffs[v1] > 0 ) title.axis = sprintf("%s (%s)",colnames(data)[v1],cutoffs[v1])
                      else title.axis = colnames(data)[v1]
                      if ( cutoffs[v2] > 0 ) title.axis = c(title.axis,sprintf("%s (%s)",colnames(data)[v2],cutoffs[v2]))
                      else title.axis = c(title.axis,colnames(data)[v2])
                      # start plot
                      # mgp: A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. 
                      # The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
                      plot(1,type='n',frame.plot=FALSE,xlim=c(xminval,xmaxval+10*binSize),axes=FALSE
                           ,ylim=c(yminval-2.5*binSize,ymaxval+5*binSize),xlab=title.axis[1]
                           ,ylab=title.axis[2],cex.lab=label.cex,cex.axis=set.cex.axes,mgp=set.mgp)
                      box(lwd=0.5,col="darkgrey")
                      
                      plot.idx = plot.idx + 1
                      
                      ### draw axis on the bottom and on the left
                      axis(side=1, at=scale,labels=label,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
                      axis(side=2, at=scale,labels=label,las=3,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
                      
                      ### add grid
                      if(checkGRID=="1") {
                        xgrid.steps=seq(0,(xmaxval),by=grid.step)
                        ygrid.steps=seq(0,(ymaxval),by=grid.step)
                        abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
                      }
                      
                      
                      ncells = ncells.total = nrow(tdata)
                      # q1 Quadrant unten links
                      # q2 Quadrant unten rechts
                      # q3 Quadrant oben rechts
                      # q4 Quadrant oben links
                      if ( cutoffs[v1] > 0 & cutoffs[v2] > 0 ) {
                        
                        ### count cells in quadrant
                        tdata.q1 = tdata[which( tdata[,1]<cutoffs[v1] &  tdata[,2]<cutoffs[v2] ),3]
                        tdata.q2 = tdata[which( tdata[,1]>=cutoffs[v1] &  tdata[,2]<cutoffs[v2] ),3]
                        tdata.q3 = tdata[which( tdata[,1]>=cutoffs[v1] &  tdata[,2]>=cutoffs[v2] ),3]
                        tdata.q4 = tdata[which( tdata[,1]<cutoffs[v1] &  tdata[,2]>=cutoffs[v2] ),3]
                        
                        ### q[x].total [ink=black]
                        ### percentage of cells in quadrant to total cells 
                        ### or in MSI(+): percentage of cells in quadrant to total positive cells
                        this$q1.total = abs(100 * length( tdata.q1 ) / ncells)
                        this$q2.total = abs(100 * length( tdata.q2 ) / ncells)
                        this$q3.total = abs(100 * length( tdata.q3 ) / ncells)
                        this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
                        
                        if (cutoffs[v3] > 0 ) {
                          ### number of cells which are producing cells in feature C
                          ncells = nrow(tdata[which(tdata[,3]> cutoffs[v3]),])
                          
                          ### q[x].prodcells [ink=red]
                          ### percentage of cells which are positive for feature C in quadrant to total quadrant cells
                          this$q1.prodcells = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[v3])]) /length(tdata.q1)
                          if (is.nan(this$q1.prodcells)) this$q1.prodcells = 0
                          this$q2.prodcells = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[v3])]) /length(tdata.q2)
                          if (is.nan(this$q2.prodcells)) this$q2.prodcells = 0
                          this$q3.prodcells = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[v3])]) /length(tdata.q3)
                          if (is.nan(this$q3.prodcells)) this$q3.prodcells = 0
                          this$q4.prodcells = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[v3])]) /length(tdata.q4)
                          if (is.nan(this$q4.prodcells)) this$q4.prodcells = 0
                          
                          ### only do MSI plots on producing cells only
                          if ( checkCALC == "MSI(+)" ) {
                            ### cut all cells which are not producing cells
                            tdata.plus = tdata[which(tdata[,3]> cutoffs[3]),]
                            ncells = nrow(tdata.plus)
                            
                            tdata.q1 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
                            tdata.q2 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
                            tdata.q3 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
                            tdata.q4 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
                            
                            ### q[x].total [ink=blue]
                            ### in MSI(+): percentage of cells in quadrant to total positive cells
                            this$q1.total = abs(100 * length( which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
                            this$q2.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
                            this$q3.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ) ) / ncells)
                            this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
                          }
                          
                          ### q[x].prodcellsplus [ink=green]
                          ### percentage of cells which are positive for feature C to total cells
                          this$q1.prodcellsplus = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[v3])]) / ncells.total
                          if (is.nan(this$q1.prodcellsplus)) this$q1.prodcellsplus = 0
                          this$q2.prodcellsplus = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[v3])]) / ncells.total
                          if (is.nan(this$q2.prodcellsplus)) this$q2.prodcellsplus = 0
                          this$q3.prodcellsplus = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[v3])]) / ncells.total
                          if (is.nan(this$q3.prodcellsplus)) this$q3.prodcellsplus = 0
                          this$q4.prodcellsplus = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[v3])]) / ncells.total
                          if (is.nan(this$q4.prodcellsplus)) this$q4.prodcellsplus = 0
                        }
                      } 
                      
                      
                      tdata.zero = tdata[apply(tdata[, -1], MARGIN = 1, function(row) {all(row > 0)}), ]
                      ncells.zero=nrow(tdata.zero)
                      
                      ### calc quadrants with only positive values
                      # q1 Quadrant unten links
                      # q2 Quadrant unten rechts
                      # q3 Quadrant oben rechts
                      # q4 Quadrant oben links
                      if ( cutoffs[v1] > 0 & cutoffs[v2] > 0 ) {
                        q1.zero = abs(100 * length( which( tdata.zero[,1]<cutoffs[v1] &  tdata.zero[,2]<cutoffs[v2] ) ) / ncells.zero)      
                        q2.zero = abs(100 * length( which( tdata.zero[,1]>=cutoffs[v1] &  tdata.zero[,2]<cutoffs[v2] ) ) / ncells.zero)                    
                        q3.zero = abs(100 * length( which( tdata.zero[,1]>=cutoffs[v1] &  tdata.zero[,2]>=cutoffs[v2] ) ) / ncells.zero)                    
                        q4.zero = abs(100 - q1.zero - q2.zero - q3.zero)
                      }
                      
                      #this$bintriplot(tdata,cutoffs[c(v1,v2,v3)],set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp
                      #	,binSize=binSize,mincells=mincount, quadrants.color = quadrants.col)
                      if (checkCALC == "density") {
                        this$bintriplot(data=tdata,cutoffs=cutoffs[c(v1,v2,v3)],density=TRUE,binSize=binSize,mincells=mincount,overview=TRUE
                                        ,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
                      } else if (checkCALC == "MSI(+)") {
                        this$bintriplot(data=tdata.plus,cutoffs=cutoffs[c(v1,v2,v3)],binSize=binSize,mincells=mincount,,overview=TRUE
                                        ,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp,quadrants.color = "blue", data.origin=tdata)
                      } else {
                        this$bintriplot(data=tdata,cutoffs=cutoffs[c(v1,v2,v3)],binSize=binSize,mincells=mincount,overview=TRUE
                                        ,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
                      }
                      
                      
                      if (this$OverviewGate) {this$ncell.perc = this$ncell.sel/this$origin.ncells*100}
                      
                      # add title for single plot
                      if ( checkCALC == "freq" | grepl("MSI",checkCALC) ) {
                        if ( this$OverviewGate ) firstLine=sprintf("%s(%0.1f%%): %s,%s/cof=%s",displayfile,this$ncell.perc,checkCALC,v3,this$current.cofactor)
                        else firstLine=sprintf("%s: %s,cutoff=%s/cof=%s",displayfile,checkCALC,cutoffs[v3],this$current.cofactor)
                      } else {
                        if ( this$OverviewGate ) firstLine=sprintf("%s(%0.1f%%): %s/cof=%s",displayfile,this$ncell.perc,checkCALC,this$current.cofactor)
                        else firstLine=sprintf("%s: %s/cof=%s",displayfile,checkCALC,this$current.cofactor)
                      }
                      title(main=firstLine,line=3.2,cex.main=1.0*label.cex,adj=0)
                      
                      secondLine=sprintf("%s-%s(%0.1f%%); min/max=%s/%s; binSize=%s,#bins=%s",ncells,ncells.zero,(ncells.zero/ncells*100),mincount,this$maxcells,binSize,this$bincount)
                      title(main=secondLine,line=2.4,cex.main=0.7*label.cex,adj=0)
                      
                      thirdLine=""
                      if ( cutoffs[v1] > 0 & cutoffs[v2] > 0 ) {
                        thirdLine=sprintf("Q1=%0.1f/%0.1f Q2=%0.1f/%0.1f Q3=%0.1f/%0.1f Q4=%0.1f/%0.1f",q1.zero,this$q1.total,q2.zero,this$q2.total,q3.zero,this$q3.total,q4.zero,this$q4.total)
                      }
                      title(main=thirdLine,line=1.6,cex.main=0.7*label.cex,adj=0)
                      
                      if (checkGATED == "1" | grepl("temp",file)) {
                        if ( length(this$coords.info)>2 ) {
                          fourthLine=sprintf("%s; %s",this$coords.info[1],this$coords.info[2])
                          fifthLine=sprintf("%s; %s",this$coords.info[3],this$coords.info[4])
                        } else {
                          fourthLine=sprintf("%s",paste(this$coords.info,collapse="; "))
                          fifthLine=""
                        }
                        title(main=fourthLine,line=0.9,cex.main=0.6*label.cex,adj=0)
                        title(main=fifthLine,line=0.1,cex.main=0.6*label.cex,adj=0)
                      }
                      
                      ### update progress bar
                      if ( (plot.idx != 0) & ((plot.idx-1) %% (max.nhorizplots*max.nvertiplots) == 0) ) {
                        page.idx = page.idx + 1
                        info <- paste("Creating page ",page.idx,"/",max.npages,sep="")
                        if (this$working) printf("%s with plot# %s",info,plot.idx)
                        setTkProgressBar(this$pb, page.idx, label=info)
                      }
                      
                      #add title for page in 3D-Overview tab
                      if (current.page != page.idx) {
                        mtext(title, outer = TRUE, cex = 1.5,line=1.0,pos=2,xpd=TRUE)
                        if (length(this$coords.info>0)) mtext(sprintf("(%s)",paste(this$coords.info,collapse="; ")),outer=TRUE,cex = 0.8)
                        #printf("Print on Page %s",page.idx)
                        current.page = page.idx
                      }
                    }
                    
                    ### plot empty plots
                    if ( (plot.idx != 0) & (res.nhorizplots > 0) ) {
                      for ( i in 1:res.nhorizplots ) {
                        plot.new(); 
                        plot.idx = plot.idx+1
                      } 
                    }
                    
                  }
                }
                
                ### plot histograms
                this$plotHistograms(plotter="triover",pdf=TRUE)
                
                tkmessageBox(title = "Output of overview",
                             message = paste ("PDF created in folder ",this$saveinFolder, sep=""))
                
              }
        )
      } else {
        toPDF(file=sprintf("%s_triploTOverview_%shist_%s.pdf",displayfile,len.colvec,this$version), 
              path=this$saveinFolder, 
              title=sprintf("Histograms of %s",displayfile), 
              
              width=3.21*max.nhorizplots,
              height=3.5*max.nvertiplots,
              pointsize=11,
              {
                # plot histograms
                this$plotHistograms(plotter="triover",pdf=TRUE)
                
                # set progress bar
                #setTkProgressBar(this$pb,2,label=paste("Creating page 2/",max.npages,sep=""))
                
                # plot density plot
                #this$plotDensities(plotter="triover",pdf=TRUE)
              }
        )
        
        for ( v1 in 1:len ) {
          printf("New PDF: #%s",v1)
          plot.idx = 0
          page.idx = 0
          
          ### update progress bar
          #info <- paste("Creating page ",page.idx+1,"/",max.npages,sep="")
          #print(info)
          
          #setTkProgressBar(this$pb, v1+1, label=info)
          toPDF(file=sprintf("%s_triploTOverview_%s_%s%s_%s.pdf",displayfile,colnames(data)[v1],len.colvec-1,checkCALC,this$version), 
                #toPDF(file=sprintf("%s_triploTOverview_%s%s_%s_%s.pdf",displayfile,len.colvec,checkCALC,colnames(data)[v1],this$version), 
                path=this$saveinFolder, 
                title=sprintf("triploTOverview of %s(%s)",displayfile,colnames(data)[v1]), 
                
                width=3.21*max.nhorizplots,
                height=3.5*max.nvertiplots,
                pointsize=11,
                {
                  
                  ### new
                  label.cex = 1.2
                  set.cex.axes = 1
                  set.mgp = c(1.9,0.5,0)
                  par(mfrow=c(max.nvertiplots,max.nhorizplots),oma=c(0.5,1,5,1),mar=c(3,4,5,1))
                  ###
                  
                  for ( v2 in 1:len) {
                    if ( v1 == v2 ) next
                    for ( v3 in 1:len) {
                      
                      ## skip loop if one of the axes are the same
                      if ( v1 == v3 | v2 ==v3 ) next
                      plot.idx = plot.idx + 1
                      
                      ### update progress bar
                      #if (plot.idx %% (max.nhorizplots*max.nvertiplots) == 0 ) {
                      #	page.idx = page.idx + 1
                      #	info <- paste("Creating page ",page.idx,"/",max.npages,sep="")
                      #	if (this$working) printf("%s with plot# %s",info,plot.idx)
                      #	setTkProgressBar(this$pb, page.idx, label=info)
                      #}
                      
                      #else {
                      #select columns to plot
                      tdata=as.matrix(data[,c(v1,v2,v3)])
                      this$cnames = c(colnames(data)[v1],colnames(data)[v2],colnames(data)[v3])
                      
                      if ( cutoffs[v1] > 0 ) title.axis = sprintf("%s (%s)",colnames(data)[v1],cutoffs[v1])
                      else title.axis = colnames(data)[v1]
                      if ( cutoffs[v2] > 0 ) title.axis = c(title.axis,sprintf("%s (%s)",colnames(data)[v2],cutoffs[v2]))
                      else title.axis = c(title.axis,colnames(data)[v2])
                      
                      # start plot
                      # mgp: A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. 
                      # The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
                      plot(1,type='n',frame.plot=FALSE,xlim=c(xminval,xmaxval+10*binSize),axes=FALSE
                           ,ylim=c(yminval-2.5*binSize,ymaxval+5*binSize),xlab=title.axis[1]
                           ,ylab=title.axis[2],cex.lab=label.cex,cex.axis=set.cex.axes,mgp=set.mgp)
                      box(lwd=0.5,col="darkgrey")
                      
                      #add title for page in 3D-Overview tab
                      #mtext(title, outer = TRUE, cex = 1.5,line=1.3,pos=2,xpd=TRUE)
                      #if (length(this$coords.info>0)) mtext(sprintf("(%s)",paste(this$coords.info,collapse="; ")),outer=TRUE,cex = 0.8)
                      
                      ### draw axis on the bottom and on the left
                      axis(side=1, at=scale,labels=label,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
                      axis(side=2, at=scale,labels=label,las=3,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
                      
                      ### add grid
                      if(checkGRID=="1") {
                        xgrid.steps=seq(0,(xmaxval),by=grid.step)
                        ygrid.steps=seq(0,(ymaxval),by=grid.step)
                        abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
                      }
                      
                      
                      ### calc quadrants in total
                      ncells = ncells.total = nrow(tdata)
                      # q1 Quadrant unten links
                      # q2 Quadrant unten rechts
                      # q3 Quadrant oben rechts
                      # q4 Quadrant oben links
                      if ( cutoffs[v1] > 0 & cutoffs[v2] > 0 ) {
                        
                        ### count cells in quadrant
                        tdata.q1 = tdata[which( tdata[,1]<cutoffs[v1] &  tdata[,2]<cutoffs[v2] ),3]
                        tdata.q2 = tdata[which( tdata[,1]>=cutoffs[v1] &  tdata[,2]<cutoffs[v2] ),3]
                        tdata.q3 = tdata[which( tdata[,1]>=cutoffs[v1] &  tdata[,2]>=cutoffs[v2] ),3]
                        tdata.q4 = tdata[which( tdata[,1]<cutoffs[v1] &  tdata[,2]>=cutoffs[v2] ),3]
                        
                        ### q[x].total [ink=black]
                        ### percentage of cells in quadrant to total cells 
                        ### or in MSI(+): percentage of cells in quadrant to total positive cells
                        this$q1.total = abs(100 * length( tdata.q1 ) / ncells)
                        this$q2.total = abs(100 * length( tdata.q2 ) / ncells)
                        this$q3.total = abs(100 * length( tdata.q3 ) / ncells)
                        this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
                        
                        if (cutoffs[v3] > 0 ) {
                          ### number of cells which are producing cells in feature C
                          ncells = nrow(tdata[which(tdata[,3]> cutoffs[v3]),])
                          
                          ### q[x].prodcells [ink=red]
                          ### percentage of cells which are positive for feature C in quadrant to total quadrant cells
                          this$q1.prodcells = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[v3])]) /length(tdata.q1)
                          if (is.nan(this$q1.prodcells)) this$q1.prodcells = 0
                          this$q2.prodcells = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[v3])]) /length(tdata.q2)
                          if (is.nan(this$q2.prodcells)) this$q2.prodcells = 0
                          this$q3.prodcells = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[v3])]) /length(tdata.q3)
                          if (is.nan(this$q3.prodcells)) this$q3.prodcells = 0
                          this$q4.prodcells = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[v3])]) /length(tdata.q4)
                          if (is.nan(this$q4.prodcells)) this$q4.prodcells = 0
                          
                          ### only do MSI plots on producing cells only
                          if ( checkCALC == "MSI(+)" ) {
                            ### cut all cells which are not producing cells
                            tdata.plus = tdata[which(tdata[,3]> cutoffs[3]),]
                            ncells = nrow(tdata.plus)
                            
                            tdata.q1 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
                            tdata.q2 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
                            tdata.q3 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
                            tdata.q4 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
                            
                            ### q[x].total [ink=blue]
                            ### in MSI(+): percentage of cells in quadrant to total positive cells
                            this$q1.total = abs(100 * length( which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
                            this$q2.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
                            this$q3.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ) ) / ncells)
                            this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
                          }
                          
                          ### q[x].prodcellsplus [ink=green]
                          ### percentage of cells which are positive for feature C to total cells
                          this$q1.prodcellsplus = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[v3])]) / ncells.total
                          if (is.nan(this$q1.prodcellsplus)) this$q1.prodcellsplus = 0
                          this$q2.prodcellsplus = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[v3])]) / ncells.total
                          if (is.nan(this$q2.prodcellsplus)) this$q2.prodcellsplus = 0
                          this$q3.prodcellsplus = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[v3])]) / ncells.total
                          if (is.nan(this$q3.prodcellsplus)) this$q3.prodcellsplus = 0
                          this$q4.prodcellsplus = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[v3])]) / ncells.total
                          if (is.nan(this$q4.prodcellsplus)) this$q4.prodcellsplus = 0
                        }
                      } 
                      
                      tdata.zero = tdata[apply(tdata[, -1], MARGIN = 1, function(row) {all(row > 0)}), ]
                      ncells.zero=nrow(tdata.zero)
                      
                      ### calc quadrants with only positive values
                      # q1 Quadrant unten links
                      # q2 Quadrant unten rechts
                      # q3 Quadrant oben rechts
                      # q4 Quadrant oben links
                      if ( cutoffs[v1] > 0 & cutoffs[v2] > 0 ) {
                        q1.zero = abs(100 * length( which( tdata.zero[,1]<cutoffs[v1] &  tdata.zero[,2]<cutoffs[v2] ) ) / ncells.zero)      
                        q2.zero = abs(100 * length( which( tdata.zero[,1]>=cutoffs[v1] &  tdata.zero[,2]<cutoffs[v2] ) ) / ncells.zero)                    
                        q3.zero = abs(100 * length( which( tdata.zero[,1]>=cutoffs[v1] &  tdata.zero[,2]>=cutoffs[v2] ) ) / ncells.zero)                    
                        q4.zero = abs(100 - q1.zero - q2.zero - q3.zero)
                      }
                      
                      #this$bintriplot(tdata,cutoffs[c(v1,v2,v3)],set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp
                      #	,binSize=binSize,mincells=mincount, quadrants.color = quadrants.col)
                      if (checkCALC == "density") {
                        this$bintriplot(data=tdata,cutoffs=cutoffs[c(v1,v2,v3)],density=TRUE,binSize=binSize,mincells=mincount,overview=TRUE
                                        ,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
                      } else if (checkCALC == "MSI(+)") {
                        this$bintriplot(data=tdata.plus,cutoffs=cutoffs[c(v1,v2,v3)],binSize=binSize,mincells=mincount,,overview=TRUE
                                        ,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp,quadrants.color = "blue", data.origin=tdata)
                      } else {
                        this$bintriplot(data=tdata,cutoffs=cutoffs[c(v1,v2,v3)],binSize=binSize,mincells=mincount,overview=TRUE
                                        ,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
                      }
                      
                      
                      if (this$OverviewGate) {this$ncell.perc = this$ncell.sel/this$origin.ncells*100}
                      
                      # add title for single plot
                      if ( checkCALC == "freq" | grepl("MSI",checkCALC) ) {
                        if ( this$OverviewGate ) firstLine=sprintf("%s(%0.1f%%): %s,%s/cof=%s",
                                                                   displayfile,this$ncell.perc,checkCALC,v3,this$current.cofactor)
                        else firstLine=sprintf("%s: %s,cutoff=%s/cof=%s",
                                               displayfile,checkCALC,cutoffs[v3],this$current.cofactor)
                      } else {
                        if ( this$OverviewGate ) firstLine=sprintf("%s(%0.1f%%): %s/cof=%s",
                                                                   displayfile,this$ncell.perc,checkCALC,this$current.cofactor)
                        else firstLine=sprintf("%s: %s/cof=%s",
                                               displayfile,checkCALC,this$current.cofactor)
                      }
                      title(main=firstLine,line=3.2,cex.main=1.0*label.cex,adj=0)
                      
                      secondLine=sprintf("%s-%s(%0.1f%%); min/max=%s/%s; binSize=%s,#bins=%s",ncells,ncells.zero,(ncells.zero/ncells*100),mincount,this$maxcells,binSize,this$bincount)
                      title(main=secondLine,line=2.4,cex.main=0.7*label.cex,adj=0)
                      
                      thirdLine=""
                      if ( cutoffs[v1] > 0 & cutoffs[v2] > 0 ) {
                        thirdLine=sprintf("Q1=%0.1f/%0.1f Q2=%0.1f/%0.1f Q3=%0.1f/%0.1f Q4=%0.1f/%0.1f"
                                          ,q1.zero,this$q1.total,q2.zero,this$q2.total,q3.zero,this$q3.total,q4.zero,this$q4.total)
                      }
                      title(main=thirdLine,line=1.6,cex.main=0.7*label.cex,adj=0)
                      
                      if (checkGATED == "1" | grepl("temp",file)) {
                        if ( length(this$coords.info)>2 ) {
                          fourthLine=sprintf("%s; %s",this$coords.info[1],this$coords.info[2])
                          fifthLine=sprintf("%s; %s",this$coords.info[3],this$coords.info[4])
                        } else {
                          fourthLine=sprintf("%s",paste(this$coords.info,collapse="; "))
                          fifthLine=""
                        }
                        title(main=fourthLine,line=0.9,cex.main=0.6*label.cex,adj=0)
                        title(main=fifthLine,line=0.1,cex.main=0.6*label.cex,adj=0)
                      }
                    }
                    ### update progress bar
                    if ( (plot.idx != 0) & ((plot.idx-1) %% (max.nhorizplots*max.nvertiplots) == 0) ) {
                      page.idx = page.idx + 1
                      info <- paste("Creating page ",page.idx,"/",max.npages,sep="")
                      if (this$working) printf("%s with plot# %s",info,plot.idx)
                      setTkProgressBar(this$pb, page.idx, label=info)
                    }
                    
                    #add title for page in 3D-Overview tab
                    if (current.page != page.idx) {
                      mtext(title, outer = TRUE, cex = 1.5,line=1.0,pos=2,xpd=TRUE)
                      if (length(this$coords.info>0)) mtext(sprintf("(%s)",paste(this$coords.info,collapse="; ")),outer=TRUE,cex = 0.8)
                      #printf("Print on Page %s",page.idx)
                      current.page = page.idx
                    }
                    
                    ### plot emptry plots
                    printf("last plot #:%s with X:%s Y:%s; plotting %s more empty plots",plot.idx,v1,v2,res.nhorizplots)
                    if (res.nhorizplots > 0) {
                      for ( i in 1:res.nhorizplots ) {
                        plot(1,type='n',frame.plot=FALSE,axes=FALSE,ylab="")
                        plot.idx = plot.idx+1
                      }
                    }
                  }
                }
          )
        }
        
        tkmessageBox(title = "Output of overview",
                     message = paste ("Plots created in folder ",this$saveinFolder, sep=""))
      }
      
      tkconfigure(this$tt, cursor = "left_ptr")
      
      cat("\n")
      # close progress bar
      close(this$pb)
      
      #if (checkPNG == "1") printf("PNGs are saved in %s.",file.path(this$saveinFolder,"png"))
      
      print(Sys.time()-timeSTART)
      
      cat("\n>>>> Done.\n")
    }
    # }
  # }	
  this$OverviewGate=FALSE
}

fcs$dotriploTOverviewX <- function(table=NA) {
	### triploTs Overviews only with Feature A on x axis
	this=fcs
	answer = tclVar("yes")
	freqans = tclVar("no")

	### if axes ranges are not the same
	if ( !this$checkAxesRange() ) {
		tkmessageBox(title = "An error has occured!",
			message = "Please set x and y axis with the same range.", icon = "error", type = "ok")
		stop("Set x and y axis with the same range.")
	}

	if (!this$working) answer=tkmessageBox(title = "Starting..", message = "Are you sure?", icon = "info", type = "yesno")

	if (tclvalue(answer) == "yes") {
		var1 = this$checkMarker(tclvalue(tkget(this$cbvar1)))
		### if features are not in sample
		if ( length(var1)==0  ){
			tkmessageBox(title = "An error has occured!",
				message = "Check your Feature A.")
			stop("One of the features are not existent.")
		}

		var1.idx = which(this$selected.vars==var1)
		cutoff.var1 = this$checkDigits(cutoff_id=var1.idx)

		checkCALC = tclvalue(this$rbcalc)
		checkTRANS = tclvalue(this$rbtrans)
		checkGATED = tclvalue(this$cbtgateData)
		checkGRID = tclvalue(this$cbtshowGrid)
		checkTRIMMING = tclvalue(this$cbttrimming)

		if ( checkTRANS == "asinh" ) {
		  scale = this$asinh$scale
		  label = this$asinh$label
		  grid.step = this$asinh$step
		} else {
		  scale = this$biex$scale
		  label = this$biex$label
		  grid.step = this$biex$step
		} 

		quadrants.col = "black"

		len=length(this$selected.vars)

		# column and cutoff number vector
		colvec=vector()
		cutoffs=vector()
		for (i in 1:len) {
			if ( tclvalue(this$cbVal[[i]]) == "1") {
				colvec=c(colvec,i);
				cutoffs=c(cutoffs,this$checkDigits(cutoff_id=i))
			}
		}
		### if var1 is already checked in cutoff panel
		if ( any(var1.idx==colvec) ) {
			colvec.var1.idx = which(var1.idx==colvec)
			colvec = colvec[-colvec.var1.idx]
			cutoffs = cutoffs[-colvec.var1.idx]
		}
		colvec = c(var1.idx,colvec)
		cutoffs = c(cutoff.var1,cutoffs)
		len.colvec = length(colvec)
		# title length
		titlelen = length(unlist(strsplit(tclvalue(tkget(this$title,"1.0","end-1c")),"")))
		titleans = tclVar("yes")

		if ( len.colvec < 3 ) {
			tkmessageBox(title = "An error has occured!",
				message = "Please select at least three markers.", icon = "error", type = "ok")
			stop("Please select at least three markers.")
		}
		if ( titlelen <= 1 & this$working==FALSE) {
			titleans = tkmessageBox(title = "Are you sure?",
				message = "Do not forget to set page title. Continue?", icon = "info", type = "yesno")
		} else if ( checkCALC == "freq" & !is.integer(which(cutoffs==0)) & this$working==FALSE ) {
			freqans = tkmessageBox(title = "Are you sure?",
				message = "Do not forget to set cutoffs. Continue?", icon = "info", type = "yesno")
		} else if (checkCALC=="density") {
		  tkmessageBox(title = "For your information.",
		               message = "It is not meaningful to plot method \"density\".", icon = "error", type = "ok")
		  stop("Chose another statistical method.")
		}

		if ( (tclvalue(titleans) == "yes") | tclvalue(freqans) == "yes" ) {
			if (this$working) this$saveinFolder=getwd()
			else this$saveinFolder = tk_choose.dir(default = getwd(), caption = "Select directory to save the PDFs")
		

		if ( !is.na(this$saveinFolder) ) {
			tkconfigure(this$tt, cursor = "watch")

			max.nhorizplots = 6
			max.nvertiplots = 8
			res.nhorizplots = (len.colvec-2)%%max.nhorizplots
			if (res.nhorizplots != 0) res.nhorizplots = max.nhorizplots - res.nhorizplots
			max.nplots = (len.colvec-2)*(len.colvec-3)
			max.npages = ceiling( (max.nplots+res.nhorizplots*(len.colvec-1))/(max.nhorizplots*max.nvertiplots) )

			### Create Progress Bar
			this$pb <- tkProgressBar(title="triploT overview", label=paste("Creating page 1/",max.npages,sep=""),
				min = 0, max = max.npages, initial = 1, width = 300)
			# Set progress bar
			setTkProgressBar(this$pb,1,label=paste("Creating page 1/",max.npages,sep=""))
			
			timeSTART = Sys.time()
			if (!this$OverviewGate) cat("\n\n>>>> Start triploTOverview with total data: \n\n")
			else cat("\n\n>>>> Start triploTOverview with gated data: \n\n")
			if (this$working) printf("w: %s - time started",timeSTART)

			tmp.folder = file.path(this$saveinFolder,"tmp")
			if ( grepl("linux",sessionInfo()$R.version$os) ) {
				### on linux
				system(paste("mkdir -p -v ",tmp.folder,sep=""))
			} else {
				### on Windows
				dir.create(tmp.folder)
				printf("Creating folder: %s",tmp.folder)
				#system(paste("mkdir ",tmp.folder,sep=""))
				# or shell()?
			}

			file=tclvalue(tkget(this$tkchoosefile))
			if ( this$OverviewGate ) {
			displayfile = this$shortenFilename(this$plot.attr[[1]]$file.name,title=TRUE)
			} else {
					displayfile = this$shortenFilename(file,title=TRUE)
			}

			table = this$total.projects[this$selected.projectnum]
			if ( grepl("temp",file) ) {
				file.idx = 1
				table=file
			} else {
				file.idx=this$current.filetable[which(this$current.filetable[,2]==file),1]
				this$selected.filenum = file.idx
			}

			if (this$working) printf("w: table=%s file.idx=%s",table,file.idx)

			binSize=as.numeric(tkget(this$binSize))
			mincount=as.numeric(tkget(this$minCountTri))
			xminval=as.numeric(tkget(this$minvalX))
			xmaxval=as.numeric(tkget(this$maxvalX))
			yminval=as.numeric(tkget(this$minvalY))
			ymaxval=as.numeric(tkget(this$maxvalY))


			#### Do not trim and remove doublets if data is gated
			if (this$OverviewGate) {
				# if data is gated, just recall
				data=this$getData(table,file.idx,columns=colvec)
			} else if ( checkTRIMMING == "1" ) {
				this$preprocData(mode="trim")
				data=this$data[,colvec]
			} else {
				this$getFile(table,file.idx)
				data=this$data[,colvec]
			}

			len = dim(data)[2]

			printf("Fixed Feature A: %s",var1)
			printf("Columns selected (%s): %s",dim(data)[2],paste(colvec,collapse=" "))
			printf("Column names: %s",paste(colnames(data),collapse=" "))

			title=as.character(tclvalue(tkget(this$title,"1.0","end-1c")))

			
			toPDF(file=sprintf("%s_triploTOverviewX_%s_%s%s_%s.pdf",displayfile,var1,len.colvec-1,checkCALC,this$version), 
				path=tmp.folder, 
				title=sprintf("triploTOverviewX of %s",displayfile), 
				### 
				width=3.21*max.nhorizplots,
				height=3.5*max.nvertiplots,
				pointsize=11,
				###
				{
				# set progress bar
				#setTkProgressBar(this$pb,2,label=paste("Creating page 2/",max.npages+2,sep=""))
				# plot density plot
				#this$plotDensities(plotter="triover",pdf=TRUE)

				##### start triploT overview

				### label and axes font sizes
				label.cex = 1.2
				set.cex.axes = 1
				set.mgp = c(1.9,0.5,0)
				par(mfrow=c(max.nvertiplots,max.nhorizplots),oma=c(0.5,1,5,1),mar=c(3,4,5,1))
				###

				#add title for page in 3D-Overview tab
				#mtext(title, outer = TRUE, cex = 1.5,line=1.3,pos=2,xpd=TRUE)
				if (length(this$coords.info>0)) mtext(sprintf("(%s)",paste(this$coords.info,collapse="; ")),outer=TRUE,cex = 0.8)
				
				plot.idx = 0
				page.idx = 0
				current.page = 0
				for ( v2 in 2:len ){
					for ( v3 in 2:len ){

						## skip loop if one of the axes are the same
						if ( v2 ==v3 ) next

						### update progress bar
						if (plot.idx %% (max.nhorizplots*max.nvertiplots) == 0 ) {
							page.idx = page.idx + 1
							info <- paste("Creating page ",page.idx,"/",max.npages,sep="")
							print(info)
							printf("plot #:%s",plot.idx)
							setTkProgressBar(this$pb, page.idx, label=info)
						}

						#select columns to plot
						tdata=as.matrix(data[c(1,v2,v3)])
						this$cnames = colnames(tdata)

						if ( cutoff.var1 > 0 ) title.axis = sprintf("%s(%s)",var1,cutoff.var1)
						else title.axis = var1
						if ( cutoffs[2] > 0 ) title.axis = c(title.axis,sprintf("%s(%s)",colnames(data)[v2],cutoffs[v2]))
						else title.axis = c(title.axis,colnames(data)[v2])

						# start plot
						# mgp: A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. 
						# The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
						plot(1,type='n',frame.plot=FALSE,xlim=c(xminval,xmaxval+10*binSize),axes=FALSE
							,ylim=c(yminval-2.5*binSize,ymaxval+5*binSize),xlab=title.axis[1]
							,ylab=title.axis[2],cex.lab=label.cex,cex.axis=set.cex.axes,mgp=set.mgp)
						box(lwd=0.5,col="darkgrey")

						plot.idx = plot.idx + 1

						### draw axis on the bottom and on the left
						axis(side=1, at=scale,labels=label,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
						axis(side=2, at=scale,labels=label,las=3,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")

						### add grid
						if(checkGRID=="1") {
						  xgrid.steps=seq(0,(xmaxval),by=grid.step)
						  ygrid.steps=seq(0,(ymaxval),by=grid.step)
						  abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
						}


						### calc quadrants in total
						ncells = ncells.total = nrow(tdata)
						# q1 Quadrant unten links
						# q2 Quadrant unten rechts
						# q3 Quadrant oben rechts
						# q4 Quadrant oben links
						if ( cutoff.var1 > 0 & cutoffs[v2] > 0 ) {

							### count cells in quadrant
							tdata.q1 = tdata[which( tdata[,1]<cutoff.var1 &  tdata[,2]<cutoffs[v2] ),3]
							tdata.q2 = tdata[which( tdata[,1]>=cutoff.var1 &  tdata[,2]<cutoffs[v2] ),3]
							tdata.q3 = tdata[which( tdata[,1]>=cutoff.var1 &  tdata[,2]>=cutoffs[v2] ),3]
							tdata.q4 = tdata[which( tdata[,1]<cutoff.var1 &  tdata[,2]>=cutoffs[v2] ),3]

							### q[x].total [ink=black]
							### percentage of cells in quadrant to total cells 
							### or in MSI(+): percentage of cells in quadrant to total positive cells
							this$q1.total = abs(100 * length( tdata.q1 ) / ncells)
							this$q2.total = abs(100 * length( tdata.q2 ) / ncells)
							this$q3.total = abs(100 * length( tdata.q3 ) / ncells)
							this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
						
							if (cutoffs[v3] > 0 ) {
								### number of cells which are producing cells in feature C
								ncells = nrow(tdata[which(tdata[,3]> cutoffs[v3]),])

								### q[x].prodcells [ink=red]
								### percentage of cells which are positive for feature C in quadrant to total quadrant cells
								this$q1.prodcells = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[v3])]) /length(tdata.q1)
								if (is.nan(this$q1.prodcells)) this$q1.prodcells = 0
								this$q2.prodcells = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[v3])]) /length(tdata.q2)
								if (is.nan(this$q2.prodcells)) this$q2.prodcells = 0
								this$q3.prodcells = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[v3])]) /length(tdata.q3)
								if (is.nan(this$q3.prodcells)) this$q3.prodcells = 0
								this$q4.prodcells = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[v3])]) /length(tdata.q4)
								if (is.nan(this$q4.prodcells)) this$q4.prodcells = 0

								### only do MSI plots on producing cells only
								if ( checkCALC == "MSI(+)" ) {
								  ### cut all cells which are not producing cells
								  tdata.plus = tdata[which(tdata[,3]> cutoffs[3]),]
								  ncells = nrow(tdata.plus)
								  
								  tdata.q1 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
								  tdata.q2 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
								  tdata.q3 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
								  tdata.q4 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
								  
								  ### q[x].total [ink=blue]
								  ### in MSI(+): percentage of cells in quadrant to total positive cells
								  this$q1.total = abs(100 * length( which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
								  this$q2.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
								  this$q3.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ) ) / ncells)
								  this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
								}

								### q[x].prodcellsplus [ink=green]
								### percentage of cells which are positive for feature C to total cells
								this$q1.prodcellsplus = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[v3])]) / ncells.total
								if (is.nan(this$q1.prodcellsplus)) this$q1.prodcellsplus = 0
								this$q2.prodcellsplus = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[v3])]) / ncells.total
								if (is.nan(this$q2.prodcellsplus)) this$q2.prodcellsplus = 0
								this$q3.prodcellsplus = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[v3])]) / ncells.total
								if (is.nan(this$q3.prodcellsplus)) this$q3.prodcellsplus = 0
								this$q4.prodcellsplus = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[v3])]) / ncells.total
								if (is.nan(this$q4.prodcellsplus)) this$q4.prodcellsplus = 0
							}
						}

						tdata.zero = tdata[apply(tdata[, -1], MARGIN = 1, function(row) {all(row > 0)}), ]
						ncells.zero=nrow(tdata.zero)

						### calc quadrants with only positive values
						# q1 Quadrant unten links
						# q2 Quadrant unten rechts
						# q3 Quadrant oben rechts
						# q4 Quadrant oben links
						if ( cutoff.var1 > 0 & cutoffs[v2] > 0 ) {
							q1.zero = abs(100 * length( which( tdata.zero[,1]<cutoff.var1 &  tdata.zero[,2]<cutoffs[v2] ) ) / ncells.zero)      
							q2.zero = abs(100 * length( which( tdata.zero[,1]>=cutoff.var1 &  tdata.zero[,2]<cutoffs[v2] ) ) / ncells.zero)                    
							q3.zero = abs(100 * length( which( tdata.zero[,1]>=cutoff.var1 &  tdata.zero[,2]>=cutoffs[v2] ) ) / ncells.zero)                    
							q4.zero = abs(100 - q1.zero - q2.zero - q3.zero)
						}

						#this$bintriplot(tdata,c(cutoff.var1,cutoffs[c(v2,v3)]),set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp
						#	,binSize=binSize,mincells=mincount, quadrants.color = quadrants.col)
						if (checkCALC == "density") {
							this$bintriplot(data=tdata,cutoffs=c(cutoff.var1,cutoffs[c(v2,v3)]),density=TRUE,binSize=binSize,mincells=mincount,overview=TRUE
								,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
						} else if (checkCALC == "MSI(+)") {
							this$bintriplot(data=tdata.plus,cutoffs=c(cutoff.var1,cutoffs[c(v2,v3)]),binSize=binSize,mincells=mincount,overview=TRUE
								,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp,quadrants.color = "blue", data.origin=tdata)
						} else {
							this$bintriplot(data=tdata,cutoffs=c(cutoff.var1,cutoffs[c(v2,v3)]),binSize=binSize,mincells=mincount,overview=TRUE
								,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
						}

						if (this$OverviewGate) {this$ncell.perc = this$ncell.sel/this$origin.ncells*100}

						# add title for single plot
						if ( checkCALC == "freq" | grepl("MSI",checkCALC) ) {
							if ( this$OverviewGate ) firstLine=sprintf("%s(%0.1f%%): %s(%s)/cof=%s",
							                                           displayfile,this$ncell.perc,checkCALC,this$cnames[3],this$current.cofactor)
							else firstLine=sprintf("%s: %s(%s)/cof=%s",
							                       displayfile,checkCALC,this$cnames[3],this$current.cofactor)
						} else {
							if ( this$OverviewGate ) firstLine=sprintf("%s(%0.1f%%): %s/cof=%s",
							                                           displayfile,this$ncell.perc,checkCALC,this$current.cofactor)
							else firstLine=sprintf("%s: %s/cof=%s",
							                       displayfile,checkCALC,this$current.cofactor)
						}
						title(main=firstLine,line=3.2,cex.main=1.0*label.cex,adj=0)

						secondLine=sprintf("%s-%s(%0.1f%%); min/max=%s/%s; binSize=%s,#bins=%s",ncells,ncells.zero,(ncells.zero/ncells*100),mincount,this$maxcells,binSize,this$bincount)
						title(main=secondLine,line=2.4,cex.main=0.7*label.cex,adj=0)

						thirdLine=""
						if ( cutoff.var1 > 0 & cutoffs[v2] > 0 ) {
							thirdLine=sprintf("Q1=%0.1f/%0.1f Q2=%0.1f/%0.1f Q3=%0.1f/%0.1f Q4=%0.1f/%0.1f",q1.zero,this$q1.total,q2.zero,this$q2.total,q3.zero,this$q3.total,q4.zero,this$q4.total)
						}
						title(main=thirdLine,line=1.6,cex.main=0.7*label.cex,adj=0)

						if (checkGATED == "1" | grepl("temp",file)) {
							if ( length(this$coords.info)>2 ) {
								fourthLine=sprintf("%s; %s",this$coords.info[1],this$coords.info[2])
								fifthLine=sprintf("%s; %s",this$coords.info[3],this$coords.info[4])
							} else {
								fourthLine=sprintf("%s",paste(this$coords.info,collapse="; "))
								fifthLine=""
							}
							title(main=fourthLine,line=0.9,cex.main=0.6*label.cex,adj=0)
							title(main=fifthLine,line=0.1,cex.main=0.6*label.cex,adj=0)
						}
					}
					#add title for page in 3D-Overview tab
					if (current.page != page.idx) {
						mtext(title, outer = TRUE, cex = 1.5,line=1.0,pos=2,xpd=TRUE)
						if (length(this$coords.info>0)) mtext(sprintf("(%s)",paste(this$coords.info,collapse="; ")),outer=TRUE,cex = 0.8)
						#printf("Print on Page %s",page.idx)
						current.page = page.idx
					}
					if (res.nhorizplots > 0) {
						for ( i in 1:res.nhorizplots ) {
							plot(1,type='n',frame.plot=FALSE,axes=FALSE,ylab=""); 
							plot.idx = plot.idx+1;
						}
					}
				}

				### plot histograms
				this$plotHistograms(plotter="triover",pdf=TRUE)
			})
			tkconfigure(this$tt, cursor = "left_ptr")

			cat("\n")
			# close progress bar
			close(this$pb)
			if (!this$working) tkmessageBox(title = "Output of overview",
				message = paste ("Plots created in folder ",tmp.folder, sep=""))

			print(Sys.time()-timeSTART)
			cat("\n>>>> Done.\n")
		}
		}
	}	
	this$OverviewGate=FALSE
}

fcs$dotriploTOverviewXY <- function(table=NA) {
  ### triploTs Overviews only with Feature A on x axis and feature Y on y axis
  this=fcs
  answer = tclVar("yes")
  
  ### if axes ranges are not the same
  if ( !this$checkAxesRange() ) {
    tkmessageBox(title = "An error has occured!",
                 message = "Please set x and y axis with the same range.", icon = "error", type = "ok")
    stop("Set x and y axis with the same range.")
  }
  
  ### if Feature A or Y is not in sample
  var1 = this$checkMarker(tclvalue(tkget(this$cbvar1)))
  var2 = this$checkMarker(tclvalue(tkget(this$cbvar2)))
  if ( length(var1)==0 | length(var2)==0 ){
    tkmessageBox(title = "An error has occured!",
                 message = "Check your Feature A or Y.")
    stop("One of the features are not existent.")
  }
  
  var1.idx = which(this$selected.vars==var1)
  cutoff.var1 = this$checkDigits(cutoff_id=var1.idx)
  var2.idx = which(this$selected.vars==var2)
  cutoff.var2 = this$checkDigits(cutoff_id=var2.idx)
  
  file = tclvalue(tkget(this$tkchoosefile))
  displayfile = this$shortenFilename(file)
  table = this$total.projects[this$selected.projectnum]
  if ( grepl("temp",file) ) {
    file.idx = 1
    table=file
  } else {
    file.idx=this$current.filetable[which(this$current.filetable[,2]==file),1]
    this$selected.filenum = file.idx
  }
  if (this$working) printf("w: table=%s file.idx=%s",table,file.idx)
  
  binSize = as.numeric(tkget(this$binSize))
  mincount = as.numeric(tkget(this$minCountTri))
  xminval = as.numeric(tkget(this$minvalX))
  xmaxval = as.numeric(tkget(this$maxvalX))
  yminval = as.numeric(tkget(this$minvalY))
  ymaxval = as.numeric(tkget(this$maxvalY))
  
  checkCALC = tclvalue(this$rbcalc)
  checkTRANS = tclvalue(this$rbtrans)
  checkGATED = tclvalue(this$cbtgateData)
  checkGRID = tclvalue(this$cbtshowGrid)
  checkTRIMMING = tclvalue(this$cbttrimming)
  
  if ( checkTRANS == "asinh" ) {
	  scale = this$asinh$scale
	  label = this$asinh$label
	  grid.step = this$asinh$step
  } else {
	  scale = this$biex$scale
	  label = this$biex$label
	  grid.step = this$biex$step
  } 

  quadrants.col = "black"
  
  len=length(this$selected.vars)
  
  ### select checked features
  # column and cutoff vector
  colvec=vector()
  cutoffs=vector()
  for (i in 1:len) {
    if ( tclvalue(this$cbVal[[i]]) == "1") {
      colvec=c(colvec,i);
      cutoffs=c(cutoffs,this$checkDigits(cutoff_id=i))
    }
  }
  
  ### if method is freq or MSI(+), cutoff(z) needs to be setted
  if ( any(cutoffs==0) & (checkCALC=="freq" | checkCALC=="MSI(+)") ) {
    tkmessageBox(title = "An error has occured in method: MSI(+) or freq!",
                 message = "You forgot to set cutoff for some features C.", icon = "error", type = "ok")
    stop("Missing production cutoff for features C.")
  } else if (checkCALC=="density") {
    tkmessageBox(title = "For your information.",
                 message = "It is not meaningful to plot method \"density\".", icon = "error", type = "ok")
    stop("Chose another statistical method.")
  }

    ### if var1 is already checked in cutoff panel, then remove
  if ( any(var1.idx==colvec) ) {
    colvec.var1.idx = which(var1.idx==colvec)
    colvec = colvec[-colvec.var1.idx]
    cutoffs = cutoffs[-colvec.var1.idx]
  }
  ### if var2 is already checked in cutoff panel, then remove
  if ( any(var2.idx==colvec) ) {
    colvec.var2.idx = which(var2.idx==colvec)
    colvec = colvec[-colvec.var2.idx]
    cutoffs = cutoffs[-colvec.var2.idx]
  }
  colvec = c(var1.idx,var2.idx,colvec)
  cutoffs = c(cutoff.var1,cutoff.var2,cutoffs)
  len.colvec = length(colvec)
  
  # title length
  titlelen = length(unlist(strsplit(tclvalue(tkget(this$title,"1.0","end-1c")),"")))
  titleans = tclVar("yes")
  
  ### if there is no feature C to plot
  if ( len.colvec < 3 ) {
    tkmessageBox(title = "An error has occured!",
                 message = "Please check at least one feature to the left.", icon = "error", type = "ok")
    stop("Select one or more features to the left.")
  }
  
  #### new getFile if current file is not recent one
  timeSTART = Sys.time()
  
  if ( is.null(this$data) | this$current.project != table | this$current.filenum != file.idx | 
       this$current.trans != tclvalue(this$rbtrans) | this$current.cofactor != as.numeric(tclvalue(this$rbasinh))) {
    if ( this$working ) print("Time loading data:")
    this$getFile(table,file.idx)
    if ( checkTRIMMING == "1" ) this$preprocData(mode="trim")
    if ( this$working ) print(Sys.time()-timeSTART)
  }
  data=this$data[,colvec]
  
  ### if percentage is checked
  # calculate cutoffs and set check button to zero
  for ( i in 1:length(cutoffs)) {
    if ( tclvalue(this$cbcutoffperc[[colvec[i]]]) == "1" ) {
      cutoffs[colvec[i]] = this$calcCUTOFF(this$data[colvec[i]],this$checkDigits(cutoff_id=colvec[i]),colnames(this$data[colvec[i]]),cutoffs[colvec[i]])
    }
  }
  ####
  
  printf("Fixed Feature A: %s",var1)
  printf("Fixed feature Y: %s",var2)
  printf("Columns selected (%s): %s",dim(data)[2],paste(colvec,collapse=" "))
  printf("Column names: %s",paste(colnames(data),collapse=" "))
  
  
  ### plot triploTs in plot window if there are not many features to plot
  plot.ncol=as.numeric(tclvalue(this$vncol))
  plot.nrow=as.numeric(tclvalue(this$vnrow))
  
  ### label and axes font sizes
  label.cex = 1.2
  set.cex.axes = 1
  set.mgp = c(1.9,0.5,0)
  plot.idx = 0
  
  if ( (len.colvec-2) <= (plot.ncol*plot.nrow)) {
    ### plot triplots in plot window
    dev.label=paste("plotter","tri",this$plotter.tri.num,sep=".")
    if ( length(which(dev.label==names(devList()))) == 0 ) {
      this$plotter.tri.num = this$plotter.tri.num + 1
      dev.label=paste("plotter","tri",this$plotter.tri.num,sep=".")
      devNew(type="x11",title="n-triploTs",width=plot.ncol*3.4,height=plot.nrow*3.7,label=dev.label)
      # mar in points, mai in inches
      # oma adds title lines
      # order: bottom, left, top, and right
      par(mfrow=c(plot.nrow,plot.ncol),oma=c(0.5,1,2,1),mar=c(3,3,4,2))
      
      this$plot.windows = c(this$plot.windows,dev.label)
    } else {
      devSet(devList()[which(dev.label==names(devList()))])
    }
    
    
    for ( v3 in 3:len.colvec ){
      #select columns to plot
      tdata=as.matrix(data[c(1,2,v3)])
      this$cnames = colnames(tdata)
      
      ### calculate cells which were not plotted 
      cells.overmaxFI = length(which( tdata[,1]>xmaxval | tdata[,2]>ymaxval ))
      cells.underminFI = length(which( tdata[,1]<xminval | tdata[,2]<yminval ))
      cells.overmaxFI.perc = round(100 * cells.overmaxFI / (dim(tdata)[1]-cells.underminFI))
      ### warn if more then 5% productive cells (q2+q3+q4) werent plotted
      if ( cells.overmaxFI.perc >= 5 & !this$working) {
        tkmessageBox(title = "Warning!",
                     message = sprintf("Your cells exceed %s%% of your plot max ranges. You might want to increase your max ranges.",cells.overmaxFI.perc), 
                     icon = "info", type = "ok")
      }
      
      if ( cutoff.var1 > 0 ) title.axis = sprintf("%s(%s)",var1,cutoff.var1)
      else title.axis = var1
      if ( cutoff.var2 > 0 ) title.axis = c(title.axis, sprintf("%s(%s)",var2,cutoff.var2))
      else title.axis = c(title.axis, var2)
      
      # start plot
      # mgp: A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. 
      # The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
      plot(1,type='n',frame.plot=FALSE,xlim=c(xminval,xmaxval+10*binSize),axes=FALSE
           ,ylim=c(yminval-2.5*binSize,ymaxval+5*binSize),xlab=title.axis[1]
           ,ylab=title.axis[2],cex.lab=label.cex,cex.axis=set.cex.axes,mgp=set.mgp)
      box(lwd=0.5,col="darkgrey")
      
      plot.idx = plot.idx + 1
      
      
	  ### draw axis on the bottom and on the left
	  axis(side=1, at=scale,labels=label,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
	  axis(side=2, at=scale,labels=label,las=3,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")

	  ### add grid
	  if(checkGRID=="1") {
	    xgrid.steps=seq(0,(xmaxval),by=grid.step)
	    ygrid.steps=seq(0,(ymaxval),by=grid.step)
	    abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
	  }

      
      ### calc quadrants in total
      ncells = ncells.total = nrow(tdata)
      # q1 Quadrant unten links
      # q2 Quadrant unten rechts
      # q3 Quadrant oben rechts
      # q4 Quadrant oben links
      if ( cutoff.var1 > 0 & cutoff.var2 > 0 ) {
        ### count cells in quadrant
        tdata.q1 = tdata[which( tdata[,1]<cutoff.var1 &  tdata[,2]<cutoff.var2 ),3]
        tdata.q2 = tdata[which( tdata[,1]>=cutoff.var1 &  tdata[,2]<cutoff.var2 ),3]
        tdata.q3 = tdata[which( tdata[,1]>=cutoff.var1 &  tdata[,2]>=cutoff.var2 ),3]
        tdata.q4 = tdata[which( tdata[,1]<cutoff.var1 &  tdata[,2]>=cutoff.var2 ),3]
        
        ### q[x].total [ink=black]
        ### percentage of cells in quadrant to total cells 
        ### or in MSI(+): percentage of cells in quadrant to total positive cells
        this$q1.total = abs(100 * length( tdata.q1 ) / ncells)
        this$q2.total = abs(100 * length( tdata.q2 ) / ncells)
        this$q3.total = abs(100 * length( tdata.q3 ) / ncells)
        this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
        
        if (cutoffs[v3] > 0 ) {
          ### number of cells which are producing cells in feature C
          ncells = nrow(tdata[which(tdata[,3]> cutoffs[v3]),])
          
          ### q[x].prodcells [ink=red]
          ### percentage of cells which are positive for feature C in quadrant to total quadrant cells
          this$q1.prodcells = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[v3])]) /length(tdata.q1)
          if (is.nan(this$q1.prodcells)) this$q1.prodcells = 0
          this$q2.prodcells = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[v3])]) /length(tdata.q2)
          if (is.nan(this$q2.prodcells)) this$q2.prodcells = 0
          this$q3.prodcells = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[v3])]) /length(tdata.q3)
          if (is.nan(this$q3.prodcells)) this$q3.prodcells = 0
          this$q4.prodcells = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[v3])]) /length(tdata.q4)
          if (is.nan(this$q4.prodcells)) this$q4.prodcells = 0
          
          ### only do MSI plots on producing cells only
          if ( checkCALC == "MSI(+)" ) {
            ### cut all cells which are not producing cells
            tdata.plus = tdata[which(tdata[,3]> cutoffs[3]),]
            ncells = nrow(tdata.plus)
            
            tdata.q1 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
            tdata.q2 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
            tdata.q3 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
            tdata.q4 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
            
            ### q[x].total [ink=blue]
            ### in MSI(+): percentage of cells in quadrant to total positive cells
            this$q1.total = abs(100 * length( which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
            this$q2.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
            this$q3.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ) ) / ncells)
            this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
          }
          
          ### q[x].prodcellsplus [ink=green]
          ### percentage of cells which are positive for feature C to total cells
          this$q1.prodcellsplus = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[v3])]) / ncells.total
          if (is.nan(this$q1.prodcellsplus)) this$q1.prodcellsplus = 0
          this$q2.prodcellsplus = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[v3])]) / ncells.total
          if (is.nan(this$q2.prodcellsplus)) this$q2.prodcellsplus = 0
          this$q3.prodcellsplus = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[v3])]) / ncells.total
          if (is.nan(this$q3.prodcellsplus)) this$q3.prodcellsplus = 0
          this$q4.prodcellsplus = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[v3])]) / ncells.total
          if (is.nan(this$q4.prodcellsplus)) this$q4.prodcellsplus = 0
        }
      }
      
      tdata.zero = tdata[apply(tdata[, -1], MARGIN = 1, function(row) {all(row > 0)}), ]
      ncells.zero=nrow(tdata.zero)
      
      ### calc quadrants with only positive values
      # q1 Quadrant unten links
      # q2 Quadrant unten rechts
      # q3 Quadrant oben rechts
      # q4 Quadrant oben links
      if ( cutoff.var1 > 0 & cutoff.var2 > 0 ) {
        q1.zero = abs(100 * length( which( tdata.zero[,1]<cutoff.var1 &  tdata.zero[,2]<cutoff.var2 ) ) / ncells.zero)      
        q2.zero = abs(100 * length( which( tdata.zero[,1]>=cutoff.var1 &  tdata.zero[,2]<cutoff.var2 ) ) / ncells.zero)                    
        q3.zero = abs(100 * length( which( tdata.zero[,1]>=cutoff.var1 &  tdata.zero[,2]>=cutoff.var2 ) ) / ncells.zero)                    
        q4.zero = abs(100 - q1.zero - q2.zero - q3.zero)
      }
      
      if (checkCALC == "density") {
        this$bintriplot(data=tdata,cutoffs=c(cutoff.var1,cutoff.var2,cutoffs[v3]),density=TRUE,binSize=binSize,mincells=mincount,overview=TRUE
                        ,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
      } else if (checkCALC == "MSI(+)") {
        this$bintriplot(data=tdata.plus,cutoffs=c(cutoff.var1,cutoff.var2,cutoffs[v3]),binSize=binSize,mincells=mincount,overview=TRUE
                        ,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp,quadrants.color = "blue", data.origin=tdata)
      } else {
        this$bintriplot(data=tdata,cutoffs=c(cutoff.var1,cutoff.var2,cutoffs[v3]),binSize=binSize,mincells=mincount,overview=TRUE
                        ,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
      }
      if (this$OverviewGate) {this$ncell.perc = this$ncell.sel/this$origin.ncells*100}
      
      # add title for single plot
      if ( checkCALC == "freq" | grepl("MSI",checkCALC) ) {
        if ( this$OverviewGate ) firstLine=sprintf("%s(%0.1f%%): %s,%s/cof=%s",displayfile,this$ncell.perc,checkCALC,this$cnames[3],this$current.cofactor)
        else firstLine=sprintf("%s: %s,%s/cof=%s",displayfile,checkCALC,this$cnames[3],this$current.cofactor)
      } else {
        if ( this$OverviewGate ) firstLine=sprintf("%s(%0.1f%%): %s/cof=%s",displayfile,this$ncell.perc,checkCALC,this$current.cofactor)
        else firstLine=sprintf("%s: %s/cof=%s",displayfile,checkCALC,this$current.cofactor)
      }
      title(main=firstLine,line=3.2,cex.main=1.0*label.cex,adj=0)
      
      secondLine=sprintf("%s-%s(%0.1f%%); min/max=%s/%s; binSize=%s,#bins=%s",ncells,ncells.zero,(ncells.zero/ncells*100),mincount,this$maxcells,binSize,this$bincount)
      title(main=secondLine,line=2.4,cex.main=0.7*label.cex,adj=0)
      
      thirdLine=""
      if ( cutoff.var1 > 0 & cutoff.var2 > 0 ) {
        thirdLine=sprintf("Q1=%0.1f/%0.1f Q2=%0.1f/%0.1f Q3=%0.1f/%0.1f Q4=%0.1f/%0.1f",q1.zero,this$q1.total,q2.zero,this$q2.total,q3.zero,this$q3.total,q4.zero,this$q4.total)
      }
      title(main=thirdLine,line=1.6,cex.main=0.7*label.cex,adj=0)
      
      if (checkGATED == "1" | grepl("temp",file)) {
        if ( length(this$coords.info)>2 ) {
          fourthLine=sprintf("%s; %s",this$coords.info[1],this$coords.info[2])
          fifthLine=sprintf("%s; %s",this$coords.info[3],this$coords.info[4])
        } else {
          fourthLine=sprintf("%s",paste(this$coords.info,collapse="; "))
          fifthLine=""
        }
        title(main=fourthLine,line=0.9,cex.main=0.6*label.cex,adj=0)
        title(main=fifthLine,line=0.1,cex.main=0.6*label.cex,adj=0)
      }
    }
    
    ### DONE plot triploTs in plot window
    
    
  } else {
    ### plot triploTs in pdf
    if ( titlelen <= 1 & this$working==FALSE) {
      titleans = tkmessageBox(title = "Are you sure?",
                              message = "Do not forget to set page title. Continue?", icon = "info", type = "yesno")
    } 
    
  if (tclvalue(titleans) == "yes") {
    if (this$working) this$saveinFolder = getwd()
    else this$saveinFolder = tk_choose.dir(default = getwd(), caption = "Select directory to save the PDFs")
    
    if ( !is.na(this$saveinFolder) ) {
      tkconfigure(this$tt, cursor = "watch")
      
      max.nhorizplots = 6
      max.nvertiplots = ceiling((len.colvec-2)/max.nhorizplots)
      
      timeSTART = Sys.time()
      cat("\n\n>>>> Start triploTOverviewXY with total data: \n\n")
      if (this$working) printf("w: %s - time started",timeSTART)
      
      toPDF(file=sprintf("%s_triploTOverviewXY_%s_%s_%s%s_%s.pdf",displayfile,var1,var2,len.colvec-2,checkCALC,this$version), 
            path=this$saveinFolder, 
            title=sprintf("triploTOverviewXY of %s",displayfile), 
            ### 
            width=3.21*max.nhorizplots,
            height=3.5*max.nvertiplots,
            pointsize=11,
            ###
            {
              ##### start triploT-OverviewXY
              par(mfrow=c(max.nvertiplots,max.nhorizplots),oma=c(0.5,1,5,1),mar=c(3,4,5,1))
              ###
              
              #add title for page in 3D-Overview tab
              #mtext(title, outer = TRUE, cex = 1.5,line=1.3,pos=2,xpd=TRUE)
              if (length(this$coords.info>0)) mtext(sprintf("(%s)",paste(this$coords.info,collapse="; ")),outer=TRUE,cex = 0.8)
              
              page.idx = 0
              current.page = 0
              for ( v3 in 3:len.colvec ){
                #select columns to plot
                tdata=as.matrix(data[c(1,2,v3)])
                this$cnames = colnames(tdata)
                
                ### calculate cells which were not plotted 
                cells.overmaxFI = length(which( tdata[,1]>xmaxval | tdata[,2]>ymaxval ))
                cells.underminFI = length(which( tdata[,1]<xminval | tdata[,2]<yminval ))
                cells.overmaxFI.perc = round(100 * cells.overmaxFI / (dim(tdata)[1]-cells.underminFI))
                ### warn if more then 5% productive cells (q2+q3+q4) werent plotted
                if ( cells.overmaxFI.perc >= 5 & !this$working) {
                  tkmessageBox(title = "Warning!",
                               message = sprintf("Your cells exceed %s%% of your plot max ranges. You might want to increase your max ranges.",cells.overmaxFI.perc), 
                               icon = "info", type = "ok")
                }
                
                if ( cutoff.var1 > 0 ) title.axis = sprintf("%s(%s)",var1,cutoff.var1)
                else title.axis = var1
                if ( cutoff.var2 > 0 ) title.axis = c(title.axis, sprintf("%s(%s)",var2,cutoff.var2))
                else title.axis = c(title.axis, var2)
                
                # start plot
                # mgp: A numeric vector of length 3, which sets the axis label locations relative to the edge of the inner plot window. 
                # The first value represents the location the labels (i.e. xlab and ylab in plot), the second the tick-mark labels, and third the tick marks. The default is c(3, 1, 0).
                plot(1,type='n',frame.plot=FALSE,xlim=c(xminval,xmaxval+10*binSize),axes=FALSE
                     ,ylim=c(yminval-2.5*binSize,ymaxval+5*binSize),xlab=title.axis[1]
                     ,ylab=title.axis[2],cex.lab=label.cex,cex.axis=set.cex.axes,mgp=set.mgp)
                box(lwd=0.5,col="darkgrey")
                
                plot.idx = plot.idx + 1
                
				                
				### draw axis on the bottom and on the left
				axis(side=1, at=scale,labels=label,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
				axis(side=2, at=scale,labels=label,las=3,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")

				### add grid
				if(checkGRID=="1") {
				  xgrid.steps=seq(0,(xmaxval),by=grid.step)
				  ygrid.steps=seq(0,(ymaxval),by=grid.step)
				  abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
				}
                
                ### calc quadrants in total
                ncells = ncells.total = nrow(tdata)
                # q1 Quadrant unten links
                # q2 Quadrant unten rechts
                # q3 Quadrant oben rechts
                # q4 Quadrant oben links
                if ( cutoff.var1 > 0 & cutoff.var2 > 0 ) {
                  ### count cells in quadrant
                  tdata.q1 = tdata[which( tdata[,1]<cutoff.var1 &  tdata[,2]<cutoff.var2 ),3]
                  tdata.q2 = tdata[which( tdata[,1]>=cutoff.var1 &  tdata[,2]<cutoff.var2 ),3]
                  tdata.q3 = tdata[which( tdata[,1]>=cutoff.var1 &  tdata[,2]>=cutoff.var2 ),3]
                  tdata.q4 = tdata[which( tdata[,1]<cutoff.var1 &  tdata[,2]>=cutoff.var2 ),3]
                  
                  ### q[x].total [ink=black]
                  ### percentage of cells in quadrant to total cells 
                  ### or in MSI(+): percentage of cells in quadrant to total positive cells
                  this$q1.total = abs(100 * length( tdata.q1 ) / ncells)
                  this$q2.total = abs(100 * length( tdata.q2 ) / ncells)
                  this$q3.total = abs(100 * length( tdata.q3 ) / ncells)
                  this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
                  
                  if (cutoffs[v3] > 0 ) {
                    ### number of cells which are producing cells in feature C
                    ncells = nrow(tdata[which(tdata[,3]> cutoffs[v3]),])
                    
                    ### q[x].prodcells [ink=red]
                    ### percentage of cells which are positive for feature C in quadrant to total quadrant cells
                    this$q1.prodcells = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[v3])]) /length(tdata.q1)
                    if (is.nan(this$q1.prodcells)) this$q1.prodcells = 0
                    this$q2.prodcells = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[v3])]) /length(tdata.q2)
                    if (is.nan(this$q2.prodcells)) this$q2.prodcells = 0
                    this$q3.prodcells = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[v3])]) /length(tdata.q3)
                    if (is.nan(this$q3.prodcells)) this$q3.prodcells = 0
                    this$q4.prodcells = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[v3])]) /length(tdata.q4)
                    if (is.nan(this$q4.prodcells)) this$q4.prodcells = 0
                    
                    ### only do MSI plots on producing cells only
                    if ( checkCALC == "MSI(+)" ) {
                      ### cut all cells which are not producing cells
                      tdata.plus = tdata[which(tdata[,3]> cutoffs[3]),]
                      ncells = nrow(tdata.plus)
                      
                      tdata.q1 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
                      tdata.q2 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
                      tdata.q3 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
                      tdata.q4 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
                      
                      ### q[x].total [ink=blue]
                      ### in MSI(+): percentage of cells in quadrant to total positive cells
                      this$q1.total = abs(100 * length( which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
                      this$q2.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
                      this$q3.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ) ) / ncells)
                      this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
                    }
                    
                    ### q[x].prodcellsplus [ink=green]
                    ### percentage of cells which are positive for feature C to total cells
                    this$q1.prodcellsplus = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[v3])]) / ncells.total
                    if (is.nan(this$q1.prodcellsplus)) this$q1.prodcellsplus = 0
                    this$q2.prodcellsplus = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[v3])]) / ncells.total
                    if (is.nan(this$q2.prodcellsplus)) this$q2.prodcellsplus = 0
                    this$q3.prodcellsplus = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[v3])]) / ncells.total
                    if (is.nan(this$q3.prodcellsplus)) this$q3.prodcellsplus = 0
                    this$q4.prodcellsplus = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[v3])]) / ncells.total
                    if (is.nan(this$q4.prodcellsplus)) this$q4.prodcellsplus = 0
                  }
                }
                
                tdata.zero = tdata[apply(tdata[, -1], MARGIN = 1, function(row) {all(row > 0)}), ]
                ncells.zero = nrow(tdata.zero)
                
                ### calc quadrants with only positive values
                # q1 Quadrant unten links
                # q2 Quadrant unten rechts
                # q3 Quadrant oben rechts
                # q4 Quadrant oben links
                if ( cutoff.var1 > 0 & cutoff.var2 > 0 ) {
                  q1.zero = abs(100 * length( which( tdata.zero[,1]<cutoff.var1 &  tdata.zero[,2]<cutoff.var2 ) ) / ncells.zero)      
                  q2.zero = abs(100 * length( which( tdata.zero[,1]>=cutoff.var1 &  tdata.zero[,2]<cutoff.var2 ) ) / ncells.zero)                    
                  q3.zero = abs(100 * length( which( tdata.zero[,1]>=cutoff.var1 &  tdata.zero[,2]>=cutoff.var2 ) ) / ncells.zero)                    
                  q4.zero = abs(100 - q1.zero - q2.zero - q3.zero)
                }
                
                if (checkCALC == "density") {
                  this$bintriplot(data=tdata,cutoffs=c(cutoff.var1,cutoff.var2,cutoffs[v3]),density=TRUE,binSize=binSize,mincells=mincount,overview=TRUE
                                  ,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
                } else if (checkCALC == "MSI(+)") {
                  this$bintriplot(data=tdata.plus,cutoffs=c(cutoff.var1,cutoff.var2,cutoffs[v3]),binSize=binSize,mincells=mincount,overview=TRUE
                                  ,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp,quadrants.color = "blue", data.origin=tdata)
                } else {
                  this$bintriplot(data=tdata,cutoffs=c(cutoff.var1,cutoff.var2,cutoffs[v3]),binSize=binSize,mincells=mincount,overview=TRUE
                                  ,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
                }
                if (this$OverviewGate) {this$ncell.perc = this$ncell.sel/this$origin.ncells*100}
                
                # add title for single plot
                if ( checkCALC == "freq" | grepl("MSI",checkCALC) ) {
                  if ( this$OverviewGate ) firstLine=sprintf("%s(%0.1f%%): %s(%s)/cof=%s",displayfile,this$ncell.perc,checkCALC,this$cnames[3],this$current.cofactor)
                  else firstLine=sprintf("%s: %s(%s)/cof=%s",displayfile,checkCALC,this$cnames[3],this$current.cofactor)
                } else {
                  if ( this$OverviewGate ) firstLine=sprintf("%s(%0.1f%%): %s/cof=%s",displayfile,this$ncell.perc,checkCALC,this$current.cofactor)
                  else firstLine=sprintf("%s: %s/cof=%s",displayfile,checkCALC,this$current.cofactor)
                }
                title(main=firstLine,line=3.2,cex.main=1.0*label.cex,adj=0)
                
                secondLine=sprintf("%s-%s(%0.1f%%); min/max=%s/%s; binSize=%s,#bins=%s",ncells.total,ncells.zero,(ncells.zero/ncells*100),mincount,this$maxcells,binSize,this$bincount)
                title(main=secondLine,line=2.4,cex.main=0.7*label.cex,adj=0)
                
                thirdLine=""
                if ( cutoff.var1 > 0 & cutoff.var2 > 0 ) {
                  thirdLine=sprintf("Q1=%0.1f/%0.1f Q2=%0.1f/%0.1f Q3=%0.1f/%0.1f Q4=%0.1f/%0.1f",q1.zero,this$q1.total,q2.zero,this$q2.total,q3.zero,this$q3.total,q4.zero,this$q4.total)
                }
                title(main=thirdLine,line=1.6,cex.main=0.7*label.cex,adj=0)
                
                if (checkGATED == "1" | grepl("temp",file)) {
                  if ( length(this$coords.info)>2 ) {
                    fourthLine=sprintf("%s; %s",this$coords.info[1],this$coords.info[2])
                    fifthLine=sprintf("%s; %s",this$coords.info[3],this$coords.info[4])
                  } else {
                    fourthLine=sprintf("%s",paste(this$coords.info,collapse="; "))
                    fifthLine=""
                  }
                  title(main=fourthLine,line=0.9,cex.main=0.6*label.cex,adj=0)
                  title(main=fifthLine,line=0.1,cex.main=0.6*label.cex,adj=0)
                }
              }
              #add title for page in 3D-Overview tab
              if (current.page != page.idx) {
                mtext(title, outer = TRUE, cex = 1.5,line=1.0,pos=2,xpd=TRUE)
                if (length(this$coords.info>0)) mtext(sprintf("(%s)",paste(this$coords.info,collapse="; ")),outer=TRUE,cex = 0.8)
                #printf("Print on Page %s",page.idx)
                current.page = page.idx
              }
        })
        tkconfigure(this$tt, cursor = "left_ptr")
        
        cat("\n")
        if (!this$working) tkmessageBox(title = "Output of overview",
                                        message = paste ("Plots created in folder ",this$saveinFolder, sep=""))
        
        print(Sys.time()-timeSTART)
        cat("\n>>>> Done.\n")
        this$OverviewGate=FALSE
      }
    }
  }	
}

fcs$bintriplot <- function(
	data,
	cutoffs,
	set.cex=1.1,
	set.cex.axes=1.0,
	set.mgp=c(1.5, 0.3, 0),
	binSize=0.2,
	mincells=10,
	density=FALSE,
	png=FALSE,
	bg=FALSE,
	overview=FALSE,
	checkCALC=NA,
	quadrants.color="black",
	data.origin=NA,
	file=NA) 
  {
	this=fcs

	prodcells.color="red"
	prodpluscells.color="chartreuse4"
	
	if (is.na(file)) file=tclvalue(tkget(this$tkchoosefile))
	
	displayfile = this$shortenFilename(file)
	
	metadatafile = sprintf("%s_PRItrivis_metadata.csv",this$current.project)
	if(!file.exists(metadatafile)) {
	  header = c("date","sample", "cofactor","calc"
	             ,"feat.A", "feat.B", "feat.C", "absRange.C"
	             ,"cutoff.A","cutoff.B","cutoff.C"
	             ,"q1.total","q2.total","q3.total","q4.total"
	             ,"q1.prodcells","q2.prodcells","q3.prodcells","q4.prodcells"
	             ,"q1.prodcellsplus","q2.prodcellsplus","q3.prodcellsplus","q4.prodcellsplus")
	  write.table(t(header), metadatafile, sep = "\t", row.names=F, col.names=F)
	}
	current.date = format(Sys.Date(),"%y%m%d")
	  
	data = as.matrix(data)
	if (!is.na(data.origin)) data.origin = as.matrix(data.origin)

	# axes range
	xmin.val=as.numeric(tkget(this$minvalX))
	xmax.val=as.numeric(tkget(this$maxvalX))
	ymin.val=as.numeric(tkget(this$minvalY))
	ymax.val=as.numeric(tkget(this$maxvalY))
	min.MSI=as.double(tclvalue(this$vminMSI))
	max.MSI=as.double(tclvalue(this$vmaxMSI))
	
	if (max.MSI < min.MSI) {
  	tmp = min.MSI 
  	min.MSI = max.MSI
  	max.MSI = tmp
	}
	
	# checkbutton options
	checkDYNRANGE = tclvalue(this$cbtdynRange)
	checkTRANS = tclvalue(this$rbtrans)
	if (checkTRANS =="") checkTRANS = tclvalue(this$rbtrans) = "asinh"
	checkCALC = tclvalue(this$rbcalc)
	if (checkCALC == "") checkCALC = tclvalue(this$rbcalc) = "MSI"
	if (checkCALC == "density") density=TRUE
	checkPERCENTAGE = tclvalue(this$cbtshowPercentage)
	if ( checkPERCENTAGE == "1" ) {
		this$plot.percentage = TRUE
	} else {
		this$plot.percentage = FALSE
	}
	checkSHOWMINBIN = tclvalue(this$cbtshowMinBins)
	
	if (this$working) printf("w: do bintriplot density=%s checkCALC=%s",density,checkCALC)
	
	if (checkTRANS=="asinh") scipen = this$asinh$scipen
	else scipen = this$biex$scipen
	
	options(scipen=scipen)

	# legend from blue to red
	if (bg) cols=rep("gray",12)
	else cols=this$col.rainbow
	# colors for bins <mincount
	cols.pale = this$col.rainbow.pale

	# legend title
	tmp = unlist(strsplit(colnames(data)[3],"\\."))
	if ( length(tmp) > 1 & all(lengths(tmp)>1) ) {
		if ( cutoffs[3]>0 ) legend.title = sprintf("%s(%s)",tmp[1],cutoffs[3]) 
		else legend.title = tmp[1]
	} else {
		if ( cutoffs[3]>0 ) legend.title = sprintf("%s(%s)",colnames(data)[3],cutoffs[3]) 
		else legend.title = colnames(data)[3]
	}

	# boolean for only grey plot 
	# if MSI(+) mode and there are no colorful bins to display
	grey.label = TRUE
	my.LENGTHS = FALSE
			
	# set negative values of z-axis (colnum=3) to zero
	if (!density) data[which(data[,3]<0),3] = 0
	
	ncells=nrow(data)

	if ( ncells > 0 ) {

	### tdata = cut cells which lie in plot area
	#tdata2 = data[-c(which(data[,1]<xmin.val),which(data[,2]<xmin.val)), ]
	tdata = data[which(data[,1]>=xmin.val & data[,2]>=ymin.val), ]

	### construct bins
	fX=cut(tdata[,1],breaks=seq(xmin.val,xmax.val,by=binSize),include.lowest=TRUE,dig.lab=5)
	fY=cut(tdata[,2],breaks=seq(ymin.val,ymax.val,by=binSize),include.lowest=TRUE,dig.lab=5)
	tab=table(fX,fY)

	colnames(tab)=seq(ymin.val,ymax.val-binSize,by=binSize)
	rownames(tab)=seq(xmin.val,xmax.val-binSize,by=binSize)
	fXY=as.factor(paste(fX,fY))

	### construct bins for MSI(+)
	if (checkCALC == "MSI(+)") {
		fX.origin=cut(data.origin[,1],breaks=seq(xmin.val,xmax.val,by=binSize),include.lowest=TRUE,dig.lab=5)
		fY.origin=cut(data.origin[,2],breaks=seq(ymin.val,ymax.val,by=binSize),include.lowest=TRUE,dig.lab=5)
		tab.origin=table(fX.origin,fY.origin)
		colnames(tab.origin)=seq(ymin.val,ymax.val-binSize,by=binSize)
		rownames(tab.origin)=seq(xmin.val,xmax.val-binSize,by=binSize)
	} else {
		tab.origin = tab
	}

	if (density) {
		# number of cells in bin
		my.calc=aggregate(tdata[,3],by=list(fXY),length)
	} else {
		# get means/median/freq
		if ( grepl("MSI",checkCALC) ) {
			my.calc=aggregate(tdata[,3],by=list(fXY),mean)
		#else if ( checkCALC == "medianFI" )  my.calc=aggregate(tdata[,3],by=list(fXY),median)
		} else if ( checkCALC == "SD" ) {
			my.calc=aggregate(tdata[,3],by=list(fXY),sd)
			cols=this$col.blackwhite
		} else if ( checkCALC == "SEM" ) {
			my.calc=aggregate(tdata[,3],by=list(fXY),function(x) {
				SEM = sd(x)/sqrt(length(x))
				# if normally distributed, 95,4 % of the cells should lie inside the interval mean +/- SEM
				#interval_min = mean(x) - SEM
				#interval_max = mean(x) + SEM 
				SEM
			})
			cols=this$col.blackwhite
		} else if ( checkCALC == "RSEM" ) {
			my.calc=aggregate(tdata[,3],by=list(fXY),function(x) {
				RSEM = sd(x)/sqrt(length(x))
				RSEM/mean(x)*100
			})
		} else if ( checkCALC == "freq" ) {
			my.calc = aggregate(tdata[,3],by=list(fXY),function(x) {
				y= round( 100 * length(which(x >= cutoffs[3])) / length(x))
				return(y)
			})
		} 
	}
	### http://www.allgemeinmedizin.med.uni-goettingen.de/de/media/2008_Koschack_Standardabweichung_Standardfehler.pdf
	##  Standardabweichung (SD) 
	# – ist eine Aussage über die Streuung der erhobenen Werte in einer Stichprobe
	# – ist nur wenig durch die Grösse der Stichprobe beeinflussbar
	# – hängt von der biologischen Variabilität ab
	# – ist ein beschreibendes Mass
	## Standardfehler (SEM)
	# – ist eine Aussage über die „Genauigkeit“ des Mittelwerts in einer Stichprobe
	# – hängt von der Messgenauigkeit ab
	# – ist ein statistisches Mass
	# – steht in direktem Verhältnis zur Grösse der Stichprobe

	### https://www.graphpad.com/guides/prism/6/statistics/index.htm?stat_semandsdnotsame.htm
	# It is easy to be confused about the difference between the standard deviation (SD) and the standard error of the mean (SEM). Here are the key differences:
	# •   The SD quantifies scatter — how much the values vary from one another.
	# •   The SEM quantifies how precisely you know the true mean of the population. It takes into account both the value of the SD and the sample size.
	# •   Both SD and SEM are in the same units -- the units of the data.
	# •   The SEM, by definition, is always smaller than the SD.
	# •   The SEM gets smaller as your samples get larger. This makes sense, because the mean of a large sample is likely to be closer to the true population mean than is the mean of a small sample. With a huge sample, you'll know the value of the mean with a lot of precision even if the data are very scattered.
	# •   The SD does not change predictably as you acquire more data. The SD you compute from a sample is the best possible estimate of the SD of the overall population. As you collect more data, you'll assess the SD of the population with more precision. But you can't predict whether the SD from a larger sample will be bigger or smaller than the SD from a small sample. (This is not strictly true. It is the variance -- the SD squared -- that doesn't change predictably, but the change in SD is trivial and much much smaller than the change in the SEM.)
	# Note that standard errors can be computed for almost any parameter you compute from data, not just the mean. The phrase "the standard error" is a bit ambiguous. The points above refer only to the standard error of the mean.


	my.lengths=aggregate(tdata[,3],by=list(fXY),length)
	
	my.LENGTHS = any(my.lengths[,2] >= mincells)
	}
  absRange = 0
	### if there are bins to display
	if ( my.LENGTHS ) {
		# there are bins to plot, so set grey.label to FALSE
		grey.label = FALSE

		my.calc=cbind(my.calc,ncells=my.lengths$x)

		if (FALSE) {if (checkCALC == "MEAN_SEM" & !density) {
			this$my.test = aggregate(tdata[,3],by=list(fXY),function(x) {
				SEM = sd(x)/sqrt(length(x))

				if ( SEM >= 0.5 ) {
					return(1)
				} else { 
					return (0) 
				}
				# if normally distributed, 95,4 % of the cells should lie inside the interval mean +/- 2*SEM
			})
		}
		}
		
		# delete rows with NA!
		# my.calc.NA = my.calc[-grep('NA',my.calc$Group.1),]

		### get steps for legend and plot in mode 'freq'
		decim=1
		#range=""
		col.minmax="black"
		if (checkCALC == "freq" & !density) {
			#min.legend=ymin.val + 3*binSize
			#max.legend=max.range=(ymax.val-ymin.val)/2 #- (ymax.val-ymin.val)/10
			#max.legend = max.range = ymin.val + 3*binSize + diff(c(ymin.val,ymax.val))/2
			max.range = ymin.val + 3*binSize + diff(c(ymin.val,ymax.val))/2
			
			label.steps=seq(0,100,by=10)
			
			# bin color factor
			my.calc.fac=cut(my.calc$x,breaks=seq(0,100,by=10),labels=1:10,include.lowest=TRUE)
			levels(my.calc.fac)=c(0,levels(my.calc.fac),11,12)

			decim=0
		} else if (checkCALC == "RSEM" & !density) {
			#min.legend=ymin.val
			#max.legend=max.range=(ymax.val-ymin.val)/3
			min.legend = ymin.val + 3*binSize
			max.legend = max.range = ymin.val + 3*binSize + diff(c(ymin.val,ymax.val))/3
		
			step=round(diff(range(max.legend,min.legend))/6,1)
			steps=seq(min.legend,max.legend,by=step)
			label.steps=seq(25,55,by=5)
			
			label.steps[7]=""
			label.steps[6]=">=50"
			label.steps[2]="<=25"
			label.steps[1]=""

			# bin color factor
			my.calc.fac=cut(my.calc$x,breaks=seq(0,50,by=5),labels=1:10,include.lowest=TRUE)
			levels(my.calc.fac)=c(0,levels(my.calc.fac),11,12)
			this$my.calc.fac2 = my.calc.fac
			for ( i in 1:length(my.calc.fac) ){
				if ( !is.na(my.calc$x[i]) ) {
					if (my.calc$x[i]>=50 & my.calc$ncells[i]>=mincells) {
						my.calc.fac[i] = 11
					}
				}
			}
			cols = this$col.blackred
			decim=0
		} else {
			# get steps for legend and plot
			if (density) {
				idx=which(my.calc$ncells>=mincells)
				#idx=idx[grep("NA",my.calc[idx,'Group.1'])]
				min.range=floor(min(my.calc[idx,'x'])*10)/10
				max.range=max(tab)
				#max.range = max(my.lengths[,2])
				decim=0
				#printf("maxrange=%s",max.range)

				if ( max.range < 200 ) col.minmax="red"
			} else if ( checkDYNRANGE=="1" ) {
				idx=which(my.calc$ncells>=mincells)
				#idx=idx[grep("NA",my.calc[idx,'Group.1'])]
				min.range=floor(min(my.calc[idx,'x'])*10)/10
				max.range=ceiling(max(my.calc[idx,'x'])*10)/10

				# if dynamic range is too small
				if ( grepl("MSI",checkCALC) & diff(c(min.range,max.range)) <= 0.5 ) col.minmax="red"
			} else {
				min.range=min.MSI
				max.range=max.MSI
				#range=sprintf("Manual MSI range: %0.1f-%0.1f",min.range,max.range)
			}

			# get steps
			step=round(diff(range(max.range,min.range))/10,2) 
			steps=seq(min.range,max.range,by=step)
			label.steps=steps[1:11]
			if (density | checkTRANS=="biex") {
				if(max.range>500) label.steps = round(label.steps,-2)
				else label.steps = round(label.steps)
				if (density) label.steps[1] = min.range
			} else if ( checkDYNRANGE != "1" ) {
				if (label.steps[1]!="0") label.steps[1]=sprintf("<=%s",label.steps[1])
				label.steps[11]=sprintf(">=%s",label.steps[11])
			} else {
				label.steps = round(label.steps,1)
			}

			# bin color factor
			my.calc.fac=cut(my.calc$x,breaks=steps,labels=2:11,include.lowest=TRUE)

			levels(my.calc.fac)=c(0,levels(my.calc.fac),12)
			# if x < min.range
			my.calc.fac[which(my.calc$x<steps[1] & my.calc$ncells>=mincells)]=0
			# if x > max.range
			my.calc.fac[which(my.calc$x>steps[11] & my.calc$ncells>=mincells)]=11
		}
		my.calc=cbind(my.calc,fac=as.numeric(my.calc.fac)+1)
		
		### get absolute bin range of all bins with mincells
		absRange = round(diff(range(my.calc$x[which(my.calc$ncells>=mincells)])),3)
		
		#this$my.calc.fac = my.calc.fac
		this$my.calc = my.calc

		this$bincount = 0
		this$maxcells = 0
		
		##### plot bins
		for (x in rownames(tab)) {
			for (y in colnames(tab)) {
				if ( tab[x,y]>=mincells ) {
					brackets.open = c("(","[")
		      brackets.idx.x = brackets.idx.y = 1
		      if (x==0) brackets.idx.x = 2
		      if (y==0) brackets.idx.y = 2
		      
		      fact = as.factor(paste(brackets.open[brackets.idx.x],x,',',as.numeric(x)+binSize,'] ',
		                             brackets.open[brackets.idx.y],y,',',as.numeric(y)+binSize,']',sep=''))
					idx=which(as.character(fact)==as.character(my.calc$Group.1))
					
					rect(x,y,as.numeric(x)+binSize,as.numeric(y)+binSize,col=cols[my.calc[idx,'fac']],border=NA)
					
					this$bincount = this$bincount + 1                

					if (tab[x,y]>this$maxcells) this$maxcells=tab[x,y]            
				} else if ( checkCALC == "MSI(+)" & tab.origin[x,y] >= mincells) {
					rect(x,y,as.numeric(x)+binSize,as.numeric(y)+binSize,col="gray",border=NA)
				} else if ( tab[x,y]>0 & checkSHOWMINBIN=="1" ){
				  rect(x,y,as.numeric(x)+binSize,as.numeric(y)+binSize,col=cols.pale[my.calc[idx,'fac']],border=NA)
				}
			}
		}

		### space calculation
		rect.size = diff(c(xmin.val,xmax.val))/25
		rect.size = 0.2

		### add production line
		this$addProdline(cutoffs)


		###### legend plot + label
		#if (density | checkCALC == "SEM" | checkCALC == "SD") {
		if ( checkCALC != "RSEM" | density) {
			#min.legend=ymin.val
			#max.legend=(ymax.val-ymin.val)/2 
			min.legend = ymin.val + 3*binSize
			max.legend = ymin.val + 3*binSize + diff(c(ymin.val,ymax.val))/2 

			step=round(diff(range(max.legend,min.legend))/11,1)
			steps=seq(min.legend,max.legend+step,by=step)
			set.cex =1.2
		}

		##### if its not background and not a png minimalistic picture, then print legend and label
		if ( (!bg | density) & !png & this$plot.percentage) {
			if ( cutoffs[1] > 0 & cutoffs[2] > 0 & cutoffs[3] > 0 & !density) {
				#if (this$working) print("w: add quadrant and product percentages on triplot")
				### quadrant left lower
				text(par()$usr[1]-0.01*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.03*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q1.total),col=quadrants.color,cex=1.00*set.cex,pos=4,xpd=TRUE)
				text(par()$usr[1]-0.01*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.08*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q1.prodcells),col=prodcells.color,cex=1.00*set.cex,pos=4,xpd=TRUE)
				text(par()$usr[1]-0.01*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.13*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q1.prodcellsplus),col=prodpluscells.color,cex=1.00*set.cex,pos=4,xpd=TRUE)

				### quadrant right lower
				text(par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.03*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q2.total),col=quadrants.color,cex=1.00*set.cex,pos=2,xpd=TRUE)
				text(par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.08*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q2.prodcells),col=prodcells.color,cex=1.00*set.cex,pos=2,xpd=TRUE)
				text(par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.13*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q2.prodcellsplus),col=prodpluscells.color,cex=1.00*set.cex,pos=2,xpd=TRUE)

				### quadrant right upper
				text(par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.04*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q3.total),col=quadrants.color,cex=1.00*set.cex,pos=2,xpd=TRUE)
				text(par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.09*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q3.prodcells),col=prodcells.color,cex=1.00*set.cex,pos=2,xpd=TRUE)
				text(par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.14*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q3.prodcellsplus),col=prodpluscells.color,cex=1.00*set.cex,pos=2,xpd=TRUE)

				### quadrant left upper
				text(par()$usr[1]-0.01*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.04*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q4.total),col=quadrants.color,cex=1.00*set.cex,pos=4,xpd=TRUE)
				text(par()$usr[1]-0.01*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.09*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q4.prodcells),col=prodcells.color,cex=1.00*set.cex,pos=4,xpd=TRUE)
				text(par()$usr[1]-0.01*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.14*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q4.prodcellsplus),col=prodpluscells.color,cex=1.00*set.cex,pos=4,xpd=TRUE)
			  
				
				quadrant.percs <- cbind(current.date,displayfile,this$current.cofactor,checkCALC
				                        ,colnames(data)[1],colnames(data)[2],colnames(data)[3],absRange
				                        ,cutoffs[1],cutoffs[2],cutoffs[3]
				                        ,round(this$q1.total,1),round(this$q2.total,1),round(this$q3.total,1),round(this$q4.total,1)
				                        ,round(this$q1.prodcells,1),round(this$q2.prodcells,1),round(this$q3.prodcells,1),round(this$q4.prodcells,1)
				                        ,round(this$q1.prodcellsplus,1),round(this$q2.prodcellsplus,1),round(this$q3.prodcellsplus,1),round(this$q4.prodcellsplus,1))
				write.table(quadrant.percs, sprintf("%s_PRItri_metadata.csv",this$current.project), sep = "\t", col.names = F, row.names = F, append = T)
				printf("Quadrant percs written in %s",sprintf("%s_PRItri_metadata.csv",this$current.project))
				
			} else if ( cutoffs[1] > 0 & cutoffs[2] ) {
				text(par()$usr[1]-0.01*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.03*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q1.total),col=quadrants.color,cex=1.00*set.cex,pos=4,xpd=TRUE)
				text(par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.03*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q2.total),col=quadrants.color,cex=1.00*set.cex,pos=2,xpd=TRUE)
				text(par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.04*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q3.total),col=quadrants.color,cex=1.00*set.cex,pos=2,xpd=TRUE)
				text(par()$usr[1]-0.01*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.04*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q4.total),col=quadrants.color,cex=1.00*set.cex,pos=4,xpd=TRUE)
				
				quadrant.percs <- cbind(current.date,displayfile,this$current.cofactor,checkCALC
				                        ,colnames(data)[1],colnames(data)[2],colnames(data)[3],absRange
				                        ,cutoffs[1],cutoffs[2],"NA"
				                        ,round(this$q1.total,1),round(this$q2.total,1),round(this$q3.total,1),round(this$q4.total,1))
				write.table(quadrant.percs, sprintf("%s_PRItri_metadata.csv",this$current.project), sep = "\t", col.names = F, row.names = F, append = T)
				printf("Quadrant percs written in %s",sprintf("%s_PRItri_metadata.csv",this$current.project))
			} else {
			  ### no cutoffs are set
			  ### write in csv at least the absRange
			  quadrant.percs <- cbind(current.date,displayfile,this$current.cofactor,checkCALC
			                          ,colnames(data)[1],colnames(data)[2],colnames(data)[3],absRange)
			  write.table(quadrant.percs, sprintf("%s_PRItri_metadata.csv",this$current.project), sep = "\t", col.names = F, row.names = F, append = T)
			}

			
			checkLEGEND = tclvalue(this$cbtshowLegend)
			this$step=step
			this$label.steps=label.steps
			space = 0.08*(par()$usr[4]-par()$usr[3])

			### not for history png
			#if (!png) {
			if (checkLEGEND=="1") {
			  ##### legend title
			  if (density) legend.title = "# cells "
			  if (!bg) text(par()$usr[2]+0.02*(par()$usr[2]-par()$usr[1]),steps[10]+4.7*rect.size+this$legend.space+space,label=legend.title,cex=0.85*set.cex,pos=2)
			  
				label.pos.x = par()$usr[2]-0.12*(par()$usr[2]-par()$usr[1])
				
				for (i in 1:11 ) {
					## print legend rectangles
				  if (i<11) {
  				  rect(xleft = par()$usr[2]-0.13*(par()$usr[2]-par()$usr[1]),
  				       ybottom = steps[i] + space,
  				       xright = par()$usr[2]-0.105*(par()$usr[2]-par()$usr[1]),
  				       ytop = steps[i] + space + step,
  				       col=cols[i+1],border=NA,pos=2)
				  }
				  
				  if ( checkDYNRANGE != "1" & (i==1 | i==11) ) {
				    displaylabel=label.steps[i]
				  } else {
				    displaylabel = sprintf("%.1f",as.numeric(label.steps[i]))
				  }
				    
					if (checkCALC == "RSEM" & !density) {
						text(x=label.pos.x,y=steps[i]+space,label=displaylabel,col=col.minmax,cex=0.65*set.cex,pos=4)
					} else if ( i==1 | i==11) {
						text(x=label.pos.x,y=steps[i]+space,label=displaylabel,col=col.minmax,cex=0.65*set.cex,pos=4)
					} else if ( i==6 & step>=(0.65*binSize) ) {
						text(x=label.pos.x,y=steps[i]+space,label=displaylabel,col=col.minmax,cex=0.65*set.cex,pos=4)
					} else if ( step>=(1.8*binSize) ) { 
						if (i==3) {
							display.label = sprintf("%.1f",(as.numeric(label.steps[i]) + (step/2)) )
							text(label.pos.x,steps[i]+space+0.5*step,label=display.label,
							     col=col.minmax,cex=0.65*set.cex,pos=4)
						}
						if (i==9) {
							display.label = sprintf("%.1f", (as.numeric(label.steps[i]) - (step/2)) )
							text(label.pos.x,steps[i]+space-0.5*step,label=display.label,
							     col=col.minmax,cex=0.65*set.cex,pos=4)
						}
					} 
				}
			}

		}

	}

	# if there are no bins to display in MSI(+) mode
	if (grey.label) {
		text(par()$usr[2],ymin.val+10*binSize,label=legend.title,cex=0.85*set.cex,pos=2)
		
		### quadrant left lower
		text(par()$usr[1]-0.01*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.03*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q1.total),col=quadrants.color,cex=1.00*set.cex,pos=4,xpd=TRUE)
		text(par()$usr[1]-0.01*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.08*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q1.prodcells),col=prodcells.color,cex=1.00*set.cex,pos=4,xpd=TRUE)
		
		### quadrant right lower
		text(par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.03*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q2.total),col=quadrants.color,cex=1.00*set.cex,pos=2,xpd=TRUE)
		text(par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.08*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q2.prodcells),col=prodcells.color,cex=1.00*set.cex,pos=2,xpd=TRUE)
		
		### quadrant right upper
		text(par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.04*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q3.total),col=quadrants.color,cex=1.00*set.cex,pos=2,xpd=TRUE)
		text(par()$usr[2]+0.01*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.09*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q3.prodcells),col=prodcells.color,cex=1.00*set.cex,pos=2,xpd=TRUE)
		
		### quadrant left upper
		text(par()$usr[1]-0.01*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.04*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q4.total),col=quadrants.color,cex=1.00*set.cex,pos=4,xpd=TRUE)
		text(par()$usr[1]-0.01*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.09*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",this$q4.prodcells),col=prodcells.color,cex=1.00*set.cex,pos=4,xpd=TRUE)
		
		quadrant.percs <- cbind(current.date,displayfile,this$current.cofactor,checkCALC
		                        ,colnames(data)[1],colnames(data)[2],colnames(data)[3],absRange
		                        ,cutoffs[1],cutoffs[2],cutoffs[3]
		                        ,round(this$q1.total,1),round(this$q2.total,1),round(this$q3.total,1),round(this$q4.total,1)
		                        ,round(this$q1.prodcells,1),round(this$q2.prodcells,1),round(this$q3.prodcells,1),round(this$q4.prodcells,1))
		write.table(quadrant.percs, sprintf("%s_PRItri_metadata.csv",this$current.project), sep = "\t", col.names = F, row.names = F, append = T)
		#printf("Quadrant percs written in %s",sprintf("%s_PRItri_metadata.csv",this$current.project))
		
		if (!bg & !png & !overview) tkmessageBox(title = "Insufficient cell count.",
			message = sprintf("No bins to display (ncells=%s). Lower cutoff(z).",ncells), icon = "info", type = "ok")
	}

	# in mode triploT: if autorect was selected
	if (tclvalue(this$cbtautoRect) == "1") {
		this$addRectInfo(setcex=set.cex)
		this$addCellInfo(setcex=set.cex)
	} 
}























fcs$autoSetCutoff <- function() {
	this=fcs

	if(this$working) printf("Calculating cutoffs start.")
	
	table = this$selected.project
	file=tclvalue(tkget(this$tkchoosefile))
	file.idx=this$current.filetable[which(this$current.filetable[,2]==file),1]
	this$selected.filenum = file.idx

	if ( this$current.project != table | this$current.filenum != file.idx | 
	     length(this$data)==0  | this$current.cofactor != as.numeric(tclvalue(this$rbasinh))) {
	  this$getFile(table,file.idx)
	}
	
	checkTRIMMING = tclvalue(this$cbttrimming)
	checkDOUBLETS = tclvalue(this$cbtremovedub)

	vars= this$selected.vars
	len = length(this$selected.vars)

	# column and cutoff number vector
	colvec=vector()
	for (i in 1:len) {
		if ( tclvalue(this$cbVal[[i]]) == "1") colvec=c(colvec,i)
	}
	len.colvec = length(colvec)

	### if trimming or remove doublets desired
	if ( checkDOUBLETS =="1" ) this$preprocData()
  	if ( checkTRIMMING =="1" ) this$preprocData(mode="trim")

	### select columns to plot
	data = this$data[colvec]

	### open plot window
	dev.label="histogram.autogate"
	if ( length(which(dev.label==names(devList()))) != 0 ) {
		devSet(devList()[which(dev.label==names(devList()))])
	} else {
		if ( len.colvec>36 ) {
			devNew(type="x11",title="Histograms for cutoffs",width=7*2.8,height=ceiling(len.colvec/7)*2.0,label=dev.label)
			par(mfrow=c(ceiling(len.colvec/7),7),oma=c(1.5,1,1.5,1),mar=c(2,2.4,1.5,1),mgp=c(1.3,0.4,0))
		} else if ( len.colvec>25 ) {
			devNew(type="x11",title="Histograms for cutoffs",width=6*2.5,height=ceiling(len.colvec/6)*2.0,label=dev.label)
			par(mfrow=c(ceiling(len.colvec/6),6),oma=c(1.5,1,1.5,1),mar=c(2,2.4,1.5,1),mgp=c(1.3,0.4,0))
		} else if ( len.colvec>16 ) {
			devNew(type="x11",title="Histograms for cutoffs",width=5*3,height=ceiling(len.colvec/5)*2.3,label=dev.label)
			par(mfrow=c(ceiling(len.colvec/5),5),oma=c(1.5,1,1.5,1),mar=c(2,2.4,1.5,1),mgp=c(1.3,0.4,0))
		} else if ( len.colvec>9 ) {
			devNew(type="x11",title="Histograms for cutoffs",width=4*4,height=ceiling(len.colvec/4)*2.5,label=dev.label)
			par(mfrow=c(ceiling(len.colvec/4),4),oma=c(1.5,1,1.5,1),mar=c(2,2.4,1.5,1),mgp=c(1.3,0.4,0))
		} else {
			devNew(type="x11",title="Histograms for cutoffs",width=3*4,height=ceiling(len.colvec/3)*3,label=dev.label)
			par(mfrow=c(ceiling(len.colvec/3),3),oma=c(1.5,1,1.5,1),mar=c(2,2.4,1.5,1),mgp=c(1.3,0.4,0))
		}
		this$plot.windows = c(this$plot.windows,dev.label)
	}
	
	this$cutoff_list = list()
	xmaxval = vector()
	xminval = vector()
	for ( i in 1:len.colvec ) {

		tdata=data[,i]

		d=density(tdata)
		x=d$x
		y=d$y

		### look only at higher levels of y-axis
		# index where y starts at position greater than 1/20th of max(d$y)
		start.y = 1/10*max(y)
		id_y = which(y>start.y)
		### skip compensation problems/unbinded positions
		# index where x is positive
		start.x = 1/10*max(x)
		id_x = which(x>start.x)
		# move 5 positions earlier
		xminval[i] = x[intersect(id_y,id_x)[1]-5]
		xmaxval[i] = round(max(tdata))

		
		stepsize= 6:3
		for ( step in 1:length(stepsize) ) {
			incline_list = list()
			incline_x = vector()
			idx=1
			skip=FALSE
			cutoff=0
			for ( range in 1:(length(y)-stepsize[step]) ) {
				tmp = min(y[range:(range+stepsize[step])])
				diff = diff(c(y[range],y[range+stepsize[step]]))

				# if min(y) in range is more than 1/20th of peak and if x>2.5
				if (tmp > start.y & x[range]>2.5) {
					### incline_list
					# [1] step id
					# [2] x-value from step
					# [3] incline difference in range (x[range],x[range+step])
					### incline_x:  x-value from step
					incline_list[[idx]]=range
					incline_list[[idx]][2]=round(x[range],3)
					incline_list[[idx]][3]=diff
					incline_x[idx]=x[range]
					idx = idx + 1
				}
			}

			if (idx==1) {
				for ( range in 1:(length(y)-stepsize[step]) ) {
					tmp = min(y[range:(range+stepsize[step])])
					diff = diff(c(y[range],y[range+stepsize[step]]))

					# if min(y) in range is more than 1/20th of peak and if x>2.5
					if (tmp > start.y & x[range]>start.x) {
						### incline_list
						# [1] step id
						# [2] x-value from step
						# [3] incline in range (x[range],x[range+step])
						### incline_x:  x-value from step
						incline_list[[idx]]=round(range,2)
						incline_list[[idx]][2]=round(x[range],3)
						incline_list[[idx]][3]=diff
						incline_x[idx]=x[range]
						idx = idx + 1
					}
				}
			}

			if (idx==1) skip=TRUE


			if (!skip) {
				idx.minima = idx.maxima = idx.shoulder = vector()
				for ( inc in 1: (length(incline_x)-1) ) {
					#printf("inc=%s list[3]=%s",inc,incline_list[[inc]][3])
				
					### check if there is at least one minima
				  # difference at position inc should be -, and at position inc+1 should be +
					if (incline_list[[inc]][3] < 0 & incline_list[[inc+1]][3] > 0) {
						#printf("-turn inc[%s]=%0.5f inc[%s]=%0.5f x=%s",inc,incline_list[[inc]][3],inc+1,incline_list[[inc+1]][3],incline_list[[inc]][2])
						idx.minima=c(idx.minima,inc)
					}
					
					### check if there is at least one maxima
				  # difference at position inc should be +, and at position inc+1 should be -
					if (incline_list[[inc]][3] > 0 & incline_list[[inc+1]][3] < 0) {
						#printf("+turn inc[%s]=%0.5f inc[%s]=%0.5f x=%s",inc,incline_list[[inc]][3],inc+1,incline_list[[inc+1]][3],incline_list[[inc]][2])
						idx.maxima=c(idx.maxima,inc)
					}
					
					### check if there is a shoulder
					if ( abs(incline_list[[inc]][3]) < 0.01) {
						#printf("shoulder inc[%s]=%0.5f x=%s",inc,incline_list[[inc]][3],incline_list[[inc]][2])
						idx.shoulder=c(idx.shoulder,inc)
					}
					
				}
				this$inc=incline_x
				this$inc_list=incline_list
				
				idx=vector()
				### substract all maxima indices
				# remove maxima indices near minima
				minus.minima=intersect(idx.maxima,idx.minima)
				if (length(minus.minima) > 0) {
					for ( j in (minus.minima-3):(minus.minima+3)) {
						idx = c(idx,which(idx.minima==j))
					}
					idx.minima = idx.minima[-idx]
				}
				# remove shoulder indices near minima x+-0.3
				minus.shoulder=intersect(idx.maxima,idx.shoulder)
				if (length(minus.shoulder) > 0) {
					
					x.max=incline_x[minus.shoulder]
					for (k in idx.shoulder) {
						#print(incline_x[k])
						if (abs(x.max-incline_x[k])<=0.3) idx=c(idx,which(idx.shoulder==k))
						# if (abs(x.max-incline_x[k])<=0.3) printf("k=%s x[k]=%s x.max=%s abs=%s",k,incline_x[k],x.max,abs(x.max-incline_x[k]))
						#idx = c(idx, which(incline_x[k]x.max-incline_x[k])<=0.25))
					}
					#if (this$working) printf("w: idx.shoulder=%s",paste(idx.shoulder,collapse=" "))
					#if (this$working) printf("W. x.max=%s -idx=%s",x.max,paste(idx,collapse=" "))
					idx.shoulder = idx.shoulder[-idx]
				}

				#printf("length minima=%s",length(idx.minima))
				col.idx = colvec[i]
				if ( length(idx.minima)>0 ) {
				  # if there is a minima
					tmp=round(length(idx.minima)/2)
					if (tmp==0) tmp=idx.minima[1]
					else tmp=idx.minima[round(length(idx.minima)/2)]
					tclvalue(this$vcutoffs[[col.idx]]) = incline_list[[tmp]][2]
					cutoff = as.numeric(tclvalue(this$vcutoffs[[col.idx]]))
					if (this$working) printf("w: tmp=%s cutoff[%s]=%s (min)",tmp,vars[colvec[i]],cutoff)
				} else if ( length(idx.shoulder)>0 ) {
				  # if there is a shoulder
					tmp=round(length(idx.shoulder)/2)
					if (tmp==0) tmp=idx.shoulder[1]
					else tmp = idx.shoulder[round(length(idx.shoulder)/2)]
					if (this$working) printf("w: tmp=%s idx.shoulder=%s",tmp,paste(idx.shoulder,collapse=" "))
					tclvalue(this$vcutoffs[[col.idx]]) = incline_list[[tmp]][2]
					cutoff = as.numeric(tclvalue(this$vcutoffs[[col.idx]]))
					if (this$working) printf("w: cutoff[%s]=%s (shoulder)",vars[colvec[i]],cutoff)
				} else {
				  # if there is nothing found, than do 20% quantile
					tclvalue(this$vcutoffs[[col.idx]]) = this$calcCUTOFF(tdata,20,vars[i],col.idx)
					cutoff = as.numeric(tclvalue(this$vcutoffs[[col.idx]]))
					if (this$working) printf("w: cutoff[%s]=%s (max20)",vars[colvec[i]],cutoff)
				}

			}

			if (step == 1) {
				this$cutoff_list[[i]] = cutoff
			} else {
				this$cutoff_list[[i]][step] = cutoff
			}
		}
	}

	this$xminval = xminval
	this$xmaxval = xmaxval
	this$selected.cutoffs = vector() 
	for ( i in 1:len.colvec ) {
		tdata=data[,i]
		d=density(tdata)

		cutoff = round(median(this$cutoff_list[[i]]),1)
		col.idx = colvec[i]

		tclvalue(this$vcutoffs[[col.idx]]) = cutoff
		this$selected.cutoffs[colvec[i]] = cutoff
  
		### check for inifinity
		if (is.infinite(xminval[i]) | is.na(xminval[i])) xminval[i] = 0
		if (is.infinite(xmaxval[i]) | is.na(xminval[i])) xmaxval[i] = 8
		
		
		### plot histogram
		#plot(d,main=paste(colnames(data)[i]),xlim=c(xminval[i],xmaxval[i]),xlab="")
		#polygon(d, col="palegreen2", border="palegreen2")
		hist(tdata,freq=F,breaks=50,main=paste(colnames(data)[i]),xlab="",
			col = "palegreen2", lwd=0.5)
		lines(density(tdata), lwd = 2, col = "brown")

		# vertical lines
		for ( x in ceiling(xminval[i]) : floor(xmaxval[i]) ) {
			abline(v=x,col="darkgrey")
		}
		if ( cutoff != 0 ) {
			abline(v=cutoff,lty=2,lwd=2)
			tperc = round( length(tdata[which(tdata>=cutoff)]) / length(tdata)*100,1)
			# display percentage of production
			text(cutoff, max(d$y)-0.05*max(d$y),label=sprintf("%0.1f%%",tperc),cex=1,pos=4,xpd=TRUE)
			# display asinh-value of cutoff
			text(cutoff, min(d$y)+0.05*max(d$y),label=sprintf("[%s]",cutoff),cex=1,pos=4,xpd=TRUE)
		}
	}
	this$new.cutoffs[[this$selected.filenum]] = this$selected.cutoffs
	if (this$working) printf("w: current cutoffs:%s",paste(this$selected.cutoffs,collapse=" "))
}

fcs$interactGate <- function() {
	this=fcs
  printf("w: do interactGate")
	
	table = this$selected.project
	file=tclvalue(tkget(this$tkchoosefile))
	file.idx=this$current.filetable[which(this$current.filetable[,2]==file),1]
	this$selected.filenum = file.idx
	
	if ( is.null(this$data) | this$current.project != table | this$current.filenum != file.idx |
	     this$current.cofactor != as.numeric(tclvalue(this$rbasinh)) ) {
	  this$getFile(table,file.idx)
	}
	
	checkTRIMMING = tclvalue(this$cbttrimming)
	checkDOUBLETS = tclvalue(this$cbtremovedub)

	vars= this$selected.vars
	len = length(this$selected.vars)

	# column and cutoff number vector
	cutoffs=vector()
	colvec=vector()
	for (i in 1:len) {
		if ( tclvalue(this$cbVal[[i]]) == "1") {
			cutoffs=c(cutoffs,this$selected.cutoffs[i])
			colvec=c(colvec,i);
		}
	}
	len.colvec = length(colvec)

	### if trimming or remove doublets desired
	if ( checkDOUBLETS =="1" ) this$preprocData()
  	if ( checkTRIMMING =="1" ) this$preprocData(mode="trim")

	### select columns to plot
	data = this$data[colvec]
	colnames(data)=colnames(this$data[colvec])
	print(colnames(data))

	### open plot window
	dev.label="histogram.interactivegate"
	if ( length(which(dev.label==names(devList()))) != 0 ) {
		devSet(devList()[which(dev.label==names(devList()))])
	} else {
		devNew(type="x11",title="Histograms for cutoff setting",width=6,height=5,label=dev.label)
		par(oma=c(1.5,1,1.5,1),mar=c(2,2.4,1.5,1),mgp=c(1.3,0.4,0))
		this$plot.windows = c(this$plot.windows,dev.label)
	}

	for ( i in 1:len.colvec ) {
		tdata=data[,i]
		d=density(tdata)
		x=d$x
		y=d$y

		### index where y is greater than 1/35th of max(d$y)
		id_y = which((y>(1/35*max(y))))
		### index where x is positive
		id_x = which(x>0)

		xminval = x[intersect(id_y,id_x)[1]-5]
		xmaxval = round(max(tdata))

		### check for inifinity
		if (is.infinite(xminval) | is.na(xminval)) xminval = 0
		if (is.infinite(xmaxval) | is.na(xminval)) xmaxval = 8
		
		
		### plot histogram
		#plot(d,main=paste(colnames(data)[i]),xlim=c(xminval,xmaxval),xlab="")
		#polygon(d, col="palegreen2", border="palegreen2")
		hist(tdata,freq=F,breaks=50,main=paste(colnames(data)[i]),xlab="",
			col = "palegreen2", lwd=0.5)
		lines(density(tdata), lwd = 2, col = "brown")
		# vertical grid lines
		for ( x in ceiling(xminval) : floor(xmaxval) ) {
			abline(v=x,col="darkgrey")
		}
		# plot cutoffs if already entered
		if (cutoffs[i] != 0) {
			abline(v=cutoffs[i],lty=2,lwd=2)
			tperc = round( length(tdata[which(tdata>=cutoffs[i])]) / length(tdata)*100,1)
			# display percentage of production
			text(cutoffs[i], max(d$y)-0.05*max(d$y),label=sprintf("%0.1f%%",tperc),cex=1,pos=4,xpd=TRUE)
			# display asinh-value of cutoff
			text(cutoffs[i], min(d$y)+0.05*max(d$y),label=sprintf("[%s]",cutoffs[i]),cex=1,pos=4,xpd=TRUE)
		}

		loc=lapply(locator(1),round,1)
		abline(v=loc$x[1],col="darkgrey")
		printf("cutoff(%s) =%s", colnames(data)[i], loc$x[1])

		tclvalue(this$vcutoffs[[colvec[i]]]) = loc$x[1]
		this$selected.cutoffs[colvec[i]] = loc$x[1]
	}

	devOff(devList()[which(dev.label==names(devList()))])
	

	this$new.cutoffs[[this$selected.filenum]] = this$selected.cutoffs
	if (this$working) printf("w: current cutoffs:%s",paste(this$selected.cutoffs,collapse=" "))
}

fcs$normalizeFiles <- function() {
	this=fcs

	if(this$working) printf("Normalization start.")
}

fcs$calcCUTOFF <- function (data,cutoff,var,i) {
	this = fcs

	data = as.matrix(data)

	### calculate cutoff        
	data.sort=sort(data,decreasing=T)
	cutoff.new=round(data.sort[length(data)*cutoff/100],3)

	if (this$working) printf("w: do calcCUTOFF: cutoff.old=%s cutoff.new=%s var=%s idx=%s",cutoff,cutoff.new,var,i)

	### set cutoff in di- and triploT
	setcutoff =  eval(parse(text=paste0("this$vcutoffs[[",i,"]]")))
	tclvalue(setcutoff) = as.character(round(cutoff.new,3))
	setcbtperc = eval(parse(text=paste0("this$cbcutoffperc[[",i,"]]")))
	tclvalue(setcbtperc) = "0"

	### also set cutoff in overview tabs
	idx = which(this$selected.vars == var)
	setoverviewcutoff = eval(parse(text=paste("this$vcutoffs[[",idx,"]]",sep="")))
	tclvalue(setoverviewcutoff) = as.character(round(cutoff.new,3))

	return(cutoff.new)
}

fcs$getVariables <- function (staintable=NA,table=NA,index=NA) {
	this=fcs

	if (is.na(staintable)) staintable=this$current.staintable
	
	if (is.na(index)) {
		if (exists("current.file",envir=fcs)) {
			file=tclvalue(tkget(this$tkchoosefile))
		} else {
			file=this$current.filenames[1]
		}
		index=this$current.filetable[which(this$current.filetable[,2]==file),1]
	}
	if (this$working) printf("w: do getVariables from project #%s=%s: fileidx=%s",this$selected.projectnum,this$total.projects[this$selected.projectnum],index)

	vars=staintable[which(staintable[,1]==index),4]
	
	if (this$working) printf("w: vars(%s) old=%s",length(this$selected.vars),paste(this$selected.vars,collapse=" "))
	if (this$working) printf("w: vars(%s) new=%s",length(vars),paste(vars,collapse=" "))
	if (this$working) printf("w: Done getVariables")
  
	this$selected.vars = vars
	vars
}

fcs$getCutoffs <- function (fileindex=NA,new=FALSE) {
	this=fcs
	
	if (this$working) printf("w: do getCutoffs")

	if (exists("current.file",envir=fcs)) {
		file=tclvalue(tkget(this$tkchoosefile))
	} else {    
		file=this$current.filenames[1]
	}

	if (is.na(fileindex)) {
		fileindex=this$current.filetable[which(this$current.filetable[,2]==file),1]
	}

	if (new) {
		cutoffs=this$new.cutoffs[[fileindex]]
	} else {
		cutoffs = this$saved.cutoffs[[fileindex]]
	}

	cutoffs
}

fcs$getData <- function (table,fileidx, columns=NA, stain=NA, cofactor = NA) {
	this=fcs
	# get data from database
	data=dbGetQuery(this$conn,paste("SELECT * FROM ",table," WHERE file_ID == '",fileidx,"'",sep=""))
	
	# and ignore columns with NAs
	col.NA=NA
	for (i in 1:ncol(data)) {
		if (any(is.na(data[,i])))  col.NA=c(col.NA,i)
	}
	if (!is.na(col.NA)) data=data[,-col.NA]

	# change column names
	column.names=colnames(data)
	if (!is.na(columns)) {
		data=data[,columns+1]
	}

	if ( !is.na(stain) & !is.na(columns) ) { colnames(data) = stain[columns]}
	if ( is.na(stain) & !is.na(columns) ) { colnames(data) = this$selected.vars[columns] }
	if ( !is.na(stain) & is.na(columns) ) { colnames(data) = c("file_ID",stain) }
	if ( is.na(stain) & is.na(columns) ) { colnames(data) = c("file_ID",this$selected.vars) }
	
	# set asinh cofactor to 1 if not set
	cofactor = as.numeric(tclvalue(fcs$rbasinh))
	if ( is.na(cofactor) ) cofactor = 1
	
	if ( !is.na(columns) ) {
		data=asinh(data/cofactor)
	} else {
		data=asinh(data[,(2:dim(data)[2]/cofactor)])              
	}

	this$current.project = table 
	this$current.filenum = fileidx
	this$current.vars = this$selected.vars
	this$current.cutoffs = this$selected.cutoffs
	this$current.checks = this$selected.checks

	this$data.trimmed = FALSE
	
	printf("w: do getData(%s) from table='%s' with fileidx=%s and asinh cofactor=%s",nrow(data),table,fileidx,cofactor)
	
	data
}

fcs$getRectData <- function (table,file.idx,vars=NA) {
	this=fcs

	if (is.na(vars)) vars=c(this$checkMarker(tclvalue(tkget(this$cbvar1))),this$checkMarker(tclvalue(tkget(this$cbvar2))))

	checkGATED = tclvalue(this$cbtgateData)
	
	# if manual rect was selected
	if (tclvalue(this$cbtmanRect) == "1") {
		this$coords=list(
			x=c(as.double(tclvalue(this$vx1)),as.double(tclvalue(this$vx2))),
			y=c(as.double(tclvalue(this$vy1)),as.double(tclvalue(this$vy2)))
		)
	}
	
	x1 = this$coords$x[1]
	x2 = this$coords$x[2]
	y1 = this$coords$y[1]
	y2 = this$coords$y[2]
	
	#### new: gate data from this$data
	tmp.data=this$data[which(this$data[,vars[1]]>x1 & this$data[,vars[1]]<=x2),]
	tmp.data=tmp.data[which(tmp.data[,vars[2]]>y1 & tmp.data[,vars[2]]<=y2),]

	if ( checkGATED == "1" ) this$data=tmp.data
	else this$rectdata = tmp.data
	if (this$working) printf("w: do getRectData: dim(rectdata):%s from table='%s' with fileidx=%s",paste(dim(tmp.data),collapse="x"),table,file.idx)
	####

	this$current.project = table 
	this$current.filenum = file.idx
	this$current.vars = this$selected.vars
	this$current.cutoffs = this$selected.cutoffs
	this$current.checks = this$selected.checks
}

fcs$getFile <- function (table,fileidx,stain=NA,cofactor=NA) {
	this=fcs
	
	### checkbutton options
	checkTRANS = tclvalue(this$rbtrans)
	if (checkTRANS =="") checkTRANS = tclvalue(this$rbtrans) = "asinh"

	# get table from database
	data=dbGetQuery(this$conn,paste0("SELECT * FROM ",table," WHERE file_ID == '",fileidx,"'"))

	# cut column "file_ID" and ignore columns with NAs
  data$file_ID=NULL
	col.NA=NULL
	for (i in 1:ncol(data)) {
		if (any(is.na(data[,i]))) col.NA=c(col.NA,i)
	}
	if (!is.null(col.NA)) data=data[,-col.NA]

	this$current.project = table
	this$current.filenum = fileidx
	this$current.vars = this$selected.vars
	this$current.cutoffs = this$selected.cutoffs
	this$current.checks = this$selected.checks
	
	# set asinh cofactor to 1 if not set
	cofactor = as.numeric(tclvalue(fcs$rbasinh))
	if ( is.na(cofactor) ) cofactor = 1
	
	if ( checkTRANS == "asinh" ) {
		data = asinh(data/cofactor) 

		this$current.trans = "asinh"

		this$asinh$range = c(
			tclvalue(this$vminvalX),
			tclvalue(this$vmaxvalX),
			tclvalue(this$vminvalY),
			tclvalue(this$vmaxvalY)
			)
		this$asinh$binsize = as.numeric(tkget(this$binSize))
		this$asinh$mincount = as.numeric(tkget(this$minCountTri))

		this$data.trimmed = FALSE

		this$current.cofactor = cofactor
	}

	if (FALSE) { #} else {
		this$getIndexRange()

		data = as.data.frame(this$biex$call(data))

		this$current.trans = "biex"

		this$biex$range = c(
			tclvalue(this$vminvalX),
			tclvalue(this$vmaxvalX),
			tclvalue(this$vminvalY),
			tclvalue(this$vmaxvalY)
			)
		this$biex$binsize = as.numeric(tkget(this$binSize))
		this$biex$mincount = as.numeric(tkget(this$minCountTri))

		this$current.biex = list(
			x=c(tclvalue(this$vminvalX),tclvalue(this$vmaxvalX)),
			y=c(tclvalue(this$vminvalY),tclvalue(this$vmaxvalY)),
			binSize=as.numeric(tkget(this$binSize)),
			mincount=as.numeric(tkget(this$minCountTri))
		)
	}

	if ( !is.na(stain) ) { 
			colnames(data) = stain 
	} else { 
		colnames(data) = this$selected.vars 
	}

	if (this$working) print(head(data,2))
	printf("w: do getFile: %s from table='%s' with fileidx=%s and asinh cofactor=%s",nrow(data),table,fileidx,cofactor)

	this$data = data
	### set check button "Gate data" to zero
	# since you load a new file
	#this$cbtgateData=tclVar("0")
	#tkconfigure(this$cbtgateData,onvalue=FALSE,offvalue=TRUE)
}

fcs$getDFtable <- function (table) {
	this=fcs

	table.df=dbGetQuery(this$conn,paste("SELECT * FROM ", table))

	return(table.df)
}

fcs$addRectInfo <- function (setcex=1.0) {
	this=fcs
	xmin.FI=as.double(tclvalue(this$vminvalX))
	ymax.FI=as.double(tclvalue(this$vmaxvalY))

	rect(this$coords$x[1],this$coords$y[1],this$coords$x[2],this$coords$y[2])
	text(xmin.FI,(ymax.FI+1.0),sprintf("x=[%.1f,%.1f]; y=[%.1f,%.1f]",this$coords$x[1],this$coords$x[2],this$coords$y[1],this$coords$y[2]),cex=0.8*setcex,adj=0)
	this$plot.attr[[1]]$RectInfo=TRUE
}

fcs$addCellInfo <- function (setcex=1.0) {
	this=fcs
	xmin.FI=as.double(tclvalue(this$vminvalX))
	ymin.FI=as.double(tclvalue(this$vminvalY))

	this$ncell.perc=round(this$ncell.sel/this$origin.ncells*100,2)
	text(xmin.FI,(ymin.FI-0.5),paste("#all",this$origin.ncells,"#sel",this$ncell.sel,"%",this$ncell.perc),cex=0.8*setcex,adj=0)
	this$plot.attr[[1]]$CellInfo=TRUE
}

fcs$doRect <- function () {
	this=fcs

	this$rect.lwd=this$rect.lwd+0.5

	# if manual rect was not selected
	if (tclvalue(this$cbtmanRect) == "0") {

		dev.label=paste("plotter","tri",this$plotter.tri.num,sep=".")
		devSet(devList()[tail(which(dev.label==names(devList())),n=1)])
		
		this$coords=lapply(locator(2),round,1)

		# start from left bottom corner to top upper
		if (this$coords$x[2]<this$coords$x[1]) {
			tmp = this$coords$x[2]
			this$coords$x[2] = this$coords$x[1]
			this$coords$x[1] = tmp
		}
		if (this$coords$y[2]<this$coords$y[1]) {
			tmp = this$coords$y[2]
			this$coords$y[2] = this$coords$y[1]
			this$coords$y[1] = tmp
		}

		# set coords in manual rect entries
		tclvalue(this$vx1) = as.character(this$coords$x[1])
		tclvalue(this$vx2) = as.character(this$coords$x[2])
		tclvalue(this$vy1) = as.character(this$coords$y[1])
		tclvalue(this$vy2) = as.character(this$coords$y[2])

	}
	
	# plot rectangle
	rect(this$coords$x[1],this$coords$y[1], this$coords$x[2],this$coords$y[2],lwd=this$rect.lwd)

	#file=tclvalue(tkget(this$tkchoosefile))
	#this$getInfo(file)
	
	tkraise(this$tt)
}

fcs$getInfo <- function (file,vars) {
	# get information: vars, rectangle, num cells selected (%)
	this=fcs
	checkGATED = tclvalue(this$cbtgateData)
	
	if (this$working) printf("w: do getInfo file=%s GATED=%s temp.num=%s",file,checkGATED,this$temp.num)

	# file = tclvalue(tkget(this$tkchoosefile))
	# if gate, stay at temporary data
	if (( checkGATED == "1") & this$temp.num > 0 & grepl("^temp\\d+",file)) {
			table = file
			file.idx = 0
	} else {
			table = this$selected.project
			file.idx = this$current.filetable[which(this$current.filetable[,2]==file),1]
	}
	#if (this$working) printf("saved vars for gating: %s",vars)

	if ( this$current.project != table | this$current.filenum != file.idx ) this$getRectData(table,file.idx,vars)

	# refresh info about selected cells and percentage
	this$ncell.sel = nrow(this$data)
	this$ncell.perc = round(this$ncell.sel/this$origin.ncells*100,2)
	tkconfigure(this$ncell.sel.gui,text=as.character(this$ncell.sel))
	tkconfigure(this$ncell.perc.gui,text=as.character(this$ncell.perc))
	tkconfigure(this$ncell.sel.gui.di,text=as.character(this$ncell.sel))
	tkconfigure(this$ncell.perc.gui.di,text=as.character(this$ncell.perc))

	file=tclvalue(tkget(this$tkchoosefile))
	f=unlist(strsplit(file,"\\, |\\,|\\; |\\;|\\. | "))
	f=f[length(f)]

	# also refresh rect attr from "dataframes1/2"
	this$rect.attr = list(
		table = this$current.project,
		file.name = f,
		origin.ncells = this$origin.ncells,
		sel.ncells = this$ncell.sel,
		binSize = as.numeric(tclvalue(tkget(this$binSize))),
		mincount = as.numeric(tclvalue(tkget(this$minCountTri))),
		vars = vars,
		coords=c(this$coords$x[1],this$coords$x[2],this$coords$y[1],this$coords$y[2])
	)
	#this$refreshDFtree(this$treetable)
}

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

fcs$addTitle <- function(mode="") {
	this=fcs

	if ( length( dev.cur()==devList()[grep("di",names(devList()))])>0 )  set.cex=0.9
	else set.cex=1.3

	if (mode=="date") {
		date = gsub("-","",Sys.Date())
		title(main=date,outer=T,line=1,cex.main=set.cex,adj=1)
	} else {
		text=as.character(tclvalue(tkget(this$title,"1.0","end-1c")))
	
		if ( !is.null(dev.list()) ) {
			# get title withouth new line "\n" with "end-1c"
			title(main=text,outer=T,line=1,cex.main=set.cex)
		} else {
			tkmessageBox(title = "An error has occured!",
				message = "There is no plot window to print your title.", icon = "error", type = "ok")
			stop("There is no plot window to print your title.")
		}
	}
}

fcs$newPlotWindow <- function () {
	this=fcs

	this$plotter.tri.num = this$plotter.tri.num + 1
	nrow=as.numeric(tclvalue(this$vnrow))
	ncol=as.numeric(tclvalue(this$vncol))
	dev.label=paste("plotter","tri",this$plotter.tri.num,sep=".")

	#x11(width=ncol*4,height=nrow*4)
	devNew(type="x11",title="n-triploTs",width=ncol*3.4,height=nrow*3.7,label=dev.label)
	# mar in points, mai in inches
	# oma adds title lines
	# order: bottom, left, top, and right
	par(mfrow=c(nrow,ncol),oma=c(0.5,1,2,1),mar=c(3,3,4,2))
	#this$plotWindowExists = TRUE
	#this$dev.plot = dev.cur()
	#this$dev.list = dev.list()
	this$plot.windows = c(this$plot.windows,dev.label)
}

fcs$preprocData <- function(mode="removDub") {
	this=fcs

	table = this$selected.project
	file=tclvalue(tkget(this$tkchoosefile))
	file.idx=this$current.filetable[which(this$current.filetable[,2]==file),1]

	### if this$data doesn't exists or current file index is not selected file index
	if ( !exists("data",env=fcs) | is.null(this$data) | this$current.project != table | this$current.filenum!=file.idx |
	     this$current.cofactor != as.numeric(tclvalue(this$rbasinh)) ) {
		this$getFile(table,file.idx)
	}

	if ( this$working ) printf("w: do preprocData: project=%s fileindex=%s",table, file.idx)
	data = data.frame(this$data,check.names=T)

	if ( mode == "removDub" ) {
		### if "Gate data" is checked, break
		checkGATED = tclvalue(this$cbtgateData)
		if ( checkGATED == "1" ) {
			tkmessageBox(title = "An error has occured!",
				message = "Doublets removal only on ungated data. Please uncheck \"Gate data\".", icon = "error", type = "ok")
			stop("Please uncheck \"Gate data\".")
		}

		### if FSC-A and FSC-H are available
		## doublets can be calculated
		if ( length(data$FSC.A)>0 & length(data$FSC.H)>0 ) {
			fsc.median = median(data$FSC.A/data$FSC.H)
			fsc.stdev2 = 2*sd(data$FSC.A/data$FSC.H) 
			
			flag = mapply(this$flagDoublets,data$FSC.A,data$FSC.H,fsc.median,fsc.stdev2)
			
			data = data[-(which(flag==TRUE)),]
			
			if (this$working) {
				printf("Doublets removed: -%s cells", length(which(flag==TRUE)))
			} else {
				tkmessageBox(title = "Doublets removed.",
					message = sprintf("Doublets removed: -%s cells.", length(which(flag==TRUE))), icon = "info", type = "ok")
			}

			this$ncell.sel = nrow(data)
			this$ncell.perc=round(this$ncell.sel/this$origin.ncells*100,2)
			tkconfigure(this$ncell.sel.gui,text=as.character(this$ncell.sel))
			tkconfigure(this$ncell.perc.gui,text=as.character(this$ncell.perc))
		} else {
			print("No FSC-A and FSC-H to estimate doublets.")
		}
	} else if ( mode == "trim" ) {
	  if (this$data.trimmed){
	    print("Data already trimmed.")
	  } else{
	    printf("Original data dimension: %s ",paste(dim(data),collapse=" x "))
	    cells = 0
	    for (t in 1:ncol(data)) {
	      data = data[order(data[,t]),]
	      trim = floor(nrow(data)*(1-this$trim.num))
	      cells = cells + 2*(nrow(data)-trim) -1
	      data = data[(nrow(data)-trim):trim,]

	      # trim.idx = which(data[,t] > quantile(data[,t],c(1-this$trim.num)))
      	#   data = data[-trim.idx,]
	    }
	    printf("Data trimmed: -%s cells", cells)
	    this$data.trimmed = TRUE
	  }
		
	}

	this$data = data
	colnames(this$data) = this$selected.vars
}

fcs$flagDoublets <- function(fsc.a,fsc.h,median,sd2) {
	flag = ifelse ( (fsc.a/fsc.h) > (median+sd2), TRUE, FALSE)

	flag
}

fcs$calcPartition <- function(data) {
	this=fcs

	data = as.matrix(data)

	min=floor(10*min(data))/10
	max=ceiling(10*max(data))/10
	range=max-min
	partition.size=round(as.numeric(tkget(this$partition.size))/100*range,1)

	cnames=colnames(data)
	fX=cut(data[,1],breaks=seq(min,max,by=partition.size),include.lowest=TRUE)
	fY=cut(data[,2],breaks=seq(min,max,by=partition.size),include.lowest=TRUE)
	tab=table(fX,fY)
	colnames(tab)=seq(min,max-partition.size,by=partition.size)
	rownames(tab)=seq(min,max-partition.size,by=partition.size)

	# get means
	fXY=as.factor(paste(fX,fY))
	my.calc=aggregate(data[,3],by=list(fXY),mean)
	my.lengths=aggregate(data[,3],by=list(fXY),length)
	my.calc=cbind(my.calc,ncells=my.lengths$x)
	
	# binplot
	for (x in colnames(tab)) {
		for (y in rownames(tab)) {
			#if (tab[x,y]>mincells) {
			fact=as.factor(paste('(',x,',',as.numeric(x)+partition.size,'] ','(',y,',',as.numeric(y)+partition.size,']',sep=''))
			idx=which(as.character(fact)==as.character(my.calc$Group.1))
			rect(x,y,as.numeric(x)+partition.size,as.numeric(y)+partition.size,col=cols[my.calc[idx,'fac']],border=NA)
			
			this$bincount = this$bincount + 1

			#}
		}
	}
}

fcs$clickTab <- function(...) {
	this=fcs

	tab.old=this$current.tab
	tab.new=tclvalue(tkindex(this$panelTabs,"current"))

	if (this$working) printf("w: do clickTab tab.old=%s tab.new=%s",tab.old, tab.new)

	#if ( tab.new=="2" ) this$refreshTabHistory()
	#else this$refreshMarker(tab=TRUE)
	if ( tab.new!="2" ) this$refreshMarker(tab=TRUE)


	this$current.tab = tab.new
}

fcs$clearHistory <- function(){
	this=fcs
	
	imgTemplate=tclVar()
	tkimage.create("photo",imgTemplate,file=this$template.file)

	i = 0
	for (frame in strsplit(tclvalue(tkwinfo("children",this$frm))," ")[[1]]) {
		for (item in strsplit(tclvalue(tkwinfo("children",frame)), " ")[[1]]) {
			if(tclvalue(tcl(item, "cget", "-image"))!="") {
				tcl(item,"configure","-image",imgTemplate)
			}
		}
	}
	this$plot.num = 0
	this$images = list()
	this$plot.attr = list()
}

fcs$checkAxesRange <- function(){
	this=fcs

	xmin.val=as.numeric(tkget(this$minvalX))
	xmax.val=as.numeric(tkget(this$maxvalX))
	ymin.val=as.numeric(tkget(this$minvalY))
	ymax.val=as.numeric(tkget(this$maxvalY))

	if ( round(diff(c(xmin.val,xmax.val)),1) == round(diff(c(ymin.val,ymax.val)),1) ) {
	  return(TRUE)
	}	else {
	  return(FALSE)
	}
}

fcs$checkDigits <- function(num=NA,cutoff_id=NA) {
  this=fcs 
  	
	if ( !is.na(cutoff_id) ) {
		num=tclvalue(this$vcutoffs[[cutoff_id[1]]])
	}
	
	### only accepts numbers with points
	# deletes characters
	num = gsub("[a-z]","",num)
	# converts comma to point
	num = gsub(",",".",num)
	# deletes starting zeros
	num = gsub("^0+","",num)

	### if num is empty at the end, return 0
	if (num == "" | is.na(num)) num=0

	### changes value in cutoff array
	if ( !is.na(cutoff_id) ) tclvalue(this$vcutoffs[[cutoff_id[1]]]) = num

	num = round(as.numeric(num),1)
	
	### changes number in GUI
	if ( !is.na(cutoff_id) ) {
	this$selected.cutoffs[cutoff_id[1]] = num
	tclvalue(this$vcutoffs[[cutoff_id[1]]]) = as.character(num)
	}
	### returns number
	num
}

fcs$checkMarker <- function(var) {
	this=fcs

	printf("w: do checkMarker: var=%s",var)

	if ( !any(var==this$selected.vars) ) {

		# removes all punctuations for comparison
		selected.vars = gsub("[^[:alnum:]]","",toupper(this$selected.vars))
		var.tmp = gsub("[^[:alnum:]]","",toupper(var))
	
		var = this$selected.vars[which(selected.vars==var.tmp)]
		# get the var where name stops at the end
		if ( length(var)==0 ) var = this$selected.vars[grep(sprintf("%s$",var.tmp),selected.vars)]
		# if there is no, get first name which is similar
		if ( length(var)==0 ) var = this$selected.vars[grep(var.tmp,selected.vars)][1]

	# } else {
	  # grep(sprintf("\\<%s\\>",var),this$selected.vars)
		# var = this$selected.vars[grep(var,this$selected.vars)][1]
	}
	var.idx = which(this$selected.vars==var)

	if ( length(var)==0 ) print("w: Done checkMarker: NO feature found!")
	else printf("w: Done checkMarker: var #%s=%s",var.idx,var)
  
	
	var
}

fcs$checkDoubleMarkers <- function() {
	this=fcs

	len = length(this$selected.vars)

	printf("w: do checkDoubleMarkers: vars(%s)=%s",len,paste(this$selected.vars,collapse=" "))
	
	vars.unique=unique(this$selected.vars)
	len.unique = length(vars.unique)

	if (len != len.unique) {

		vars = make.unique(this$selected.vars)

		printf("w: Done checkDoubleMarkers: double feature names!")
	} else {
		printf("w: Done checkDoubleMarkers")
	}

	vars
}

fcs$clickinsertDF <- function(widget,...) {
		# mouse click prints dataframe to treetable
		this=fcs
		
		shade = c("none","gray")
		tkdelete(this$treetable,tkchildren(this$treetable,""))
		
		# get data frame
		env=...[1]
		df=eval(parse(text=paste("this$",env,sep="")))

		# insert df header
		tkconfigure(this$treetable,columns=colnames(df))
		for ( i in 0:(ncol(df)-1)){
			tcl(this$treetable,"heading",i,text=colnames(df)[i+1])
		}
		
		# insert df content
		for (i in 1:nrow(df)) {
			tkinsert(this$treetable,"","end",values=as.character(df[i,]),tag=shade[i%%2+1])
		}
		
		tktag.configure(this$treetable,"gray",background="gray95")
}

fcs$debug <- function() {
	this=fcs
	this$ttdebug=tktoplevel()
	
	# set window title
	tkwm.title(this$ttdebug,"Debug Functions")
	#set window size
	tkwm.geometry(this$ttdebug,"900x600")
	
	this$GUIdebug(this$ttdebug)
}

fcs$clickinsert <- function(widget,...) {
	# mouse click prints function to text field
	this=fcs
	
	child=...[1]
	parent=tkparent(this$treewidget,child)
	
	# first delete all text in text field
	tkdelete(this$txt,"1.0","end")
	
	# get function text
	if (as.character(parent) != "root") {

		# then print [child]$[function]=function
		tkinsert(this$txt,"end",paste(parent,"$",child,"=",sep=""))
		
		# get function
		geter=get(child,eval(parse(text=as.character(parent))))

		# now print function
		tkinsert(this$txt,"end", paste(deparse(geter), collapse="\n"), collapse="\n\n")
	} else {
		tkinsert(this$txt,"end","This is our parent environment.")
	}
}

fcs$printall <- function () {
	this=fcs
	
	# read script, separate each line
	script.header=scan("../YH_3D_plotter.r", character(0), sep ="\n",quiet=TRUE,nlines=53) 

	# get only header and tail
	script.header=noquote(script.header)
	if (this$working) printf("w: %s",script.header)
	
	# first delete all text in text field
	tkdelete(this$txt,"1.0","end")

	# insert header
	tkinsert(this$txt,"end",paste(script.header,collapse="\n"))
	tkinsert(this$txt,"end","\n\n")
		
	for (i in 1:length(this$parent.env)) {
		children=ls(eval(parse(text=paste(this$parent.env[i]))))

		for (j in 1:length(children)) {
			if (class(eval(parse(text=paste("this$",children[j],sep=""))))[1] == "function") {
				
				parent=tkparent(this$treewidget,children[j])
				
				# then print [children]$[function]=function
				tkinsert(this$txt,"end",paste(parent,"$",children[j],"=",sep=""))

				# get function
				geter=get(children[j],eval(parse(text=as.character(parent))))
				
				# now print function
				tkinsert(this$txt,"end", paste(deparse(geter), collapse="\n"), collapse="\n\n")
			}
		}
	}
}

fcs$refreshCode <- function(){
	this=fcs
	printf("w: do refreshCode :: destroy GUI and redo")
	# 1.0 means first row first character
	code <- tkget(this$txt, "1.0", "end")
	eval(parse(text=paste(as.character(tclvalue(code)),collapse="")))
	tkdestroy(this$tt)
	this$GUImain()
}

fcs$shortenFilename <- function(name,title=FALSE) {
	### shortens file name
	name = sub(".csv$","",name)
	name = sub(".fcs$","",name)
	name = sub("^temp\\d\\d_","",name)
	name = sub("^export_","",name)
	name = sub("^(\\w\\w)_*","\\1\\2",name)
	name = sub("Stain_*|stain_|Stain*|stain*","S\\1",name)
	name = sub("Panel|panel|Panel*|panel*","P\\1",name)
	name = sub("Patient*|patient*|Patient_*|patient_*","Pat\\1",name)
	name = sub("Surface*|surface*|Surface_*|surface_*","Surf\\1",name)
	name = sub("Stimulated|stimulated","Stim",name)
	name = sub("Gruppe|gruppe|group","Grp",name)
	name = sub("Day|day","d",name)
	name = sub("Antibodies|Antibody|antibodies|antibody","Ab",name)
	name = sub("Untreated|untreated|untreat","untr",name)
	name = sub("PMA_Iono|PMA_iono|PMAIono|PMAiono|pmaiono|pmaIono|pma_iono","P+I",name)
	name = sub("(\\dh)_*","\\1",name)
	
	name = sub("(\\d{4})_(\\d{2})_(\\d{2})","\\1\\2\\3",name,perl=TRUE)

	if (!title) {
		### if file name starts with 0-2 letters and follows with 6-8 digits
		# (short or long version of date) YY/MM/DD or YYYY/MM/DD)
		# then print only with last 3 digits
		name = sub("^\\w\\w\\d{3,5}(\\d)(\\d)(\\d)_","\\1\\2\\3_",name,perl=TRUE)
	}

	name
}

fcs$shortenMarkername <- function(markerNames) {
	### shortens marker names
	# take only first part of name if there is a "." in the name
	# and capitalize it
	unlist(lapply(markerNames,function(x) {
		toupper(strsplit(x,"[.]")[[1]][1])
		}))

	### NEW VERSION
	# allows dots in name
	unlist(lapply(markerNames,function(x) {
        len = length(strsplit(x,"[.]")[[1]])
        y = toupper(strsplit(x,"[.]")[[1]][1:(len-1)])
        paste(y, collapse=".")
        }))
}

fcs$saveWindow <- function (type="diploT") {
	this=fcs
	
	displayfile=tclvalue(tkget(this$tkchoosefile))
	ncol=as.numeric(tclvalue(this$vncol))
	nrow=as.numeric(tclvalue(this$vnrow))

	displayfile = this$shortenFilename(displayfile,title=TRUE)
	
	# if path didnt save yet
	if (!exists("lastpath",where=eval(parse(text="fcs")))){
		this$lastpath = getwd()
	}
  
	active.window = names(which(devList()==dev.cur()))
	var1.di = this$checkMarker(tclvalue(tkget(this$cbvar1di)))
	if (grepl("histogram", active.window)) {
	  plot.type = "histograms"
	  plot.title = sprintf("histograms of %s",displayfile)
	} else if (grepl("digraph", active.window)) {
	  plot.type = sprintf("diGraphs_%s",var1.di)
	  plot.title = sprintf("diGraphs of %s (%s)",var1.di,displayfile)
	} else if (grepl("di.", active.window)) {
	  plot.type = sprintf("diploTs_%s",var1.di)
	  plot.title = sprintf("diploTs of %s (%s)",var1.di,displayfile)
	} else {
	  plot.type = "triploTs" 
	  plot.title = sprintf("triploTs of %s",displayfile)
	}
	file = tkgetSaveFile(defaultextension=".pdf", initialdir= this$lastpath,
	                     initialfile=sprintf("%s_%s_%s.pdf",displayfile,plot.type,this$version),
	                     title = plot.title)
	printf("save file: %s",tclvalue(file))
	### remove file name to get file path
	filepath = unlist(strsplit(as.character(file),"/"))
	filepath = filepath[-length(filepath)]
	filepath = paste(filepath,collapse="/")
	this$lastpath = filepath
	
	### save pdf in original size
	if (type == "diploT") {
		dev.label=paste("plotter","di",this$plotter.di.num,sep=".")
		if (length(which(dev.label==names(devList())))!=0) devSet(devList()[which(dev.label==names(devList()))])
		dev.copy2pdf(file=tclvalue(file),
			title = plot.title)
		#devOff()
	} else {
		### triploT windows
		# does not matter how many rows/columns displayed
		dev.copy2pdf(file=tclvalue(file),
			pointsize=12*(ncol/4),width=3.66*ncol,height=4*nrow,
			title = plot.title)
	}

	print(paste("PDF saved:",file))
	
	this$lastpath = filepath
}

fcs$savePNG <- function (tdata=NULL,single=FALSE) {
	this=fcs

	if (this$working) printf("w do savePNG: single=%s", single)

	xminval = as.numeric(tkget(this$minvalX))
	xmaxval = as.numeric(tkget(this$maxvalX))
	yminval = as.numeric(tkget(this$minvalY))
	ymaxval = as.numeric(tkget(this$maxvalY))

	binSize = as.numeric(tkget(this$binSize))
	mincount = as.numeric(tkget(this$minCountTri))

	checkCALC = tclvalue(this$rbcalc)
	checkTRANS = tclvalue(this$rbtrans)
	checkGATED = tclvalue(this$cbtgateData)
	checkGRID = tclvalue(this$cbtshowGrid)

	prodcells.col = "black"
	quadrants.col = "black"
	
	# if path didnt save yet
	if (!exists("lastpath",where=eval(parse(text="fcs")))){
		this$lastpath = getwd()
	}

	if (single==TRUE) {

		displayfile = tclvalue(tkget(this$tkchoosefile))
		displayfile = this$shortenFilename(displayfile)

		v1=tclvalue(tkget(this$cbvar1))
		v2=tclvalue(tkget(this$cbvar2))
		v3=tclvalue(tkget(this$cbvar3))
		vars = c(v1,v2,v3)

		png.name = file.path(sprintf("%s_%s_%s_%s_%s.png",displayfile,paste(vars,collapse="_"),checkTRANS,checkCALC,this$version))
		file <- tkgetSaveFile(initialdir= this$lastpath, defaultextension=".png", initialfile=png.name)

		cutoff_idx = which(this$selected.vars==v1)
		cutoff_idy = which(this$selected.vars==v2)
		cutoff_idz = which(this$selected.vars==v3)
		cutoffs = c(cutoff_idx,cutoff_idy,cutoff_idz)

		tdata = as.matrix(this$data[vars])

		### if percentage is checked
		# calculate cutoffs and set check button to zero
		for ( i in 1:length(cutoffs)) {
			if ( tclvalue(this$cbcutoffperc[[cutoffs[i]]]) == "1" ) {
				cutoffs[i] = this$calcCUTOFF(tdata[vars[i]],this$checkDigits(cutoff_id=cutoffs[i]),vars[i],cutoffs[i])
			} else {
				cutoffs[i] = this$checkDigits(cutoff_id=cutoffs[i])
			}
		}
	} else {
		vars = colnames(tdata)

		png.folder = file.path(this$lastpath,"png")

		if ( grepl("linux",sessionInfo()$R.version$os) ) {
			system(paste("mkdir -p -v ",png.folder,sep=""))
		} else {
			system(paste("mkdir ",png.folder,sep=""))
		}
		
		if ( this$OverviewGate ) {
			displayfile = this$plot.attr[[1]]$file.name
		} else {
			displayfile=tclvalue(tkget(this$tkchoosefile))
		}
		displayfile = this$shortenFilename(displayfile)

		file = file.path(png.folder,sprintf("%s_%s_%s_%s_%s.png",displayfile,paste(vars,collapse="_"),checkTRANS,checkCALC,this$version))
	}

	### remove file name to get file path
	filepath = unlist(strsplit(as.character(file),"/"))
	filename = filepath[length(filepath)]
	filename = sub(".png$","",filename)

	filepath = filepath[-length(filepath)]
	filepath = paste(filepath,collapse="/")
	this$lastpath = filepath
	
	### png plot
	toPNG(name=filename,path=filepath,width=800,height=850,bg="white", {
		par(oma=c(2,2,2,1),mar=c(4.5,4.5,5,1))
		#png.cur=dev.cur()

		#### plot start
		set.cex=1.4
		set.cex.axes=1.4
		set.mgp=c(1.9, 0.6, 0)
		if ( cutoffs[1] > 0 ) title.axis = sprintf("%s (%s)",vars[1],cutoffs[1])
		else title.axis = vars[1]
		if ( cutoffs[2] > 0 ) title.axis = c(title.axis,sprintf("%s (%s)",vars[2],cutoffs[2]))
		else title.axis = c(title.axis,vars[2])

		plot(1,type='n',frame.plot=FALSE,xlim=c(xminval,xmaxval+10*binSize),axes=FALSE,
				ylim=c(yminval-2.5*binSize,ymaxval+5*binSize),xlab=title.axis[1],ylab=title.axis[2],cex.lab=set.cex,cex.axis=0.5*set.cex.axes,mgp=set.mgp)
		box(lwd=0.8,col="darkgrey")

		if ( checkTRANS == "asinh") {
			scale = this$asinh$scale
			label = this$asinh$label
			grid.step = this$asinh$step
		} else {
			scale = this$biex$scale
			label = this$biex$label
			grid.step = this$biex$step
		} 

		### draw axis on the bottom and on the left
		axis(side=1, at=scale,labels=label,las=1,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")
		axis(side=2, at=scale,labels=label,las=3,cex.axis=set.cex.axes,mgp=set.mgp,col="darkgrey")

		### add grid
		if(checkGRID=="1") {
			xgrid.steps=seq(0,(xmaxval),by=grid.step)
			ygrid.steps=seq(0,(ymaxval),by=grid.step)
			abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
		}
		

		### calc quadrants in total
		ncells = ncells.total = nrow(tdata)
		# q1 Quadrant unten links
		# q2 Quadrant unten rechts
		# q3 Quadrant oben rechts
		# q4 Quadrant oben links
		if ( cutoffs[1] > 0 & cutoffs[2] > 0 ) {

			### count cells in quadrant
			tdata.q1 = tdata[which( tdata[,1]<cutoffs[1] &  tdata[,2]<cutoffs[2] ),3]
			tdata.q2 = tdata[which( tdata[,1]>=cutoffs[1] &  tdata[,2]<cutoffs[2] ),3]
			tdata.q3 = tdata[which( tdata[,1]>=cutoffs[1] &  tdata[,2]>=cutoffs[2] ),3]
			tdata.q4 = tdata[which( tdata[,1]<cutoffs[1] &  tdata[,2]>=cutoffs[2] ),3]

			### q[x].total [ink=black]
			### percentage of cells in quadrant to total cells 
			### or in MSI(+): percentage of cells in quadrant to total positive cells
			this$q1.total = abs(100 * length( tdata.q1 ) / ncells)
			this$q2.total = abs(100 * length( tdata.q2 ) / ncells)
			this$q3.total = abs(100 * length( tdata.q3 ) / ncells)
			this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
		
			if (cutoffs[3] > 0 ) {
				### number of cells which are producing cells in feature C
				ncells = nrow(tdata[which(tdata[,3]> cutoffs[3]),])

				### q[x].prodcells [ink=red]
				### percentage of cells which are positive for feature C in quadrant to total quadrant cells
				this$q1.prodcells = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[3])]) /length(tdata.q1)
				if (is.nan(this$q1.prodcells)) this$q1.prodcells = 0
				this$q2.prodcells = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[3])]) /length(tdata.q2)
				if (is.nan(this$q2.prodcells)) this$q2.prodcells = 0
				this$q3.prodcells = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[3])]) /length(tdata.q3)
				if (is.nan(this$q3.prodcells)) this$q3.prodcells = 0
				this$q4.prodcells = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[3])]) /length(tdata.q4)
				if (is.nan(this$q4.prodcells)) this$q4.prodcells = 0

				### only do MSI plots on producing cells only
				if ( checkCALC == "MSI(+)" ) {
					### cut all cells which are not producing cells
					tdata.plus = tdata[which(tdata[,3]> cutoffs[3]),]
					ncells = nrow(tdata.plus)

					tdata.q1 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
					tdata.q2 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ),3]
					tdata.q3 = tdata.plus[which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]
					tdata.q4 = tdata.plus[which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ),3]

					### q[x].total [ink=blue]
					### in MSI(+): percentage of cells in quadrant to total positive cells
					this$q1.total = abs(100 * length( which( tdata.plus[,1]<cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
					this$q2.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]<cutoffs[2] ) ) / ncells)
					this$q3.total = abs(100 * length( which( tdata.plus[,1]>=cutoffs[1] &  tdata.plus[,2]>=cutoffs[2] ) ) / ncells)
					this$q4.total = abs(100 - this$q1.total - this$q2.total - this$q3.total)
				}

				### q[x].prodcellsplus [ink=green]
				### percentage of cells which are positive for feature C to total cells
				this$q1.prodcellsplus = 100 * length(tdata.q1[which(tdata.q1>=cutoffs[3])]) / ncells.total
				if (is.nan(this$q1.prodcellsplus)) this$q1.prodcellsplus = 0
				this$q2.prodcellsplus = 100 * length(tdata.q2[which(tdata.q2>=cutoffs[3])]) / ncells.total
				if (is.nan(this$q2.prodcellsplus)) this$q2.prodcellsplus = 0
				this$q3.prodcellsplus = 100 * length(tdata.q3[which(tdata.q3>=cutoffs[3])]) / ncells.total
				if (is.nan(this$q3.prodcellsplus)) this$q3.prodcellsplus = 0
				this$q4.prodcellsplus = 100 * length(tdata.q4[which(tdata.q4>=cutoffs[3])]) / ncells.total
				if (is.nan(this$q4.prodcellsplus)) this$q4.prodcellsplus = 0
			}
		}

		tdata.zero = tdata[apply(tdata[, -1], MARGIN = 1, function(row) {all(row > 0)}), ]
		ncells.zero=nrow(tdata.zero)

		### calc quadrants with only positive values
		# q1 Quadrant unten links
		# q2 Quadrant unten rechts
		# q3 Quadrant oben rechts
		# q4 Quadrant oben links
		if ( cutoffs[1] > 0 & cutoffs[2] > 0 ) {
			q1.zero = abs(100 * length( which( tdata.zero[,1]<cutoffs[1] &  tdata.zero[,2]<cutoffs[2] ) ) / ncells.zero)      
			q2.zero = abs(100 * length( which( tdata.zero[,1]>=cutoffs[1] &  tdata.zero[,2]<cutoffs[2] ) ) / ncells.zero)                    
			q3.zero = abs(100 * length( which( tdata.zero[,1]>=cutoffs[1] &  tdata.zero[,2]>=cutoffs[2] ) ) / ncells.zero)                    
			q4.zero = abs(100 - q1.zero - q2.zero - q3.zero)
		}

		#this$bintriplot(tdata,cutoffs,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp,binSize=binSize,mincells=mincount, quadrants.color = quadrants.col)
		if (checkCALC == "density") {
			this$bintriplot(data=tdata,cutoffs=cutoffs,density=TRUE,binSize=binSize,mincells=mincount
				,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
		} else if (checkCALC == "MSI(+)") {
			this$bintriplot(data=tdata.plus,cutoffs=cutoffs,binSize=binSize,mincells=mincount,
				,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp,quadrants.color = "blue", data.origin=tdata)
		} else {
			this$bintriplot(data=tdata,cutoffs=cutoffs,binSize=binSize,mincells=mincount
				,set.cex=label.cex,set.cex.axes=set.cex.axes,set.mgp=set.mgp)
		}

		if (checkGATED == "1") firstLine=sprintf("%s(%0.1f%%): %s/cof=%s",displayfile,this$ncell.perc,checkCALC,this$current.cofactor)
		else firstLine=sprintf("%s: %s/cof=%s",displayfile,checkCALC,this$current.cofactor)
		title(main=firstLine,line=3.2,cex.main=1.4,adj=0)

		#if (checkCALC == "freq" | checkCALC == "MSI(+)") secondLine=sprintf("cells(min/max)=%s/%s; cutoff=%s",mincount,this$maxcells, cutoffs[3])
		if (checkCALC == "freq" | grepl("MSI",checkCALC)) secondLine=sprintf("cells(min/max)=%s/%s;%s",mincount,this$maxcells,v3)
		else secondLine=sprintf("cells(min/max)=%s/%s",mincount,this$maxcells)
		title(main=secondLine,line=2.3,cex.main=1,adj=0)
												
		thirdLine=sprintf("%s-%s(%0.1f%%); binSize=%s,#bins=%s",ncells,ncells.zero,(ncells.zero/ncells*100),binSize,this$bincount)
		title(main=thirdLine,line=1.5,cex.main=1,adj=0)


								
		if (checkGATED == "1" | grepl("temp",displayfile)) {
			if ( length(this$coords.info)>2 ) {
				fourthLine=sprintf("%s; %s",this$coords.info[1],this$coords.info[2])
				fifthLine=sprintf("%s; %s",this$coords.info[3],this$coords.info[4])
			} else {
				fourthLine=sprintf("%s",paste(this$coords.info,collapse="; "))
				fifthLine=""
			}
			title(main=fourthLine,line=0.8,cex.main=0.9,adj=0)
			title(main=fifthLine,line=0.1,cex.main=0.9,adj=0)
		}

	}
	)
	#dev.off(png.cur)

	if (single) {
		print(paste("File saved:",file))
		this$lastpng = file
	} 
}

fcs$saveDF <- function(n){
	this=fcs
	
	attr = this$rect.attr

	# if path didnt save yet
	if (!exists("lastpath",where=eval(parse(text="fcs")))){
		this$lastpath = getwd()
	}

	file = tkgetSaveFile(defaultextension=".csv", initialdir= this$lastpath, initialfile=attr$file.name)

	### remove file name to get file path
	filepath = unlist(strsplit(as.character(file),"/"))
	filepath = filepath[-length(filepath)]
	filepath = paste(filepath,collapse="/")
	this$lastpath = filepath
	

	comments=vector(length=)
	comments[1] = paste("# project:",attr$table)
	comments[2] = paste("# file name:",attr$file.name)
	comments[3] = paste("# selected columns:",attr$vars)
	comments[4] = sprintf("# rectangle at: x=[%.1f,%.1f]; y=[%.1f,%.1f]",attr$coords[1],attr$coords[2],attr$coords[3],attr$coords[4])
	comments[5] = paste("# bin size:",attr$binSize)
	comments[6] = paste("# min(cells) per bin:",attr$mincout)
	comments[7] = paste("# original cell count:",attr$origin.ncells)
	comments[8] = paste("# selected cell count:",attr$sel.ncells)
	comments[9] = paste("# selected cell count(%):",round((attr$sel.ncells/attr$origin.ncells*100),2))

	# insert df content
	if(!tclvalue(file)=="") {
		chn = tcl("open", file, "w")

		# insert comments
		for (i in 1:length(comments)){
			tcl("puts",chn,comments[i])
			print(comments[i])
		}
		
		# insert dataframe
		df <- eval(parse(text=paste("this$dataframes",n, sep="")))
		tmp = colnames(df)[1]

		for ( i in 2:ncol(df)){
			tmp = paste(tmp,colnames(df)[i],sep=",")
		}
		print(tmp)
		tcl("puts", chn, tmp)

		for (j in 1:nrow(df)) {
			tmp =df[j,1]
			for (k in 2:ncol(df)){
				tmp =  paste(tmp,df[j,k],sep=",")
			}
			tcl("puts", chn,tmp)
		}
		
		tcl("close", chn)
		print(paste("Table Saved: ",tclvalue(file)))
	
	}
}

fcs$saveCode <- function(){
	this=fcs
	
	# if didn't save yet
	#if (!exists("lastfile",where=eval(parse(text=tkparent(this$treewidget,"lastfile"))))){
	#if (!exists("lastfile",where=eval(parse(text="fcs")))){
	#	this$lastfile = paste(getwd(),"/3D_plotter_",this$version,".r",sep="")
	#}

	# if path didnt save yet
	if (!exists("lastpath",where=eval(parse(text="fcs")))){
		this$lastpath = getwd()
	}
			
	file = tkgetSaveFile(defaultextension=".r", initialdir= this$lastpath, initialfile=paste("/3D_plotter_",this$version,".r",sep=""))
	
	### remove file name to get file path
	filepath = unlist(strsplit(as.character(file),"/"))
	filepath = filepath[-length(filepath)]
	filepath = paste(filepath,collapse="/")
	this$lastpath = filepath
	

	if(!tclvalue(file)=="") {
		chn = tcl("open", file, "w")
		tcl("puts", chn, tkget(this$txt,"0.0","end"))
		tcl("close", chn)
	}
}

fcs$saveCutoffsToDB <- function() {
	this=fcs

	if (this$working) {
		printf("w: saveCutoffsToDB: cutoffs saved for project=%s",this$selected.project)
		printf("w: %s",paste(unlist(this$new.cutoffs),collapse=" "))
	}
	new.staintable=cbind(this$current.staintable[,1:4],unlist(this$new.cutoffs[!is.na(this$new.cutoffs)]),unlist(this$new.checks[!is.na(this$new.checks)]))
	colnames(new.staintable)[5] = "file_savedCutoffs"
	colnames(new.staintable)[6] = "file_savedChecks"
	if (this$working) print(head(new.staintable))
	#idx.stain=grep(paste0(this$current.project,this$markerid_name),dbListTables(this$conn))
	#idx.stain=which(dbListTables(this$conn)==paste0(this$current.project,this$markerid_name))

	### get name of stain table
	name.stain=paste0(this$selected.project,this$markerid_name)
	### remove and rewrite stain table
	dbRemoveTable(this$conn,name.stain)
	dbWriteTable(this$conn,name.stain,new.staintable)

	printf("Cutoffs saved in %s",name.stain)
}

fcs$saveCutoffsToAll <- function() {
  this=fcs
  
  answer = tkmessageBox(title = "Saving cutoffs and checks to ALL files in project", message = "Are you sure?", icon = "info", type = "yesno")
  if (tclvalue(answer) == "yes") {
    this$selected.checks = vector()
    this$selected.cutoffs = vector()
    for ( i in 1:length(this$selected.vars)) {
      this$selected.checks[i] = tclvalue(this$cbVal[[i]])
      this$selected.cutoffs[i] = as.numeric(tclvalue(this$vcutoffs[[i]]))
    }
    
    if (this$working) {
      printf("w: saveCutoffsToall: cutoffs/checks saved for ALL files in project=%s",this$selected.project)
      printf("w: %s",paste(unlist(this$selected.cutoffs),collapse=" "))
    }
    # this$current.staintable = this$getDFtable(paste0(this$selected.project,this$markerid_name))
    # this$selected.filenum = this$current.filetable[which(this$current.filetable[,2]==this$current.filenames[1]),1]
    
    len.staintable = length(this$selected.cutoffs)*length(this$current.filenames)
    new.staintable=cbind(this$current.staintable[1:len.staintable,1:4]
                         ,rep(unlist(this$selected.cutoffs),length(this$current.filenames))
                         ,rep(unlist(this$selected.checks),length(this$current.filenames))
                    )
    
    colnames(new.staintable)[5] = "file_savedCutoffs"
    colnames(new.staintable)[6] = "file_savedChecks"
    if (this$working) print(head(new.staintable))
    #idx.stain=grep(paste0(this$current.project,this$markerid_name),dbListTables(this$conn))
    #idx.stain=which(dbListTables(this$conn)==paste0(this$current.project,this$markerid_name))
    
    ### get name of stain table
    name.stain=paste0(this$selected.project,this$markerid_name)
    ### remove and rewrite stain table
    dbRemoveTable(this$conn,name.stain)
    dbWriteTable(this$conn,name.stain,new.staintable)
    
    printf("Cutoffs saved in %s",name.stain)
    
    
  }
}

fcs$changeFI.range <- function (mode) {
	this=fcs

	if ( mode == 1 ) {
		tclvalue(this$rbtrans) = "asinh"
		tclvalue(this$vminvalX) = this$asinh$range[1]
		tclvalue(this$vmaxvalX) = this$asinh$range[2]
		tclvalue(this$vminvalY) = this$asinh$range[3]
		tclvalue(this$vmaxvalY) = this$asinh$range[4]

		## change in triploT
		tkset(this$binSize,this$asinh$binsize)
		tkset(this$minCountTri,this$asinh$mincount)

		## change in tetraploT
		#tkset(this$binSizetetra,this$current.asinh$binSize)
		#tkset(this$minCountTritetra,this$current.asinh$mincount)

		## change in triploT Overview
		#tkset(this$binSizeOverview,this$current.asinh$binSize)
		#tkset(this$minCountTriOverview,this$current.asinh$mincount)

		## change in diploT
		tkset(this$binSizedi,this$asinh$binsize)
		tkset(this$minCountDi,this$asinh$mincount)
	} else if (mode == 2) {
		tclvalue(this$rbtrans) = "biex"
		tclvalue(this$vminvalX) = this$asinh$range[1]
		tclvalue(this$vmaxvalX) = this$asinh$range[2]
		tclvalue(this$vminvalY) = this$asinh$range[3]
		tclvalue(this$vmaxvalY) = this$asinh$range[4]

		## change in triploT
		tkset(this$binSize,this$asinh$binsize)
		tkset(this$minCountTri,this$asinh$mincount)
		tkconfigure(this$binSize,values=this$binSizes.biex)
		tkset(this$binSize,this$current.biex$binSize)
		tkset(this$minCountTri,this$current.biex$mincount)

		## change in tetraploT
		#tkconfigure(this$binSizetetra,values=this$binSizes.biex)
		#tkset(this$binSizetetra,this$current.biex$binSize)
		#tkset(this$minCountTritetra,this$current.biex$mincount)

		## change in triploT Overview
		#tkconfigure(this$binSizeOverview,values=this$binSizes.biex)
		#tkset(this$binSizeOverview,this$current.biex$binSize)
		#tkset(this$minCountTriOverview,this$current.biex$mincount)

		## change in diploT
		tkconfigure(this$binSizedi,values=this$binSizes.biex)
		tkset(this$binSizedi,this$current.biex$binSize)
		tkset(this$minCountDi,this$current.biex$mincount)
	}
}

fcs$connectDb<- function(fname) {
	this=fcs
	this$conn <- dbConnect(SQLite(), dbname = fname)
	print(paste("Database opened:",fname))
}

fcs$disconnectDb <- function () {
	this=fcs
	dbDisconnect(this$conn)
	print(paste("Database closed:",file.path(this$db.path,this$db.name)))
}

fcs$mytkopentree <- function (widget,node,...) {

	tkitemconfigure(widget,node,image=imageBookOpen)
}

fcs$mytkclosetree <- function (widget,node,...) { 

	tkitemconfigure(widget,node,image=imageBook)
}

fcs$start <- function () {
	this=fcs
	if (file.exists(file.path(fcs$db.path,fcs$db.name))) {
		fcs$connectDb(file.path(fcs$db.path,fcs$db.name))
	} else {
		file <- tclvalue(tkgetOpenFile(initialdir=fcs$db.path,defaultextension="sqlite3"))
		if (file.exists(file)) {
			this$connectDb(file)
			this$db.name=file
		}
	}
	this$plotWindowExists=FALSE

	this$GUImain()
}

fcs$exit <- function() {
	this=fcs
	continue = FALSE
	saveCutoffs = FALSE
	answer=tclVar("no")

	### save current cutoffs
	this$refreshCutoffs(current=TRUE)

	### test if cutoffs had been changed
	new = unlist(this$new.cutoffs)
	saved = unlist(this$saved.cutoffs)
	both <- new[!is.na(new)] == saved[!is.na(saved)]

	new2 = unlist(this$new.checks)
	saved2 = unlist(this$saved.checks)
	both2 <- new2[!is.na(new2)] == saved2[!is.na(saved2)]

	### if cutoffs did not change
	# if length(new) == length(saved)
	# if both == TRUE and cutoffs were already saved in DB
	#if ( length(new) == length(saved) ) {
	#if ( ((all(both) & length(saved)>0) | all(new==0))  & ((all(both2) & length(saved2)>0) | all(new2==0) ) )  {
	if ( ((all(both) & length(saved)>0) & length(both)>0) & ((all(both2) & length(saved2)>0) & length(both2)>0) )  {
		### if cutoffs did not change
		if (this$working==FALSE) {
			answer = tkmessageBox(title = "Quitting..", message = "Are you sure?", icon = "info", type = "yesno")
			if ( tclvalue(answer) == "yes" ) continue = TRUE
		}
	} else {
		### if cutoffs have been changed
		answer=tkmessageBox(title = "Quitting..", message = "Save cutoffs?", icon = "info", type = "yesnocancel")
		if ( tclvalue(answer) == "yes" | tclvalue(answer) == "no" ) continue = TRUE
		if ( tclvalue(answer) == "yes" ) saveCutoffs = TRUE
	}


	if (this$working) {continue=TRUE; tclvalue(answer)="yes"; printf("w: do exit saveCutoffs=%s continue=%s answer=%s",saveCutoffs,continue,tclvalue(answer))}

	if ( continue ) {

		### save workspace
		# feature diplot
		#tab=tclvalue(tkindex(this$panelTabs,"current"))
		var1 = this$checkMarker(tclvalue(tkget(this$cbvar1di)))
		param$currVarDi1 = which(this$selected.vars==var1)
		var2 = this$checkMarker(tclvalue(tkget(this$cbvar2di)))
		param$currVarDi2 = which(this$selected.vars==var2)
		printf("diploT :: var1=%s var2=%s",param$currVarDi1,param$currVarDi2)
		mincountDi=as.numeric(tclvalue(tkget(this$minCountDi)))
		param$minCountDiPos = which(this$min.counts.diploT==mincountDi)
		binsizeDi = as.numeric(tclvalue(tkget(this$binSizedi)))
		param$binSizesDiPos = which(this$binSizes==binsizeDi)
		# feature triplot
		var1 = this$checkMarker(tclvalue(tkget(this$cbvar1)))
		param$currVarTri1 = which(this$selected.vars==var1)
		var2 = this$checkMarker(tclvalue(tkget(this$cbvar2)))
		param$currVarTri2 = which(this$selected.vars==var2)
		var3 = this$checkMarker(tclvalue(tkget(this$cbvar3)))
		param$currVarTri3 = which(this$selected.vars==var3)
		printf("triploT :: var1=%s var2=%s var3=%s",param$currVarTri1,param$currVarTri2,param$currVarTri3)
		mincount=as.numeric(tclvalue(tkget(this$minCountTri)))
		param$minCountTriPos = which(this$min.counts==mincount)
		binsize = as.numeric(tclvalue(tkget(this$binSize)))
		param$binSizesTriPos = which(this$binSizes==binsize)
		# general
		param$cofactor = as.character(this$rbasinh)
		# FI range
		minX = tclvalue(this$vminvalX)
		param$minvalX = minX
		maxX = tclvalue(this$vmaxvalX)
		param$maxvalX = maxX
		minY = tclvalue(this$vminvalY)
		param$minvalY = minY
		maxY = tclvalue(this$vmaxvalY)
		param$maxvalY= maxY
		minMSI = tclvalue(this$vminMSI)
		param$minMSI = minMSI
		maxMSI = tclvalue(this$vmaxMSI)
		param$maxMSI = maxMSI
		# radio and check buttons
		checkCALC = tclvalue(this$rbcalc)
		param$rbcalc = checkCALC
		checkTRANS = tclvalue(this$rbtrans)
		param$rbtrans = checkTRANS
		checkfeatA = tclvalue(this$cbtfeatA)
		param$cbtfeatA = checkfeatA
		checkfeatB = tclvalue(this$cbtfeatB)
		param$cbtfeatB = checkfeatB
		checkfeatC = tclvalue(this$cbtfeatC)
		param$cbtfeatC = checkfeatC
		checkADDFILE = tclvalue(this$cbtaddFilename)
		param$cbtaddFilename = checkADDFILE
		checkADDDATE = tclvalue(this$cbtaddDate)
		param$cbtaddDate = checkADDDATE
		checkDISPLAYA = tclvalue(this$cbtdisplayA)
		param$cbtdisplayA = checkDISPLAYA
		checkGRID = tclvalue(this$cbtshowGrid)
		param$cbtshowGrid = checkGRID
		checkPERCENTAGE = tclvalue(this$cbtshowPercentage)
		param$cbtshowPercentage = checkPERCENTAGE
		param$cbttrimming = tclvalue(this$cbttrimming)
		checkMANRECT = tclvalue(this$cbtmanRect)
		param$cbtmanRect = checkMANRECT
		checkDYNRANGE = tclvalue(this$cbtdynRange)
		param$cbtdynRange = checkDYNRANGE
		checkLEGEND = tclvalue(this$cbtshowLegend)
		param$cbtshowLegend = checkLEGEND
		checkSHOWMINBIN = tclvalue(this$cbtshowMinBins)
		param$cbtshowMinBins = checkSHOWMINBIN
		checkGATEDATA = tclvalue(this$cbtgateData)
		param$cbtgateData = checkGATEDATA
		checkAUTORECT = tclvalue(this$cbtautoRect)
		param$cbtautoRect = checkAUTORECT
		# graphics
		graphROW = tclvalue(this$vnrow)
		param$nrow = graphROW
		graphCOL = tclvalue(this$vncol)
		param$ncol = graphCOL
		save(param, file = "myPRIparam.rda")
		###

		tkdestroy(this$tt)
		if ( class(this$pb) != "NULL" ) close(this$pb)

		if (exists("ttdebug",where=eval(parse(text="fcs")))) {tkdestroy(this$ttdebug)}
		if (exists("ttpreproc",where=eval(parse(text="fcs")))) {tkdestroy(this$ttpreproc)}
		if (this$ttprojectsopen) {tkdestroy(this$ttprojects)}
		graphics.off()

		if(file.exists(this$png.file)) system(paste("rm ",this$png.file,sep=""))

		temp.table.idx = grep("^temp\\d+",dbListTables(this$conn))
		len = length(temp.table.idx)
		if ( len > 0 ) {
			for ( i in len:1 ) {
				table.name = dbListTables(this$conn)[temp.table.idx[i]]
				dbRemoveTable(this$conn,table.name)
				print(paste0("Removing temporary table '",table.name,"'"))
			}
		}

		### save cutoffs into MarkerInfo table
		if ( saveCutoffs ) this$saveCutoffsToDB()
		this$disconnectDb()
	}
}

fcs$exitandstart <- function () {
	this=fcs

	this$exit()
	source("YH_3D_plotter.r")
	load(file = "myPRIparam.rda")
}

if (length(strsplit(fcs$db.name,"")[[1]]) > 0) {
	fcs$connectDb(file.path(fcs$db.path,fcs$db.name))
} else {
	file <- tclvalue(tkgetOpenFile(initialdir=fcs$db.path,defaultextension="sqlite3"))
	fcs$connectDb(file)
	fcs$db.name=file
}

fcs$GUImain()
