#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019
# Goood2018

### Basal representatives
# UPN12/15 low risk
# UPN10/22 high ristk

rm(list = ls())



# initiate ------------------------------------------------
cofactor = 0.2
trimming = FALSE
project.name = "BCR"
project.name.long ="BCR-Crosslink"
subgroup = "Validation"
subgroup = "Training"
coverage = "full"
### set paths
folder.path = "/scratch/drfz/Good2018/glmnet"
setwd(folder.path)


# load libraries ------------------------------------------
# xlsx            :   Library for Excel reading and creating
# RSQLite         :   interact with database
# dplyr           :   faster binding of columns bind_cols() and rows bind_rows()
# reshape2        :   melt()
# gpplot2         :   for plotting density, box plots
# limma           :   plotMDS()
# ggrepel         :   plot MDS plot
# RColorBrewer    :   nice color collections
libraries = c("xlsx","RSQLite","dplyr","reshape2","ggplot2","limma","ggrepel","RColorBrewer")
lapply(libraries,require, character.only = TRUE)
### load functions
source("/scratch/drfz/Good2018/PRI_funs.R")



# load data base -----------------------------------------
db.path=file.path("","data","databases")
db.name="RB_20191002_Good2018.sqlite3"

fcs$connectDb(file.path(db.path,db.name))
print(head(dbListTables(fcs$conn)))


printf("Creat dataframe from %s::%s....",subgroup,project.name)

### metatable file
# cohort_full=readxl::read_excel("/scratch/drfz/Good2018/patient_cohort.xlsx",1,col_names=TRUE,progress=T)
cohort_full=read.csv("/scratch/drfz/Good2018/patient_cohort.csv")
cohort = cohort_full[,c(1,5,8,11,16,15)]

sub.set = cohort[which(cohort$Cohort==subgroup),1:5]
# valid.set = cohort[which(cohort$Cohort=="Validation"),1:5]
# test.set = cohort[which(cohort$Cohort=="NA"),1:5]

### select project
project.idx = which(dbListTables(fcs$conn)==project.name)
fileID = fcs$getDFtable(paste0(project.name,"_fileIdentity"))
stainID = fcs$getDFtable(paste0(project.name,"_markerIdentity"))
### and according file name
data.table.name = sprintf("%s/Rdata/%s_%s_%s_df_cof%s",folder.path,project.name,subgroup,coverage,cofactor)
patient.not.found.name = sprintf("/scratch/drfz/Good2018/%s_%s_patientNOTFOUND.txt",project.name,subgroup)
new.sub.set.name = sprintf("/scratch/drfz/Good2018/%s_%s_patient_cohort.csv",project.name,subgroup)

# go through training set and load single files and merge into one temp.data.all -----------------
### activate to create the two table structures
if (TRUE) {
  label.name = vector()
  ncells.trimmed.total = 0
  temp.data.all = data.man.all = data.frame()
  pat.not.found = new.sub.set.id = vector()
  for ( i in 1:nrow(sub.set)) {
    pat.id = sub.set[i,1]
    file.name = paste0(pat.id,"_",project.name.long,".fcs")
    file.idx = fileID$file_ID[which(fileID$filename==file.name)]
    if (length(file.idx)==0) {
      file.name = paste0(pat.id,"_",tolower(project.name.long),".fcs") # file names vary..
      file.idx = fileID$file_ID[which(fileID$filename==file.name)]
    }
    
    ### if still not found, collect the patient IDs and write them out.
    if (length(file.idx)==0) {
      pat.not.found = c(pat.not.found,as.vector(pat.id))
      printf("Could not find patient %s",file.name)
      next;
    }
    new.sub.set.id = c(new.sub.set.id,as.vector(pat.id))
    ### get data with cofactor
    temp.data = fcs$getData(table = project.name,
                            fileidx = file.idx,
                            cofactor=cofactor)
    
    ### set negative values to zero
    temp.data[temp.data<0] <- 0
    
    ### get short names
    colnames(temp.data) = stainID$shortname[which(stainID$file_ID==file.idx)]
    
    ### trim outliers
    if (trimming) {
      trim.size = 0.0005
      ncells.trimmed = 0
      for (t in 1:ncol(temp.data)) {
        trim.idx = which(temp.data[,t] > quantile(temp.data[,t],c(1-trim.size)))
        temp.data = temp.data[-trim.idx,]
        
        ncells.trimmed = ncells.trimmed + length(trim.idx)
      }
      ncells.trimmed.total = ncells.trimmed.total + ncells.trimmed
      printf("Data trimmed: -%s cells", ncells.trimmed)
    }
    
    ### create one table structure with 4 columns only
    data.man.idx = rep(file.idx,nrow(temp.data))
    data.man = data.frame(file_id=data.man.idx,temp.data)
    data.man = melt(data.man, id.var = "file_id", value.name = "expression", variable.name = "marker")
    ### match the condition "Relapse Status"
    matchmeta = match(pat.id, sub.set$Patient.ID)
    data.man$condition = sub.set$Relapse.Status[matchmeta]
    
    if ( i==1 ) {
      data.man.all = data.man
    }  else {
      ### merge() TAKES A LOT LONGER THAN bind_rows() (>25times)
      # data.man.all <- merge(data.man.all,data.man,all=T)
      data.man.all = bind_rows(data.man.all,data.man)
      # print(dim(data.man.all))
    }
    
    ### create second table with regular structure
    # concenate temp.data into one matrix: temp.data.all
    # generate sample IDs corresponding to each cell in the 'data' matrix
    col.fileid = rep(unlist(pat.id),nrow(temp.data))
    temp.data = bind_cols(file_id=col.fileid,temp.data)
    if ( i==1 ) {
      temp.data.all = temp.data
    }  else {
      ### merge() TAKES A LOT LONGER THAN bind_rows() (>25times)
      # temp.data.all <- merge(temp.data.all,temp.data,all=T)
      temp.data.all <- bind_rows(temp.data.all,temp.data)
      # print(dim(temp.data.all))
    }
  }
  printf("Total cells trimmed: %s",ncells.trimmed.total)
  print("Done loading and manipulating data for Day3.")
  print(dim(temp.data.all))
  # Validation: 470928     40
  # Training: 1943896      40
  saveRDS(temp.data.all, file = paste0(data.table.name,".rds"))
  printf("Saved temp.data.all in %s",paste0(data.table.name,".rds"))
  
  print(dim(data.man.all))
  # Validation: 18366192        4
  # Training:  75811944        4
  saveRDS(data.man.all, file = paste0(data.table.name,"_man.rds"))
  printf("Saved data.man.all in %s",paste0(data.table.name,"_man.rds"))
  
  ### patient ID not found
  write.csv(paste(subgroup,paste(pat.not.found)),patient.not.found.name)
  
  ### save existing patient IDs in separate file name
  new.sub.set = sub.set[which(sub.set$Patient.ID==new.sub.set.id),]
  write.csv(new.sub.set,new.sub.set.name)
  sub.set = new.sub.set
} else {
  if (file.exists(new.sub.set.name)) sub.set = read.csv(new.sub.set.name)
  
  ### load temp.data.all
  temp.data.all = readRDS( file = paste0(data.table.name,".rds"))
  print(dim(temp.data.all))
  # Validation: 470928     40
  # Training: 3697393      40
  
  # load matrix data.all
  data.man.all = readRDS(file = paste0(data.table.name,"_man.rds"))
  print(dim(data.man.all))
  # Validation: 18366192        4
  # Training:  144198327         4
}
dbDisconnect(fcs$conn)

#dunkelrot="#b2182b",hellrot="#f46d43"
#dunkelgrün="#006837",hellgrün="#66bd63"

if (TRUE) {
### darkgreen, red, darkred
sub.set$Relapse.Status = factor(sub.set$Relapse.Status)
color_conditions = c("#006837","#f46d43","#b2182b")
names(color_conditions) = levels(sub.set$Relapse.Status)

# plot all densitites (IT IS HUGE!!!) ---------------------------------------
ggp = ggplot(data.man.all, aes(x = expression, color = condition,
                               group = file_id)) +
  geom_density() +
  facet_wrap( ~ marker, nrow = 10, scales = "free") +
  theme_bw() +
  theme( text = element_text(size = 16)) +      
  scale_color_manual(values = color_conditions)

### get plot dimension
n_panels <- length(unique(ggplot_build(ggp)$data[[1]]$PANEL))
plot.dims = wrap_dims(n_panels)
### save density as pdf
ggsave(paste0(data.table.name,"_density.pdf"),width=3*plot.dims[1],height=3*plot.dims[2])


# plot get median marker exprs per sample ----------------------------------
expr_median_sample_tbl = data.man.all[,c(1,2,3)] %>% group_by(file_id,marker) %>% summarize_all(funs(median))
expr_median_sample_tbl = as.data.frame(expr_median_sample_tbl)
# and create matrix with markers in rows and samples in columns
median.mat = matrix(
  expr_median_sample_tbl[,3]
  , nrow = length(unique(expr_median_sample_tbl$marker))
  , ncol = nrow(sub.set)
)
colnames(median.mat) = sub.set$Patient.ID
rownames(median.mat) = unique(expr_median_sample_tbl$marker)


# MDS plot samples -----------------------------------------------------------
mds = plotMDS(median.mat,plot=F)

ggdf.mds <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                       sample_id = colnames(median.mat))
matchmeta = match(ggdf.mds$sample_id, sub.set$Patient.ID)
ggdf.mds$condition <- sub.set$Relapse.Status[matchmeta]

ggmds = ggplot(ggdf.mds, aes(x = MDS1, y = MDS2, color = condition)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = sample_id)) +
  theme_bw() +
  theme(text = element_text(size=14)) +
  scale_color_manual(values = color_conditions) +
  coord_fixed()
ggsave(paste0(data.table.name,"_MDS.pdf"), width=15, height=8)


# mds plot markers --------------------------------------------------------
mds_t = plotMDS(t(median.mat),plot=F)
ggdf.mds_t <- data.frame(MDS1 = mds_t$x, MDS2 = mds_t$y,
                         marker = rownames(median.mat))

color_rainbow = colorRampPalette(brewer.pal(9,"Set1"))(39)
names(color_rainbow) = rownames(median.mat)
ggmds_t = ggplot(ggdf.mds_t, aes(x = MDS1, y = MDS2, color = rownames(median.mat))) +
  geom_point(size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = marker)) +
  theme_bw() +
  theme_gray(base_size = 14) +
  theme(text = element_text(size=11)) +
  scale_color_manual(values = color_rainbow) 
ggsave(paste0(data.table.name,"_MDSsamples.pdf"), width=15, height=8)


# Non-redundancy score ----------------------------------------------------
### calculate the score
data.NRS = list()
it = 1
for ( i in 1:nrow(sub.set) ) {
  data.tmp = temp.data.all[which(temp.data.all[,1]==as.vector(sub.set$Patient.ID[i])),(2:ncol(temp.data.all))]
  data.NRS[[it]] = NRS(data.tmp)
  it = it +1
}
data.nrs = t(as.data.frame((data.NRS)))
rownames(data.nrs) = sub.set$Patient.ID
nrs.mean = colMeans(data.nrs)

### plot all NRS for ordered markers 
nrs.sorted = names(sort(nrs.mean, decreasing = TRUE))
data.nrs.long = data.frame(data.nrs)
data.nrs.long$sample_id = rownames(data.nrs.long)
colnames(data.nrs.long) = c(colnames(data.tmp),"sample_id")
ggdf.nrs = melt(data.nrs.long, id.var = "sample_id",
                value.name = "nrs", variable.name = "marker")
ggdf.nrs$marker = factor(ggdf.nrs$marker, levels=nrs.sorted)
matchmeta = match(ggdf.nrs$sample_id, sub.set$Patient.ID)
ggdf.nrs$condition <- sub.set$Relapse.Status[matchmeta]

gg_nrs = ggplot(ggdf.nrs, aes(x = marker, y = nrs)) +
  geom_point(aes(color = condition), alpha = 0.9,
             position = position_jitter(width = 0.3, height = 0)) +
  geom_boxplot(outlier.color = NA, fill = NA) +
  stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = color_conditions)
ggsave(paste0(data.table.name,"_NRSfull.pdf"), width=10, height=6)

### plot highest 20 NRS for ordered markers 
nrs.sorted.short = nrs.sorted[1:20]
data.nrs.short = data.frame(data.nrs)
short.idx = which( colnames(data.nrs.short) %in% nrs.sorted.short )
data.nrs.short = data.nrs.short[,short.idx]
data.nrs.short$sample_id = rownames(data.nrs.short)

ggdf.nrs = melt(data.nrs.short, id.var = "sample_id",
                value.name = "nrs", variable.name = "marker")
ggdf.nrs$marker = factor(ggdf.nrs$marker, levels=nrs.sorted)
matchmeta = match(ggdf.nrs$sample_id, sub.set$Patient.ID)
ggdf.nrs$condition <- sub.set$Relapse.Status[matchmeta]

gg_nrs = ggplot(ggdf.nrs, aes(x = marker, y = nrs)) +
  geom_point(aes(color = condition), alpha = 0.9,
             position = position_jitter(width = 0.3, height = 0)) +
  geom_boxplot(outlier.color = NA, fill = NA) +
  stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = color_conditions)
ggsave(paste0(data.table.name,"_NRSshort.pdf"), width=10, height=6)



# plot cells bar plot ------------------------------------------------------------
cells_table = table(temp.data.all$file_id)
if (any(cells_table==0)) cells_table=cells_table[-which(cells_table==0)]
ggdf.cellsbar = data.frame(sample_id = names(cells_table), cell_counts = as.numeric(cells_table) )
matchmeta = match(ggdf.cellsbar$sample_id, sub.set$Patient.ID)
ggdf.cellsbar$condition <- sub.set$Relapse.Status[matchmeta]

ggplot(ggdf.cellsbar, aes(x = sample_id, y = cell_counts, fill = condition)) +
  geom_bar(stat = "identity") +
  geom_text(text = element_text(size=11),aes(label = cell_counts), hjust=0.5, vjust=-0.5, size = 2.5) +
  theme_bw() +
  theme(text = element_text(size=11), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = color_conditions, drop = FALSE) +
  scale_x_discrete(drop = FALSE)
ggsave(paste0(data.table.name,"_barplot.pdf"),width=12,height=5)

}
