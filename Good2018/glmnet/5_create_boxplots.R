#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019
# Goood2018

### Basal representatives
# UPN12/15 low risk
# UPN10/22 high ristk

rm(list = ls())



############### initiate
cofactor = 0.2
set.alpha = 0.01
stat.info = "freq_green"
# stat.info = "absRange"
rmse_threshold = 0.5
coverage = "func"
# coverage = "full"
iterations = 500
date = 191118
current.date = format(Sys.Date(),"%y%m%d")
top = 20

### set paths
folder.path = "/scratch/drfz/Good2018/glmnet"
setwd(folder.path)
project.name = "Basal"


### directory of variables
sub.path = file.path(folder.path, sprintf("%s_FR_%s_%s_a%s_RMSE%s",project.name,coverage,stat.info,set.alpha,rmse_threshold))
file.name = sprintf("%s/%s_signif_quadrants_it%s_table.csv",sub.path,current.date,iterations)
plot.file = sprintf("%s/%s_boxplots_top%s.pdf",sub.path,stat.info,top)

################

### read data
marker.comb = read.table(file.name,sep="\t",header=T)[1:top,2]
# remove NA element, if less than top highest variables is chosen
if (any(is.na(marker.comb))) marker.comb = marker.comb[-which(is.na(marker.comb))]
Training.table.name = sprintf("%s/Rdata/%s_Training_%s_quadrant_%s_cof%s.rds",folder.path,project.name,coverage,stat.info,cofactor)
Validation.table.name = sprintf("%s/Rdata/%s_Validation_%s_quadrant_%s_cof%s.rds",folder.path,project.name,coverage,stat.info,cofactor)

### load libraries
# xlsx            :   Library for Excel reading and creating
# ggpubr          :   https://www.r-bloggers.com/add-p-values-and-significance-levels-to-ggplots/
# reshape2        :   function melt()
libraries = c("xlsx","ggpubr","reshape2")
lapply(libraries,require, character.only = TRUE)
### load functions
source("/scratch/drfz/Good2018/PRI_funs.R")

### metatable file
cohort_full=readxl::read_excel("/scratch/drfz/Good2018/patient_cohort.xlsx",1,col_names=TRUE,progress=T)
cohort = cohort_full[,c(1,5,8,11,16,15)]
train.set = cohort[which(cohort$Cohort=="Training"),1:5]
valid.set = cohort[which(cohort$Cohort=="Validation"),1:5]
trainval.set = bind_rows(train.set,valid.set)

### read table
df.Training = readRDS(Training.table.name)
df.Validation = readRDS(Validation.table.name)
### combine both tables
df.total = bind_rows(df.Training,df.Validation)
sample.size = ncol(df.total)
rownames(df.total) = c(rownames(df.Training),rownames(df.Validation))

## column index of predicted variable in dataset
typeColNum = 1

### condition 1/0
condition = trainval.set$`Relapse Status`[which(trainval.set$`Patient ID` %in% rownames(df.total))]
### change condition "Yes2" / "Yes" to 1
condition[condition %in% c("Yes","Yes2")] = 1
### change condition "No" to 0
condition[condition %in% c("No")] = 0
condition = as.numeric(condition)

df.total= bind_cols(as.data.frame(condition),df.total)
names(df.total)[typeColNum] = "condition"
df.total = as.matrix(df.total)

if (stat.info =="absRange") {
  ### convert NaN/+-Inf to -1
  df.total[is.nan(df.total) | is.infinite(df.total)] <- -1
  ### convert NAs to sample group mean
  for ( i in 2:ncol(df.total)) {
    NA.idx = which(is.na(df.total[,i]))
    
    for (j in NA.idx) {
      tmp = df.total[which(df.total[j,1]==df.total[,1]),i]
      tmp = tmp[-which(is.na(tmp))]
      tmp = round(sample(seq(mean(tmp)-sd(tmp),mean(tmp)+sd(tmp),by=0.01),1),2)
      
      df.total[j,i] = tmp
    }
  }
}

df.total = as.data.frame(df.total)
quad.idx = which(colnames(df.total) %in% as.character(marker.comb))
df.predict = df.total[order(rownames(df.total)),c(1,quad.idx)]



### manipulate dataframe
dat.m = melt(df.predict, id.vars="condition")


#"#b2182b","#006837"
#"#f46d43","#66bd63"
### save boxplot as pdf
plot.pred = ggboxplot(dat.m, x="condition",y="value",ylim = c(min(dat.m$value),max(dat.m$value)),
                      #title=sprintf("Box plots of top 10 sections with \"%s\", p=%s",checkCALC,sample.size),
                      color = "condition", palette=c("#b2182b","#006837"), add = "jitter",
                      facet.by = "variable",ncol = round(ncol(df.predict)/3),
                      xlab="",width=0.3
                      #bxp.errorbar=TRUE,bxp.errorbar.width=0.2,
                      # label.rectangle=TRUE
) +
  # facet_wrap(~variable,scales="free") +
  theme(text = element_text(size=12),
        axis.title.x.top = element_text(colour="grey20",
                                        size=10,angle=0,hjust=.5,vjust=0,face="plain"),
        axis.title.x.bottom = element_text(colour="black",size=12),
        plot.margin = unit(c(0.5,1,0.3,1), "cm")) +
  stat_compare_means(method="t.test",label = "p.signif", label.x = 1.5, label.y = 0)
ggsave(plot.file,plot.pred,width=2.5*round(ncol(df.predict)/3),height=round(nrow(df.predict)/round(ncol(df.predict)/3)))

printf("Written in %s",plot.file)







