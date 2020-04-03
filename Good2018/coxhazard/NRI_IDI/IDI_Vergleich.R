#lade eine DF mit allen nötigen Informationen ( MDR,ROm-Criteria, DDPR Status, Pri-Status....)
df <- read.csv2("C:/Users/eric-/Dropbox/DRFZ PROJEKT/RDS/IDI.csv")

#zum erstellen der SurvObj in diesem R Paket
outcome <- data.frame(df$Survival.Time..Day.,df$Relapse.Status)
colnames(outcome) <- c("TIME","STATUS")
outcome$STATUS <- as.character(outcome$STATUS)
outcome$STATUS[which(outcome$STATUS == "Yes")] <- 1
outcome$STATUS[which(outcome$STATUS == "No")] <- 0
outcome$STATUS <- as.numeric(outcome$STATUS)

#Anpassen der Kriterien an das Paket

#MRD
MRD <- df$MRD.Risk
MRD<-as.character(MRD)
MRD[which(MRD== "Intermediate")] <- 1
MRD[which(MRD== "Standard")] <- 0
MRD[which(MRD== "High")] <- 1
MRD<-as.numeric(MRD)
MRD_nans <- which(is.na(MRD))

#ROM
ROM <- df$NCI.Rome.Risk
ROM <- as.character(ROM)
ROM[which(ROM == "Standard")] <- 0
ROM[which(ROM == "High")] <- 1
ROM<-as.numeric(ROM)
ROM_nans <- which(is.na(ROM))

#DDPR
DDPR <- df$DDPR.Risk

#pri
pri <- df$predictors_Pri
pri<- as.numeric(pri)


#5 Jahre wie in Pri
t0=365*5

#Vergleich der Modell 
#MDR +Pri
idi.MRD <- IDI.INF(outcome[-MRD_nans,],as.matrix(MRD[-MRD_nans]),as.matrix(data_frame(MRD,pri)[-MRD_nans,]),t0)
#Rom Criteria +Pri
idi.ROM <- IDI.INF(outcome[-ROM_nans,],as.matrix(ROM[-ROM_nans]),as.matrix(data_frame(ROM,pri)[-ROM_nans,]),t0)
#DDPR MOdell +Pri
idi.DDPR <- IDI.INF(outcome,as.matrix(DDPR),as.matrix(data_frame(DDPR,pri)),t0)

#Ergebnis output 
print("MDR + PRI")
IDI.INF.OUT(idi.MRD)
print("ROM + PRI")
IDI.INF.OUT(idi.ROM)
print("DDPR +PRI")
IDI.INF.OUT(idi.DDPR)
