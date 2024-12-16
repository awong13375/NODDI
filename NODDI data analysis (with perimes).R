library(irr)
library(stringr)
library(blandr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggpubr)
library(ggeasy)
library(readxl)
library(FSA)
library(MASS)
library(Hmisc)
library(brant)
library(generalhoslem)
library(DescTools)
library(aod)
library(pROC)

setwd("G:/My Drive/Desktop files/McGill DR/NODDI")

# Load data ------------------------------------
df_SAH = read_excel("NODDI Dataset.xlsx", sheet="SAH")
df_PM = read_excel("NODDI Dataset.xlsx", sheet="Perimes")

# remove missing values
df_SAH = subset(df_SAH, df_SAH$DCI!="NA")
df_PM = subset(df_PM, df_PM$DCI!="NA")


# Create column to categorize outcome (0 = nothing, 1 = radiologic vasospasm, 2 = DCI)
df_SAH['Outcome'] = c(0)
df_PM['Outcome'] = c(3)

df_SAH[df_SAH$`Radiologic vasospasm` == 1, 'Outcome'] = 1
df_SAH[df_SAH$DCI == 1, 'Outcome'] = 2

df_SAH = rbind(df_SAH, df_PM)

# Table 1 values ------------------------------------
outcomes = c(3, 0, 1, 2)

for (i in outcomes){
  col = c()
  data = subset(df_SAH, df_SAH$Outcome == i)
  
  col = append(col, nrow(data))
  ## Age ------------------------------------
  col = append(col, paste(mean(data$Age, na.rm=TRUE), "(", round(sd(data$Age, na.rm=TRUE), digit = 1), ")"))
  # Gender ------------------------------------
  col=append(col, paste(nrow(subset(data, data$Gender==0)), "(", round(nrow(subset(data, data$Gender==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$Gender==1)), "(", round(nrow(subset(data, data$Gender==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$Gender))), "(", round(nrow(subset(data, is.na(data$Gender)))/nrow(data)*100, digit = 1), ")"))
  # Ethnicity ------------------------------------
  col=append(col, paste(nrow(subset(data, data$Ethnicity=="Caucasian")), "(", round(nrow(subset(data, data$Ethnicity=="Caucasian"))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$Ethnicity!="Caucasian" & data$Ethnicity!="NA")), "(", round(nrow(subset(data, data$Ethnicity!="Caucasian" & data$Ethnicity!="NA"))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$Ethnicity))), "(", round(nrow(subset(data, is.na(data$Ethnicity)))/nrow(data)*100, digit = 1), ")"))
  # Smoke ------------------------------------
  col=append(col, paste(nrow(subset(data, data$Smoke==0)), "(", round(nrow(subset(data, data$Smoke==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$Smoke==1)), "(", round(nrow(subset(data, data$Smoke==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$Smoke))), "(", round(nrow(subset(data, is.na(data$Smoke)))/nrow(data)*100, digit = 1), ")"))
  # HTN ------------------------------------
  col=append(col, paste(nrow(subset(data, data$HTN==0)), "(", round(nrow(subset(data, data$HTN==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$HTN==1)), "(", round(nrow(subset(data, data$HTN==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$HTN))), "(", round(nrow(subset(data, is.na(data$HTN)))/nrow(data)*100, digit = 1), ")"))
  # EtOH ------------------------------------
  col=append(col, paste(nrow(subset(data, data$EtOH==0)), "(", round(nrow(subset(data, data$EtOH==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$EtOH==1)), "(", round(nrow(subset(data, data$EtOH==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$EtOH))), "(", round(nrow(subset(data, is.na(data$EtOH)))/nrow(data)*100, digit = 1), ")"))
  # DLP ------------------------------------
  col=append(col, paste(nrow(subset(data, data$DLP==0)), "(", round(nrow(subset(data, data$DLP==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$DLP==1)), "(", round(nrow(subset(data, data$DLP==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$DLP))), "(", round(nrow(subset(data, is.na(data$DLP)))/nrow(data)*100, digit = 1), ")"))
  # Diabetes ------------------------------------
  col=append(col, paste(nrow(subset(data, data$Diabetes==0)), "(", round(nrow(subset(data, data$Diabetes==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$Diabetes==1)), "(", round(nrow(subset(data, data$Diabetes==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$Diabetes))), "(", round(nrow(subset(data, is.na(data$Diabetes)))/nrow(data)*100, digit = 1), ")"))
  # FmHx ------------------------------------
  col=append(col, paste(nrow(subset(data, data$FmHx==0)), "(", round(nrow(subset(data, data$FmHx==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$FmHx==1)), "(", round(nrow(subset(data, data$FmHx==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$FmHx))), "(", round(nrow(subset(data, is.na(data$FmHx)))/nrow(data)*100, digit = 1), ")"))
  # Handedness ------------------------------------
  col=append(col, paste(nrow(subset(data, data$Handedness==0)), "(", round(nrow(subset(data, data$Handedness==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$Handedness==1)), "(", round(nrow(subset(data, data$Handedness==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$Handedness))), "(", round(nrow(subset(data, is.na(data$Handedness)))/nrow(data)*100, digit = 1), ")"))
  # Prior stroke ------------------------------------
  col=append(col, paste(nrow(subset(data, data$`Prior stroke`==0)), "(", round(nrow(subset(data, data$`Prior stroke`==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Prior stroke`==1)), "(", round(nrow(subset(data, data$`Prior stroke`==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$`Prior stroke`))), "(", round(nrow(subset(data, is.na(data$`Prior stroke`)))/nrow(data)*100, digit = 1), ")"))
  # Prior aneurysm repair ------------------------------------
  col=append(col, paste(nrow(subset(data, data$`Prior aneurysm repair`==0)), "(", round(nrow(subset(data, data$`Prior aneurysm repair`==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Prior aneurysm repair`==1)), "(", round(nrow(subset(data, data$`Prior aneurysm repair`==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$`Prior aneurysm repair`))), "(", round(nrow(subset(data, is.na(data$`Prior aneurysm repair`)))/nrow(data)*100, digit = 1), ")"))
  # CTD ------------------------------------
  col=append(col, paste(nrow(subset(data, data$`Connective tissue disease/PCKD`==0)), "(", round(nrow(subset(data, data$`Connective tissue disease/PCKD`==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Connective tissue disease/PCKD`==1)), "(", round(nrow(subset(data, data$`Connective tissue disease/PCKD`==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$`Connective tissue disease/PCKD`))), "(", round(nrow(subset(data, is.na(data$`Connective tissue disease/PCKD`)))/nrow(data)*100, digit = 1), ")"))
  # A_P circulation ------------------------------------
  col=append(col, paste(nrow(subset(data, data$`A_P circulation`==0)), "(", round(nrow(subset(data, data$`A_P circulation`==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`A_P circulation`==1)), "(", round(nrow(subset(data, data$`A_P circulation`==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`A_P circulation`==2)), "(", round(nrow(subset(data, data$`A_P circulation`==2))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$`A_P circulation`))), "(", round(nrow(subset(data, is.na(data$`A_P circulation`)))/nrow(data)*100, digit = 1), ")"))
  # Coil/Clip ------------------------------------
  col=append(col, paste(nrow(subset(data, data$`Coil/Clip`==0)), "(", round(nrow(subset(data, data$`Coil/Clip`==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Coil/Clip`==1)), "(", round(nrow(subset(data, data$`Coil/Clip`==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Coil/Clip`==2)), "(", round(nrow(subset(data, data$`Coil/Clip`==2))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$`Coil/Clip`))), "(", round(nrow(subset(data, is.na(data$`Coil/Clip`)))/nrow(data)*100, digit = 1), ")"))
  # SEBES ------------------------------------
  col=append(col, paste(nrow(subset(data, data$SEBES==0)), "(", round(nrow(subset(data, data$SEBES==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$SEBES==1)), "(", round(nrow(subset(data, data$SEBES==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$SEBES==2)), "(", round(nrow(subset(data, data$SEBES==2))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$SEBES==3)), "(", round(nrow(subset(data, data$SEBES==3))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$SEBES==4)), "(", round(nrow(subset(data, data$SEBES==4))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$SEBES))), "(", round(nrow(subset(data, is.na(data$SEBES)))/nrow(data)*100, digit = 1), ")"))
  # mFisher ------------------------------------
  col=append(col, paste(nrow(subset(data, data$mFisher==1)), "(", round(nrow(subset(data, data$mFisher==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$mFisher==2)), "(", round(nrow(subset(data, data$mFisher==2))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$mFisher==3)), "(", round(nrow(subset(data, data$mFisher==3))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$mFisher==4)), "(", round(nrow(subset(data, data$mFisher==4))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$mFisher))), "(", round(nrow(subset(data, is.na(data$mFisher)))/nrow(data)*100, digit = 1), ")"))
  # WFNS ------------------------------------
  col=append(col, paste(nrow(subset(data, data$WFNS==1)), "(", round(nrow(subset(data, data$WFNS==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$WFNS==2)), "(", round(nrow(subset(data, data$WFNS==2))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$WFNS==3)), "(", round(nrow(subset(data, data$WFNS==3))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$WFNS==4)), "(", round(nrow(subset(data, data$WFNS==4))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$WFNS==5)), "(", round(nrow(subset(data, data$WFNS==5))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$WFNS))), "(", round(nrow(subset(data, is.na(data$WFNS)))/nrow(data)*100, digit = 1), ")"))
  # HH ------------------------------------
  col=append(col, paste(nrow(subset(data, data$`Hunt and Hess`==1)), "(", round(nrow(subset(data, data$`Hunt and Hess`==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Hunt and Hess`==2)), "(", round(nrow(subset(data, data$`Hunt and Hess`==2))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Hunt and Hess`==3)), "(", round(nrow(subset(data, data$`Hunt and Hess`==3))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Hunt and Hess`==4)), "(", round(nrow(subset(data, data$`Hunt and Hess`==4))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Hunt and Hess`==5)), "(", round(nrow(subset(data, data$`Hunt and Hess`==5))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$`Hunt and Hess`))), "(", round(nrow(subset(data, is.na(data$`Hunt and Hess`)))/nrow(data)*100, digit = 1), ")"))
  # Milrinone ------------------------------------
  col=append(col, paste(nrow(subset(data, data$Milrinone==0)), "(", round(nrow(subset(data, data$Milrinone==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$Milrinone==1)), "(", round(nrow(subset(data, data$Milrinone==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$Milrinone))), "(", round(nrow(subset(data, is.na(data$Milrinone)))/nrow(data)*100, digit = 1), ")"))
  # EVD ------------------------------------
  col=append(col, paste(nrow(subset(data, data$EVD==0)), "(", round(nrow(subset(data, data$EVD==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$EVD==1)), "(", round(nrow(subset(data, data$EVD==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$EVD))), "(", round(nrow(subset(data, is.na(data$EVD)))/nrow(data)*100, digit = 1), ")"))
  
  if (i == 3){
    tab1 = as.data.frame(col)
  } else {
    tab1 = cbind(tab1, col)
  }
  
}
colnames(tab1) = c("PM","No complications", "Radiologic Vasospasm", "DCI")
rownames(tab1) = c("n","mean Age",
                   "M","F","Gender NA",
                   "Caucasian","Non-caucasian","Ethnicity NA",
                   "No smoke","Smoke","Smoke NA",
                   "No HTN","HTN","HTN NA",
                   "No EtOH","EtOH","EtOH NA",
                   "No DLP","DLP","DLP NA",
                   "No Diabetes","Diabetes","Diabetes NA",
                   "No FmHx","FmHx","FmHx NA",
                   "R Hand","L Hand","Handedness NA",
                   "No stroke","Stroke","Stroke NA",
                   "No repair","Repair","Repair NA",
                   "No CTD","CTD","CTD NA",
                   "Ant Circ","Post Circ","No aneurysm","Circ NA",
                   "Coil","Clip","NoTx","TxNA",
                   "SEBES0","SEBES1","SEBES2","SEBES3","SEBES4","SEBESNA",
                   "mFish1","mFish2","mFish3","mFish4","mFishNA",
                   "WFNS1","WFNS2","WFNS3","WFNS4","WFNS5","WFNSNA",
                   "HH1","HH2","HH3","HH4","HH5","HHNA",
                   "No Milrinone","Milrinone","Milrinone NA",
                   "No EVD","EVD","EVD NA"
)

write.csv(tab1, "PM Table 1 values.csv")

# Table 1 HCP values ------------------------------------
HCP = c(0, 1)

for (i in HCP){
  col = c()
  data = subset(df_SAH, df_SAH$HCP == i)
  
  col = append(col, nrow(data))
  ## Age ------------------------------------
  col = append(col, paste(mean(data$Age, na.rm=TRUE), "(", round(sd(data$Age, na.rm=TRUE), digit = 1), ")"))
  # Gender ------------------------------------
  col=append(col, paste(nrow(subset(data, data$Gender==0)), "(", round(nrow(subset(data, data$Gender==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$Gender==1)), "(", round(nrow(subset(data, data$Gender==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$Gender))), "(", round(nrow(subset(data, is.na(data$Gender)))/nrow(data)*100, digit = 1), ")"))
  # Ethnicity ------------------------------------
  col=append(col, paste(nrow(subset(data, data$Ethnicity=="Caucasian")), "(", round(nrow(subset(data, data$Ethnicity=="Caucasian"))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$Ethnicity!="Caucasian" & data$Ethnicity!="NA")), "(", round(nrow(subset(data, data$Ethnicity!="Caucasian" & data$Ethnicity!="NA"))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$Ethnicity))), "(", round(nrow(subset(data, is.na(data$Ethnicity)))/nrow(data)*100, digit = 1), ")"))
  # Smoke ------------------------------------
  col=append(col, paste(nrow(subset(data, data$Smoke==0)), "(", round(nrow(subset(data, data$Smoke==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$Smoke==1)), "(", round(nrow(subset(data, data$Smoke==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$Smoke))), "(", round(nrow(subset(data, is.na(data$Smoke)))/nrow(data)*100, digit = 1), ")"))
  # HTN ------------------------------------
  col=append(col, paste(nrow(subset(data, data$HTN==0)), "(", round(nrow(subset(data, data$HTN==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$HTN==1)), "(", round(nrow(subset(data, data$HTN==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$HTN))), "(", round(nrow(subset(data, is.na(data$HTN)))/nrow(data)*100, digit = 1), ")"))
  # EtOH ------------------------------------
  col=append(col, paste(nrow(subset(data, data$EtOH==0)), "(", round(nrow(subset(data, data$EtOH==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$EtOH==1)), "(", round(nrow(subset(data, data$EtOH==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$EtOH))), "(", round(nrow(subset(data, is.na(data$EtOH)))/nrow(data)*100, digit = 1), ")"))
  # DLP ------------------------------------
  col=append(col, paste(nrow(subset(data, data$DLP==0)), "(", round(nrow(subset(data, data$DLP==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$DLP==1)), "(", round(nrow(subset(data, data$DLP==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$DLP))), "(", round(nrow(subset(data, is.na(data$DLP)))/nrow(data)*100, digit = 1), ")"))
  # Diabetes ------------------------------------
  col=append(col, paste(nrow(subset(data, data$Diabetes==0)), "(", round(nrow(subset(data, data$Diabetes==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$Diabetes==1)), "(", round(nrow(subset(data, data$Diabetes==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$Diabetes))), "(", round(nrow(subset(data, is.na(data$Diabetes)))/nrow(data)*100, digit = 1), ")"))
  # FmHx ------------------------------------
  col=append(col, paste(nrow(subset(data, data$FmHx==0)), "(", round(nrow(subset(data, data$FmHx==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$FmHx==1)), "(", round(nrow(subset(data, data$FmHx==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$FmHx))), "(", round(nrow(subset(data, is.na(data$FmHx)))/nrow(data)*100, digit = 1), ")"))
  # Handedness ------------------------------------
  col=append(col, paste(nrow(subset(data, data$Handedness==0)), "(", round(nrow(subset(data, data$Handedness==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$Handedness==1)), "(", round(nrow(subset(data, data$Handedness==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$Handedness))), "(", round(nrow(subset(data, is.na(data$Handedness)))/nrow(data)*100, digit = 1), ")"))
  # Prior stroke ------------------------------------
  col=append(col, paste(nrow(subset(data, data$`Prior stroke`==0)), "(", round(nrow(subset(data, data$`Prior stroke`==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Prior stroke`==1)), "(", round(nrow(subset(data, data$`Prior stroke`==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$`Prior stroke`))), "(", round(nrow(subset(data, is.na(data$`Prior stroke`)))/nrow(data)*100, digit = 1), ")"))
  # Prior aneurysm repair ------------------------------------
  col=append(col, paste(nrow(subset(data, data$`Prior aneurysm repair`==0)), "(", round(nrow(subset(data, data$`Prior aneurysm repair`==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Prior aneurysm repair`==1)), "(", round(nrow(subset(data, data$`Prior aneurysm repair`==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$`Prior aneurysm repair`))), "(", round(nrow(subset(data, is.na(data$`Prior aneurysm repair`)))/nrow(data)*100, digit = 1), ")"))
  # CTD ------------------------------------
  col=append(col, paste(nrow(subset(data, data$`Connective tissue disease/PCKD`==0)), "(", round(nrow(subset(data, data$`Connective tissue disease/PCKD`==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Connective tissue disease/PCKD`==1)), "(", round(nrow(subset(data, data$`Connective tissue disease/PCKD`==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$`Connective tissue disease/PCKD`))), "(", round(nrow(subset(data, is.na(data$`Connective tissue disease/PCKD`)))/nrow(data)*100, digit = 1), ")"))
  # A_P circulation ------------------------------------
  col=append(col, paste(nrow(subset(data, data$`A_P circulation`==0)), "(", round(nrow(subset(data, data$`A_P circulation`==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`A_P circulation`==1)), "(", round(nrow(subset(data, data$`A_P circulation`==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`A_P circulation`==2)), "(", round(nrow(subset(data, data$`A_P circulation`==2))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$`A_P circulation`))), "(", round(nrow(subset(data, is.na(data$`A_P circulation`)))/nrow(data)*100, digit = 1), ")"))
  # Coil/Clip ------------------------------------
  col=append(col, paste(nrow(subset(data, data$`Coil/Clip`==0)), "(", round(nrow(subset(data, data$`Coil/Clip`==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Coil/Clip`==1)), "(", round(nrow(subset(data, data$`Coil/Clip`==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Coil/Clip`==2)), "(", round(nrow(subset(data, data$`Coil/Clip`==2))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$`Coil/Clip`))), "(", round(nrow(subset(data, is.na(data$`Coil/Clip`)))/nrow(data)*100, digit = 1), ")"))
  # SEBES ------------------------------------
  col=append(col, paste(nrow(subset(data, data$SEBES==0)), "(", round(nrow(subset(data, data$SEBES==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$SEBES==1)), "(", round(nrow(subset(data, data$SEBES==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$SEBES==2)), "(", round(nrow(subset(data, data$SEBES==2))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$SEBES==3)), "(", round(nrow(subset(data, data$SEBES==3))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$SEBES==4)), "(", round(nrow(subset(data, data$SEBES==4))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$SEBES))), "(", round(nrow(subset(data, is.na(data$SEBES)))/nrow(data)*100, digit = 1), ")"))
  # mFisher ------------------------------------
  col=append(col, paste(nrow(subset(data, data$mFisher==1)), "(", round(nrow(subset(data, data$mFisher==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$mFisher==2)), "(", round(nrow(subset(data, data$mFisher==2))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$mFisher==3)), "(", round(nrow(subset(data, data$mFisher==3))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$mFisher==4)), "(", round(nrow(subset(data, data$mFisher==4))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$mFisher))), "(", round(nrow(subset(data, is.na(data$mFisher)))/nrow(data)*100, digit = 1), ")"))
  # WFNS ------------------------------------
  col=append(col, paste(nrow(subset(data, data$WFNS==1)), "(", round(nrow(subset(data, data$WFNS==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$WFNS==2)), "(", round(nrow(subset(data, data$WFNS==2))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$WFNS==3)), "(", round(nrow(subset(data, data$WFNS==3))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$WFNS==4)), "(", round(nrow(subset(data, data$WFNS==4))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$WFNS==5)), "(", round(nrow(subset(data, data$WFNS==5))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$WFNS))), "(", round(nrow(subset(data, is.na(data$WFNS)))/nrow(data)*100, digit = 1), ")"))
  # HH ------------------------------------
  col=append(col, paste(nrow(subset(data, data$`Hunt and Hess`==1)), "(", round(nrow(subset(data, data$`Hunt and Hess`==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Hunt and Hess`==2)), "(", round(nrow(subset(data, data$`Hunt and Hess`==2))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Hunt and Hess`==3)), "(", round(nrow(subset(data, data$`Hunt and Hess`==3))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Hunt and Hess`==4)), "(", round(nrow(subset(data, data$`Hunt and Hess`==4))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$`Hunt and Hess`==5)), "(", round(nrow(subset(data, data$`Hunt and Hess`==5))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$`Hunt and Hess`))), "(", round(nrow(subset(data, is.na(data$`Hunt and Hess`)))/nrow(data)*100, digit = 1), ")"))
  # Milrinone ------------------------------------
  col=append(col, paste(nrow(subset(data, data$Milrinone==0)), "(", round(nrow(subset(data, data$Milrinone==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$Milrinone==1)), "(", round(nrow(subset(data, data$Milrinone==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$Milrinone))), "(", round(nrow(subset(data, is.na(data$Milrinone)))/nrow(data)*100, digit = 1), ")"))
  # EVD ------------------------------------
  col=append(col, paste(nrow(subset(data, data$EVD==0)), "(", round(nrow(subset(data, data$EVD==0))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, data$EVD==1)), "(", round(nrow(subset(data, data$EVD==1))/nrow(data)*100, digit = 1), ")"))
  col=append(col, paste(nrow(subset(data, is.na(data$EVD))), "(", round(nrow(subset(data, is.na(data$EVD)))/nrow(data)*100, digit = 1), ")"))
  
  if (i == 0){
    tab1 = as.data.frame(col)
  } else {
    tab1 = cbind(tab1, col)
  }
  
}
colnames(tab1) = c("No HCP", "HCP")
rownames(tab1) = c("n","mean Age",
                   "M","F","Gender NA",
                   "Caucasian","Non-caucasian","Ethnicity NA",
                   "No smoke","Smoke","Smoke NA",
                   "No HTN","HTN","HTN NA",
                   "No EtOH","EtOH","EtOH NA",
                   "No DLP","DLP","DLP NA",
                   "No Diabetes","Diabetes","Diabetes NA",
                   "No FmHx","FmHx","FmHx NA",
                   "R Hand","L Hand","Handedness NA",
                   "No stroke","Stroke","Stroke NA",
                   "No repair","Repair","Repair NA",
                   "No CTD","CTD","CTD NA",
                   "Ant Circ","Post Circ","No aneurysm","Circ NA",
                   "Coil","Clip","NoTx","TxNA",
                   "SEBES0","SEBES1","SEBES2","SEBES3","SEBES4","SEBESNA",
                   "mFish1","mFish2","mFish3","mFish4","mFishNA",
                   "WFNS1","WFNS2","WFNS3","WFNS4","WFNS5","WFNSNA",
                   "HH1","HH2","HH3","HH4","HH5","HHNA",
                   "No Milrinone","Milrinone","Milrinone NA",
                   "No EVD","EVD","EVD NA"
)

write.csv(tab1, "Table 1 HCP values.csv")

# Table 1 stats ------------------------------------

tab1row = c()

## Age ------------------------------------
df_SAH_PM = subset(df_SAH, df_SAH$Outcome == 3)
df_SAH_N = subset(df_SAH, df_SAH$Outcome == 0)
df_SAH_VS = subset(df_SAH, df_SAH$Outcome == 1)
df_SAH_DCI = subset(df_SAH, df_SAH$Outcome == 2)

p_PM = shapiro.test(subset(df_SAH, df_SAH$Outcome == 3)$Age)$p.value
p_N = shapiro.test(subset(df_SAH, df_SAH$Outcome == 0)$Age)$p.value
p_VS = shapiro.test(subset(df_SAH, df_SAH$Outcome == 1)$Age)$p.value
p_DCI = shapiro.test(subset(df_SAH, df_SAH$Outcome == 2)$Age)$p.value

if (p_PM < 0.05 |p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
  tab1row = append(tab1row, kruskal.test(Age ~ Outcome, data = df_SAH)$p.value)
  tab1row = append(tab1row, "kruskal")
} else {
  tab1row = append(tab1row, oneway.test(Age ~ Outcome, data = df_SAH)$p.value)
  tab1row = append(tab1row, "ANOVA")
}

tab1data = t(as.data.frame(tab1row))

## Gender ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$Gender == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$Gender == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$Gender == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$Gender == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$Gender == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$Gender == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$Gender == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$Gender == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Ethnicity ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$Ethnicity == "Caucasian")), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & (df_SAH$Ethnicity != "Caucasian" & df_SAH$Ethnicity != "NA")))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$Ethnicity == "Caucasian")), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & (df_SAH$Ethnicity != "Caucasian" & df_SAH$Ethnicity != "NA")))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$Ethnicity == "Caucasian")), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$Ethnicity != "Caucasian" & df_SAH$Ethnicity != "NA"))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$Ethnicity == "Caucasian")), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$Ethnicity != "Caucasian" & df_SAH$Ethnicity != "NA")))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Smoke ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$Smoke == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$Smoke == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$Smoke == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$Smoke == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$Smoke == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$Smoke == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$Smoke == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$Smoke == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## HTN ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$HTN == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$HTN == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$HTN == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$HTN == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$HTN == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$HTN == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$HTN == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$HTN == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## EtOH ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$EtOH == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$EtOH == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$EtOH == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$EtOH == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$EtOH == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$EtOH == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$EtOH == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$EtOH == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## DLP ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$DLP == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$DLP == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$DLP == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$DLP == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$DLP == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$DLP == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$DLP == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$DLP == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Diabetes ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$Diabetes == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$Diabetes == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$Diabetes == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$Diabetes == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$Diabetes == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$Diabetes == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$Diabetes == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$Diabetes == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## FmHx ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$FmHx == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$FmHx == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$FmHx == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$FmHx == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$FmHx == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$FmHx == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$FmHx == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$FmHx == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Handedness ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$Handedness == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$Handedness == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$Handedness == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$Handedness == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$Handedness == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$Handedness == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$Handedness == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$Handedness == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Prior stroke ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`Prior stroke` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`Prior stroke` == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`Prior stroke` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`Prior stroke` == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`Prior stroke` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`Prior stroke` == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`Prior stroke` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`Prior stroke` == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Prior aneurysm repair ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`Prior aneurysm repair` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`Prior aneurysm repair` == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`Prior aneurysm repair` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`Prior aneurysm repair` == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`Prior aneurysm repair` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`Prior aneurysm repair` == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`Prior aneurysm repair` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`Prior aneurysm repair` == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## CTD ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`Connective tissue disease/PCKD` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`Connective tissue disease/PCKD` == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`Connective tissue disease/PCKD` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`Connective tissue disease/PCKD` == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`Connective tissue disease/PCKD` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`Connective tissue disease/PCKD` == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`Connective tissue disease/PCKD` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`Connective tissue disease/PCKD` == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## A_P circulation ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`A_P circulation` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`A_P circulation` == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`A_P circulation` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`A_P circulation` == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`A_P circulation` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`A_P circulation` == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`A_P circulation` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`A_P circulation` == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)
## Coil/Clip ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`Coil/Clip` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`Coil/Clip` == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`Coil/Clip` == 2))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`Coil/Clip` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`Coil/Clip` == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`Coil/Clip` == 2))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`Coil/Clip` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`Coil/Clip` == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`Coil/Clip` == 2))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`Coil/Clip` == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`Coil/Clip` == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`Coil/Clip` == 2)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## SEBES ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$SEBES == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$SEBES == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$SEBES == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$SEBES == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$SEBES == 4))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$SEBES == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$SEBES == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$SEBES == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$SEBES == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$SEBES == 4))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$SEBES == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$SEBES == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$SEBES == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$SEBES == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$SEBES == 4))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$SEBES == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$SEBES == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$SEBES == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$SEBES == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$SEBES == 4)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## mFisher ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$mFisher == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$mFisher == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$mFisher == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$mFisher == 4))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$mFisher == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$mFisher == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$mFisher == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$mFisher == 4))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$mFisher == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$mFisher == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$mFisher == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$mFisher == 4))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$mFisher == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$mFisher == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$mFisher == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$mFisher == 4)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## WFNS ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$WFNS == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$WFNS == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$WFNS == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$WFNS == 4)),
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$WFNS == 5))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$WFNS == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$WFNS == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$WFNS == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$WFNS == 4)),
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$WFNS == 5))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$WFNS == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$WFNS == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$WFNS == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$WFNS == 4)),
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$WFNS == 5))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$WFNS == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$WFNS == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$WFNS == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$WFNS == 4)),
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$WFNS == 5)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Hunt and Hess ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`Hunt and Hess` == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`Hunt and Hess` == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`Hunt and Hess` == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`Hunt and Hess` == 4)),
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$`Hunt and Hess` == 5))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`Hunt and Hess` == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`Hunt and Hess` == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`Hunt and Hess` == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`Hunt and Hess` == 4)),
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$`Hunt and Hess` == 5))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`Hunt and Hess` == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`Hunt and Hess` == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`Hunt and Hess` == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`Hunt and Hess` == 4)),
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$`Hunt and Hess` == 5))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`Hunt and Hess` == 1)),
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`Hunt and Hess` == 2)),
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`Hunt and Hess` == 3)),
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`Hunt and Hess` == 4)),
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$`Hunt and Hess` == 5)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Milrinone ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$Milrinone == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$Milrinone == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$Milrinone == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$Milrinone == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$Milrinone == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$Milrinone == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$Milrinone == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$Milrinone == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## EVD ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$EVD == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 3 & df_SAH$EVD == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$EVD == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 0 & df_SAH$EVD == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$EVD == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 1 & df_SAH$EVD == 1))),
    c(nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$EVD == 0)), 
      nrow(subset(df_SAH, df_SAH$Outcome == 2 & df_SAH$EVD == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Generate table ------------------------------------
colnames(tab1row) = c("p-value","stat test")
rownames(tab1row) = c("Age","Gender","Ethnicity","Smoke","HTN","EtOH","DLP","Diabetes","FmHx",
                      "Handedness","Prior stroke","Prior aneurysm repair","Connective tissue disease/PCKD","A_P circulation","Coil/Clip",
                      "SEBES","mFisher","WFNS","Hunt and Hess","Milrinone","EVD")

write.csv(tab1row, "PM Table 1 stats.csv")


# Table 1 HCP stats ------------------------------------

tab1row = c()

## Age ------------------------------------
df_SAH_noHCP = subset(df_SAH, df_SAH$HCP == 0)
df_SAH_HCP = subset(df_SAH, df_SAH$HCP == 1)

p_noHCP = shapiro.test(subset(df_SAH, df_SAH$HCP == 0)$Age)$p.value
p_HCP = shapiro.test(subset(df_SAH, df_SAH$HCP == 1)$Age)$p.value

if (p_noHCP < 0.05 | p_HCP < 0.05){
  tab1row = append(tab1row, kruskal.test(Age ~ HCP, data = df_SAH)$p.value)
  tab1row = append(tab1row, "kruskal")
} else {
  tab1row = append(tab1row, oneway.test(Age ~ HCP, data = df_SAH)$p.value)
  tab1row = append(tab1row, "ANOVA")
}

tab1data = t(as.data.frame(tab1row))

## Gender ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$Gender == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$Gender == 1))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$Gender == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$Gender == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Ethnicity ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$Ethnicity == "Caucasian")), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & (df_SAH$Ethnicity != "Caucasian" & df_SAH$Ethnicity != "NA")))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$Ethnicity == "Caucasian")), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$Ethnicity != "Caucasian" & df_SAH$Ethnicity != "NA")))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Smoke ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$Smoke == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$Smoke == 1))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$Smoke == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$Smoke == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## HTN ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$HTN == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$HTN == 1))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$HTN == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$HTN == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## EtOH ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$EtOH == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$EtOH == 1))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$EtOH == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$EtOH == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## DLP ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$DLP == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$DLP == 1))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$DLP == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$DLP == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Diabetes ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$Diabetes == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$Diabetes == 1))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$Diabetes == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$Diabetes == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## FmHx ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$FmHx == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$FmHx == 1))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$FmHx == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$FmHx == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Handedness ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$Handedness == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$Handedness == 1))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$Handedness == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$Handedness == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Prior stroke ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`Prior stroke` == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`Prior stroke` == 1))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`Prior stroke` == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`Prior stroke` == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Prior aneurysm repair ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`Prior aneurysm repair` == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`Prior aneurysm repair` == 1))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`Prior aneurysm repair` == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`Prior aneurysm repair` == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## CTD ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`Connective tissue disease/PCKD` == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`Connective tissue disease/PCKD` == 1))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`Connective tissue disease/PCKD` == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`Connective tissue disease/PCKD` == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## A_P circulation ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`A_P circulation` == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`A_P circulation` == 1))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`A_P circulation` == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`A_P circulation` == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)
## Coil/Clip ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`Coil/Clip` == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`Coil/Clip` == 1)),
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`Coil/Clip` == 2))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`Coil/Clip` == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`Coil/Clip` == 1)),
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`Coil/Clip` == 2)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## SEBES ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$SEBES == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$SEBES == 1)),
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$SEBES == 2)),
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$SEBES == 3)),
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$SEBES == 4))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$SEBES == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$SEBES == 1)),
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$SEBES == 2)),
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$SEBES == 3)),
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$SEBES == 4)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## mFisher ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$mFisher == 1)),
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$mFisher == 2)),
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$mFisher == 3)),
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$mFisher == 4))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$mFisher == 1)),
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$mFisher == 2)),
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$mFisher == 3)),
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$mFisher == 4)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## WFNS ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$WFNS == 1)),
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$WFNS == 2)),
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$WFNS == 3)),
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$WFNS == 4)),
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$WFNS == 5))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$WFNS == 1)),
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$WFNS == 2)),
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$WFNS == 3)),
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$WFNS == 4)),
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$WFNS == 5)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Hunt and Hess ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`Hunt and Hess` == 1)),
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`Hunt and Hess` == 2)),
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`Hunt and Hess` == 3)),
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`Hunt and Hess` == 4)),
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$`Hunt and Hess` == 5))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`Hunt and Hess` == 1)),
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`Hunt and Hess` == 2)),
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`Hunt and Hess` == 3)),
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`Hunt and Hess` == 4)),
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$`Hunt and Hess` == 5)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Milrinone ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$Milrinone == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$Milrinone == 1))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$Milrinone == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$Milrinone == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## EVD ------------------------------------
M = as.table(
  cbind(
    c(nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$EVD == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 0 & df_SAH$EVD == 1))),
    c(nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$EVD == 0)), 
      nrow(subset(df_SAH, df_SAH$HCP == 1 & df_SAH$EVD == 1)))
  )
)
tab1row = rbind(tab1row, tryCatch({
  c(chisq.test(M)$p.value, "Chi")
}, warning = function(w) {
  c(fisher.test(M)$p.value, "Fisher")
})
)

## Generate table ------------------------------------
colnames(tab1row) = c("p-value","stat test")
rownames(tab1row) = c("Age","Gender","Ethnicity","Smoke","HTN","EtOH","DLP","Diabetes","FmHx",
                      "Handedness","Prior stroke","Prior aneurysm repair","Connective tissue disease/PCKD","A_P circulation","Coil/Clip",
                      "SEBES","mFisher","WFNS","Hunt and Hess","Milrinone","EVD")

write.csv(tab1row, "Table 1 HCP stats.csv")





# Table 2 values ------------------------------------

outcomes = c(3, 0, 1, 2)

## DCE ------------------------------------

for (i in outcomes){
  col = c()
  data = subset(df_SAH, df_SAH$Outcome == i)
  
  col = append(col, nrow(data))
  
  col = append(col, paste(mean(as.numeric(data$DCE_Ktrans), na.rm=TRUE), "(", sd(as.numeric(data$DCE_Ktrans), na.rm=TRUE), ")"))
  
  col = append(col, paste(mean(as.numeric(data$Vp), na.rm=TRUE), "(", sd(as.numeric(data$Vp), na.rm=TRUE), ")"))
  
  if (i == 3){
    tab2_DCE = as.data.frame(col)
  } else {
    tab2_DCE = cbind(tab2_DCE, col)
  }
  
}
colnames(tab2_DCE) = c("PM", "No complications", "Radiologic vasospasm","DCI")
rownames(tab2_DCE) = c("n","Ktrans", "Vp")



## DTI_ALPS ------------------------------------

for (i in outcomes){
  col = c()
  data = subset(df_SAH, df_SAH$Outcome == i)
  
  col = append(col, paste(mean(as.numeric(data$DTI_ALPS_left_manual), na.rm=TRUE), "(", sd(as.numeric(data$DTI_ALPS_left_manual), na.rm=TRUE), ")"))
  
  col = append(col, paste(mean(as.numeric(data$DTI_ALPS_right_manual), na.rm=TRUE), "(", sd(as.numeric(data$DTI_ALPS_right_manual), na.rm=TRUE), ")"))
  
  col = append(col, paste(mean(as.numeric(data$DTI_ALPS_dom_manual), na.rm=TRUE), "(", sd(as.numeric(data$DTI_ALPS_dom_manual), na.rm=TRUE), ")"))
  
  col = append(col, paste(mean(as.numeric(data$DTI_ALPS_nondom_manual), na.rm=TRUE),"(", sd(as.numeric(data$DTI_ALPS_nondom_manual), na.rm=TRUE), ")"))
  
  if (i == 3){
    tab2_DTIALPS = as.data.frame(col)
  } else {
    tab2_DTIALPS = cbind(tab2_DTIALPS, col)
  }
  
}
colnames(tab2_DTIALPS) = c("PM", "No complications", "Radiologic vasospasm","DCI")
rownames(tab2_DTIALPS) = c("DTIALPS_left_manual",
                           "DTIALPS_right_manual",
                           "DTIALPS_dom_manual",
                           "DTIALPS_nondom_manual")

tab2 = rbind(tab2_DCE, tab2_DTIALPS)

## NODDI ------------------------------------

# remove non NODDI columns
df_SAH_NODDI = df_SAH[,54:ncol(df_SAH)]

# extract mean columns
df_SAH_NODDI = df_SAH_NODDI[,seq(2, ncol(df_SAH_NODDI),3)]

# add outcomes column
df_SAH_NODDI = cbind(df_SAH_NODDI, df_SAH$Outcome)
colnames(df_SAH_NODDI)[colnames(df_SAH_NODDI) == 'df_SAH$Outcome'] <- 'Outcome'


n = c(colnames(df_SAH_NODDI[seq(1, ncol(df_SAH_NODDI)-1)]))

outcomes = c(3, 0, 1, 2)

for (i in outcomes){
  col = c()
  
  data = subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == i)
  
  for (k in seq(1,ncol(df_SAH_NODDI)-1)){
    col = append(col, paste(mean(as.numeric(data[,k]), na.rm=TRUE),"(",sd(as.numeric(data[,k]), na.rm=TRUE),")"))
  }
  
  if (i == 3){
    tab2_NODDI = as.data.frame(col)
  } else {
    tab2_NODDI = cbind(tab2_NODDI, col)
  }
  
}

rownames(tab2_NODDI) = as.vector(colnames(df_SAH_NODDI)[seq(1,ncol(df_SAH_NODDI)-1)])
colnames(tab2_NODDI) =c("PM", "No complications", "Radiologic vasospasm", "DCI")
tab2 = rbind(tab2, tab2_NODDI)

write.csv(tab2, "PM Table 2 values.csv")











# Univariate analysis ------------------------------------

outcomes = c(3, 0, 1, 2)

## DCE ------------------------------------

df_SAH = as.data.frame(lapply(df_SAH, as.numeric))

df_SAH_PM = subset(df_SAH, df_SAH$Outcome == 3)
df_SAH_N = subset(df_SAH, df_SAH$Outcome == 0)
df_SAH_VS = subset(df_SAH, df_SAH$Outcome == 1)
df_SAH_DCI = subset(df_SAH, df_SAH$Outcome == 2)

p_PM = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 3)$DCE_Ktrans))$p.value
p_N = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 0)$DCE_Ktrans))$p.value
p_VS = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 1)$DCE_Ktrans))$p.value
p_DCI = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 2)$DCE_Ktrans))$p.value

tab2row = c()

if (p_PM < 0.05 | p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('DCE_Ktrans ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('DCE_Ktrans ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DCE = t(as.data.frame(tab2row))

tab2row = c()

p_PM = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 3)$Vp))$p.value
p_N = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 0)$Vp))$p.value
p_VS = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 1)$Vp))$p.value
p_DCI = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 2)$Vp))$p.value


if (p_PM < 0.05 | p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('Vp ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('Vp ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DCE = rbind(tab2stat_DCE, tab2row)

rownames(tab2stat_DCE) = c("Ktrans","Vp")
colnames(tab2stat_DCE) = c("p-value","stat")



## DTI_ALPS ------------------------------------

df_SAH = as.data.frame(lapply(df_SAH, as.numeric))

df_SAH_PM = subset(df_SAH, df_SAH$Outcome == 3)
df_SAH_N = subset(df_SAH, df_SAH$Outcome == 0)
df_SAH_VS = subset(df_SAH, df_SAH$Outcome == 1)
df_SAH_DCI = subset(df_SAH, df_SAH$Outcome == 2)

p_PM = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 3)$DTI_ALPS_left_manual))$p.value
p_N = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 0)$DTI_ALPS_left_manual))$p.value
p_VS = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 1)$DTI_ALPS_left_manual))$p.value
p_DCI = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 2)$DTI_ALPS_left_manual))$p.value

tab2row = c()

if (p_PM < 0.05 | p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('DTI_ALPS_left_manual ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('DTI_ALPS_left_manual ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DTI = t(as.data.frame(tab2row))

p_PM = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 3)$DTI_ALPS_right_manual))$p.value
p_N = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 0)$DTI_ALPS_right_manual))$p.value
p_VS = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 1)$DTI_ALPS_right_manual))$p.value
p_DCI = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 2)$DTI_ALPS_right_manual))$p.value

tab2row = c()

if (p_PM < 0.05 | p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('DTI_ALPS_right_manual ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('DTI_ALPS_right_manual ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DTI = rbind(tab2stat_DTI, tab2row)

p_PM = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 3)$DTI_ALPS_dom_manual))$p.value
p_N = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 0)$DTI_ALPS_dom_manual))$p.value
p_VS = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 1)$DTI_ALPS_dom_manual))$p.value
p_DCI = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 2)$DTI_ALPS_dom_manual))$p.value

tab2row = c()

if (p_PM < 0.05 | p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('DTI_ALPS_dom_manual ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('DTI_ALPS_dom_manual ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DTI = rbind(tab2stat_DTI, tab2row)

p_PM = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 3)$DTI_ALPS_nondom_manual))$p.value
p_N = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 0)$DTI_ALPS_nondom_manual))$p.value
p_VS = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 1)$DTI_ALPS_nondom_manual))$p.value
p_DCI = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 2)$DTI_ALPS_nondom_manual))$p.value

tab2row = c()

if (p_PM < 0.05 | p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('DTI_ALPS_nondom_manual ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('DTI_ALPS_nondom_manual ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DTI = rbind(tab2stat_DTI, tab2row)

rownames(tab2stat_DTI) = c("DTI_ALPS_left_manual","DTI_ALPS_right_manual","DTI_ALPS_dom_manual","DTI_ALPS_nondom_manual")
colnames(tab2stat_DTI) = c("p-value","stat")

## NODDI ------------------------------------

df_SAH_NODDI_PM = subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == 3)
df_SAH_NODDI_N = subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == 0)
df_SAH_NODDI_VS = subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == 1)
df_SAH_NODDI_DCI = subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == 2)

df_SAH_NODDI = as.data.frame(lapply(df_SAH_NODDI,as.numeric))

for (k in seq(1,ncol(df_SAH_NODDI)-1)){
  p_PM = shapiro.test(as.numeric(subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == 3)[,k]))$p.value
  p_N = shapiro.test(as.numeric(subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == 0)[,k]))$p.value
  p_VS = shapiro.test(as.numeric(subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == 1)[,k]))$p.value
  p_DCI = shapiro.test(as.numeric(subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == 2)[,k]))$p.value
  
  tab2row = c()
  
  if (p_PM < 0.05 |p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
    tab2row = append(tab2row, kruskal.test(as.formula(paste(colnames(df_SAH_NODDI)[k],'~ Outcome')), data = df_SAH_NODDI)$p.value)
    tab2row = append(tab2row, "kruskal")
  } else {
    tab2row = append(tab2row, oneway.test(as.formula(paste(colnames(df_SAH_NODDI)[k],'~ Outcome')), data = df_SAH_NODDI)$p.value)
    tab2row = append(tab2row, "ANOVA")
  }
  
  if (k == 1){
    tab2stat_NODDI = t(as.data.frame(tab2row))
  } else {
    tab2stat_NODDI = rbind(tab2stat_NODDI, tab2row)
  }
}

rownames(tab2stat_NODDI) = as.vector(colnames(df_SAH_NODDI)[seq(1,ncol(df_SAH_NODDI)-1)])
colnames(tab2stat_NODDI) = c("p-value","stat")
tab2stat_NODDI = as.data.frame(tab2stat_NODDI)




tab2stat = rbind(tab2stat_DCE, tab2stat_DTI, tab2stat_NODDI)

tab2stat_VOI = subset(tab2stat, tab2stat$'p-value'<0.05)

write.csv(tab2stat, "PM Univariate analysis.csv")
write.csv(tab2stat_VOI, "PM Univariate analysis variables of interest.csv")


# Post hoc tests for variables of interest ------------------------------------
VOI_tukey = subset(tab2stat_VOI, tab2stat_VOI$stat == "ANOVA")
VOI_tukey = rownames(VOI_tukey)
VOI_dunn = subset(tab2stat_VOI, tab2stat_VOI$stat == "kruskal")
VOI_dunn = rownames(VOI_dunn)

## Tukey HSD ------------------------------------

for (i in seq(1:length(VOI_tukey))){
  Tukey = TukeyHSD(aov(as.formula(paste(VOI_tukey[i], "~ factor(Outcome)")), data = df_SAH))[1]
  
  if (i == 1){
    Tukey_list = Tukey
  } else {
    Tukey_list = append(Tukey_list, Tukey)
  }
  
}

names(Tukey_list) = VOI_tukey

capture.output(Tukey_list, file = "PM TukeyHSD.txt")

## Dunn test ------------------------------------

for (i in seq(1:length(VOI_dunn))){
  dunn = dunnTest(as.formula(paste(VOI_dunn[i], "~ factor(Outcome)")), data = df_SAH)[2]
  
  if (i == 1){
    dunn_list = dunn
  } else {
    dunn_list = append(dunn_list, dunn)
  }
  
}

names(dunn_list) = VOI_dunn

capture.output(dunn_list, file = "PM dunn.txt")


# Box plots ------------------------------------
## DTI ALPS ------------------------------------

df_SAH$outcome_name = c(0)
df_SAH[df_SAH$Outcome == 3, 'outcome_name'] = "PM"
df_SAH[df_SAH$Outcome == 0, 'outcome_name'] = "No Complications"
df_SAH[df_SAH$Outcome == 1, 'outcome_name'] = "Vasospasm"
df_SAH[df_SAH$Outcome == 2, 'outcome_name'] = "DCI"

df_SAH$outcome_name = factor(df_SAH$outcome_name, levels = c("PM", "No Complications","Vasospasm","DCI"))

plot_DTIALPS_right = ggplot(df_SAH, aes(x=outcome_name, y=DTI_ALPS_right_manual))
plot_DTIALPS_right = plot_DTIALPS_right + geom_boxplot() + theme_classic() + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)

plot_DTIALPS_dom = ggplot(df_SAH, aes(x=outcome_name, y=DTI_ALPS_dom_manual))
plot_DTIALPS_dom = plot_DTIALPS_dom + geom_boxplot() + theme_classic() + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)


plot_DTIALPS_left = ggplot(df_SAH, aes(x=outcome_name, y=DTI_ALPS_left_manual))
plot_DTIALPS_left = plot_DTIALPS_left + geom_boxplot() + theme_classic() + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)

plot_DTIALPS_nondom = ggplot(df_SAH, aes(x=outcome_name, y=DTI_ALPS_nondom_manual))
plot_DTIALPS_nondom = plot_DTIALPS_nondom + geom_boxplot() + theme_classic() + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)

plot_grid(plot_DTIALPS_right, plot_DTIALPS_dom, plot_DTIALPS_left, plot_DTIALPS_nondom, 
          labels = c("Manual_method"),
          ncol = 2, nrow = 2)

## DCE ktrans ------------------------------------

plot_DCE_ktrans = ggplot(df_SAH, aes(x=outcome_name, y=DCE_Ktrans))
plot_DCE_ktrans = plot_DCE_ktrans + geom_boxplot() + theme_classic() + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_y_continuous(limits = c(0, 0.000001)) +
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)

plot_DCE_Vp = ggplot(df_SAH, aes(x=outcome_name, y=Vp))
plot_DCE_Vp = plot_DCE_Vp + geom_boxplot() + theme_classic() + theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)

plot_grid(plot_DCE_ktrans, plot_DCE_Vp, 
          labels = c("DCE"),
          ncol = 2, nrow = 1)


# Multivariable analyses ------------------------------------

## Crude models  ------------------------------------



### DTI_ALPS Right ------------------------------------

model_dtialps_right = polr(outcome_name ~ DTI_ALPS_right_manual, data = df_SAH, Hess=TRUE)
ctable <- coef(summary(model_dtialps_right))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
ctable <- cbind(ctable, "p value" = p)
ci <- confint(model_dtialps_right)
model_dtialps_right
exp(coef(model_dtialps_right))
exp(rbind(OR = coef(model_dtialps_right), ci))



### DTI_ALPS non-dom ------------------------------------

model_dtialps_nondom = polr(outcome_name ~ DTI_ALPS_nondom_manual, data = df_SAH, Hess=TRUE)
ctable <- coef(summary(model_dtialps_nondom))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
ctable <- cbind(ctable, "p value" = p)
ci <- confint(model_dtialps_nondom)
model_dtialps_nondom
exp(coef(model_dtialps_nondom))
exp(rbind(OR = coef(model_dtialps_nondom), ci))


## Adjusted models ------------------------------------

### DTI_ALPS Right ------------------------------------

model_dtialps_right = polr(outcome_name ~ DTI_ALPS_right_manual + Ethnicity, data = df_SAH, Hess=TRUE)
ctable <- coef(summary(model_dtialps_right))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
ctable <- cbind(ctable, "p value" = p)
ci <- confint(model_dtialps_right)
model_dtialps_right
exp(coef(model_dtialps_right))
exp(rbind(OR = coef(model_dtialps_right), ci))









# Table 2 HCP values ------------------------------------

HCP = c(0, 1)

## DCE ------------------------------------

for (i in HCP){
  col = c()
  data = subset(df_SAH, df_SAH$HCP == i)
  
  col = append(col, nrow(data))
  
  col = append(col, paste(mean(as.numeric(data$DCE_Ktrans), na.rm=TRUE), "(", sd(as.numeric(data$DCE_Ktrans), na.rm=TRUE), ")"))
  
  if (i == 0){
    tab2_DCE = as.data.frame(col)
  } else {
    tab2_DCE = cbind(tab2_DCE, col)
  }
  
}
colnames(tab2_DCE) = c("No HCP", "HCP")
rownames(tab2_DCE) = c("n","Ktrans")

## DTI_ALPS ------------------------------------

for (i in HCP){
  col = c()
  data = subset(df_SAH, df_SAH$HCP == i)
  
  col = append(col, paste(mean(as.numeric(data$DTI_ALPS_left), na.rm=TRUE), "(", sd(as.numeric(data$DTI_ALPS_left), na.rm=TRUE), ")"))
  
  col = append(col, paste(mean(as.numeric(data$DTI_ALPS_right), na.rm=TRUE), "(", sd(as.numeric(data$DTI_ALPS_right), na.rm=TRUE), ")"))
  
  col = append(col, paste(mean(as.numeric(data$DTI_ALPS_dom), na.rm=TRUE), "(", sd(as.numeric(data$DTI_ALPS_dom), na.rm=TRUE), ")"))
  
  col = append(col, paste(mean(as.numeric(data$DTI_ALPS_nondom), na.rm=TRUE),"(", sd(as.numeric(data$DTI_ALPS_nondom), na.rm=TRUE), ")"))
  
  if (i == 0){
    tab2_DTIALPS = as.data.frame(col)
  } else {
    tab2_DTIALPS = cbind(tab2_DTIALPS, col)
  }
  
}
colnames(tab2_DTIALPS) = c("No HCP", "HCP")
rownames(tab2_DTIALPS) = c("DTIALPS_left",
                           "DTIALPS_right",
                           "DTIALPS_dom",
                           "DTIALPS_nondom")

tab2 = rbind(tab2_DCE, tab2_DTIALPS)

## NODDI ------------------------------------

# remove non NODDI columns
df_SAH_NODDI = df_SAH[,50:ncol(df_SAH)]

# extract mean columns
df_SAH_NODDI = df_SAH_NODDI[,seq(2, ncol(df_SAH_NODDI),3)]

# add HCP column
df_SAH_NODDI = cbind(df_SAH_NODDI, df_SAH$HCP)
colnames(df_SAH_NODDI)[colnames(df_SAH_NODDI) == 'df_SAH$HCP'] <- 'HCP'


n = c(colnames(df_SAH_NODDI[seq(1, ncol(df_SAH_NODDI)-1)]))

HCP = c(0, 1)

for (i in HCP){
  col = c()
  
  data = subset(df_SAH_NODDI, df_SAH_NODDI$HCP == i)
  
  for (k in seq(1,ncol(df_SAH_NODDI)-1)){
    col = append(col, paste(mean(as.numeric(data[,k]), na.rm=TRUE),"(",sd(as.numeric(data[,k]), na.rm=TRUE),")"))
  }
  
  if (i == 0){
    tab2_NODDI = as.data.frame(col)
  } else {
    tab2_NODDI = cbind(tab2_NODDI, col)
  }
  
}

rownames(tab2_NODDI) = as.vector(colnames(df_SAH_NODDI)[seq(1,ncol(df_SAH_NODDI)-1)])
colnames(tab2_NODDI) =c("No HCP", "HCP")
tab2 = rbind(tab2, tab2_NODDI)

write.csv(tab2, "Table 2 HCP values.csv")


# HCP Univariate analysis ------------------------------------

HCP = c(0, 1)

## DCE ------------------------------------

df_SAH = as.data.frame(lapply(df_SAH, as.numeric))

df_SAH_noHCP = subset(df_SAH, df_SAH$HCP == 0)
df_SAH_HCP = subset(df_SAH, df_SAH$HCP == 1)

p_noHCP = shapiro.test(as.numeric(subset(df_SAH, df_SAH$HCP == 0)$DCE_Ktrans))$p.value
p_HCP = shapiro.test(as.numeric(subset(df_SAH, df_SAH$HCP == 1)$DCE_Ktrans))$p.value

tab2row = c()

if (p_noHCP < 0.05 | p_HCP < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('DCE_Ktrans ~ HCP')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('DCE_Ktrans ~ HCP')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DCE = t(as.data.frame(tab2row))

rownames(tab2stat_DCE) = c("Ktrans")
colnames(tab2stat_DCE) = c("p-value","stat")


## DTI_ALPS ------------------------------------

df_SAH = as.data.frame(lapply(df_SAH, as.numeric))

df_SAH_noHCP = subset(df_SAH, df_SAH$HCP == 0)
df_SAH_HCP = subset(df_SAH, df_SAH$HCP == 1)

p_noHCP = shapiro.test(as.numeric(subset(df_SAH, df_SAH$HCP == 0)$DTI_ALPS_left))$p.value
p_HCP = shapiro.test(as.numeric(subset(df_SAH, df_SAH$HCP == 1)$DTI_ALPS_left))$p.value

tab2row = c()

if (p_noHCP < 0.05 | p_HCP < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('DTI_ALPS_left ~ HCP')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('DTI_ALPS_left ~ HCP')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DTI = t(as.data.frame(tab2row))

p_noHCP = shapiro.test(as.numeric(subset(df_SAH, df_SAH$HCP == 0)$DTI_ALPS_right))$p.value
p_HCP = shapiro.test(as.numeric(subset(df_SAH, df_SAH$HCP == 1)$DTI_ALPS_right))$p.value

tab2row = c()

if (p_noHCP < 0.05 | p_HCP < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('DTI_ALPS_right ~ HCP')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('DTI_ALPS_right ~ HCP')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DTI = rbind(tab2stat_DTI, tab2row)

p_noHCP = shapiro.test(as.numeric(subset(df_SAH, df_SAH$HCP == 0)$DTI_ALPS_dom))$p.value
p_HCP = shapiro.test(as.numeric(subset(df_SAH, df_SAH$HCP == 1)$DTI_ALPS_dom))$p.value

tab2row = c()

if (p_noHCP < 0.05 | p_HCP < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('DTI_ALPS_dom ~ HCP')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('DTI_ALPS_dom ~ HCP')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DTI = rbind(tab2stat_DTI, tab2row)

p_noHCP = shapiro.test(as.numeric(subset(df_SAH, df_SAH$HCP == 0)$DTI_ALPS_nondom))$p.value
p_HCP = shapiro.test(as.numeric(subset(df_SAH, df_SAH$HCP == 1)$DTI_ALPS_nondom))$p.value

tab2row = c()

if (p_noHCP < 0.05 | p_HCP < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('DTI_ALPS_nondom ~ HCP')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('DTI_ALPS_nondom ~ HCP')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DTI = rbind(tab2stat_DTI, tab2row)

rownames(tab2stat_DTI) = c("DTI_ALPS_left","DTI_ALPS_right","DTI_ALPS_dom","DTI_ALPS_nondom")
colnames(tab2stat_DTI) = c("p-value","stat")

## NODDI ------------------------------------

df_SAH_NODDI_noHCP = subset(df_SAH_NODDI, df_SAH_NODDI$HCP == 0)
df_SAH_NODDI_HCP = subset(df_SAH_NODDI, df_SAH_NODDI$HCP == 1)

df_SAH_NODDI = as.data.frame(lapply(df_SAH_NODDI,as.numeric))

for (k in seq(1,ncol(df_SAH_NODDI)-1)){
  p_noHCP = shapiro.test(as.numeric(subset(df_SAH_NODDI, df_SAH_NODDI$HCP == 0)[,k]))$p.value
  p_HCP = shapiro.test(as.numeric(subset(df_SAH_NODDI, df_SAH_NODDI$HCP == 1)[,k]))$p.value
  
  tab2row = c()
  
  if (p_noHCP < 0.05 | p_HCP < 0.05){
    tab2row = append(tab2row, kruskal.test(as.formula(paste(colnames(df_SAH_NODDI)[k],'~ HCP')), data = df_SAH_NODDI)$p.value)
    tab2row = append(tab2row, "kruskal")
  } else {
    tab2row = append(tab2row, oneway.test(as.formula(paste(colnames(df_SAH_NODDI)[k],'~ HCP')), data = df_SAH_NODDI)$p.value)
    tab2row = append(tab2row, "ANOVA")
  }
  
  if (k == 1){
    tab2stat_NODDI = t(as.data.frame(tab2row))
  } else {
    tab2stat_NODDI = rbind(tab2stat_NODDI, tab2row)
  }
}

rownames(tab2stat_NODDI) = as.vector(colnames(df_SAH_NODDI)[seq(1,ncol(df_SAH_NODDI)-1)])
colnames(tab2stat_NODDI) = c("p-value","stat")
tab2stat_NODDI = as.data.frame(tab2stat_NODDI)

tab2stat = rbind(tab2stat_DCE, tab2stat_DTI, tab2stat_NODDI)

tab2stat_VOI = subset(tab2stat, tab2stat$'p-value'<0.05)

write.csv(tab2stat, "Univariate analysis HCP.csv")
write.csv(tab2stat_VOI, "Univariate analysis variables of interest HCP.csv")


# HCP Post hoc tests for variables of interest ------------------------------------
VOI_tukey = subset(tab2stat_VOI, tab2stat_VOI$stat == "ANOVA")
VOI_tukey = rownames(VOI_tukey)
VOI_dunn = subset(tab2stat_VOI, tab2stat_VOI$stat == "kruskal")
VOI_dunn = rownames(VOI_dunn)

## Tukey HSD ------------------------------------

for (i in seq(1:length(VOI_tukey))){
  Tukey = TukeyHSD(aov(as.formula(paste(VOI_tukey[i], "~ factor(HCP)")), data = df_SAH))[1]
  
  if (i == 1){
    Tukey_list = Tukey
  } else {
    Tukey_list = append(Tukey_list, Tukey)
  }
  
}

names(Tukey_list) = VOI_tukey

capture.output(Tukey_list, file = "TukeyHSD_HCP.txt")

## Dunn test ------------------------------------

for (i in seq(1:length(VOI_dunn))){
  dunn = dunnTest(as.formula(paste(VOI_dunn[i], "~ factor(HCP)")), data = df_SAH)[2]
  
  if (i == 1){
    dunn_list = dunn
  } else {
    dunn_list = append(dunn_list, dunn)
  }
  
}

names(dunn_list) = VOI_dunn

capture.output(dunn_list, file = "dunn_HCP.txt")


# HCP Box plots ------------------------------------
## DTI ALPS ------------------------------------

df_SAH$HCP_name = c(0)
df_SAH[df_SAH$HCP == 0, 'HCP_name'] = "No HCP"
df_SAH[df_SAH$HCP == 1, 'HCP_name'] = "HCP"

df_SAH$HCP_name = factor(df_SAH$HCP_name, levels = c("No HCP","HCP"))

plot_DTIALPS_right = ggplot(df_SAH, aes(x=HCP_name, y=DTI_ALPS_right_manual))
plot_DTIALPS_right = plot_DTIALPS_right + geom_boxplot() + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(method = "anova", label.y = 2.2) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)

plot_DTIALPS_dom = ggplot(df_SAH, aes(x=HCP_name, y=DTI_ALPS_dom_manual))
plot_DTIALPS_dom = plot_DTIALPS_dom + geom_boxplot() + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(method = "anova", label.y = 2.2) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)


plot_DTIALPS_left = ggplot(df_SAH, aes(x=HCP_name, y=DTI_ALPS_left_manual))
plot_DTIALPS_left = plot_DTIALPS_left + geom_boxplot() + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(method = "anova", label.y = 1.85) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)

plot_DTIALPS_nondom = ggplot(df_SAH, aes(x=HCP_name, y=DTI_ALPS_nondom_manual))
plot_DTIALPS_nondom = plot_DTIALPS_nondom + geom_boxplot() + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(method = "anova", label.y = 1.85) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)

plot_grid(plot_DTIALPS_right, plot_DTIALPS_dom, plot_DTIALPS_left, plot_DTIALPS_nondom, 
          labels = c("DTI_ALPS"),
          ncol = 2, nrow = 2)

# Multivariable analyses ------------------------------------

## Crude models  ------------------------------------



### DTI_ALPS Left ------------------------------------

model_dtialps_left = polr(outcome_name ~ DTI_ALPS_left, data = df_SAH, Hess=TRUE)
ctable <- coef(summary(model_dtialps_left))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
ctable <- cbind(ctable, "p value" = p)
ci <- confint(model_dtialps_left)
model_dtialps_left
exp(coef(model_dtialps_left))
exp(rbind(OR = coef(model_dtialps_left), ci))



### DTI_ALPS non-dom ------------------------------------

model_dtialps_nondom = polr(outcome_name ~ DTI_ALPS_nondom, data = df_SAH, Hess=TRUE)
ctable <- coef(summary(model_dtialps_nondom))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
ctable <- cbind(ctable, "p value" = p)
ci <- confint(model_dtialps_nondom)
model_dtialps_nondom
exp(coef(model_dtialps_nondom))
exp(rbind(OR = coef(model_dtialps_nondom), ci))


## Adjusted models ------------------------------------

### DTI_ALPS Left ------------------------------------

model_dtialps_left = polr(outcome_name ~ DTI_ALPS_left + Milrinone, data = df_SAH, Hess=TRUE)
ctable <- coef(summary(model_dtialps_left))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
ctable <- cbind(ctable, "p value" = p)
ci <- confint(model_dtialps_left)
model_dtialps_left
exp(coef(model_dtialps_left))
exp(rbind(OR = coef(model_dtialps_left), ci))









