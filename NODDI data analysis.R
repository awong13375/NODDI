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
df_PM['Outcome'] = c(0)

df_SAH[df_SAH$`Radiologic vasospasm` == 1, 'Outcome'] = 1
df_SAH[df_SAH$DCI == 1, 'Outcome'] = 2

df_PM[df_PM$`Radiologic vasospasm` == 1, 'Outcome'] = 1
df_PM[df_PM$DCI == 1, 'Outcome'] = 2


# Table 1 values ------------------------------------
outcomes = c(0, 1, 2)

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

  if (i == 0){
    tab1 = as.data.frame(col)
  } else {
    tab1 = cbind(tab1, col)
  }
  
}
colnames(tab1) = c("No complications", "Radiologic Vasospasm", "DCI")
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

write.csv(tab1, "Table 1 values.csv")

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
df_SAH_N = subset(df_SAH, df_SAH$Outcome == 0)
df_SAH_VS = subset(df_SAH, df_SAH$Outcome == 1)
df_SAH_DCI = subset(df_SAH, df_SAH$Outcome == 2)

p_N = shapiro.test(subset(df_SAH, df_SAH$Outcome == 0)$Age)$p.value
p_VS = shapiro.test(subset(df_SAH, df_SAH$Outcome == 1)$Age)$p.value
p_DCI = shapiro.test(subset(df_SAH, df_SAH$Outcome == 2)$Age)$p.value

if (p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
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

write.csv(tab1row, "Table 1 stats.csv")


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

outcomes = c(0, 1, 2)

## DCE ------------------------------------

for (i in outcomes){
  col = c()
  data = subset(df_SAH, df_SAH$Outcome == i)
  
  col = append(col, nrow(data))
  
  col = append(col, paste(mean(as.numeric(data$DCE_Ktrans), na.rm=TRUE), "(", sd(as.numeric(data$DCE_Ktrans), na.rm=TRUE), ")"))

  col = append(col, paste(mean(as.numeric(data$Vp), na.rm=TRUE), "(", sd(as.numeric(data$Vp), na.rm=TRUE), ")"))
  
  if (i == 0){
    tab2_DCE = as.data.frame(col)
  } else {
    tab2_DCE = cbind(tab2_DCE, col)
  }
  
}
colnames(tab2_DCE) = c("No complications", "Radiologic vasospasm","DCI")
rownames(tab2_DCE) = c("n","Ktrans", "Vp")



## DTI_ALPS ------------------------------------

for (i in outcomes){
  col = c()
  data = subset(df_SAH, df_SAH$Outcome == i)
  
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
colnames(tab2_DTIALPS) = c("No complications", "Radiologic vasospasm","DCI")
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

# add outcomes column
df_SAH_NODDI = cbind(df_SAH_NODDI, df_SAH$Outcome)
colnames(df_SAH_NODDI)[colnames(df_SAH_NODDI) == 'df_SAH$Outcome'] <- 'Outcome'


n = c(colnames(df_SAH_NODDI[seq(1, ncol(df_SAH_NODDI)-1)]))

outcomes = c(0, 1, 2)

for (i in outcomes){
  col = c()
  
  data = subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == i)
  
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
colnames(tab2_NODDI) =c("No complications", "Radiologic vasospasm", "DCI")
tab2 = rbind(tab2, tab2_NODDI)

write.csv(tab2, "Table 2 values.csv")











# Univariate analysis ------------------------------------

outcomes = c(0, 1, 2)

## DCE ------------------------------------

df_SAH = as.data.frame(lapply(df_SAH, as.numeric))

df_SAH_N = subset(df_SAH, df_SAH$Outcome == 0)
df_SAH_VS = subset(df_SAH, df_SAH$Outcome == 1)
df_SAH_DCI = subset(df_SAH, df_SAH$Outcome == 2)

p_N = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 0)$DCE_Ktrans))$p.value
p_VS = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 1)$DCE_Ktrans))$p.value
p_DCI = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 2)$DCE_Ktrans))$p.value

tab2row = c()

if (p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('DCE_Ktrans ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('DCE_Ktrans ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DCE = t(as.data.frame(tab2row))

tab2row = c()

p_N = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 0)$Vp))$p.value
p_VS = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 1)$Vp))$p.value
p_DCI = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 2)$Vp))$p.value


if (p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
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


df_SAH_N = subset(df_SAH, df_SAH$Outcome == 0)
df_SAH_VS = subset(df_SAH, df_SAH$Outcome == 1)
df_SAH_DCI = subset(df_SAH, df_SAH$Outcome == 2)

p_N = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 0)$DTI_ALPS_left))$p.value
p_VS = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 1)$DTI_ALPS_left))$p.value
p_DCI = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 2)$DTI_ALPS_left))$p.value

tab2row = c()

if (p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('DTI_ALPS_left ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('DTI_ALPS_left ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DTI = t(as.data.frame(tab2row))

p_N = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 0)$DTI_ALPS_right))$p.value
p_VS = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 1)$DTI_ALPS_right))$p.value
p_DCI = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 2)$DTI_ALPS_right))$p.value

tab2row = c()

if (p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('DTI_ALPS_right ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('DTI_ALPS_right ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DTI = rbind(tab2stat_DTI, tab2row)

p_N = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 0)$DTI_ALPS_dom))$p.value
p_VS = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 1)$DTI_ALPS_dom))$p.value
p_DCI = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 2)$DTI_ALPS_dom))$p.value

tab2row = c()

if (p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('DTI_ALPS_dom ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('DTI_ALPS_dom ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DTI = rbind(tab2stat_DTI, tab2row)

p_N = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 0)$DTI_ALPS_nondom))$p.value
p_VS = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 1)$DTI_ALPS_nondom))$p.value
p_DCI = shapiro.test(as.numeric(subset(df_SAH, df_SAH$Outcome == 2)$DTI_ALPS_nondom))$p.value

tab2row = c()

if (p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
  tab2row = append(tab2row, kruskal.test(as.formula(paste('DTI_ALPS_nondom ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "kruskal")
} else {
  tab2row = append(tab2row, oneway.test(as.formula(paste('DTI_ALPS_nondom ~ Outcome')), data = df_SAH)$p.value)
  tab2row = append(tab2row, "ANOVA")
}

tab2stat_DTI = rbind(tab2stat_DTI, tab2row)

rownames(tab2stat_DTI) = c("DTI_ALPS_left","DTI_ALPS_right","DTI_ALPS_dom","DTI_ALPS_nondom")
colnames(tab2stat_DTI) = c("p-value","stat")

## NODDI ------------------------------------

df_SAH_NODDI_N = subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == 0)
df_SAH_NODDI_VS = subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == 1)
df_SAH_NODDI_DCI = subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == 2)

df_SAH_NODDI = as.data.frame(lapply(df_SAH_NODDI,as.numeric))

for (k in seq(1,ncol(df_SAH_NODDI)-1)){
  p_N = shapiro.test(as.numeric(subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == 0)[,k]))$p.value
  p_VS = shapiro.test(as.numeric(subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == 1)[,k]))$p.value
  p_DCI = shapiro.test(as.numeric(subset(df_SAH_NODDI, df_SAH_NODDI$Outcome == 2)[,k]))$p.value
  
  tab2row = c()
  
  if (p_N < 0.05 | p_VS < 0.05 | p_DCI < 0.05){
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

write.csv(tab2stat, "Univariate analysis.csv")
write.csv(tab2stat_VOI, "Univariate analysis variables of interest.csv")


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

capture.output(Tukey_list, file = "TukeyHSD.txt")

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

capture.output(dunn_list, file = "dunn.txt")


# Box plots ------------------------------------
## DTI ALPS ------------------------------------

df_SAH$outcome_name = c(0)
df_SAH[df_SAH$Outcome == 0, 'outcome_name'] = "No Complications"
df_SAH[df_SAH$Outcome == 1, 'outcome_name'] = "Vasospasm"
df_SAH[df_SAH$Outcome == 2, 'outcome_name'] = "DCI"

df_SAH$outcome_name = factor(df_SAH$outcome_name, levels = c("No Complications","Vasospasm","DCI"))

plot_DTIALPS_right = ggplot(df_SAH, aes(x=outcome_name, y=DTI_ALPS_right_manual))
plot_DTIALPS_right = plot_DTIALPS_right + geom_boxplot() + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(method = "anova", label.y = 2.2) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)

plot_DTIALPS_dom = ggplot(df_SAH, aes(x=outcome_name, y=DTI_ALPS_dom))
plot_DTIALPS_dom = plot_DTIALPS_dom + geom_boxplot() + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(method = "anova", label.y = 2.2) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)


plot_DTIALPS_left = ggplot(df_SAH, aes(x=outcome_name, y=DTI_ALPS_left_manual))
plot_DTIALPS_left = plot_DTIALPS_left + geom_boxplot() + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(method = "anova", label.y = 1.85) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)

plot_DTIALPS_nondom = ggplot(df_SAH, aes(x=outcome_name, y=DTI_ALPS_nondom))
plot_DTIALPS_nondom = plot_DTIALPS_nondom + geom_boxplot() + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(method = "anova", label.y = 1.85) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)

plot_grid(plot_DTIALPS_right, plot_DTIALPS_dom, plot_DTIALPS_left, plot_DTIALPS_nondom, 
          labels = c("Automatic_method"),
          ncol = 2, nrow = 2)

## DCE ktrans ------------------------------------

plot_DCE = ggplot(df_SAH, aes(x=outcome_name, y=Vp))
plot_DCE = plot_DCE + geom_boxplot() + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(method = "anova", label.y = 0.17) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)



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

model_dtialps_left = polr(outcome_name ~ DTI_ALPS_left + Ethnicity, data = df_SAH, Hess=TRUE)
ctable <- coef(summary(model_dtialps_left))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
ctable <- cbind(ctable, "p value" = p)
ci <- confint(model_dtialps_left)
model_dtialps_left
exp(coef(model_dtialps_left))
exp(rbind(OR = coef(model_dtialps_left), ci))









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

plot_DTIALPS_right = ggplot(df_SAH, aes(x=HCP_name, y=DTI_ALPS_right))
plot_DTIALPS_right = plot_DTIALPS_right + geom_boxplot() + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(method = "anova", label.y = 2.2) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)

plot_DTIALPS_dom = ggplot(df_SAH, aes(x=HCP_name, y=DTI_ALPS_dom))
plot_DTIALPS_dom = plot_DTIALPS_dom + geom_boxplot() + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(method = "anova", label.y = 2.2) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)


plot_DTIALPS_left = ggplot(df_SAH, aes(x=HCP_name, y=DTI_ALPS_left))
plot_DTIALPS_left = plot_DTIALPS_left + geom_boxplot() + theme_classic() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  stat_compare_means(method = "anova", label.y = 1.85) + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.", hide.ns = TRUE)

plot_DTIALPS_nondom = ggplot(df_SAH, aes(x=HCP_name, y=DTI_ALPS_nondom))
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











#load data
mni_NDI_data = read.csv('MNI_NDI_raw_intensities.csv')
mni_FWF_data = read.csv('MNI_FWF_raw_intensities.csv')
mni_ODI_data = read.csv('MNI_ODI_raw_intensities.csv')

HO_cort_NDI_data = read.csv('HO_cort_NDI_raw_intensities.csv')
HO_cort_FWF_data = read.csv('HO_cort_FWF_raw_intensities.csv')
HO_cort_ODI_data = read.csv('HO_cort_ODI_raw_intensities.csv')

HO_subcort_NDI_data = read.csv('HO_subcort_NDI_raw_intensities.csv')
HO_subcort_FWF_data = read.csv('HO_subcort_FWF_raw_intensities.csv')
HO_subcort_ODI_data = read.csv('HO_subcort_ODI_raw_intensities.csv')

JHU_label_NDI_data = read.csv('JHU_label_NDI_raw_intensities.csv')
JHU_label_FWF_data = read.csv('JHU_label_FWF_raw_intensities.csv')
JHU_label_ODI_data = read.csv('JHU_label_ODI_raw_intensities.csv')

JHU_tract_NDI_data = read.csv('JHU_tract_NDI_raw_intensities.csv')
JHU_tract_FWF_data = read.csv('JHU_tract_FWF_raw_intensities.csv')
JHU_tract_ODI_data = read.csv('JHU_tract_ODI_raw_intensities.csv')

mni_NDI_data <- mutate_all(mni_NDI_data, function(x) as.numeric(as.character(x)))
mni_FWF_data <- mutate_all(mni_FWF_data, function(x) as.numeric(as.character(x)))
mni_ODI_data <- mutate_all(mni_ODI_data, function(x) as.numeric(as.character(x)))

HO_cort_NDI_data <- mutate_all(HO_cort_NDI_data, function(x) as.numeric(as.character(x)))
HO_cort_FWF_data <- mutate_all(HO_cort_FWF_data, function(x) as.numeric(as.character(x)))
HO_cort_ODI_data <- mutate_all(HO_cort_ODI_data, function(x) as.numeric(as.character(x)))

HO_subcort_NDI_data <- mutate_all(HO_subcort_NDI_data, function(x) as.numeric(as.character(x)))
HO_subcort_FWF_data <- mutate_all(HO_subcort_FWF_data, function(x) as.numeric(as.character(x)))
HO_subcort_ODI_data <- mutate_all(HO_subcort_ODI_data, function(x) as.numeric(as.character(x)))

JHU_label_NDI_data <- mutate_all(JHU_label_NDI_data, function(x) as.numeric(as.character(x)))
JHU_label_FWF_data <- mutate_all(JHU_label_FWF_data, function(x) as.numeric(as.character(x)))
JHU_label_ODI_data <- mutate_all(JHU_label_ODI_data, function(x) as.numeric(as.character(x)))

JHU_tract_NDI_data <- mutate_all(JHU_tract_NDI_data, function(x) as.numeric(as.character(x)))
JHU_tract_FWF_data <- mutate_all(JHU_tract_FWF_data, function(x) as.numeric(as.character(x)))
JHU_tract_ODI_data <- mutate_all(JHU_tract_ODI_data, function(x) as.numeric(as.character(x)))


mni_NDI_data[mni_NDI_data == 0] <- NA
mni_NDI_data <- mni_NDI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

mni_FWF_data[mni_FWF_data == 0] <- NA
mni_FWF_data <- mni_FWF_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

mni_ODI_data[mni_ODI_data == 0] <- NA
mni_ODI_data <- mni_ODI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))


HO_cort_NDI_data[HO_cort_NDI_data == 0] <- NA
HO_cort_NDI_data <- HO_cort_NDI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

HO_cort_FWF_data[HO_cort_FWF_data == 0] <- NA
HO_cort_FWF_data <- HO_cort_FWF_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

HO_cort_ODI_data[HO_cort_ODI_data == 0] <- NA
HO_cort_ODI_data <- HO_cort_ODI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))


HO_subcort_NDI_data[HO_subcort_NDI_data == 0] <- NA
HO_subcort_NDI_data <- HO_subcort_NDI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

HO_subcort_FWF_data[HO_subcort_FWF_data == 0] <- NA
HO_subcort_FWF_data <- HO_subcort_FWF_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

HO_subcort_ODI_data[HO_subcort_ODI_data == 0] <- NA
HO_subcort_ODI_data <- HO_subcort_ODI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))


JHU_label_NDI_data[JHU_label_NDI_data == 0] <- NA
JHU_label_NDI_data <- JHU_label_NDI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

JHU_label_FWF_data[JHU_label_FWF_data == 0] <- NA
JHU_label_FWF_data <- JHU_label_FWF_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

JHU_label_ODI_data[JHU_label_ODI_data == 0] <- NA
JHU_label_ODI_data <- JHU_label_ODI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))


JHU_tract_NDI_data[JHU_tract_NDI_data == 0] <- NA
JHU_tract_NDI_data <- JHU_tract_NDI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

JHU_tract_FWF_data[JHU_tract_FWF_data == 0] <- NA
JHU_tract_FWF_data <- JHU_tract_FWF_data %>% mutate_all(~ifelse(is.nan(.), NA, .))

JHU_tract_ODI_data[JHU_tract_ODI_data == 0] <- NA
JHU_tract_ODI_data <- JHU_tract_ODI_data %>% mutate_all(~ifelse(is.nan(.), NA, .))


# MNI

col_list = c("Caudate","Cerebellum","Frontal_Lobe","Insula","Occipital_Lobe","Parietal_Lobe","Temporal_Lobe","Thalamus")
for (n in 1:length(col_list)){
    col = col_list[n]
    print (paste("Max:", col, "=", max(mni_NDI_data[col], na.rm=T)))
    print (paste("Min:", col, "=", min(mni_NDI_data[col], na.rm=T)))
    
  
}
for (n in 1:length(col_list)){
  col = col_list[n]
  print (paste("Max:", col, "=", max(mni_FWF_data[col], na.rm=T)))
  print (paste("Min:", col, "=", min(mni_FWF_data[col], na.rm=T)))
  
  
}
for (n in 1:length(col_list)){
  col = col_list[n]
  print (paste("Max:", col, "=", max(mni_ODI_data[col], na.rm=T)))
  print (paste("Min:", col, "=", min(mni_ODI_data[col], na.rm=T)))
  
  
}


#create plots
#NDI
Caudate = gghistogram(mni_NDI_data$Caudate,
                    main = "Caudate",
                    xlab = "Intensity"
                   ) + ggeasy::easy_center_title()
 
Cerebellum = gghistogram(mni_NDI_data$Cerebellum,
                       main = "Cerebellum",
                       xlab = "Intensity"
                      ) + ggeasy::easy_center_title()

Frontal_Lobe = gghistogram(mni_NDI_data$Frontal_Lobe,
                         main = "Frontal_Lobe",
                         xlab = "Intensity"
                        ) + ggeasy::easy_center_title()

Insula = gghistogram(mni_NDI_data$Insula,
                   main = "Insula",
                   xlab = "Intensity"
                  ) + ggeasy::easy_center_title()

Occipital_Lobe = gghistogram(mni_NDI_data$Occipital_Lobe,
                           main = "Occipital_Lobe",
                           xlab = "Intensity"
                          ) + ggeasy::easy_center_title()

Parietal_Lobe = gghistogram(mni_NDI_data$Parietal_Lobe,
                          main = "Parietal_Lobe",
                          xlab = "Intensity"
                         ) + ggeasy::easy_center_title()

Temporal_Lobe = gghistogram(mni_NDI_data$Temporal_Lobe,
                          main = "Temporal_Lobe",
                          xlab = "Intensity"
                         ) + ggeasy::easy_center_title()

Thalamus = gghistogram(mni_NDI_data$Thalamus,
                     main = "Thalamus",
                     xlab = "Intensity"
                    ) + ggeasy::easy_center_title()


plot_grid(Caudate, Cerebellum, Frontal_Lobe, Insula, Occipital_Lobe, Parietal_Lobe, Temporal_Lobe, Thalamus, 
          labels = c("MNI_NDI"),
          ncol = 4, nrow = 2)

#FWF
Caudate = gghistogram(mni_FWF_data$Caudate,
                    main = "caudate",
                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebellum = gghistogram(mni_FWF_data$Cerebellum,
                       main = "Cerebellum",
                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Lobe = gghistogram(mni_FWF_data$Frontal_Lobe,
                         main = "Frontal_Lobe",
                         xlab = "Intensit y"
) + ggeasy::easy_center_title()

Insula = gghistogram(mni_FWF_data$Insula,
                   main = "Insula",
                   xlab = "Intensity",
                   xlim = c(0,0.5)
) + ggeasy::easy_center_title()

Occipital_Lobe = gghistogram(mni_FWF_data$Occipital_Lobe,
                           main = "Occipital_Lobe",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Parietal_Lobe = gghistogram(mni_FWF_data$Parietal_Lobe,
                          main = "Parietal_Lobe",
                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Lobe = gghistogram(mni_FWF_data$Temporal_Lobe,
                          main = "Temporal_Lobe",
                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Thalamus = gghistogram(mni_FWF_data$Thalamus,
                     main = "Thalamus",
                     xlab = "Intensity"
) + ggeasy::easy_center_title()


plot_grid(Caudate, Cerebellum, Frontal_Lobe, Insula, Occipital_Lobe, Parietal_Lobe, Temporal_Lobe, Thalamus, 
          labels = c("MNI_FWF"),
          ncol = 4, nrow = 2)

#ODI
Caudate = gghistogram(mni_ODI_data$Caudate,
                    main = "caudate",
                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebellum = gghistogram(mni_ODI_data$Cerebellum,
                       main = "Cerebellum",
                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Lobe = gghistogram(mni_ODI_data$Frontal_Lobe,
                         main = "Frontal_Lobe",
                         xlab = "Intensit y"
) + ggeasy::easy_center_title()

Insula = gghistogram(mni_ODI_data$Insula,
                   main = "Insula",
                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Lobe = gghistogram(mni_ODI_data$Occipital_Lobe,
                           main = "Occipital_Lobe",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Parietal_Lobe = gghistogram(mni_ODI_data$Parietal_Lobe,
                          main = "Parietal_Lobe",
                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Lobe = gghistogram(mni_ODI_data$Temporal_Lobe,
                          main = "Temporal_Lobe",
                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Thalamus = gghistogram(mni_ODI_data$Thalamus,
                     main = "Thalamus",
                     xlab = "Intensity",
                     xlim = c(0,0.8)
) + ggeasy::easy_center_title()


plot_grid(Caudate, Cerebellum, Frontal_Lobe, Insula, Occipital_Lobe, Parietal_Lobe, Temporal_Lobe, Thalamus, 
          labels = c("MNI_ODI"),
          ncol = 4, nrow = 2)

#HO Cort 

col_list = c("Frontal_Pole","Parahippocampal_Gyrus_posterior_division","Parahippocampal_Gyrus_anterior_division","Frontal_Orbital_Cortex","Cuneal_Cortex","Precuneous_Cortex",
             "Cingulate_Gyrus_posterior_division","Cingulate_Gyrus_anterior_division","Paracingulate_Gyrus","Subcallosal_Cortex","Juxtapositional_Lobule_Cortex_formerly_Supplementa","Frontal_Medial_Cortex",
             "Intracalcarine_Cortex","Lateral_Occipital_Cortex_inferior_division","Lateral_Occipital_Cortex_superior_division","Angular_Gyrus","Inferior_Temporal_Gyrus_temporooccipital_part",
             "Inferior_Temporal_Gyrus_posterior_division","Inferior_Temporal_Gyrus_anterior_division","Middle_Temporal_Gyrus_temporooccipital_part","Middle_Temporal_Gyrus_posterior_division",
             "Middle_Temporal_Gyrus_anterior_division","Superior_Temporal_Gyrus_posterior_division","Superior_Temporal_Gyrus_anterior_division","Temporal_Pole","Precentral_Gyrus","Inferior_Frontal_Gyrus_pars_opercularis",
             "Inferior_Frontal_Gyrus_pars_triangularis","Middle_Frontal_Gyrus","Superior_Frontal_Gyrus","Insular_Cortex","Occipital_Pole","Supracalcarine_Cortex","Planum_Temporale","Heschls_Gyrus_includes_H1_and_H2",
             "Planum_Polare","Parietal_Operculum_Cortex","Central_Opercular_Cortex","Frontal_Operculum_Cortex","Occipital_Fusiform_Gyrus","Temporal_Occipital_Fusiform_Cortex","Temporal_Fusiform_Cortex_posterior_division",
             "Temporal_Fusiform_Cortex_anterior_division","Lingual_Gyrus")
#NDI
Frontal_Pole = gghistogram(HO_cort_NDI_data$Frontal_Pole,
                      main = "Frontal_Pole",
                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Parahippocampal_Gyrus_posterior_division = gghistogram(HO_cort_NDI_data$Parahippocampal_Gyrus_posterior_division,
                         main = "Parahippocampal_Gyrus_posterior_division",
                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Parahippocampal_Gyrus_anterior_division = gghistogram(HO_cort_NDI_data$Parahippocampal_Gyrus_anterior_division,
                           main = "Parahippocampal_Gyrus_anterior_division",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Orbital_Cortex = gghistogram(HO_cort_NDI_data$Frontal_Orbital_Cortex,
                     main = "Frontal_Orbital_Cortex",
                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cuneal_Cortex = gghistogram(HO_cort_NDI_data$Cuneal_Cortex,
                             main = "Cuneal_Cortex",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Precuneous_Cortex = gghistogram(HO_cort_NDI_data$Precuneous_Cortex,
                            main = "Precuneous_Cortex",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulate_Gyrus_posterior_division = gghistogram(HO_cort_NDI_data$Cingulate_Gyrus_posterior_division,
                            main = "Cingulate_Gyrus_posterior_division",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulate_Gyrus_anterior_division = gghistogram(HO_cort_NDI_data$Cingulate_Gyrus_anterior_division,
                       main = "Cingulate_Gyrus_anterior_division",
                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Paracingulate_Gyrus = gghistogram(HO_cort_NDI_data$Paracingulate_Gyrus,
                                                main = "Paracingulate_Gyrus",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Subcallosal_Cortex = gghistogram(HO_cort_NDI_data$Subcallosal_Cortex,
                                                main = "Subcallosal_Cortex",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Juxtapositional_Lobule_Cortex_formerly_Supplementa = gghistogram(HO_cort_NDI_data$Juxtapositional_Lobule_Cortex_formerly_Supplementa,
                                                main = "Juxtapositional_Lobule_Cortex_formerly_Supplementa",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Medial_Cortex = gghistogram(HO_cort_NDI_data$Frontal_Medial_Cortex,
                                                main = "Frontal_Medial_Cortex",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Intracalcarine_Cortex = gghistogram(HO_cort_NDI_data$Intracalcarine_Cortex,
                                                main = "Intracalcarine_Cortex",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Lateral_Occipital_Cortex_inferior_division = gghistogram(HO_cort_NDI_data$Lateral_Occipital_Cortex_inferior_division,
                                                main = "Lateral_Occipital_Cortex_inferior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Lateral_Occipital_Cortex_superior_division = gghistogram(HO_cort_NDI_data$Lateral_Occipital_Cortex_superior_division,
                                                main = "Lateral_Occipital_Cortex_superior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Angular_Gyrus = gghistogram(HO_cort_NDI_data$Angular_Gyrus,
                                                main = "Angular_Gyrus",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Supramarginal_Gyrus_posterior_division = gghistogram(HO_cort_NDI_data$Supramarginal_Gyrus_posterior_division,
                                                main = "Supramarginal_Gyrus_posterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Supramarginal_Gyrus_anterior_division = gghistogram(HO_cort_NDI_data$Supramarginal_Gyrus_anterior_division,
                                                main = "Supramarginal_Gyrus_anterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Postcentral_Gyrus = gghistogram(HO_cort_NDI_data$Postcentral_Gyrus,
                                                main = "Postcentral_Gyrus",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_temporooccipital_part = gghistogram(HO_cort_NDI_data$Inferior_Temporal_Gyrus_temporooccipital_part,
                                                main = "Inferior_Temporal_Gyrus_temporooccipital_part",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_NDI_data$Inferior_Temporal_Gyrus_posterior_division,
                                                main = "Inferior_Temporal_Gyrus_posterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_NDI_data$Inferior_Temporal_Gyrus_anterior_division,
                                                main = "Inferior_Temporal_Gyrus_anterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_temporooccipital_part = gghistogram(HO_cort_NDI_data$Middle_Temporal_Gyrus_temporooccipital_part,
                                                main = "Middle_Temporal_Gyrus_temporooccipital_part",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_NDI_data$Middle_Temporal_Gyrus_posterior_division,
                                                main = "Cingulate_Gyrus_anterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_NDI_data$Middle_Temporal_Gyrus_anterior_division,
                                                main = "Middle_Temporal_Gyrus_anterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_NDI_data$Superior_Temporal_Gyrus_posterior_division,
                                                main = "Superior_Temporal_Gyrus_posterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_NDI_data$Superior_Temporal_Gyrus_anterior_division,
                                                main = "Superior_Temporal_Gyrus_anterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Pole = gghistogram(HO_cort_NDI_data$Temporal_Pole,
                                                main = "Temporal_Pole",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Precentral_Gyrus = gghistogram(HO_cort_NDI_data$Precentral_Gyrus,
                                                main = "Precentral_Gyrus",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Frontal_Gyrus_pars_opercularis = gghistogram(HO_cort_NDI_data$Inferior_Frontal_Gyrus_pars_opercularis,
                                                main = "Inferior_Frontal_Gyrus_pars_opercularis",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Frontal_Gyrus_pars_triangularis = gghistogram(HO_cort_NDI_data$Inferior_Frontal_Gyrus_pars_triangularis,
                                                main = "Inferior_Frontal_Gyrus_pars_triangularis",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Frontal_Gyrus = gghistogram(HO_cort_NDI_data$Middle_Frontal_Gyrus,
                                                main = "Middle_Frontal_Gyrus",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Frontal_Gyrus = gghistogram(HO_cort_NDI_data$Superior_Frontal_Gyrus,
                                                main = "Superior_Frontal_Gyrus",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Insular_Cortex = gghistogram(HO_cort_NDI_data$Insular_Cortex,
                                     main = "Insular_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Pole = gghistogram(HO_cort_NDI_data$Occipital_Pole,
                                     main = "Occipital_Pole",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Supracalcarine_Cortex = gghistogram(HO_cort_NDI_data$Supracalcarine_Cortex,
                                     main = "Supracalcarine_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Planum_Temporale = gghistogram(HO_cort_NDI_data$Planum_Temporale,
                                     main = "Planum_Temporale",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Heschls_Gyrus_includes_H1_and_H2 = gghistogram(HO_cort_NDI_data$Heschls_Gyrus_includes_H1_and_H2,
                                     main = "Heschls_Gyrus_includes_H1_and_H2",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Planum_Polare = gghistogram(HO_cort_NDI_data$Planum_Polare,
                                     main = "Planum_Polare",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Parietal_Operculum_Cortex = gghistogram(HO_cort_NDI_data$Parietal_Operculum_Cortex,
                                     main = "Parietal_Operculum_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Central_Opercular_Cortex = gghistogram(HO_cort_NDI_data$Central_Opercular_Cortex,
                                     main = "Central_Opercular_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Operculum_Cortex = gghistogram(HO_cort_NDI_data$Frontal_Operculum_Cortex,
                                     main = "Frontal_Operculum_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Fusiform_Gyrus = gghistogram(HO_cort_NDI_data$Occipital_Fusiform_Gyrus,
                                     main = "Occipital_Fusiform_Gyrus",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Occipital_Fusiform_Cortex = gghistogram(HO_cort_NDI_data$Temporal_Occipital_Fusiform_Cortex,
                                     main = "Temporal_Occipital_Fusiform_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Fusiform_Cortex_posterior_division = gghistogram(HO_cort_NDI_data$Temporal_Fusiform_Cortex_posterior_division,
                                     main = "Temporal_Fusiform_Cortex_posterior_division",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Fusiform_Cortex_anterior_division = gghistogram(HO_cort_NDI_data$Temporal_Fusiform_Cortex_anterior_division,
                                     main = "Temporal_Fusiform_Cortex_anterior_division",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Lingual_Gyrus = gghistogram(HO_cort_NDI_data$Lingual_Gyrus,
                                     main = "Lingual_Gyrus",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Frontal_Pole,Parahippocampal_Gyrus_posterior_division,Parahippocampal_Gyrus_anterior_division,Frontal_Orbital_Cortex,Cuneal_Cortex,Precuneous_Cortex,
          Cingulate_Gyrus_posterior_division,Cingulate_Gyrus_anterior_division,Paracingulate_Gyrus,Subcallosal_Cortex,Juxtapositional_Lobule_Cortex_formerly_Supplementa,Frontal_Medial_Cortex,
          Intracalcarine_Cortex,Lateral_Occipital_Cortex_inferior_division,Lateral_Occipital_Cortex_superior_division,Angular_Gyrus,Inferior_Temporal_Gyrus_temporooccipital_part,
          Inferior_Temporal_Gyrus_posterior_division,Inferior_Temporal_Gyrus_anterior_division,Middle_Temporal_Gyrus_temporooccipital_part,Middle_Temporal_Gyrus_posterior_division,
          Middle_Temporal_Gyrus_anterior_division,Superior_Temporal_Gyrus_posterior_division,Superior_Temporal_Gyrus_anterior_division,Temporal_Pole,Precentral_Gyrus,Inferior_Frontal_Gyrus_pars_opercularis,
          Inferior_Frontal_Gyrus_pars_triangularis,Middle_Frontal_Gyrus,Superior_Frontal_Gyrus,Insular_Cortex,Occipital_Pole,Supracalcarine_Cortex,Planum_Temporale,Heschls_Gyrus_includes_H1_and_H2,
          Planum_Polare,Parietal_Operculum_Cortex,Central_Opercular_Cortex,Frontal_Operculum_Cortex,Occipital_Fusiform_Gyrus,Temporal_Occipital_Fusiform_Cortex,Temporal_Fusiform_Cortex_posterior_division,
          Temporal_Fusiform_Cortex_anterior_division,Lingual_Gyrus, 
          labels = c("HO_cort_NDI"),
          ncol = 4, nrow = 12)


#FWF
Frontal_Pole = gghistogram(HO_cort_FWF_data$Frontal_Pole,
                           main = "Frontal_Pole",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Parahippocampal_Gyrus_posterior_division = gghistogram(HO_cort_FWF_data$Parahippocampal_Gyrus_posterior_division,
                                                       main = "Parahippocampal_Gyrus_posterior_division",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Parahippocampal_Gyrus_anterior_division = gghistogram(HO_cort_FWF_data$Parahippocampal_Gyrus_anterior_division,
                                                      main = "Parahippocampal_Gyrus_anterior_division",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Orbital_Cortex = gghistogram(HO_cort_FWF_data$Frontal_Orbital_Cortex,
                                     main = "Frontal_Orbital_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cuneal_Cortex = gghistogram(HO_cort_FWF_data$Cuneal_Cortex,
                            main = "Cuneal_Cortex",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Precuneous_Cortex = gghistogram(HO_cort_FWF_data$Precuneous_Cortex,
                                main = "Precuneous_Cortex",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulate_Gyrus_posterior_division = gghistogram(HO_cort_FWF_data$Cingulate_Gyrus_posterior_division,
                                                 main = "Cingulate_Gyrus_posterior_division",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulate_Gyrus_anterior_division = gghistogram(HO_cort_FWF_data$Cingulate_Gyrus_anterior_division,
                                                main = "Cingulate_Gyrus_anterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Paracingulate_Gyrus = gghistogram(HO_cort_FWF_data$Paracingulate_Gyrus,
                                  main = "Paracingulate_Gyrus",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Subcallosal_Cortex = gghistogram(HO_cort_FWF_data$Subcallosal_Cortex,
                                 main = "Subcallosal_Cortex",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Juxtapositional_Lobule_Cortex_formerly_Supplementa = gghistogram(HO_cort_FWF_data$Juxtapositional_Lobule_Cortex_formerly_Supplementa,
                                                                 main = "Juxtapositional_Lobule_Cortex_formerly_Supplementa",
                                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Medial_Cortex = gghistogram(HO_cort_FWF_data$Frontal_Medial_Cortex,
                                    main = "Frontal_Medial_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Intracalcarine_Cortex = gghistogram(HO_cort_FWF_data$Intracalcarine_Cortex,
                                    main = "Intracalcarine_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Lateral_Occipital_Cortex_inferior_division = gghistogram(HO_cort_FWF_data$Lateral_Occipital_Cortex_inferior_division,
                                                         main = "Lateral_Occipital_Cortex_inferior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Lateral_Occipital_Cortex_superior_division = gghistogram(HO_cort_FWF_data$Lateral_Occipital_Cortex_superior_division,
                                                         main = "Lateral_Occipital_Cortex_superior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Angular_Gyrus = gghistogram(HO_cort_FWF_data$Angular_Gyrus,
                            main = "Angular_Gyrus",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Supramarginal_Gyrus_posterior_division = gghistogram(HO_cort_FWF_data$Supramarginal_Gyrus_posterior_division,
                                                     main = "Supramarginal_Gyrus_posterior_division",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Supramarginal_Gyrus_anterior_division = gghistogram(HO_cort_FWF_data$Supramarginal_Gyrus_anterior_division,
                                                    main = "Supramarginal_Gyrus_anterior_division",
                                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Postcentral_Gyrus = gghistogram(HO_cort_FWF_data$Postcentral_Gyrus,
                                main = "Postcentral_Gyrus",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_temporooccipital_part = gghistogram(HO_cort_FWF_data$Inferior_Temporal_Gyrus_temporooccipital_part,
                                                            main = "Inferior_Temporal_Gyrus_temporooccipital_part",
                                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_FWF_data$Inferior_Temporal_Gyrus_posterior_division,
                                                         main = "Inferior_Temporal_Gyrus_posterior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_FWF_data$Inferior_Temporal_Gyrus_anterior_division,
                                                        main = "Inferior_Temporal_Gyrus_anterior_division",
                                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_temporooccipital_part = gghistogram(HO_cort_FWF_data$Middle_Temporal_Gyrus_temporooccipital_part,
                                                          main = "Middle_Temporal_Gyrus_temporooccipital_part",
                                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_FWF_data$Middle_Temporal_Gyrus_posterior_division,
                                                       main = "Cingulate_Gyrus_anterior_division",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_FWF_data$Middle_Temporal_Gyrus_anterior_division,
                                                      main = "Middle_Temporal_Gyrus_anterior_division",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_FWF_data$Superior_Temporal_Gyrus_posterior_division,
                                                         main = "Superior_Temporal_Gyrus_posterior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_FWF_data$Superior_Temporal_Gyrus_anterior_division,
                                                        main = "Superior_Temporal_Gyrus_anterior_division",
                                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Pole = gghistogram(HO_cort_FWF_data$Temporal_Pole,
                            main = "Temporal_Pole",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Precentral_Gyrus = gghistogram(HO_cort_FWF_data$Precentral_Gyrus,
                               main = "Precentral_Gyrus",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Frontal_Gyrus_pars_opercularis = gghistogram(HO_cort_FWF_data$Inferior_Frontal_Gyrus_pars_opercularis,
                                                      main = "Inferior_Frontal_Gyrus_pars_opercularis",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Frontal_Gyrus_pars_triangularis = gghistogram(HO_cort_FWF_data$Inferior_Frontal_Gyrus_pars_triangularis,
                                                       main = "Inferior_Frontal_Gyrus_pars_triangularis",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Frontal_Gyrus = gghistogram(HO_cort_FWF_data$Middle_Frontal_Gyrus,
                                   main = "Middle_Frontal_Gyrus",
                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Frontal_Gyrus = gghistogram(HO_cort_FWF_data$Superior_Frontal_Gyrus,
                                     main = "Superior_Frontal_Gyrus",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Insular_Cortex = gghistogram(HO_cort_FWF_data$Insular_Cortex,
                             main = "Insular_Cortex",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Pole = gghistogram(HO_cort_FWF_data$Occipital_Pole,
                             main = "Occipital_Pole",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Supracalcarine_Cortex = gghistogram(HO_cort_FWF_data$Supracalcarine_Cortex,
                                    main = "Supracalcarine_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Planum_Temporale = gghistogram(HO_cort_FWF_data$Planum_Temporale,
                               main = "Planum_Temporale",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Heschls_Gyrus_includes_H1_and_H2 = gghistogram(HO_cort_FWF_data$Heschls_Gyrus_includes_H1_and_H2,
                                               main = "Heschls_Gyrus_includes_H1_and_H2",
                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Planum_Polare = gghistogram(HO_cort_FWF_data$Planum_Polare,
                            main = "Planum_Polare",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Parietal_Operculum_Cortex = gghistogram(HO_cort_FWF_data$Parietal_Operculum_Cortex,
                                        main = "Parietal_Operculum_Cortex",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Central_Opercular_Cortex = gghistogram(HO_cort_FWF_data$Central_Opercular_Cortex,
                                       main = "Central_Opercular_Cortex",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Operculum_Cortex = gghistogram(HO_cort_FWF_data$Frontal_Operculum_Cortex,
                                       main = "Frontal_Operculum_Cortex",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Fusiform_Gyrus = gghistogram(HO_cort_FWF_data$Occipital_Fusiform_Gyrus,
                                       main = "Occipital_Fusiform_Gyrus",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Occipital_Fusiform_Cortex = gghistogram(HO_cort_FWF_data$Temporal_Occipital_Fusiform_Cortex,
                                                 main = "Temporal_Occipital_Fusiform_Cortex",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Fusiform_Cortex_posterior_division = gghistogram(HO_cort_FWF_data$Temporal_Fusiform_Cortex_posterior_division,
                                                          main = "Temporal_Fusiform_Cortex_posterior_division",
                                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Fusiform_Cortex_anterior_division = gghistogram(HO_cort_FWF_data$Temporal_Fusiform_Cortex_anterior_division,
                                                         main = "Temporal_Fusiform_Cortex_anterior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Lingual_Gyrus = gghistogram(HO_cort_FWF_data$Lingual_Gyrus,
                            main = "Lingual_Gyrus",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Frontal_Pole,Parahippocampal_Gyrus_posterior_division,Parahippocampal_Gyrus_anterior_division,Frontal_Orbital_Cortex,Cuneal_Cortex,Precuneous_Cortex,
          Cingulate_Gyrus_posterior_division,Cingulate_Gyrus_anterior_division,Paracingulate_Gyrus,Subcallosal_Cortex,Juxtapositional_Lobule_Cortex_formerly_Supplementa,Frontal_Medial_Cortex,
          Intracalcarine_Cortex,Lateral_Occipital_Cortex_inferior_division,Lateral_Occipital_Cortex_superior_division,Angular_Gyrus,Inferior_Temporal_Gyrus_temporooccipital_part,
          Inferior_Temporal_Gyrus_posterior_division,Inferior_Temporal_Gyrus_anterior_division,Middle_Temporal_Gyrus_temporooccipital_part,Middle_Temporal_Gyrus_posterior_division,
          Middle_Temporal_Gyrus_anterior_division,Superior_Temporal_Gyrus_posterior_division,Superior_Temporal_Gyrus_anterior_division,Temporal_Pole,Precentral_Gyrus,Inferior_Frontal_Gyrus_pars_opercularis,
          Inferior_Frontal_Gyrus_pars_triangularis,Middle_Frontal_Gyrus,Superior_Frontal_Gyrus,Insular_Cortex,Occipital_Pole,Supracalcarine_Cortex,Planum_Temporale,Heschls_Gyrus_includes_H1_and_H2,
          Planum_Polare,Parietal_Operculum_Cortex,Central_Opercular_Cortex,Frontal_Operculum_Cortex,Occipital_Fusiform_Gyrus,Temporal_Occipital_Fusiform_Cortex,Temporal_Fusiform_Cortex_posterior_division,
          Temporal_Fusiform_Cortex_anterior_division,Lingual_Gyrus, 
          labels = c("HO_cort_FWF"),
          ncol = 4, nrow = 12)



#ODI

Frontal_Pole = gghistogram(HO_cort_ODI_data$Frontal_Pole,
                           main = "Frontal_Pole",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Parahippocampal_Gyrus_posterior_division = gghistogram(HO_cort_ODI_data$Parahippocampal_Gyrus_posterior_division,
                                                       main = "Parahippocampal_Gyrus_posterior_division",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Parahippocampal_Gyrus_anterior_division = gghistogram(HO_cort_ODI_data$Parahippocampal_Gyrus_anterior_division,
                                                      main = "Parahippocampal_Gyrus_anterior_division",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Orbital_Cortex = gghistogram(HO_cort_ODI_data$Frontal_Orbital_Cortex,
                                     main = "Frontal_Orbital_Cortex",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cuneal_Cortex = gghistogram(HO_cort_ODI_data$Cuneal_Cortex,
                            main = "Cuneal_Cortex",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Precuneous_Cortex = gghistogram(HO_cort_ODI_data$Precuneous_Cortex,
                                main = "Precuneous_Cortex",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulate_Gyrus_posterior_division = gghistogram(HO_cort_ODI_data$Cingulate_Gyrus_posterior_division,
                                                 main = "Cingulate_Gyrus_posterior_division",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulate_Gyrus_anterior_division = gghistogram(HO_cort_ODI_data$Cingulate_Gyrus_anterior_division,
                                                main = "Cingulate_Gyrus_anterior_division",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Paracingulate_Gyrus = gghistogram(HO_cort_ODI_data$Paracingulate_Gyrus,
                                  main = "Paracingulate_Gyrus",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Subcallosal_Cortex = gghistogram(HO_cort_ODI_data$Subcallosal_Cortex,
                                 main = "Subcallosal_Cortex",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Juxtapositional_Lobule_Cortex_formerly_Supplementa = gghistogram(HO_cort_ODI_data$Juxtapositional_Lobule_Cortex_formerly_Supplementa,
                                                                 main = "Juxtapositional_Lobule_Cortex_formerly_Supplementa",
                                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Medial_Cortex = gghistogram(HO_cort_ODI_data$Frontal_Medial_Cortex,
                                    main = "Frontal_Medial_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Intracalcarine_Cortex = gghistogram(HO_cort_ODI_data$Intracalcarine_Cortex,
                                    main = "Intracalcarine_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Lateral_Occipital_Cortex_inferior_division = gghistogram(HO_cort_ODI_data$Lateral_Occipital_Cortex_inferior_division,
                                                         main = "Lateral_Occipital_Cortex_inferior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Lateral_Occipital_Cortex_superior_division = gghistogram(HO_cort_ODI_data$Lateral_Occipital_Cortex_superior_division,
                                                         main = "Lateral_Occipital_Cortex_superior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Angular_Gyrus = gghistogram(HO_cort_ODI_data$Angular_Gyrus,
                            main = "Angular_Gyrus",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Supramarginal_Gyrus_posterior_division = gghistogram(HO_cort_ODI_data$Supramarginal_Gyrus_posterior_division,
                                                     main = "Supramarginal_Gyrus_posterior_division",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Supramarginal_Gyrus_anterior_division = gghistogram(HO_cort_ODI_data$Supramarginal_Gyrus_anterior_division,
                                                    main = "Supramarginal_Gyrus_anterior_division",
                                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Postcentral_Gyrus = gghistogram(HO_cort_ODI_data$Postcentral_Gyrus,
                                main = "Postcentral_Gyrus",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_temporooccipital_part = gghistogram(HO_cort_ODI_data$Inferior_Temporal_Gyrus_temporooccipital_part,
                                                            main = "Inferior_Temporal_Gyrus_temporooccipital_part",
                                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_ODI_data$Inferior_Temporal_Gyrus_posterior_division,
                                                         main = "Inferior_Temporal_Gyrus_posterior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_ODI_data$Inferior_Temporal_Gyrus_anterior_division,
                                                        main = "Inferior_Temporal_Gyrus_anterior_division",
                                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_temporooccipital_part = gghistogram(HO_cort_ODI_data$Middle_Temporal_Gyrus_temporooccipital_part,
                                                          main = "Middle_Temporal_Gyrus_temporooccipital_part",
                                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_ODI_data$Middle_Temporal_Gyrus_posterior_division,
                                                       main = "Cingulate_Gyrus_anterior_division",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_ODI_data$Middle_Temporal_Gyrus_anterior_division,
                                                      main = "Middle_Temporal_Gyrus_anterior_division",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Temporal_Gyrus_posterior_division = gghistogram(HO_cort_ODI_data$Superior_Temporal_Gyrus_posterior_division,
                                                         main = "Superior_Temporal_Gyrus_posterior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Temporal_Gyrus_anterior_division = gghistogram(HO_cort_ODI_data$Superior_Temporal_Gyrus_anterior_division,
                                                        main = "Superior_Temporal_Gyrus_anterior_division",
                                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Pole = gghistogram(HO_cort_ODI_data$Temporal_Pole,
                            main = "Temporal_Pole",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Precentral_Gyrus = gghistogram(HO_cort_ODI_data$Precentral_Gyrus,
                               main = "Precentral_Gyrus",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Frontal_Gyrus_pars_opercularis = gghistogram(HO_cort_ODI_data$Inferior_Frontal_Gyrus_pars_opercularis,
                                                      main = "Inferior_Frontal_Gyrus_pars_opercularis",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_Frontal_Gyrus_pars_triangularis = gghistogram(HO_cort_ODI_data$Inferior_Frontal_Gyrus_pars_triangularis,
                                                       main = "Inferior_Frontal_Gyrus_pars_triangularis",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_Frontal_Gyrus = gghistogram(HO_cort_ODI_data$Middle_Frontal_Gyrus,
                                   main = "Middle_Frontal_Gyrus",
                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_Frontal_Gyrus = gghistogram(HO_cort_ODI_data$Superior_Frontal_Gyrus,
                                     main = "Superior_Frontal_Gyrus",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Insular_Cortex = gghistogram(HO_cort_ODI_data$Insular_Cortex,
                             main = "Insular_Cortex",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Pole = gghistogram(HO_cort_ODI_data$Occipital_Pole,
                             main = "Occipital_Pole",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Supracalcarine_Cortex = gghistogram(HO_cort_ODI_data$Supracalcarine_Cortex,
                                    main = "Supracalcarine_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Planum_Temporale = gghistogram(HO_cort_ODI_data$Planum_Temporale,
                               main = "Planum_Temporale",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Heschls_Gyrus_includes_H1_and_H2 = gghistogram(HO_cort_ODI_data$Heschls_Gyrus_includes_H1_and_H2,
                                               main = "Heschls_Gyrus_includes_H1_and_H2",
                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Planum_Polare = gghistogram(HO_cort_ODI_data$Planum_Polare,
                            main = "Planum_Polare",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Parietal_Operculum_Cortex = gghistogram(HO_cort_ODI_data$Parietal_Operculum_Cortex,
                                        main = "Parietal_Operculum_Cortex",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Central_Opercular_Cortex = gghistogram(HO_cort_ODI_data$Central_Opercular_Cortex,
                                       main = "Central_Opercular_Cortex",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Frontal_Operculum_Cortex = gghistogram(HO_cort_ODI_data$Frontal_Operculum_Cortex,
                                       main = "Frontal_Operculum_Cortex",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Occipital_Fusiform_Gyrus = gghistogram(HO_cort_ODI_data$Occipital_Fusiform_Gyrus,
                                       main = "Occipital_Fusiform_Gyrus",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Occipital_Fusiform_Cortex = gghistogram(HO_cort_ODI_data$Temporal_Occipital_Fusiform_Cortex,
                                                 main = "Temporal_Occipital_Fusiform_Cortex",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Fusiform_Cortex_posterior_division = gghistogram(HO_cort_ODI_data$Temporal_Fusiform_Cortex_posterior_division,
                                                          main = "Temporal_Fusiform_Cortex_posterior_division",
                                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Temporal_Fusiform_Cortex_anterior_division = gghistogram(HO_cort_ODI_data$Temporal_Fusiform_Cortex_anterior_division,
                                                         main = "Temporal_Fusiform_Cortex_anterior_division",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Lingual_Gyrus = gghistogram(HO_cort_ODI_data$Lingual_Gyrus,
                            main = "Lingual_Gyrus",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Frontal_Pole,Parahippocampal_Gyrus_posterior_division,Parahippocampal_Gyrus_anterior_division,Frontal_Orbital_Cortex,Cuneal_Cortex,Precuneous_Cortex,
          Cingulate_Gyrus_posterior_division,Cingulate_Gyrus_anterior_division,Paracingulate_Gyrus,Subcallosal_Cortex,Juxtapositional_Lobule_Cortex_formerly_Supplementa,Frontal_Medial_Cortex,
          Intracalcarine_Cortex,Lateral_Occipital_Cortex_inferior_division,Lateral_Occipital_Cortex_superior_division,Angular_Gyrus,Inferior_Temporal_Gyrus_temporooccipital_part,
          Inferior_Temporal_Gyrus_posterior_division,Inferior_Temporal_Gyrus_anterior_division,Middle_Temporal_Gyrus_temporooccipital_part,Middle_Temporal_Gyrus_posterior_division,
          Middle_Temporal_Gyrus_anterior_division,Superior_Temporal_Gyrus_posterior_division,Superior_Temporal_Gyrus_anterior_division,Temporal_Pole,Precentral_Gyrus,Inferior_Frontal_Gyrus_pars_opercularis,
          Inferior_Frontal_Gyrus_pars_triangularis,Middle_Frontal_Gyrus,Superior_Frontal_Gyrus,Insular_Cortex,Occipital_Pole,Supracalcarine_Cortex,Planum_Temporale,Heschls_Gyrus_includes_H1_and_H2,
          Planum_Polare,Parietal_Operculum_Cortex,Central_Opercular_Cortex,Frontal_Operculum_Cortex,Occipital_Fusiform_Gyrus,Temporal_Occipital_Fusiform_Cortex,Temporal_Fusiform_Cortex_posterior_division,
          Temporal_Fusiform_Cortex_anterior_division,Lingual_Gyrus, 
          labels = c("HO_cort_ODI"),
          ncol = 4, nrow = 12)



# HO subcort

col_list = c("Left_Cerebral_White_Matter","Left_Cerebral_Cortex","Left_Lateral_Ventricle","Left_Thalamus","Left_Caudate","Left_Putamen","Left_Pallidum","Brain.Stem","Left_Hippocampus","Left_Amygdala",
             "Left_Accumbens","Right_Cerebral_White_Matter","Right_Cerebral_Cortex","Right_Lateral_Ventricle","Right_Thalamus","Right_Caudate","Right_Putamen","Right_Pallidum","Right_Hippocampus",
             "Right_Amygdala","Right_Accumbens")
#NDI
Left_Cerebral_White_Matter = gghistogram(HO_subcort_NDI_data$Left_Cerebral_White_Matter,
                           main = "Left_Cerebral_White_Matter",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Cerebral_Cortex = gghistogram(HO_subcort_NDI_data$Left_Cerebral_Cortex,
                                                       main = "Left_Cerebral_Cortex",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Lateral_Ventricle = gghistogram(HO_subcort_NDI_data$Left_Lateral_Ventricle,
                                                      main = "Left_Lateral_Ventricle",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Thalamus = gghistogram(HO_subcort_NDI_data$Left_Thalamus,
                                     main = "Left_Thalamus",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Caudate = gghistogram(HO_subcort_NDI_data$Left_Caudate,
                            main = "Left_Caudate",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Putamen = gghistogram(HO_subcort_NDI_data$Left_Putamen,
                                main = "Left_Putamen",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Pallidum = gghistogram(HO_subcort_NDI_data$Left_Pallidum,
                                                 main = "Left_Pallidum",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Brain.Stem = gghistogram(HO_subcort_NDI_data$Brain.Stem,
                                                main = "Brain.Stem",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Hippocampus = gghistogram(HO_subcort_NDI_data$Left_Hippocampus,
                                  main = "Left_Hippocampus",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Amygdala = gghistogram(HO_subcort_NDI_data$Left_Amygdala,
                                 main = "Left_Amygdala",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Accumbens = gghistogram(HO_subcort_NDI_data$Left_Accumbens,
                                                                 main = "Left_Accumbens",
                                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Cerebral_White_Matter = gghistogram(HO_subcort_NDI_data$Right_Cerebral_White_Matter,
                                    main = "Right_Cerebral_White_Matter",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Cerebral_Cortex = gghistogram(HO_subcort_NDI_data$Right_Cerebral_Cortex,
                                    main = "Right_Cerebral_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Lateral_Ventricle = gghistogram(HO_subcort_NDI_data$Right_Lateral_Ventricle,
                                                         main = "Right_Lateral_Ventricle",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Thalamus = gghistogram(HO_subcort_NDI_data$Right_Thalamus,
                                                         main = "Right_Thalamus",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Caudate = gghistogram(HO_subcort_NDI_data$Right_Caudate,
                            main = "Right_Caudate",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Putamen = gghistogram(HO_subcort_NDI_data$Right_Putamen,
                                                     main = "Right_Putamen",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Pallidum = gghistogram(HO_subcort_NDI_data$Right_Pallidum,
                                                    main = "Right_Pallidum",
                                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Hippocampus = gghistogram(HO_subcort_NDI_data$Right_Hippocampus,
                                main = "Right_Hippocampus",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Amygdala = gghistogram(HO_subcort_NDI_data$Right_Amygdala,
                                                            main = "Right_Amygdala",
                                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Accumbens = gghistogram(HO_subcort_NDI_data$Right_Accumbens,
                                                         main = "Right_Accumbens",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Left_Cerebral_White_Matter,Left_Cerebral_Cortex,Left_Lateral_Ventricle,Left_Thalamus,Left_Caudate,Left_Putamen,Left_Pallidum,Brain.Stem,Left_Hippocampus,Left_Amygdala,
          Left_Accumbens,Right_Cerebral_White_Matter,Right_Cerebral_Cortex,Right_Lateral_Ventricle,Right_Thalamus,Right_Caudate,Right_Putamen,Right_Pallidum,Right_Hippocampus,
          Right_Amygdala,Right_Accumbens, 
          labels = c("HO_subcort_NDI"),
          ncol = 3, nrow = 7)

#FWF
Left_Cerebral_White_Matter = gghistogram(HO_subcort_FWF_data$Left_Cerebral_White_Matter,
                                         main = "Left_Cerebral_White_Matter",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Cerebral_Cortex = gghistogram(HO_subcort_FWF_data$Left_Cerebral_Cortex,
                                   main = "Left_Cerebral_Cortex",
                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Lateral_Ventricle = gghistogram(HO_subcort_FWF_data$Left_Lateral_Ventricle,
                                     main = "Left_Lateral_Ventricle",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Thalamus = gghistogram(HO_subcort_FWF_data$Left_Thalamus,
                            main = "Left_Thalamus",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Caudate = gghistogram(HO_subcort_FWF_data$Left_Caudate,
                           main = "Left_Caudate",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Putamen = gghistogram(HO_subcort_FWF_data$Left_Putamen,
                           main = "Left_Putamen",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Pallidum = gghistogram(HO_subcort_FWF_data$Left_Pallidum,
                            main = "Left_Pallidum",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Brain.Stem = gghistogram(HO_subcort_FWF_data$Brain.Stem,
                         main = "Brain.Stem",
                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Hippocampus = gghistogram(HO_subcort_FWF_data$Left_Hippocampus,
                               main = "Left_Hippocampus",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Amygdala = gghistogram(HO_subcort_FWF_data$Left_Amygdala,
                            main = "Left_Amygdala",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Accumbens = gghistogram(HO_subcort_FWF_data$Left_Accumbens,
                             main = "Left_Accumbens",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Cerebral_White_Matter = gghistogram(HO_subcort_FWF_data$Right_Cerebral_White_Matter,
                                          main = "Right_Cerebral_White_Matter",
                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Cerebral_Cortex = gghistogram(HO_subcort_FWF_data$Right_Cerebral_Cortex,
                                    main = "Right_Cerebral_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Lateral_Ventricle = gghistogram(HO_subcort_FWF_data$Right_Lateral_Ventricle,
                                      main = "Right_Lateral_Ventricle",
                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Thalamus = gghistogram(HO_subcort_FWF_data$Right_Thalamus,
                             main = "Right_Thalamus",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Caudate = gghistogram(HO_subcort_FWF_data$Right_Caudate,
                            main = "Right_Caudate",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Putamen = gghistogram(HO_subcort_FWF_data$Right_Putamen,
                            main = "Right_Putamen",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Pallidum = gghistogram(HO_subcort_FWF_data$Right_Pallidum,
                             main = "Right_Pallidum",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Hippocampus = gghistogram(HO_subcort_FWF_data$Right_Hippocampus,
                                main = "Right_Hippocampus",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Amygdala = gghistogram(HO_subcort_FWF_data$Right_Amygdala,
                             main = "Right_Amygdala",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Accumbens = gghistogram(HO_subcort_FWF_data$Right_Accumbens,
                              main = "Right_Accumbens",
                              xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Left_Cerebral_White_Matter,Left_Cerebral_Cortex,Left_Lateral_Ventricle,Left_Thalamus,Left_Caudate,Left_Putamen,Left_Pallidum,Brain.Stem,Left_Hippocampus,Left_Amygdala,
          Left_Accumbens,Right_Cerebral_White_Matter,Right_Cerebral_Cortex,Right_Lateral_Ventricle,Right_Thalamus,Right_Caudate,Right_Putamen,Right_Pallidum,Right_Hippocampus,
          Right_Amygdala,Right_Accumbens, 
          labels = c("HO_subcort_FWF"),
          ncol = 3, nrow = 7)


#ODI
Left_Cerebral_White_Matter = gghistogram(HO_subcort_ODI_data$Left_Cerebral_White_Matter,
                                         main = "Left_Cerebral_White_Matter",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Cerebral_Cortex = gghistogram(HO_subcort_ODI_data$Left_Cerebral_Cortex,
                                   main = "Left_Cerebral_Cortex",
                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Lateral_Ventricle = gghistogram(HO_subcort_ODI_data$Left_Lateral_Ventricle,
                                     main = "Left_Lateral_Ventricle",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Thalamus = gghistogram(HO_subcort_ODI_data$Left_Thalamus,
                            main = "Left_Thalamus",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Caudate = gghistogram(HO_subcort_ODI_data$Left_Caudate,
                           main = "Left_Caudate",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Putamen = gghistogram(HO_subcort_ODI_data$Left_Putamen,
                           main = "Left_Putamen",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Pallidum = gghistogram(HO_subcort_ODI_data$Left_Pallidum,
                            main = "Left_Pallidum",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Brain.Stem = gghistogram(HO_subcort_ODI_data$Brain.Stem,
                         main = "Brain.Stem",
                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Hippocampus = gghistogram(HO_subcort_ODI_data$Left_Hippocampus,
                               main = "Left_Hippocampus",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Amygdala = gghistogram(HO_subcort_ODI_data$Left_Amygdala,
                            main = "Left_Amygdala",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Left_Accumbens = gghistogram(HO_subcort_ODI_data$Left_Accumbens,
                             main = "Left_Accumbens",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Cerebral_White_Matter = gghistogram(HO_subcort_ODI_data$Right_Cerebral_White_Matter,
                                          main = "Right_Cerebral_White_Matter",
                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Cerebral_Cortex = gghistogram(HO_subcort_ODI_data$Right_Cerebral_Cortex,
                                    main = "Right_Cerebral_Cortex",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Lateral_Ventricle = gghistogram(HO_subcort_ODI_data$Right_Lateral_Ventricle,
                                      main = "Right_Lateral_Ventricle",
                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Thalamus = gghistogram(HO_subcort_ODI_data$Right_Thalamus,
                             main = "Right_Thalamus",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Caudate = gghistogram(HO_subcort_ODI_data$Right_Caudate,
                            main = "Right_Caudate",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Putamen = gghistogram(HO_subcort_ODI_data$Right_Putamen,
                            main = "Right_Putamen",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Pallidum = gghistogram(HO_subcort_ODI_data$Right_Pallidum,
                             main = "Right_Pallidum",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Hippocampus = gghistogram(HO_subcort_ODI_data$Right_Hippocampus,
                                main = "Right_Hippocampus",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Amygdala = gghistogram(HO_subcort_ODI_data$Right_Amygdala,
                             main = "Right_Amygdala",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Right_Accumbens = gghistogram(HO_subcort_ODI_data$Right_Accumbens,
                              main = "Right_Accumbens",
                              xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Left_Cerebral_White_Matter,Left_Cerebral_Cortex,Left_Lateral_Ventricle,Left_Thalamus,Left_Caudate,Left_Putamen,Left_Pallidum,Brain.Stem,Left_Hippocampus,Left_Amygdala,
          Left_Accumbens,Right_Cerebral_White_Matter,Right_Cerebral_Cortex,Right_Lateral_Ventricle,Right_Thalamus,Right_Caudate,Right_Putamen,Right_Pallidum,Right_Hippocampus,
          Right_Amygdala,Right_Accumbens, 
          labels = c("HO_subcort_ODI"),
          ncol = 3, nrow = 7)

#JHU labels

col_list = c("Anterior_corona_radiata_L","Anterior_corona_radiata_R","Superior_corona_radiata_R","Superior_corona_radiata_L","Retrolenticular_part_of_internal_capsule_L","Retrolenticular_part_of_internal_capsule_R",
             "Posterior_limb_of_internal_capsule_L","Posterior_limb_of_internal_capsule_R","Anterior_limb_of_internal_capsule_L","Anterior_limb_of_internal_capsule_R","Cerebral_peduncle_L","Cerebral_peduncle_R",
             "Superior_cerebellar_peduncle_L","Superior_cerebellar_peduncle_R","Inferior_cerebellar_peduncle_L","Inferior_cerebellar_peduncle_R","Medial_lemniscus_L","Medial_lemniscus_R","Corticospinal_tract_L",
             "Corticospinal_tract_R","Fornix_column_and_body_of_fornix","Splenium_of_corpus_callosum","Body_of_corpus_callosum","Genu_of_corpus_callosum","Pontine_crossing_tract_a_part_of_MCP","Middle_cerebellar_peduncle",
             "Tapetum_L","Tapetum_R","Uncinate_fasciculus_L","Uncinate_fasciculus_R","Inferior_fronto.occipital_fasciculus_L","Inferior_fronto.occipital_fasciculus_R","Superior_fronto.occipital_fasciculus_L",
             "Superior_fronto.occipital_fasciculus_R","Superior_longitudinal_fasciculus_L","Superior_longitudinal_fasciculus_R","Fornix_cres_L","Fornix_cres_R","Cingulum_hippocampus_L","Cingulum_hippocampus_R",
             "Cingulum_cingulate_gyrus_L","Cingulum_cingulate_gyrus_R","External_capsule_L","External_capsule_R","Sagittal_stratum_R","Sagittal_stratum_L","Posterior_thalamic_radiation_include_optic_radiation_L",
             "Posterior_thalamic_radiation_include_optic_radiation_R","Posterior_corona_radiata_L","Posterior_corona_radiata_R")

#NDI

Anterior_corona_radiata_L = gghistogram(JHU_label_NDI_data$Anterior_corona_radiata_L,
                           main = "Anterior_corona_radiata_L",
                           xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_corona_radiata_R = gghistogram(JHU_label_NDI_data$Anterior_corona_radiata_R,
                                                       main = "Anterior_corona_radiata_R",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_corona_radiata_R = gghistogram(JHU_label_NDI_data$Superior_corona_radiata_R,
                                                      main = "Superior_corona_radiata_R",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_corona_radiata_L = gghistogram(JHU_label_NDI_data$Superior_corona_radiata_L,
                                     main = "Superior_corona_radiata_L",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Retrolenticular_part_of_internal_capsule_L = gghistogram(JHU_label_NDI_data$Retrolenticular_part_of_internal_capsule_L,
                            main = "Retrolenticular_part_of_internal_capsule_L",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Retrolenticular_part_of_internal_capsule_R = gghistogram(JHU_label_NDI_data$Retrolenticular_part_of_internal_capsule_R,
                                main = "Retrolenticular_part_of_internal_capsule_R",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_limb_of_internal_capsule_L = gghistogram(JHU_label_NDI_data$Posterior_limb_of_internal_capsule_L,
                                                 main = "Posterior_limb_of_internal_capsule_L",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_limb_of_internal_capsule_R = gghistogram(JHU_label_NDI_data$Posterior_limb_of_internal_capsule_R,
                                                main = "Posterior_limb_of_internal_capsule_R",
                                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_limb_of_internal_capsule_L = gghistogram(JHU_label_NDI_data$Anterior_limb_of_internal_capsule_L,
                                  main = "Anterior_limb_of_internal_capsule_L",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_limb_of_internal_capsule_R = gghistogram(JHU_label_NDI_data$Anterior_limb_of_internal_capsule_R,
                                 main = "Anterior_limb_of_internal_capsule_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebral_peduncle_L = gghistogram(JHU_label_NDI_data$Cerebral_peduncle_L,
                                                                 main = "Cerebral_peduncle_L",
                                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebral_peduncle_R = gghistogram(JHU_label_NDI_data$Cerebral_peduncle_R,
                                    main = "Cerebral_peduncle_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_cerebellar_peduncle_L = gghistogram(JHU_label_NDI_data$Superior_cerebellar_peduncle_L,
                                    main = "Superior_cerebellar_peduncle_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_cerebellar_peduncle_R = gghistogram(JHU_label_NDI_data$Superior_cerebellar_peduncle_R,
                                                         main = "Superior_cerebellar_peduncle_R",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_cerebellar_peduncle_L = gghistogram(JHU_label_NDI_data$Inferior_cerebellar_peduncle_L,
                                                         main = "Inferior_cerebellar_peduncle_L",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_cerebellar_peduncle_R = gghistogram(JHU_label_NDI_data$Inferior_cerebellar_peduncle_R,
                            main = "Inferior_cerebellar_peduncle_R",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Medial_lemniscus_L = gghistogram(JHU_label_NDI_data$Medial_lemniscus_L,
                                                     main = "Medial_lemniscus_L",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Medial_lemniscus_R = gghistogram(JHU_label_NDI_data$Medial_lemniscus_R,
                                                    main = "Medial_lemniscus_R",
                                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_L = gghistogram(JHU_label_NDI_data$Corticospinal_tract_L,
                                main = "Corticospinal_tract_L",
                                xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_R = gghistogram(JHU_label_NDI_data$Corticospinal_tract_R,
                                                            main = "Corticospinal_tract_R",
                                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_column_and_body_of_fornix = gghistogram(JHU_label_NDI_data$Fornix_column_and_body_of_fornix,
                                                         main = "Fornix_column_and_body_of_fornix",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Splenium_of_corpus_callosum = gghistogram(JHU_label_NDI_data$Splenium_of_corpus_callosum,
                                                        main = "Splenium_of_corpus_callosum",
                                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Body_of_corpus_callosum = gghistogram(JHU_label_NDI_data$Body_of_corpus_callosum,
                                                          main = "Body_of_corpus_callosum",
                                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Genu_of_corpus_callosum = gghistogram(JHU_label_NDI_data$Genu_of_corpus_callosum,
                                                       main = "Genu_of_corpus_callosum",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Pontine_crossing_tract_a_part_of_MCP = gghistogram(JHU_label_NDI_data$Pontine_crossing_tract_a_part_of_MCP,
                                                      main = "Pontine_crossing_tract_a_part_of_MCP",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_cerebellar_peduncle = gghistogram(JHU_label_NDI_data$Middle_cerebellar_peduncle,
                                                         main = "Middle_cerebellar_peduncle",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Tapetum_L = gghistogram(JHU_label_NDI_data$Tapetum_L,
                                                        main = "Tapetum_L",
                                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Tapetum_R = gghistogram(JHU_label_NDI_data$Tapetum_R,
                            main = "Tapetum_R",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_L = gghistogram(JHU_label_NDI_data$Uncinate_fasciculus_L,
                               main = "Uncinate_fasciculus_L",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_R = gghistogram(JHU_label_NDI_data$Uncinate_fasciculus_R,
                                                      main = "Uncinate_fasciculus_R",
                                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_L = gghistogram(JHU_label_NDI_data$Inferior_fronto.occipital_fasciculus_L,
                                                       main = "Inferior_fronto.occipital_fasciculus_L",
                                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_R = gghistogram(JHU_label_NDI_data$Inferior_fronto.occipital_fasciculus_R,
                                   main = "Inferior_fronto.occipital_fasciculus_R",
                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_fronto.occipital_fasciculus_L = gghistogram(JHU_label_NDI_data$Superior_fronto.occipital_fasciculus_L,
                                     main = "Superior_fronto.occipital_fasciculus_L",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_fronto.occipital_fasciculus_R = gghistogram(JHU_label_NDI_data$Superior_fronto.occipital_fasciculus_R,
                             main = "Superior_fronto.occipital_fasciculus_R",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_L = gghistogram(JHU_label_NDI_data$Superior_longitudinal_fasciculus_L,
                             main = "Superior_longitudinal_fasciculus_L",
                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_R = gghistogram(JHU_label_NDI_data$Superior_longitudinal_fasciculus_R,
                                    main = "Superior_longitudinal_fasciculus_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_cres_L = gghistogram(JHU_label_NDI_data$Fornix_cres_L,
                               main = "Fornix_cres_L",
                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_cres_R = gghistogram(JHU_label_NDI_data$Fornix_cres_R,
                                               main = "Fornix_cres_R",
                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_L = gghistogram(JHU_label_NDI_data$Cingulum_hippocampus_L,
                            main = "Cingulum_hippocampus_L",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_R = gghistogram(JHU_label_NDI_data$Cingulum_hippocampus_R,
                                        main = "Cingulum_hippocampus_R",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_L = gghistogram(JHU_label_NDI_data$Cingulum_cingulate_gyrus_L,
                                       main = "Cingulum_cingulate_gyrus_L",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_R = gghistogram(JHU_label_NDI_data$Cingulum_cingulate_gyrus_R,
                                       main = "Cingulum_cingulate_gyrus_R",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

External_capsule_L = gghistogram(JHU_label_NDI_data$External_capsule_L,
                                       main = "External_capsule_L",
                                       xlab = "Intensity"
) + ggeasy::easy_center_title()

External_capsule_R = gghistogram(JHU_label_NDI_data$External_capsule_R,
                                                 main = "External_capsule_R",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Sagittal_stratum_R = gghistogram(JHU_label_NDI_data$Sagittal_stratum_R,
                                                          main = "Sagittal_stratum_R",
                                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Sagittal_stratum_L = gghistogram(JHU_label_NDI_data$Sagittal_stratum_L,
                                                         main = "Sagittal_stratum_L",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_thalamic_radiation_include_optic_radiation_L = gghistogram(JHU_label_NDI_data$Posterior_thalamic_radiation_include_optic_radiation_L,
                            main = "Posterior_thalamic_radiation_include_optic_radiation_L",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_thalamic_radiation_include_optic_radiation_R = gghistogram(JHU_label_NDI_data$Posterior_thalamic_radiation_include_optic_radiation_R,
                            main = "Posterior_thalamic_radiation_include_optic_radiation_R",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()


Posterior_corona_radiata_L = gghistogram(JHU_label_NDI_data$Posterior_corona_radiata_L,
                            main = "Posterior_corona_radiata_L",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_corona_radiata_R = gghistogram(JHU_label_NDI_data$Posterior_corona_radiata_R,
                                         main = "Posterior_corona_radiata_R",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()


plot_grid(Frontal_Pole,Parahippocampal_Gyrus_posterior_division,Parahippocampal_Gyrus_anterior_division,Frontal_Orbital_Cortex,Cuneal_Cortex,Precuneous_Cortex,
          Cingulate_Gyrus_posterior_division,Cingulate_Gyrus_anterior_division,Paracingulate_Gyrus,Subcallosal_Cortex,Juxtapositional_Lobule_Cortex_formerly_Supplementa,Frontal_Medial_Cortex,
          Intracalcarine_Cortex,Lateral_Occipital_Cortex_inferior_division,Lateral_Occipital_Cortex_superior_division,Angular_Gyrus,Inferior_Temporal_Gyrus_temporooccipital_part,
          Inferior_Temporal_Gyrus_posterior_division,Inferior_Temporal_Gyrus_anterior_division,Middle_Temporal_Gyrus_temporooccipital_part,Middle_Temporal_Gyrus_posterior_division,
          Middle_Temporal_Gyrus_anterior_division,Superior_Temporal_Gyrus_posterior_division,Superior_Temporal_Gyrus_anterior_division,Temporal_Pole,Precentral_Gyrus,Inferior_Frontal_Gyrus_pars_opercularis,
          Inferior_Frontal_Gyrus_pars_triangularis,Middle_Frontal_Gyrus,Superior_Frontal_Gyrus,Insular_Cortex,Occipital_Pole,Supracalcarine_Cortex,Planum_Temporale,Heschls_Gyrus_includes_H1_and_H2,
          Planum_Polare,Parietal_Operculum_Cortex,Central_Opercular_Cortex,Frontal_Operculum_Cortex,Occipital_Fusiform_Gyrus,Temporal_Occipital_Fusiform_Cortex,Temporal_Fusiform_Cortex_posterior_division,
          Temporal_Fusiform_Cortex_anterior_division,Lingual_Gyrus, 
          labels = c("JHU_label_NDI"),
          ncol = 5, nrow = 8)

#FWF

Anterior_corona_radiata_L = gghistogram(JHU_label_FWF_data$Anterior_corona_radiata_L,
                                        main = "Anterior_corona_radiata_L",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_corona_radiata_R = gghistogram(JHU_label_FWF_data$Anterior_corona_radiata_R,
                                        main = "Anterior_corona_radiata_R",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_corona_radiata_R = gghistogram(JHU_label_FWF_data$Superior_corona_radiata_R,
                                        main = "Superior_corona_radiata_R",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_corona_radiata_L = gghistogram(JHU_label_FWF_data$Superior_corona_radiata_L,
                                        main = "Superior_corona_radiata_L",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Retrolenticular_part_of_internal_capsule_L = gghistogram(JHU_label_FWF_data$Retrolenticular_part_of_internal_capsule_L,
                                                         main = "Retrolenticular_part_of_internal_capsule_L",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Retrolenticular_part_of_internal_capsule_R = gghistogram(JHU_label_FWF_data$Retrolenticular_part_of_internal_capsule_R,
                                                         main = "Retrolenticular_part_of_internal_capsule_R",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_limb_of_internal_capsule_L = gghistogram(JHU_label_FWF_data$Posterior_limb_of_internal_capsule_L,
                                                   main = "Posterior_limb_of_internal_capsule_L",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_limb_of_internal_capsule_R = gghistogram(JHU_label_FWF_data$Posterior_limb_of_internal_capsule_R,
                                                   main = "Posterior_limb_of_internal_capsule_R",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_limb_of_internal_capsule_L = gghistogram(JHU_label_FWF_data$Anterior_limb_of_internal_capsule_L,
                                                  main = "Anterior_limb_of_internal_capsule_L",
                                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_limb_of_internal_capsule_R = gghistogram(JHU_label_FWF_data$Anterior_limb_of_internal_capsule_R,
                                                  main = "Anterior_limb_of_internal_capsule_R",
                                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebral_peduncle_L = gghistogram(JHU_label_FWF_data$Cerebral_peduncle_L,
                                  main = "Cerebral_peduncle_L",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebral_peduncle_R = gghistogram(JHU_label_FWF_data$Cerebral_peduncle_R,
                                  main = "Cerebral_peduncle_R",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_cerebellar_peduncle_L = gghistogram(JHU_label_FWF_data$Superior_cerebellar_peduncle_L,
                                             main = "Superior_cerebellar_peduncle_L",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_cerebellar_peduncle_R = gghistogram(JHU_label_FWF_data$Superior_cerebellar_peduncle_R,
                                             main = "Superior_cerebellar_peduncle_R",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_cerebellar_peduncle_L = gghistogram(JHU_label_FWF_data$Inferior_cerebellar_peduncle_L,
                                             main = "Inferior_cerebellar_peduncle_L",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_cerebellar_peduncle_R = gghistogram(JHU_label_FWF_data$Inferior_cerebellar_peduncle_R,
                                             main = "Inferior_cerebellar_peduncle_R",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Medial_lemniscus_L = gghistogram(JHU_label_FWF_data$Medial_lemniscus_L,
                                 main = "Medial_lemniscus_L",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Medial_lemniscus_R = gghistogram(JHU_label_FWF_data$Medial_lemniscus_R,
                                 main = "Medial_lemniscus_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_L = gghistogram(JHU_label_FWF_data$Corticospinal_tract_L,
                                    main = "Corticospinal_tract_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_R = gghistogram(JHU_label_FWF_data$Corticospinal_tract_R,
                                    main = "Corticospinal_tract_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_column_and_body_of_fornix = gghistogram(JHU_label_FWF_data$Fornix_column_and_body_of_fornix,
                                               main = "Fornix_column_and_body_of_fornix",
                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Splenium_of_corpus_callosum = gghistogram(JHU_label_FWF_data$Splenium_of_corpus_callosum,
                                          main = "Splenium_of_corpus_callosum",
                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Body_of_corpus_callosum = gghistogram(JHU_label_FWF_data$Body_of_corpus_callosum,
                                      main = "Body_of_corpus_callosum",
                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Genu_of_corpus_callosum = gghistogram(JHU_label_FWF_data$Genu_of_corpus_callosum,
                                      main = "Genu_of_corpus_callosum",
                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Pontine_crossing_tract_a_part_of_MCP = gghistogram(JHU_label_FWF_data$Pontine_crossing_tract_a_part_of_MCP,
                                                   main = "Pontine_crossing_tract_a_part_of_MCP",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_cerebellar_peduncle = gghistogram(JHU_label_FWF_data$Middle_cerebellar_peduncle,
                                         main = "Middle_cerebellar_peduncle",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Tapetum_L = gghistogram(JHU_label_FWF_data$Tapetum_L,
                        main = "Tapetum_L",
                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Tapetum_R = gghistogram(JHU_label_FWF_data$Tapetum_R,
                        main = "Tapetum_R",
                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_L = gghistogram(JHU_label_FWF_data$Uncinate_fasciculus_L,
                                    main = "Uncinate_fasciculus_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_R = gghistogram(JHU_label_FWF_data$Uncinate_fasciculus_R,
                                    main = "Uncinate_fasciculus_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_L = gghistogram(JHU_label_FWF_data$Inferior_fronto.occipital_fasciculus_L,
                                                     main = "Inferior_fronto.occipital_fasciculus_L",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_R = gghistogram(JHU_label_FWF_data$Inferior_fronto.occipital_fasciculus_R,
                                                     main = "Inferior_fronto.occipital_fasciculus_R",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_fronto.occipital_fasciculus_L = gghistogram(JHU_label_FWF_data$Superior_fronto.occipital_fasciculus_L,
                                                     main = "Superior_fronto.occipital_fasciculus_L",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_fronto.occipital_fasciculus_R = gghistogram(JHU_label_FWF_data$Superior_fronto.occipital_fasciculus_R,
                                                     main = "Superior_fronto.occipital_fasciculus_R",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_L = gghistogram(JHU_label_FWF_data$Superior_longitudinal_fasciculus_L,
                                                 main = "Superior_longitudinal_fasciculus_L",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_R = gghistogram(JHU_label_FWF_data$Superior_longitudinal_fasciculus_R,
                                                 main = "Superior_longitudinal_fasciculus_R",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_cres_L = gghistogram(JHU_label_FWF_data$Fornix_cres_L,
                            main = "Fornix_cres_L",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_cres_R = gghistogram(JHU_label_FWF_data$Fornix_cres_R,
                            main = "Fornix_cres_R",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_L = gghistogram(JHU_label_FWF_data$Cingulum_hippocampus_L,
                                     main = "Cingulum_hippocampus_L",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_R = gghistogram(JHU_label_FWF_data$Cingulum_hippocampus_R,
                                     main = "Cingulum_hippocampus_R",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_L = gghistogram(JHU_label_FWF_data$Cingulum_cingulate_gyrus_L,
                                         main = "Cingulum_cingulate_gyrus_L",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_R = gghistogram(JHU_label_FWF_data$Cingulum_cingulate_gyrus_R,
                                         main = "Cingulum_cingulate_gyrus_R",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

External_capsule_L = gghistogram(JHU_label_FWF_data$External_capsule_L,
                                 main = "External_capsule_L",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

External_capsule_R = gghistogram(JHU_label_FWF_data$External_capsule_R,
                                 main = "External_capsule_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Sagittal_stratum_R = gghistogram(JHU_label_FWF_data$Sagittal_stratum_R,
                                 main = "Sagittal_stratum_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Sagittal_stratum_L = gghistogram(JHU_label_FWF_data$Sagittal_stratum_L,
                                 main = "Sagittal_stratum_L",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_thalamic_radiation_include_optic_radiation_L = gghistogram(JHU_label_FWF_data$Posterior_thalamic_radiation_include_optic_radiation_L,
                                                                     main = "Posterior_thalamic_radiation_include_optic_radiation_L",
                                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_thalamic_radiation_include_optic_radiation_R = gghistogram(JHU_label_FWF_data$Posterior_thalamic_radiation_include_optic_radiation_R,
                                                                     main = "Posterior_thalamic_radiation_include_optic_radiation_R",
                                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()


Posterior_corona_radiata_L = gghistogram(JHU_label_FWF_data$Posterior_corona_radiata_L,
                                         main = "Posterior_corona_radiata_L",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_corona_radiata_R = gghistogram(JHU_label_FWF_data$Posterior_corona_radiata_R,
                                         main = "Posterior_corona_radiata_R",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()


plot_grid(Frontal_Pole,Parahippocampal_Gyrus_posterior_division,Parahippocampal_Gyrus_anterior_division,Frontal_Orbital_Cortex,Cuneal_Cortex,Precuneous_Cortex,
          Cingulate_Gyrus_posterior_division,Cingulate_Gyrus_anterior_division,Paracingulate_Gyrus,Subcallosal_Cortex,Juxtapositional_Lobule_Cortex_formerly_Supplementa,Frontal_Medial_Cortex,
          Intracalcarine_Cortex,Lateral_Occipital_Cortex_inferior_division,Lateral_Occipital_Cortex_superior_division,Angular_Gyrus,Inferior_Temporal_Gyrus_temporooccipital_part,
          Inferior_Temporal_Gyrus_posterior_division,Inferior_Temporal_Gyrus_anterior_division,Middle_Temporal_Gyrus_temporooccipital_part,Middle_Temporal_Gyrus_posterior_division,
          Middle_Temporal_Gyrus_anterior_division,Superior_Temporal_Gyrus_posterior_division,Superior_Temporal_Gyrus_anterior_division,Temporal_Pole,Precentral_Gyrus,Inferior_Frontal_Gyrus_pars_opercularis,
          Inferior_Frontal_Gyrus_pars_triangularis,Middle_Frontal_Gyrus,Superior_Frontal_Gyrus,Insular_Cortex,Occipital_Pole,Supracalcarine_Cortex,Planum_Temporale,Heschls_Gyrus_includes_H1_and_H2,
          Planum_Polare,Parietal_Operculum_Cortex,Central_Opercular_Cortex,Frontal_Operculum_Cortex,Occipital_Fusiform_Gyrus,Temporal_Occipital_Fusiform_Cortex,Temporal_Fusiform_Cortex_posterior_division,
          Temporal_Fusiform_Cortex_anterior_division,Lingual_Gyrus, 
          labels = c("JHU_label_FWF"),
          ncol = 5, nrow = 8)


#ODI

Anterior_corona_radiata_L = gghistogram(JHU_label_ODI_data$Anterior_corona_radiata_L,
                                        main = "Anterior_corona_radiata_L",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_corona_radiata_R = gghistogram(JHU_label_ODI_data$Anterior_corona_radiata_R,
                                        main = "Anterior_corona_radiata_R",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_corona_radiata_R = gghistogram(JHU_label_ODI_data$Superior_corona_radiata_R,
                                        main = "Superior_corona_radiata_R",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_corona_radiata_L = gghistogram(JHU_label_ODI_data$Superior_corona_radiata_L,
                                        main = "Superior_corona_radiata_L",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Retrolenticular_part_of_internal_capsule_L = gghistogram(JHU_label_ODI_data$Retrolenticular_part_of_internal_capsule_L,
                                                         main = "Retrolenticular_part_of_internal_capsule_L",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Retrolenticular_part_of_internal_capsule_R = gghistogram(JHU_label_ODI_data$Retrolenticular_part_of_internal_capsule_R,
                                                         main = "Retrolenticular_part_of_internal_capsule_R",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_limb_of_internal_capsule_L = gghistogram(JHU_label_ODI_data$Posterior_limb_of_internal_capsule_L,
                                                   main = "Posterior_limb_of_internal_capsule_L",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_limb_of_internal_capsule_R = gghistogram(JHU_label_ODI_data$Posterior_limb_of_internal_capsule_R,
                                                   main = "Posterior_limb_of_internal_capsule_R",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_limb_of_internal_capsule_L = gghistogram(JHU_label_ODI_data$Anterior_limb_of_internal_capsule_L,
                                                  main = "Anterior_limb_of_internal_capsule_L",
                                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_limb_of_internal_capsule_R = gghistogram(JHU_label_ODI_data$Anterior_limb_of_internal_capsule_R,
                                                  main = "Anterior_limb_of_internal_capsule_R",
                                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebral_peduncle_L = gghistogram(JHU_label_ODI_data$Cerebral_peduncle_L,
                                  main = "Cerebral_peduncle_L",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Cerebral_peduncle_R = gghistogram(JHU_label_ODI_data$Cerebral_peduncle_R,
                                  main = "Cerebral_peduncle_R",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_cerebellar_peduncle_L = gghistogram(JHU_label_ODI_data$Superior_cerebellar_peduncle_L,
                                             main = "Superior_cerebellar_peduncle_L",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_cerebellar_peduncle_R = gghistogram(JHU_label_ODI_data$Superior_cerebellar_peduncle_R,
                                             main = "Superior_cerebellar_peduncle_R",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_cerebellar_peduncle_L = gghistogram(JHU_label_ODI_data$Inferior_cerebellar_peduncle_L,
                                             main = "Inferior_cerebellar_peduncle_L",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_cerebellar_peduncle_R = gghistogram(JHU_label_ODI_data$Inferior_cerebellar_peduncle_R,
                                             main = "Inferior_cerebellar_peduncle_R",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Medial_lemniscus_L = gghistogram(JHU_label_ODI_data$Medial_lemniscus_L,
                                 main = "Medial_lemniscus_L",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Medial_lemniscus_R = gghistogram(JHU_label_ODI_data$Medial_lemniscus_R,
                                 main = "Medial_lemniscus_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_L = gghistogram(JHU_label_ODI_data$Corticospinal_tract_L,
                                    main = "Corticospinal_tract_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_R = gghistogram(JHU_label_ODI_data$Corticospinal_tract_R,
                                    main = "Corticospinal_tract_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_column_and_body_of_fornix = gghistogram(JHU_label_ODI_data$Fornix_column_and_body_of_fornix,
                                               main = "Fornix_column_and_body_of_fornix",
                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Splenium_of_corpus_callosum = gghistogram(JHU_label_ODI_data$Splenium_of_corpus_callosum,
                                          main = "Splenium_of_corpus_callosum",
                                          xlab = "Intensity"
) + ggeasy::easy_center_title()

Body_of_corpus_callosum = gghistogram(JHU_label_ODI_data$Body_of_corpus_callosum,
                                      main = "Body_of_corpus_callosum",
                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Genu_of_corpus_callosum = gghistogram(JHU_label_ODI_data$Genu_of_corpus_callosum,
                                      main = "Genu_of_corpus_callosum",
                                      xlab = "Intensity"
) + ggeasy::easy_center_title()

Pontine_crossing_tract_a_part_of_MCP = gghistogram(JHU_label_ODI_data$Pontine_crossing_tract_a_part_of_MCP,
                                                   main = "Pontine_crossing_tract_a_part_of_MCP",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Middle_cerebellar_peduncle = gghistogram(JHU_label_ODI_data$Middle_cerebellar_peduncle,
                                         main = "Middle_cerebellar_peduncle",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Tapetum_L = gghistogram(JHU_label_ODI_data$Tapetum_L,
                        main = "Tapetum_L",
                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Tapetum_R = gghistogram(JHU_label_ODI_data$Tapetum_R,
                        main = "Tapetum_R",
                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_L = gghistogram(JHU_label_ODI_data$Uncinate_fasciculus_L,
                                    main = "Uncinate_fasciculus_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_R = gghistogram(JHU_label_ODI_data$Uncinate_fasciculus_R,
                                    main = "Uncinate_fasciculus_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_L = gghistogram(JHU_label_ODI_data$Inferior_fronto.occipital_fasciculus_L,
                                                     main = "Inferior_fronto.occipital_fasciculus_L",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_R = gghistogram(JHU_label_ODI_data$Inferior_fronto.occipital_fasciculus_R,
                                                     main = "Inferior_fronto.occipital_fasciculus_R",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_fronto.occipital_fasciculus_L = gghistogram(JHU_label_ODI_data$Superior_fronto.occipital_fasciculus_L,
                                                     main = "Superior_fronto.occipital_fasciculus_L",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_fronto.occipital_fasciculus_R = gghistogram(JHU_label_ODI_data$Superior_fronto.occipital_fasciculus_R,
                                                     main = "Superior_fronto.occipital_fasciculus_R",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_L = gghistogram(JHU_label_ODI_data$Superior_longitudinal_fasciculus_L,
                                                 main = "Superior_longitudinal_fasciculus_L",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_R = gghistogram(JHU_label_ODI_data$Superior_longitudinal_fasciculus_R,
                                                 main = "Superior_longitudinal_fasciculus_R",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_cres_L = gghistogram(JHU_label_ODI_data$Fornix_cres_L,
                            main = "Fornix_cres_L",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Fornix_cres_R = gghistogram(JHU_label_ODI_data$Fornix_cres_R,
                            main = "Fornix_cres_R",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_L = gghistogram(JHU_label_ODI_data$Cingulum_hippocampus_L,
                                     main = "Cingulum_hippocampus_L",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_R = gghistogram(JHU_label_ODI_data$Cingulum_hippocampus_R,
                                     main = "Cingulum_hippocampus_R",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_L = gghistogram(JHU_label_ODI_data$Cingulum_cingulate_gyrus_L,
                                         main = "Cingulum_cingulate_gyrus_L",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_R = gghistogram(JHU_label_ODI_data$Cingulum_cingulate_gyrus_R,
                                         main = "Cingulum_cingulate_gyrus_R",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

External_capsule_L = gghistogram(JHU_label_ODI_data$External_capsule_L,
                                 main = "External_capsule_L",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

External_capsule_R = gghistogram(JHU_label_ODI_data$External_capsule_R,
                                 main = "External_capsule_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Sagittal_stratum_R = gghistogram(JHU_label_ODI_data$Sagittal_stratum_R,
                                 main = "Sagittal_stratum_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Sagittal_stratum_L = gghistogram(JHU_label_ODI_data$Sagittal_stratum_L,
                                 main = "Sagittal_stratum_L",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_thalamic_radiation_include_optic_radiation_L = gghistogram(JHU_label_ODI_data$Posterior_thalamic_radiation_include_optic_radiation_L,
                                                                     main = "Posterior_thalamic_radiation_include_optic_radiation_L",
                                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_thalamic_radiation_include_optic_radiation_R = gghistogram(JHU_label_ODI_data$Posterior_thalamic_radiation_include_optic_radiation_R,
                                                                     main = "Posterior_thalamic_radiation_include_optic_radiation_R",
                                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()


Posterior_corona_radiata_L = gghistogram(JHU_label_ODI_data$Posterior_corona_radiata_L,
                                         main = "Posterior_corona_radiata_L",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Posterior_corona_radiata_R = gghistogram(JHU_label_ODI_data$Posterior_corona_radiata_R,
                                         main = "Posterior_corona_radiata_R",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()


plot_grid(Frontal_Pole,Parahippocampal_Gyrus_posterior_division,Parahippocampal_Gyrus_anterior_division,Frontal_Orbital_Cortex,Cuneal_Cortex,Precuneous_Cortex,
          Cingulate_Gyrus_posterior_division,Cingulate_Gyrus_anterior_division,Paracingulate_Gyrus,Subcallosal_Cortex,Juxtapositional_Lobule_Cortex_formerly_Supplementa,Frontal_Medial_Cortex,
          Intracalcarine_Cortex,Lateral_Occipital_Cortex_inferior_division,Lateral_Occipital_Cortex_superior_division,Angular_Gyrus,Inferior_Temporal_Gyrus_temporooccipital_part,
          Inferior_Temporal_Gyrus_posterior_division,Inferior_Temporal_Gyrus_anterior_division,Middle_Temporal_Gyrus_temporooccipital_part,Middle_Temporal_Gyrus_posterior_division,
          Middle_Temporal_Gyrus_anterior_division,Superior_Temporal_Gyrus_posterior_division,Superior_Temporal_Gyrus_anterior_division,Temporal_Pole,Precentral_Gyrus,Inferior_Frontal_Gyrus_pars_opercularis,
          Inferior_Frontal_Gyrus_pars_triangularis,Middle_Frontal_Gyrus,Superior_Frontal_Gyrus,Insular_Cortex,Occipital_Pole,Supracalcarine_Cortex,Planum_Temporale,Heschls_Gyrus_includes_H1_and_H2,
          Planum_Polare,Parietal_Operculum_Cortex,Central_Opercular_Cortex,Frontal_Operculum_Cortex,Occipital_Fusiform_Gyrus,Temporal_Occipital_Fusiform_Cortex,Temporal_Fusiform_Cortex_posterior_division,
          Temporal_Fusiform_Cortex_anterior_division,Lingual_Gyrus, 
          labels = c("JHU_label_ODI"),
          ncol = 5, nrow = 8)


#JHU tract

col_list = c("Superior_longitudinal_fasciculus_temporal_part_R","Superior_longitudinal_fasciculus_temporal_part_L","Uncinate_fasciculus_R","Uncinate_fasciculus_L","Superior_longitudinal_fasciculus_R",
             "Superior_longitudinal_fasciculus_L","Inferior_longitudinal_fasciculus_R","Inferior_longitudinal_fasciculus_L","Inferior_fronto.occipital_fasciculus_R","Inferior_fronto.occipital_fasciculus_L",
             "Forceps_minor","Forceps_major","Cingulum_hippocampus_R","Cingulum_hippocampus_L","Cingulum_cingulate_gyrus_R","Cingulum_cingulate_gyrus_L","Corticospinal_tract_R","Corticospinal_tract_L",
             "Anterior_thalamic_radiation_R","Anterior_thalamic_radiation_L")

#NDI

Superior_longitudinal_fasciculus_temporal_part_R = gghistogram(JHU_tract_NDI_data$Superior_longitudinal_fasciculus_temporal_part_R,
                                        main = "Superior_longitudinal_fasciculus_temporal_part_R",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_temporal_part_L = gghistogram(JHU_tract_NDI_data$Superior_longitudinal_fasciculus_temporal_part_L,
                                        main = "Superior_longitudinal_fasciculus_temporal_part_L",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_R = gghistogram(JHU_tract_NDI_data$Uncinate_fasciculus_R,
                                        main = "Uncinate_fasciculus_R",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_L = gghistogram(JHU_tract_NDI_data$Uncinate_fasciculus_L,
                                        main = "Uncinate_fasciculus_L",
                                        xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_R = gghistogram(JHU_tract_NDI_data$Superior_longitudinal_fasciculus_R,
                                                         main = "Superior_longitudinal_fasciculus_R",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_L = gghistogram(JHU_tract_NDI_data$Superior_longitudinal_fasciculus_L,
                                                         main = "Superior_longitudinal_fasciculus_L",
                                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_longitudinal_fasciculus_R = gghistogram(JHU_tract_NDI_data$Inferior_longitudinal_fasciculus_R,
                                                   main = "Inferior_longitudinal_fasciculus_R",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_longitudinal_fasciculus_L = gghistogram(JHU_tract_NDI_data$Inferior_longitudinal_fasciculus_L,
                                                   main = "Inferior_longitudinal_fasciculus_L",
                                                   xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_R = gghistogram(JHU_tract_NDI_data$Inferior_fronto.occipital_fasciculus_R,
                                                  main = "Inferior_fronto.occipital_fasciculus_R",
                                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_L = gghistogram(JHU_tract_NDI_data$Inferior_fronto.occipital_fasciculus_L,
                                                  main = "Inferior_fronto.occipital_fasciculus_L",
                                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Forceps_minor = gghistogram(JHU_tract_NDI_data$Forceps_minor,
                                  main = "Forceps_minor",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Forceps_major = gghistogram(JHU_tract_NDI_data$Forceps_major,
                                  main = "Forceps_major",
                                  xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_R = gghistogram(JHU_tract_NDI_data$Cingulum_hippocampus_R,
                                             main = "Cingulum_hippocampus_R",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_L = gghistogram(JHU_tract_NDI_data$Cingulum_hippocampus_L,
                                             main = "Cingulum_hippocampus_L",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_R = gghistogram(JHU_tract_NDI_data$Cingulum_cingulate_gyrus_R,
                                             main = "Cingulum_cingulate_gyrus_R",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_L = gghistogram(JHU_tract_NDI_data$Cingulum_cingulate_gyrus_L,
                                             main = "Cingulum_cingulate_gyrus_L",
                                             xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_R = gghistogram(JHU_tract_NDI_data$Corticospinal_tract_R,
                                 main = "Corticospinal_tract_R",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_L = gghistogram(JHU_tract_NDI_data$Corticospinal_tract_L,
                                 main = "Corticospinal_tract_L",
                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_thalamic_radiation_R = gghistogram(JHU_tract_NDI_data$Anterior_thalamic_radiation_R,
                                    main = "Anterior_thalamic_radiation_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_thalamic_radiation_L = gghistogram(JHU_tract_NDI_data$Anterior_thalamic_radiation_L,
                                    main = "Anterior_thalamic_radiation_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Superior_longitudinal_fasciculus_temporal_part_R,Superior_longitudinal_fasciculus_temporal_part_L,Uncinate_fasciculus_R,Uncinate_fasciculus_L,Superior_longitudinal_fasciculus_R,
          Superior_longitudinal_fasciculus_L,Inferior_longitudinal_fasciculus_R,Inferior_longitudinal_fasciculus_L,Inferior_fronto.occipital_fasciculus_R,Inferior_fronto.occipital_fasciculus_L,
          Forceps_minor,Forceps_major,Cingulum_hippocampus_R,Cingulum_hippocampus_L,Cingulum_cingulate_gyrus_R,Cingulum_cingulate_gyrus_L,Corticospinal_tract_R,Corticospinal_tract_L,
          Anterior_thalamic_radiation_R,Anterior_thalamic_radiation_L, 
          labels = c("JHU_tract_NDI"),
          ncol = 5, nrow = 4)



# FWF

Superior_longitudinal_fasciculus_temporal_part_R = gghistogram(JHU_tract_FWF_data$Superior_longitudinal_fasciculus_temporal_part_R,
                                                               main = "Superior_longitudinal_fasciculus_temporal_part_R",
                                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_temporal_part_L = gghistogram(JHU_tract_FWF_data$Superior_longitudinal_fasciculus_temporal_part_L,
                                                               main = "Superior_longitudinal_fasciculus_temporal_part_L",
                                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_R = gghistogram(JHU_tract_FWF_data$Uncinate_fasciculus_R,
                                    main = "Uncinate_fasciculus_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_L = gghistogram(JHU_tract_FWF_data$Uncinate_fasciculus_L,
                                    main = "Uncinate_fasciculus_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_R = gghistogram(JHU_tract_FWF_data$Superior_longitudinal_fasciculus_R,
                                                 main = "Superior_longitudinal_fasciculus_R",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_L = gghistogram(JHU_tract_FWF_data$Superior_longitudinal_fasciculus_L,
                                                 main = "Superior_longitudinal_fasciculus_L",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_longitudinal_fasciculus_R = gghistogram(JHU_tract_FWF_data$Inferior_longitudinal_fasciculus_R,
                                                 main = "Inferior_longitudinal_fasciculus_R",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_longitudinal_fasciculus_L = gghistogram(JHU_tract_FWF_data$Inferior_longitudinal_fasciculus_L,
                                                 main = "Inferior_longitudinal_fasciculus_L",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_R = gghistogram(JHU_tract_FWF_data$Inferior_fronto.occipital_fasciculus_R,
                                                     main = "Inferior_fronto.occipital_fasciculus_R",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_L = gghistogram(JHU_tract_FWF_data$Inferior_fronto.occipital_fasciculus_L,
                                                     main = "Inferior_fronto.occipital_fasciculus_L",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Forceps_minor = gghistogram(JHU_tract_FWF_data$Forceps_minor,
                            main = "Forceps_minor",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Forceps_major = gghistogram(JHU_tract_FWF_data$Forceps_major,
                            main = "Forceps_major",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_R = gghistogram(JHU_tract_FWF_data$Cingulum_hippocampus_R,
                                     main = "Cingulum_hippocampus_R",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_L = gghistogram(JHU_tract_FWF_data$Cingulum_hippocampus_L,
                                     main = "Cingulum_hippocampus_L",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_R = gghistogram(JHU_tract_FWF_data$Cingulum_cingulate_gyrus_R,
                                         main = "Cingulum_cingulate_gyrus_R",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_L = gghistogram(JHU_tract_FWF_data$Cingulum_cingulate_gyrus_L,
                                         main = "Cingulum_cingulate_gyrus_L",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_R = gghistogram(JHU_tract_FWF_data$Corticospinal_tract_R,
                                    main = "Corticospinal_tract_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_L = gghistogram(JHU_tract_FWF_data$Corticospinal_tract_L,
                                    main = "Corticospinal_tract_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_thalamic_radiation_R = gghistogram(JHU_tract_FWF_data$Anterior_thalamic_radiation_R,
                                            main = "Anterior_thalamic_radiation_R",
                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_thalamic_radiation_L = gghistogram(JHU_tract_FWF_data$Anterior_thalamic_radiation_L,
                                            main = "Anterior_thalamic_radiation_L",
                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Superior_longitudinal_fasciculus_temporal_part_R,Superior_longitudinal_fasciculus_temporal_part_L,Uncinate_fasciculus_R,Uncinate_fasciculus_L,Superior_longitudinal_fasciculus_R,
          Superior_longitudinal_fasciculus_L,Inferior_longitudinal_fasciculus_R,Inferior_longitudinal_fasciculus_L,Inferior_fronto.occipital_fasciculus_R,Inferior_fronto.occipital_fasciculus_L,
          Forceps_minor,Forceps_major,Cingulum_hippocampus_R,Cingulum_hippocampus_L,Cingulum_cingulate_gyrus_R,Cingulum_cingulate_gyrus_L,Corticospinal_tract_R,Corticospinal_tract_L,
          Anterior_thalamic_radiation_R,Anterior_thalamic_radiation_L, 
          labels = c("JHU_tract_FWF"),
          ncol = 5, nrow = 4)

# ODI

Superior_longitudinal_fasciculus_temporal_part_R = gghistogram(JHU_tract_ODI_data$Superior_longitudinal_fasciculus_temporal_part_R,
                                                               main = "Superior_longitudinal_fasciculus_temporal_part_R",
                                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_temporal_part_L = gghistogram(JHU_tract_ODI_data$Superior_longitudinal_fasciculus_temporal_part_L,
                                                               main = "Superior_longitudinal_fasciculus_temporal_part_L",
                                                               xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_R = gghistogram(JHU_tract_ODI_data$Uncinate_fasciculus_R,
                                    main = "Uncinate_fasciculus_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Uncinate_fasciculus_L = gghistogram(JHU_tract_ODI_data$Uncinate_fasciculus_L,
                                    main = "Uncinate_fasciculus_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_R = gghistogram(JHU_tract_ODI_data$Superior_longitudinal_fasciculus_R,
                                                 main = "Superior_longitudinal_fasciculus_R",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Superior_longitudinal_fasciculus_L = gghistogram(JHU_tract_ODI_data$Superior_longitudinal_fasciculus_L,
                                                 main = "Superior_longitudinal_fasciculus_L",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_longitudinal_fasciculus_R = gghistogram(JHU_tract_ODI_data$Inferior_longitudinal_fasciculus_R,
                                                 main = "Inferior_longitudinal_fasciculus_R",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_longitudinal_fasciculus_L = gghistogram(JHU_tract_ODI_data$Inferior_longitudinal_fasciculus_L,
                                                 main = "Inferior_longitudinal_fasciculus_L",
                                                 xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_R = gghistogram(JHU_tract_ODI_data$Inferior_fronto.occipital_fasciculus_R,
                                                     main = "Inferior_fronto.occipital_fasciculus_R",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Inferior_fronto.occipital_fasciculus_L = gghistogram(JHU_tract_ODI_data$Inferior_fronto.occipital_fasciculus_L,
                                                     main = "Inferior_fronto.occipital_fasciculus_L",
                                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Forceps_minor = gghistogram(JHU_tract_ODI_data$Forceps_minor,
                            main = "Forceps_minor",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Forceps_major = gghistogram(JHU_tract_ODI_data$Forceps_major,
                            main = "Forceps_major",
                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_R = gghistogram(JHU_tract_ODI_data$Cingulum_hippocampus_R,
                                     main = "Cingulum_hippocampus_R",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_hippocampus_L = gghistogram(JHU_tract_ODI_data$Cingulum_hippocampus_L,
                                     main = "Cingulum_hippocampus_L",
                                     xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_R = gghistogram(JHU_tract_ODI_data$Cingulum_cingulate_gyrus_R,
                                         main = "Cingulum_cingulate_gyrus_R",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Cingulum_cingulate_gyrus_L = gghistogram(JHU_tract_ODI_data$Cingulum_cingulate_gyrus_L,
                                         main = "Cingulum_cingulate_gyrus_L",
                                         xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_R = gghistogram(JHU_tract_ODI_data$Corticospinal_tract_R,
                                    main = "Corticospinal_tract_R",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Corticospinal_tract_L = gghistogram(JHU_tract_ODI_data$Corticospinal_tract_L,
                                    main = "Corticospinal_tract_L",
                                    xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_thalamic_radiation_R = gghistogram(JHU_tract_ODI_data$Anterior_thalamic_radiation_R,
                                            main = "Anterior_thalamic_radiation_R",
                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

Anterior_thalamic_radiation_L = gghistogram(JHU_tract_ODI_data$Anterior_thalamic_radiation_L,
                                            main = "Anterior_thalamic_radiation_L",
                                            xlab = "Intensity"
) + ggeasy::easy_center_title()

plot_grid(Superior_longitudinal_fasciculus_temporal_part_R,Superior_longitudinal_fasciculus_temporal_part_L,Uncinate_fasciculus_R,Uncinate_fasciculus_L,Superior_longitudinal_fasciculus_R,
          Superior_longitudinal_fasciculus_L,Inferior_longitudinal_fasciculus_R,Inferior_longitudinal_fasciculus_L,Inferior_fronto.occipital_fasciculus_R,Inferior_fronto.occipital_fasciculus_L,
          Forceps_minor,Forceps_major,Cingulum_hippocampus_R,Cingulum_hippocampus_L,Cingulum_cingulate_gyrus_R,Cingulum_cingulate_gyrus_L,Corticospinal_tract_R,Corticospinal_tract_L,
          Anterior_thalamic_radiation_R,Anterior_thalamic_radiation_L, 
          labels = c("JHU_tract_ODI"),
          ncol = 5, nrow = 4)




