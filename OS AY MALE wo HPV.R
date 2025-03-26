library("survival")
library("survminer")
library(readr) #import txt
library(readxl) #import excel
library(tidyverse) # data manipulation
library(limma) # DEG
library(edgeR) # DEG
library(AnnotationDbi) # Converting Entrez Gene ID to Symbol
library(org.Hs.eg.db) # Converting Entrez Gene ID to Symbol
library(ComplexHeatmap) #Making the Heatmap
library(circlize) #colours
library(RColorBrewer) # colours

TCGA_clinicaldata <- read_excel("~/Desktop/PROMOTION/DEG/MALE/UNI-:MultiVARIAT/TCGA-HNSC Master clinical data MALE.xlsx")
PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative")

PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)


#---------Age------------

PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative") 

PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)


PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$Age))

res.cox <- coxph(Surv(OS.time, OS) ~ Age, PTGS_lasso)
res.cox

summary(res.cox)


#----------tobacco-------

PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative") 

PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$Smoking))
PTGS_lasso <- PTGS_lasso %>%
  filter(Smoking == "Yes"| Smoking == "No")
res.cox <- coxph(Surv(OS.time, OS) ~ Smoking, PTGS_lasso)
res.cox

summary(res.cox)

#----------ALCOHOL----------------

PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative")

PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$Alcohol))
PTGS_lasso <- PTGS_lasso %>%
  filter(Alcohol == "Yes"| Alcohol == "No")
res.cox <- coxph(Surv(OS.time, OS) ~ Alcohol, PTGS_lasso)
res.cox

summary(res.cox)


#----------ALCOHOL----------------

PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative")

PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$Alcohol))
PTGS_lasso <- PTGS_lasso %>%
  filter(Alcohol == "Yes"| Alcohol == "No")
PTGS_lasso <- PTGS_lasso %>%
  filter(Radiation == "Yes"| Radiation == "No")
res.cox <- coxph(Surv(OS.time, OS) ~ Radiation, PTGS_lasso)
res.cox

summary(res.cox)

#----------HPV16----------------

PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative") 

PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$HPV16))
PTGS_lasso <- PTGS_lasso %>%
  filter(HPV16 == "Negative"| HPV16 == "Positive")
res.cox <- coxph(Surv(OS.time, OS) ~ HPV16, PTGS_lasso)
res.cox

summary(res.cox)

#----------Tumorsize cT----------------

PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative")

PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$cT))

res.cox <- coxph(Surv(OS.time, OS) ~ cT, PTGS_lasso)
res.cox

summary(res.cox)

#----------Tumorsize pT----------------

PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative") 

PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$pT))

res.cox <- coxph(Surv(OS.time, OS) ~ pT, PTGS_lasso)
res.cox

summary(res.cox)

#----------Lymphnodes cN----------------

PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative")

PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$cN))

res.cox <- coxph(Surv(OS.time, OS) ~ cN, PTGS_lasso)
res.cox

summary(res.cox)

#----------Lymphnodes pN----------------

PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative") 

PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$pN))

res.cox <- coxph(Surv(OS.time, OS) ~ pN, PTGS_lasso)
res.cox

summary(res.cox)


#----------Grading----------------

PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative")

PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$Grading))

res.cox <- coxph(Surv(OS.time, OS) ~ Grading, PTGS_lasso)
res.cox

summary(res.cox)



#----------margin_status----------------

PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative") 

PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$margin_status))

res.cox <- coxph(Surv(OS.time, OS) ~ margin_status, PTGS_lasso)
res.cox

summary(res.cox)

#----------ALI-------

PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative")

PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$ALI))
PTGS_lasso <- PTGS_lasso %>%
  filter(ALI == "Yes"| ALI == "No")
res.cox <- coxph(Surv(OS.time, OS) ~ ALI, PTGS_lasso)
res.cox

summary(res.cox)


#----------PNI-------

PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative") 

PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$PNI))
PTGS_lasso <- PTGS_lasso %>%
  filter(PNI == "Yes"| PNI == "No")
res.cox <- coxph(Surv(OS.time, OS) ~ PNI, PTGS_lasso)
res.cox

summary(res.cox)

#------high and low risk------------

RiskScore <- read_excel("~/Desktop/PROMOTION/DEG/MALE/DEG ohne HPV/LassoCOx/OS/AY/RISK MALE OS AY WO HPV.xlsx")

res.cox <- coxph(Surv(OS.time, OS) ~ risk, RiskScore)
res.cox

summary(res.cox)

#------OS mit Endpunkt DSS------------
PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative") 

RiskScore <- read_excel("Desktop/PROMOTION/DEG/MALE/DEG ohne HPV/LassoCOx/OS/AY/RISK MALE OS AY WO HPV.xlsx")


PTGS_lasso <- left_join(RiskScore, PTGS_lasso, by = "Patient_ID") #adding to your clinicaldata
PTGS_lasso$DSS.time <- PTGS_lasso$DSS_months
PTGS_lasso$DSS <- PTGS_lasso$DSS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$DSS_months <- NULL

res.cox <- coxph(Surv(DSS.time, DSS) ~ risk, PTGS_lasso)
res.cox

summary(res.cox)

#------OS mit Endpunkt PFI------------
PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative")

RiskScore <- read_excel("Desktop/PROMOTION/DEG/MALE/DEG ohne HPV/LassoCOx/OS/AY/RISK MALE OS AY WO HPV.xlsx")

PTGS_lasso <- left_join(RiskScore, PTGS_lasso, by = "Patient_ID") #adding to your clinicaldata
PTGS_lasso$PFI.time <- PTGS_lasso$PFI_months
PTGS_lasso$PFI <- PTGS_lasso$PFI
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$PFI_months <- NULL

res.cox <- coxph(Surv(PFI.time, PFI) ~ risk, PTGS_lasso)
res.cox

summary(res.cox)





#------------MULTIVARIAT------

PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative")


PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)

RiskScore <- RiskScore %>% 
  dplyr::select(Patient_ID, risk)

PTGS_lasso <- left_join(RiskScore, PTGS_lasso, by = "Patient_ID") #adding to your clinicaldata

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$HPV16))


res.cox <- coxph(Surv(OS.time, OS) ~   Radiation + pT  + cN + pN + margin_status  + PNI + risk, PTGS_lasso)
summary(res.cox)


#------------MULTIVARIAT------


PTGS_lasso <- TCGA_clinicaldata %>%
  filter(HPV16== "Negative")


PTGS_lasso$OS.time <- PTGS_lasso$OS_months
PTGS_lasso$OS <- PTGS_lasso$OS
#PTGS_lasso$OS_5y_event <- NULL
PTGS_lasso$OS_months <- NULL

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$OS))

PTGS_lasso <- PTGS_lasso %>%
  filter(OS == "0"| OS == "1")

PTGS_lasso $OS <- as.numeric(PTGS_lasso $OS)

RiskScore <- RiskScore %>% 
  dplyr::select(Patient_ID, risk)

PTGS_lasso <- left_join(RiskScore, PTGS_lasso, by = "Patient_ID") #adding to your clinicaldata

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$HPV16))

PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$Smoking))
PTGS_lasso <- PTGS_lasso %>%
  filter(Smoking == "Yes"| Smoking == "No")
PTGS_lasso <- PTGS_lasso %>%
  filter(!is.na(PTGS_lasso$Alcohol))
PTGS_lasso <- PTGS_lasso %>%
  filter(Alcohol == "Yes"| Alcohol == "No")

res.cox <- coxph(Surv(OS.time, OS) ~Age +Smoking + Alcohol + Radiation + cT + pT+ cN+ pN + Grading + margin_status + risk, PTGS_lasso)
summary(res.cox)