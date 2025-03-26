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
library(survival)
library(survminer)
library(dplyr)



# Füge die Clinicaldata mit dem Cluster zusammen
TCGA_clinicaldata <- read_excel("~/Desktop/PROMOTION/DEG/MALE/UNI-:MultiVARIAT/TCGA-HNSC Master clinical data MALE.xlsx")
RiskScore <- read_excel("~/Desktop/PROMOTION/DEG/MALE/Lasso cox/OS/all/new/RISK MALE OS UNBE.xlsx")

RiskScore <- RiskScore %>%
  dplyr::select(risk, Patient_ID)


TCGA_clinicaldata <- left_join(RiskScore, TCGA_clinicaldata, by = "Patient_ID") #adding to your clinicaldata

#Annotations========

sn <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$risk)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select( HPV16, Alcohol, Smoking, Gender, DSS_5y_event, PNI, ALI, Age, Subsite, pathologic_stage,cT, cN, margin_status, risk) %>% # selecting your annotations
  as.data.frame(.) #data frame



#HPV

sn <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$risk)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select( HPV16, Alcohol, Smoking, Gender, DSS_5y_event, PNI, ALI, Age, Subsite, pathologic_stage,cT, cN, margin_status, risk) %>% # selecting your annotations
  as.data.frame(.) #data frameame

sn <- sn %>%
  filter(HPV16 == "Negative"| HPV16 == "Positive")
sn.tab
sn.tab <- table(sn$HPV16,sn$risk)
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
sn.tab
#names(sn)
chisq.test(sn.tab)
#da viele werte unter 5 -> ChiQ kann inkorrekt sein, deswegen jetzt Fisher Testfür exakten p wert:
fisher.test(sn$HPV16, sn$risk)

#PNI

sn <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$risk)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select( HPV16, Alcohol, Smoking, Gender, DSS_5y_event, PNI, ALI, Age, Subsite, pathologic_stage,cT, cN, margin_status, risk) %>% # selecting your annotations
  as.data.frame(.) #data frame
sn <- sn %>%
  filter(!is.na(sn$PNI))

sn.tab <- table(sn$PNI,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

#Gender

sn <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$risk)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select( HPV16, Alcohol, Smoking, Gender, DSS_5y_event, PNI, ALI, Age, Subsite, pathologic_stage,cT, cN, margin_status, risk) %>% # selecting your annotations
  as.data.frame(.) #data frame
sn <- sn %>%
  filter(!is.na(sn$Gender))

sn.tab <- table(sn$Gender,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

#-----------SMOKING-----------

sn <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$risk)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select( HPV16, Alcohol, Smoking, Gender, DSS_5y_event, PNI, ALI, Age, Subsite, pathologic_stage,cT, cN, margin_status, risk) %>% # selecting your annotations
  as.data.frame(.) #data frame
sn <- sn %>%
  filter(!is.na(sn$Smoking))
sn <- sn %>%
  filter(Smoking == "Yes"| Smoking == "No")

sn.tab <- table(sn$Smoking,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

#--------ALCOHOL--------

sn <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$risk)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select( HPV16, Alcohol, Smoking, Gender, DSS_5y_event, PNI, ALI, Age, Subsite, pathologic_stage,cT, cN, margin_status, risk) %>% # selecting your annotations
  as.data.frame(.) #data frame
sn <- sn %>%
  filter(!is.na(sn$Alcohol))
sn <- sn %>%
  filter(Alcohol == "Yes"| Alcohol == "No")

sn.tab <- table(sn$Alcohol,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)


#--------os_5y_event---------

sn <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$risk)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select( HPV16, Alcohol, Smoking, Gender, DSS_5y_event, PNI, ALI, Age, Subsite, pathologic_stage,cT, cN, margin_status, risk) %>% # selecting your annotations
  as.data.frame(.) #data frame
sn <- sn %>%
  filter(!is.na(sn$DSS_5y_event))
sn <- sn %>%
  filter(DSS_5y_event == "0"| DSS_5y_event == "1")

sn.tab <- table(sn$DSS_5y_event,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)


#--------ALI---------

sn <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$risk)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select( HPV16, Alcohol, Smoking, Gender, DSS_5y_event, PNI, ALI, Age, Subsite, pathologic_stage,cT, cN, margin_status, risk) %>% # selecting your annotations
  as.data.frame(.) #data frame
sn <- sn %>%
  filter(!is.na(sn$ALI))
sn <- sn %>%
  filter(ALI == "Yes"| ALI == "No")

sn.tab <- table(sn$ALI,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

#---------- Age------------

sn <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$risk)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select( HPV16, Alcohol, Smoking, Gender, DSS_5y_event, PNI, ALI, Age, Subsite, pathologic_stage,cT, cN, margin_status, risk) %>% # selecting your annotations
  as.data.frame(.) #data frame
sn <- sn %>%
  filter(!is.na(sn$Age))


sn.tab <- table(sn$Age,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)


#--------SUBSITE------

sn <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$risk)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select( HPV16, Alcohol, Smoking, Gender, DSS_5y_event, PNI, ALI, Age, Subsite, pathologic_stage,cT, cN, margin_status, risk) %>% # selecting your annotations
  as.data.frame(.) #data frame
sn <- sn %>%
  filter(!is.na(sn$Subsite))


sn.tab <- table(sn$Subsite,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)


#-------cT----------

sn <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$risk)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select( HPV16, Alcohol, Smoking, Gender, DSS_5y_event, PNI, ALI, Age, Subsite, pathologic_stage,cT, cN, margin_status, risk) %>% # selecting your annotations
  as.data.frame(.) #data frame
sn <- sn %>%
  filter(!is.na(sn$cT))


sn.tab <- table(sn$cT,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

#------cN--------

sn <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$risk)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select( HPV16, Alcohol, Smoking, Gender, DSS_5y_event, PNI, ALI, Age, Subsite, pathologic_stage,cT, cN, margin_status, risk) %>% # selecting your annotations
  as.data.frame(.) #data frame
sn <- sn %>%
  filter(!is.na(sn$cN))


sn.tab <- table(sn$cN,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

#--------pathological stage

sn <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$risk)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select( HPV16, Alcohol, Smoking, Gender, DSS_5y_event, PNI, ALI, Age, Subsite, pathologic_stage,cT, cN, margin_status, risk) %>% # selecting your annotations
  as.data.frame(.) #data frame
sn <- sn %>%
  filter(!is.na(sn$pathologic_stage))


sn.tab <- table(sn$pathologic_stage,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

#--------margin_Status--------

sn <- TCGA_clinicaldata %>%
  filter(!is.na(TCGA_clinicaldata$risk)) %>% #selecting only cases where RNA-Seq/Cluster data is available
  dplyr::select( HPV16, Alcohol, Smoking, Gender, DSS_5y_event, PNI, ALI, Age, Subsite, pathologic_stage,cT, cN, margin_status, risk) %>% # selecting your annotations
  as.data.frame(.) #data frame
sn <- sn %>%
  filter(!is.na(sn$margin_status))


sn.tab <- table(sn$margin_status,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)


## Kaplan Meier
Overall_Survival <- read.xlsx("~/Desktop/Doktorarbeit/AXL/KM Survival/KM_Plot__Overall_Survival__(months).xlsx")
column_annotations <- read_xlsx("~/Desktop/m/Column_annotation1.xlsx", sheet = 2)
kmCluster_Death <- read.xlsx("~/Desktop/Doktorarbeit/AXL/KM Survival/kmCluster_Death.xlsx")

attach(kmCluster_Death)

surv_m <- survfit(formula = Surv(OS_STATUS,OS_MONTHS)~CLUSTERS, data = kmCluster_Death)


kmCluster_Death$CLUSTERS, kmCluster_Death$OS_STATUS,

kmCluster_Death = survfit(formula = Surv(CLUSTERS,OS_STATUS),
                          data = kmCluster_Death)

plot(KM)

#kmCluster <- Overall_Survival %>%
#  filter(Overall_Survival$Patients %in% column_annotations$Patients) # selecting important genes from counts 
#write.table(kmCluster, "kmCluster.txt", sep="\t", row.names=FALSE, na="")

kmCluster <- survfit( Surv(OS_months, OS) ~ CLUSTERS, data = kmCluster_Death)

abline(h=0.5, col="red")

plot(kmCluster, 
     conf.int=FALSE,
     xlab = "OS_months",ylab = "Survival Rate",
     yscale = 100,
     las=1,
     lwd=2,
     col =c("red1","blue","green","purple"))

legend("bottomleft",                                      #man kann auch Koordinaten angeben wie (18, 0.95) oder "bottomleft"
       title= "Clusters",
       legend=c("B1","B2b","A1","B2a"),
       col = c("red1","blue","green","purple"),  #man kann auch die Linien zeichnen lassen mit fill/col=c("red1", "blue"...)
       lty=1,
       lwd=2,
       bty = "",
       cex=1)

survdiff(Surv(OS_months, OS) ~ CLUSTERS,data = kmCluster_Death)

