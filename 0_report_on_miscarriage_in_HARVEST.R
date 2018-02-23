

### MISCARRIAGE GWAS META-ANALYSIS (EGCUT) 2018
### used to generate initial report on the phenotype

setwd("~/Dropbox/GIT/MOBA_MISCARRIAGE_META_EGCUT/") # /report_on_miscarriage_in_HARVEST.R

mfr = read.table("~/Biostuff/MOBA_PDB1724/Data/PDB1724_MBRN_460_v9.csv",h=T,sep=",",stringsAsFactors = F)
mfr = mfr[,c("PREG_ID_1724","SPABORT_12_5","SPABORT_23_5","MORS_ALDER","PARITET_5")]
mfr[1:10,]
          
flg = read.table("~/Biostuff/MOBA_PDB1724/sample_flag_list-20170522.txt",h=T,stringsAsFactors = F)
flg = flg[which(flg$coreOK==TRUE),]
flg = flg[which(flg$genotypesOK==TRUE),]
flg = flg[which(flg$phenotypesOK==TRUE),]
flg = flg[,c("IID","ROLE","genotypesOK","coreOK","phenotypesOK")]
head(flg)

key = read.table("~/Biostuff/MOBA_PDB1724/Linkage_PDB1724.csv",h=T,sep=",",stringsAsFactors = F)
key = key[which(key$Role=="Mother"),]
key = key[which(key$SentrixID_1 %in% flg$IID),]
key = key[,c("PREG_ID_1724","SentrixID_1")]

mrg = merge(mfr,key,by="PREG_ID_1724",all = F)

# question about menarche
q1 = read.table("~/Biostuff/MOBA_PDB1724/Data/PDB1724_Q1_v9.csv",sep=",",h=T,stringsAsFactors = F)
q1 = q1[,c("PREG_ID_1724","AA12")]
q1$AA12[which(q1$AA12>20)] = NA
q1$AA12[which(q1$AA12<5)] = NA
dat = merge(mrg,q1,by="PREG_ID_1724",all.x=T)

table(table(dat$SentrixID_1))

table(dat$SPABORT_12_5,useNA = "a")
table(dat$SPABORT_23_5,useNA = "a")


table(par=dat$PARITET_5, misc12=dat$SPABORT_12_5)
table(par=dat$PARITET_5, misc23=dat$SPABORT_23_5)


table(dat$KILDE)
# PARITET_MOR for kilde = 1 is the sum of prior live borns and stillborns as 
# reported by the mother (without specification of gestational age for stillborns). 
# For kilde = 2 and 3 PARITET_MOR is based upon the number of prior live- and stillborns after 12 weeks of gestation. 
# PARITET_MFR is based upon the sum of prior births for this mother, from 16 weeks of gestation in kilde = 1, and 
# from 12 weeks of gestation for kilde = 2 and 3.

library(dplyr)
tbl = group_by(dat,SentrixID_1) %>% 
        summarise(n=n(),
                  mx12=max(SPABORT_12_5,na.rm=T),
                  mx23=max(SPABORT_23_5,na.rm=T),
                  sum12NA=sum(is.na(SPABORT_12_5)),
                  sum23NA=sum(is.na(SPABORT_23_5)),
                  mxAge=max(MORS_ALDER,na.rm=T),
                  menarc = mean(AA12,na.rm=T),
                  mxPAR=max(PARITET_5,na.rm=T)) %>% 
        ungroup()
tbl = as.data.frame(tbl)
table(tbl$n)
table(tbl$mx12,useNA = "a")
table(tbl$mx23,useNA = "a")

table(tbl$mxAge)
hist(tbl$mxAge,breaks=100,col="grey")

### "real" controls

table(tbl$mx12,tbl$sum12NA,useNA = "a")
table(tbl$mx23,tbl$sum23NA,useNA = "a")

table((tbl$mx12==0)&(tbl$mx23==0))


table(parity=tbl$mxPAR, spabort12=tbl$mx12,useNA = "a")
table(parity=tbl$mxPAR, spabort23=tbl$mx23,useNA = "a")




# mothers who have no indication at all for spontaneous abortion
tbl$noInfo = (is.na(tbl$mx12))&(is.na(tbl$mx23))
table(tbl$noInfo)

tbl$any = tbl$mx12 + tbl$mx23
table(tbl$any)
tbl[1:100,]



#### menarche exclusion
hist(tbl$menarc,breaks=100,col="grey")
sum(tbl$menarc<9,na.rm=T)
sum(tbl$menarc>17,na.rm=T)
sum(is.na(tbl$menarc))
nrow(tbl)



##############################################
table(dat$SPABORT_12_5,useNA = "a")
table(dat$SPABORT_23_5,useNA = "a")
table(dat$PARITET_5,useNA = "a")

# SPABORT_12_5 = Previous miscarriages/stillbirths befoore 12 weeks of gestation
# SPABORT_23_5 = Previous miscarriages/stillbirths 12-23 weeks of gestation


### eliminate previous pregnancies of each mother (leave only the latest pregnancy)
dat = dat[order(dat$SentrixID_1,dat$PARITET_5,decreasing = T),]
dat = dat[which(!duplicated(dat$SentrixID_1)),]
head(dat)

table(dat$MORS_ALDER)


table(mfr$FAAR,useNA = "a")
table(mfr$MORS_ALDER,useNA = "a")
