
### MISCARRIAGE GWAS META-ANALYSIS (organizer: EGCUT)
### prepare phenotype files
### 2018-02-20. Jonas B.

setwd("~/Dropbox/GIT/MOBA_MISCARRIAGE_META_EGCUT/") # 1_prepare_phenotype_file.R

####  reproducibility section:  code-hash and time-stamp
report = system("git status",intern = T)
if ( length(grep("modified",report))>0) {
        git_hash = "--???--" # warn user that script does not represent the contents
} else {
        git_hash = system("git log --pretty=format:'%h' -n 1",intern = T)  # get recent git commit version
}
time_stamp = gsub("-","",substr(Sys.time(),1,10))



mfr = read.table("~/Biostuff/MOBA_PDB1724/Data/PDB1724_MBRN_460_v9.csv",h=T,sep=",",stringsAsFactors = F)
mfr = mfr[,c("PREG_ID_1724","SPABORT_12_5","SPABORT_23_5","MORS_ALDER","PARITET_5","FAAR","KILDE")]
#mfr[1:10,]

# estimate maternal birth year
mfr$MBY = mfr$FAAR - mfr$MORS_ALDER  # -> MatBirtyYear


flg = read.table("~/Biostuff/MOBA_PDB1724/sample_flag_list-20170522.txt",h=T,stringsAsFactors = F)
flg = flg[which(flg$coreOK==TRUE),]
flg = flg[which(flg$genotypesOK==TRUE),]
flg = flg[which(flg$phenotypesOK==TRUE),]
flg = flg[,c("IID","ROLE","genotypesOK","coreOK","phenotypesOK")]
#head(flg)

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

table(table(dat$SentrixID_1))  # multiple pregs of the same mother are present (ok)


table(dat$SPABORT_12_5,useNA = "a")
table(dat$SPABORT_23_5,useNA = "a")
table(par=dat$PARITET_5, misc12=dat$SPABORT_12_5)
table(par=dat$PARITET_5, misc23=dat$SPABORT_23_5)


library(dplyr)
tbl = group_by(dat,SentrixID_1) %>% 
        summarise(n=n(),
                  mx12=max(SPABORT_12_5,na.rm=T),
                  mx23=max(SPABORT_23_5,na.rm=T),
                  sum12NA=sum(is.na(SPABORT_12_5)),
                  sum23NA=sum(is.na(SPABORT_23_5)),
                  mxAge=max(MORS_ALDER,na.rm=T),
                  menarc = mean(AA12,na.rm=T),
                  MatBirthYear = mean(MBY,na.rm=T),
                  mxPAR=max(PARITET_5,na.rm=T)) %>% 
        ungroup()
tbl = as.data.frame(tbl)
hist(tbl$MatBirthYear,breaks=1000,col="grey")  # maternal birth year

### apply exclusion criteria
tbl = tbl[which((tbl$menarc>=9)&(tbl$menarc<=17)),]  


#####################
#####################  PHENOTYPE 1
#####################



table(type12=tbl$mx12,type23=tbl$mx23,useNA = "a")
table(tbl$mxPAR[which((tbl$mx12==0)&(tbl$mx23==0))],useNA = "a")

# define cases:
ix1 = which( (tbl$mx12==0) & (tbl$mx23 %in% c(1,2)) )
ix2 = which( (tbl$mx12==1) & ((tbl$mx23 %in% c(0,1))|(is.na(tbl$mx23))) )
ix3 = which( (tbl$mx12==2) & ((tbl$mx23  == 0)|(is.na(tbl$mx23))) )
ix4 = which( (is.na(tbl$mx12)) & (tbl$mx23 %in% c(1,2)) )
case_ix = unique(c(ix1,ix2,ix3,ix4))
length(case_ix)

# define controls:
contr_ix = which((tbl$mxPAR>0)&(tbl$mx12==0)&(tbl$mx23==0))
length(contr_ix)

cases = tbl$SentrixID_1[case_ix]
contr = tbl$SentrixID_1[contr_ix]

#####
##### prepare PHENOTYPE files
#####

##
##  AUTOSOMES
##

phe = read.table("~/Biostuff/mount_hunt/GWAS_MISCARRIAGE/momIDs_autosom.txt",h=F,stringsAsFactors = F)
colnames(phe)[1] = "ID_1"
phe$ID_2 = phe$ID_1
phe$missing = 0
phe$miscrr = NA
phe$miscrr[which(phe$ID_1 %in% cases)] = 1
phe$miscrr[which(phe$ID_1 %in% contr)] = 0
phe$seq = seq(nrow(phe))

mby = tbl[,c("SentrixID_1","MatBirthYear")]
tmp = merge(phe,mby,by.x="ID_1",by.y="SentrixID_1",all.x=T)
tmp = tmp[order(tmp$seq),]
phe_autosomes = tmp[,c("ID_1","ID_2","missing","miscrr","MatBirthYear")]

##
##  X-chromosome
##

phe = read.table("~/Biostuff/mount_hunt/GWAS_MISCARRIAGE/momIDs_xchromo.txt",h=F,stringsAsFactors = F)
colnames(phe)[1] = "ID_1"
phe$ID_2 = phe$ID_1
phe$missing = 0
phe$miscrr = NA
phe$miscrr[which(phe$ID_1 %in% cases)] = 1
phe$miscrr[which(phe$ID_1 %in% contr)] = 0
phe$seq = seq(nrow(phe))

mby = tbl[,c("SentrixID_1","MatBirthYear")]
tmp = merge(phe,mby,by.x="ID_1",by.y="SentrixID_1",all.x=T)
tmp = tmp[order(tmp$seq),]
phe_xchromosome = tmp[,c("ID_1","ID_2","missing","miscrr","MatBirthYear")]

#####
##### export pheno files
#####

file_name_A = paste("phe1_auto_nCaCo-",
                  sum(phe_autosomes$miscrr==1,na.rm = T),"-",
                  sum(phe_autosomes$miscrr==0,na.rm = T),"_",
                  git_hash,"_",time_stamp,"_JB.txt",sep="")
write.table(phe_autosomes,file_name_A,row.names = F,col.names = T,sep="\t",quote=F)  # local save

file_name_X = paste("phe1_xchr_nCaCo-",
                    sum(phe_xchromosome$miscrr==1,na.rm = T),"-",
                    sum(phe_xchromosome$miscrr==0,na.rm = T),"_",
                    git_hash,"_",time_stamp,"_JB.txt",sep="")
write.table(phe_xchromosome,file_name_X,row.names = F,col.names = T,sep="\t",quote=F)  # local save

# note, for SNPTEST a second header must be added manually (by Julius, after adding PCs) !


#####################
#####################  PHENOTYPE 2:  3 or more miscarriages
#####################

table(type12=tbl$mx12,type23=tbl$mx23,useNA = "a")
table(tbl$mxPAR[which((tbl$mx12==0)&(tbl$mx23==0))],useNA = "a")

# define cases:
ix1 = which( (tbl$mx12==0) & (tbl$mx23 %in% c(3,4)) )
ix2 = which( (tbl$mx12==1) & (tbl$mx23 %in% c(2,3,4)) )
ix3 = which( (tbl$mx12==2) & (tbl$mx23 %in% c(1,2,3,4)) )
ix4 = which(  tbl$mx12==3 )
ix5 = which(  tbl$mx12==4 )
ix6 = which( (is.na(tbl$mx12)) & (tbl$mx23 %in% c(3,4)) )
case_ix = unique(c(ix1,ix2,ix3,ix4,ix5,ix6))
length(case_ix)

# define controls:
contr_ix = which((tbl$mxPAR>0)&(tbl$mx12==0)&(tbl$mx23==0))
length(contr_ix)

cases = tbl$SentrixID_1[case_ix]
contr = tbl$SentrixID_1[contr_ix]

#####
##### prepare PHENOTYPE files
#####

##
##  AUTOSOMES
##

phe = read.table("~/Biostuff/mount_hunt/GWAS_MISCARRIAGE/momIDs_autosom.txt",h=F,stringsAsFactors = F)
colnames(phe)[1] = "ID_1"
phe$ID_2 = phe$ID_1
phe$missing = 0
phe$miscrr = NA
phe$miscrr[which(phe$ID_1 %in% cases)] = 1
phe$miscrr[which(phe$ID_1 %in% contr)] = 0
phe$seq = seq(nrow(phe))

mby = tbl[,c("SentrixID_1","MatBirthYear")]
tmp = merge(phe,mby,by.x="ID_1",by.y="SentrixID_1",all.x=T)
tmp = tmp[order(tmp$seq),]
phe_autosomes = tmp[,c("ID_1","ID_2","missing","miscrr","MatBirthYear")]

##
##  X-chromosome
##

phe = read.table("~/Biostuff/mount_hunt/GWAS_MISCARRIAGE/momIDs_xchromo.txt",h=F,stringsAsFactors = F)
colnames(phe)[1] = "ID_1"
phe$ID_2 = phe$ID_1
phe$missing = 0
phe$miscrr = NA
phe$miscrr[which(phe$ID_1 %in% cases)] = 1
phe$miscrr[which(phe$ID_1 %in% contr)] = 0
phe$seq = seq(nrow(phe))

mby = tbl[,c("SentrixID_1","MatBirthYear")]
tmp = merge(phe,mby,by.x="ID_1",by.y="SentrixID_1",all.x=T)
tmp = tmp[order(tmp$seq),]
phe_xchromosome = tmp[,c("ID_1","ID_2","missing","miscrr","MatBirthYear")]

#####
##### export pheno files
#####

file_name_A = paste("phe2_auto_nCaCo-",
                    sum(phe_autosomes$miscrr==1,na.rm = T),"-",
                    sum(phe_autosomes$miscrr==0,na.rm = T),"_",
                    git_hash,"_",time_stamp,"_JB.txt",sep="")
write.table(phe_autosomes,file_name_A,row.names = F,col.names = T,sep="\t",quote=F)  # local save

file_name_X = paste("phe2_xchr_nCaCo-",
                    sum(phe_xchromosome$miscrr==1,na.rm = T),"-",
                    sum(phe_xchromosome$miscrr==0,na.rm = T),"_",
                    git_hash,"_",time_stamp,"_JB.txt",sep="")
write.table(phe_xchromosome,file_name_X,row.names = F,col.names = T,sep="\t",quote=F)  # local save

# note, for SNPTEST a second header must be added manually (by Julius, after adding PCs) !



#####################
#####################  PHENOTYPE 4:  at least 1 miscarriage
#####################

table(type12=tbl$mx12,type23=tbl$mx23,useNA = "a")
table(tbl$mxPAR[which((tbl$mx12==0)&(tbl$mx23==0))],useNA = "a")

# define cases:
ix = which( (tbl$mx12 %in% c(1,2,3,4)) | (tbl$mx23 %in% c(1,2,3,4)) )
case_ix = unique(ix)
length(case_ix)

# define controls:
contr_ix = which((tbl$mxPAR>0)&(tbl$mx12==0)&(tbl$mx23==0))
length(contr_ix)

cases = tbl$SentrixID_1[case_ix]
contr = tbl$SentrixID_1[contr_ix]

#####
##### prepare PHENOTYPE files
#####

##
##  AUTOSOMES
##

phe = read.table("~/Biostuff/mount_hunt/GWAS_MISCARRIAGE/momIDs_autosom.txt",h=F,stringsAsFactors = F)
colnames(phe)[1] = "ID_1"
phe$ID_2 = phe$ID_1
phe$missing = 0
phe$miscrr = NA
phe$miscrr[which(phe$ID_1 %in% cases)] = 1
phe$miscrr[which(phe$ID_1 %in% contr)] = 0
phe$seq = seq(nrow(phe))

mby = tbl[,c("SentrixID_1","MatBirthYear")]
tmp = merge(phe,mby,by.x="ID_1",by.y="SentrixID_1",all.x=T)
tmp = tmp[order(tmp$seq),]
phe_autosomes = tmp[,c("ID_1","ID_2","missing","miscrr","MatBirthYear")]

##
##  X-chromosome
##

phe = read.table("~/Biostuff/mount_hunt/GWAS_MISCARRIAGE/momIDs_xchromo.txt",h=F,stringsAsFactors = F)
colnames(phe)[1] = "ID_1"
phe$ID_2 = phe$ID_1
phe$missing = 0
phe$miscrr = NA
phe$miscrr[which(phe$ID_1 %in% cases)] = 1
phe$miscrr[which(phe$ID_1 %in% contr)] = 0
phe$seq = seq(nrow(phe))

mby = tbl[,c("SentrixID_1","MatBirthYear")]
tmp = merge(phe,mby,by.x="ID_1",by.y="SentrixID_1",all.x=T)
tmp = tmp[order(tmp$seq),]
phe_xchromosome = tmp[,c("ID_1","ID_2","missing","miscrr","MatBirthYear")]

#####
##### export pheno files
#####

file_name_A = paste("phe4_auto_nCaCo-",
                    sum(phe_autosomes$miscrr==1,na.rm = T),"-",
                    sum(phe_autosomes$miscrr==0,na.rm = T),"_",
                    git_hash,"_",time_stamp,"_JB.txt",sep="")
write.table(phe_autosomes,file_name_A,row.names = F,col.names = T,sep="\t",quote=F)  # local save

file_name_X = paste("phe4_xchr_nCaCo-",
                    sum(phe_xchromosome$miscrr==1,na.rm = T),"-",
                    sum(phe_xchromosome$miscrr==0,na.rm = T),"_",
                    git_hash,"_",time_stamp,"_JB.txt",sep="")
write.table(phe_xchromosome,file_name_X,row.names = F,col.names = T,sep="\t",quote=F)  # local save

# note, for SNPTEST a second header must be added manually (by Julius, after adding PCs) !