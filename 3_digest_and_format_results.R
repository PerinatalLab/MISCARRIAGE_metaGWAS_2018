#!/usr/bin/Rscript

##  digest and format HARVEST miscarriage metaGWAS results for EGCUT (run by Julius on SNPTEST) 2018-02-22
##  Jonas B.
##  (to be run on HUNT cloud)

## file locations
bim12_file = "/media/local-disk2/jjuod/erc-genotypes-results/postqc/to_imputation/m12-ready-for-imputation.bim"
bim24_file = "/media/local-disk2/jjuod/erc-genotypes-results/postqc/to_imputation/m24-ready-for-imputation.bim"
rez_suffix = "/media/local-disk2/jjuod/other/results_meta/" #/moms_miscarr1_1.txt        
out_suffix = "./results/MOBA-HARVEST_pheno"  # final example: "MOBA-HARVEST_pheno1_20180222_jb.txt.gz"

### prepare a list of genotyped SNPs (positions)
bim12 = read.table(bim12_file,h=F,stringsAsFactors = F)
bim12$chrpos = paste(bim12$V1,bim12$V4,sep="_")
bim24 = read.table(bim24_file,h=F,stringsAsFactors = F)
bim24$chrpos = paste(bim24$V1,bim24$V4,sep="_")
# for restoring SNP IDs when they are missing in the results
bim = rbind(bim12[,c(2,7)],bim24[,c(2,7)])
bim = bim[which(!duplicated(bim$V2)),]
# list of genotyped SNP names and positions
genotyped_snps = unique(bim12$V2[which(bim12$V2 %in% bim24$V2)])  # strict (only with genotypes on both batches)
genotyped_chrpos = unique(bim12$chrpos[which(bim12$chrpos %in% bim24$chrpos)]) # strict (only with genotypes on both batches)


#length(genotyped_snps)
#length(genotyped_chrpos)
#table(bim12$V1,useNA = "a")
#table(substr(bim12$V2,1,2),useNA = "a")
rm(bim12,bim24)


for (pheno in c(1,2,4)) {   # which phenotype is being digested
        
        print(paste("  DIGESTING PHENOTYPE ",pheno,"...",sep=""))
        
        for (chr in c(1:22,"X")) {
                print(paste("digesting chr",chr,"...",sep=""))
                
                # read-in the results
                #a = read.table("~/Biostuff/MOBA_META_MISCARRIAGE/moms_miscarr1_22.txt",h=T,stringsAsFactors = F) # locally
                file_name = paste(rez_suffix,"moms_miscarr",pheno,"_",chr,".txt",sep="")
                a = read.table(file_name,h=T,stringsAsFactors = F)
                
                # special case for chromosome coding:
                if (chr=="X") { a$chromosome = 23 }
                
                # estimate effect-allele frequency (in SNPTEST - B allele frequency)
                a$bfreq = (a$cohort_1_BB*2 + a$cohort_1_AB) / ((a$cohort_1_AA + a$cohort_1_AB + a$cohort_1_BB)*2)
                
                # issing variables
                a$STRAND = "+"
                a$BUILD = "b37"
                a$IMPUTED = 1
                a$INFO_TYPE = 3
                
                a$chrpos = paste(a$chromosome,a$position,sep="_")
                
                # sum(a$rsid  %in% genotyped_snps)
                # sum(a$chrpos %in% genotyped_chrpos)
                # sum(a$info ==1)
                # 
                # ix = which(a$rsid  %in% genotyped_snps)
                # sum(a$info[ix]!=1) #hist(a$info[ix],breaks=100)
                # ix = which(a$chrpos %in% genotyped_chrpos)
                # sum(a$info[ix]!=1) #hist(a$info[ix],breaks=100)
                # ix = which( (a$chrpos %in% genotyped_chrpos)&(a$rsid %in% genotyped_snps))
                # sum(a$info[ix]!=1) #hist(a$info[ix],breaks=100)
                
                # mark those SNPs that are suspected to be genotyped
                a$IMPUTED[which( (a$chrpos %in% genotyped_chrpos)&(a$info==1)&(a$rsid %in% genotyped_snps) )] = 0
                #table(a$IMPUTED,useNA = "a")
                
                
                ## restore SNP names when those are missing
                a = a[,c("chromosome","position","rsid", "STRAND","BUILD","all_total",
                         "alleleB","alleleA","bfreq","frequentist_add_beta_1","all_OR",
                         "frequentist_add_se_1","frequentist_add_pvalue","IMPUTED","INFO_TYPE","info","chrpos")]
                
                tmp = merge(a,bim,by="chrpos",all.x=T)
                a = tmp; rm(tmp)
                ix = which((a$rsid==".")&(substr(a$V2,1,2)=="rs"))
                a$rsid[ix] = a$V2[ix]
                rm(ix)
                a = a[order(a$position),]
                
                # restore column order as it should be
                a = a[,c("chromosome","position","rsid", "STRAND","BUILD","all_total",
                         "alleleB","alleleA","bfreq","frequentist_add_beta_1","all_OR",
                         "frequentist_add_se_1","frequentist_add_pvalue","IMPUTED","INFO_TYPE","info")]
                
                # correct order of columns
                colnames(a) = c("CHR","POS","SNP","SRAND","BUILD","N","EFFECT_ALLELE","NON_EFFECT_ALLELE",
                                "EFFECT_ALLELE_FREQ","BETA","OR","SE","P","IMPUTED","INFO_TYPE","INFO")
                
                a$EFFECT_ALLELE_FREQ = format(a$EFFECT_ALLELE_FREQ,scientific = F,digits = 2) # not explicitly requested
                #a$BETA = format(a$BETA,scientific = T,trim = T,nsmall = 4)  # not explicitly requested
                #a$OR = format(a$OR,scientific = T,trim = T,nsmall = 4)  # not explicitly requested
                #a$SE = format(a$SE,scientific = T,trim = T,nsmall = 4)  # not explicitly requested
                a$P = format(a$P,scientific = T,trim = T,nsmall = 4) # requested!
                
                # recode missing values
                a$BETA[which((a$BETA=="NA")|(is.na(a$BETA)))] = "."
                a$OR[which((a$OR=="NA")|(is.na(a$OR)))] = "."
                a$SE[which((a$SE=="NA")|(is.na(a$SE)))] = "."
                a$P[which((a$P=="NA")|(is.na(a$P)))] = "."
                
                
                # save the results
                output_name = paste(out_suffix,pheno,"_20180222_jb.txt",sep="")
                if(chr == 1) {
                        write.table(a,output_name,row.names = F,col.names = T,quote=F,sep="\t",append = FALSE)        
                } else {
                        write.table(a,output_name,row.names = F,col.names = F,quote=F,sep="\t",append = TRUE)        
                }
                
                # detect problems on the go
                # print(paste("chr = ",a$CHR[1]))
                # print(paste("tail(a$POS) = ",tail(a$POS)[6]))
                # print(paste("class(a$POS) = ",class(a$POS)))
                # print("table(substr(a$POS,1,1)) : ")
                # print(table(substr(a$POS,1,1),useNA = "a"))
                # print("table(substr(a$POS,2,2)) : ")
                # print(table(substr(a$POS,2,2),useNA = "a"))
                # print(paste("max(a$POS) = ",max(a$POS)))
                
                rm(a,file_name,output_name)

        } # end of one chromosome
} # end of one phenotype


