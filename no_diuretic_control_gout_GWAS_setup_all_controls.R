# Murray Cadzow
# University of Otago
# 21 April 2016

# script for creating normal gout GWAS using all cases for each gout criteria and all controls
#control criteria:  no diuretics
#gout criteria: all, self report, self report + ULT, winnard, hospital

# run ukbiobank_gout.R first to create required affection columns
# only change affection status based on gout/control criteria, everything else is -9 and 
# plink will handle subsetting during analysis

#GWAS_all_controls/{controls_no_diuretics,controls_diuretics,controls}/{all,winnard,hosp,all_male,hosp_male,self,self_ult}/{adjusted,unadjusted}



xsan <- '/media/xsan/'

# read in fam file for sample filtered plink files
data_dir <- paste0(xsan,"/staff_groups/merrimanlab/Documents/Murray/ukbiobank_util/data/")
scratch_dir <- paste0(xsan,"/scratch/merrimanlab/murray/working_dir/UkBio/")
control_cond <- 'controls_no_diuretics' 

load(paste0(data_dir,"ukbiobank_genotyped2016-04-26.RData"))

genotyped$waist_to_height_ratio <- genotyped$f.48.0.0 / genotyped$f.50.0.0


fam_file <- read.table(paste0(scratch_dir,"chr1impv1.fam"), header=FALSE, stringsAsFactors = FALSE)
colnames(fam_file) <- c("FID","IID","PID","MID","SEX","AFF")
fam_file$sort <- as.numeric(rownames(fam_file))

# f.210000.0.0 = ethnicity
# f.31.0.0 = sex
# f.21003.0.0 = age
# f.21001.0.0 = bmi
# f.48.0.0 = waist
# f.50.0.0 = height
# f.22001.0.0 = genetic sex


#merge
new_fam_file <- merge(fam_file, genotyped[,c("f.eid", "f.21000.0.0", "f.31.0.0", "f.21003.0.0","f.21001.0.0",'waist_to_height_ratio', 'f.48.0.0', 'f.50.0.0','f.22001.0.0', "goutaff", "control", "goutwinnard","goutself", "gouthosp","gout_winnard_self","goutall","goutult","gout_self_ult", 'diuretics', colnames(genotyped)[grep('22009',colnames(genotyped))] )], by.x = "IID", by.y = "f.eid", all.x=TRUE)
#resort
new_fam_file <- new_fam_file[order(new_fam_file$sort),]
#check lengths are equal
length(fam_file[,1]) == length(new_fam_file[,1])



#reset case/control as precaution
new_fam_file$AFF <- -9
new_fam_file$SEX <- -9
#in ukbio males are coded as female = 0, males = 1
new_fam_file[new_fam_file$f.31.0.0 == 1 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 1
new_fam_file[new_fam_file$f.31.0.0 == 0 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 2

# affstat for all gout
gout_cond <- 'all'
new_fam_file[new_fam_file$goutall == 1 & !is.na(new_fam_file$goutall), "AFF"] <- 2
new_fam_file[new_fam_file$control == 1 & !is.na(new_fam_file$control) & is.na(new_fam_file$diuretics),"AFF"] <- 1

table(new_fam_file$AFF, new_fam_file$f.21000.0.0, exclude=NULL)
# blank out non-white ethnicities
new_fam_file[!(!is.na(new_fam_file$f.21000.0.0) & (new_fam_file$f.21000.0.0 == 1001 | new_fam_file$f.21000.0.0 == 1002 | new_fam_file$f.21000.0.0 == 1003)) , "AFF"] <- -9

write.table(new_fam_file[,c("FID","IID","PID","MID","SEX","AFF")], file = paste0(scratch_dir,"GWAS_all_controls/",control_cond,'/',gout_cond,'/chrallimpv1.fam_',gout_cond), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = ' ')


### all_male
#reset case/control as precaution
new_fam_file$AFF <- -9
new_fam_file$SEX <- -9
#in ukbio males are coded as female = 0, males = 1
new_fam_file[new_fam_file$f.31.0.0 == 1 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 1
new_fam_file[new_fam_file$f.31.0.0 == 0 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- NA

# affstat for all gout
gout_cond <- 'all_male'
new_fam_file[new_fam_file$goutall == 1 & !is.na(new_fam_file$goutall), "AFF"] <- 2
new_fam_file[new_fam_file$control == 1 & !is.na(new_fam_file$control) & is.na(new_fam_file$diuretics),"AFF"] <- 1

table(new_fam_file$AFF, new_fam_file$f.21000.0.0, exclude=NULL)
# blank out non-white ethnicities
new_fam_file[!(!is.na(new_fam_file$f.21000.0.0) & (new_fam_file$f.21000.0.0 == 1001 | new_fam_file$f.21000.0.0 == 1002 | new_fam_file$f.21000.0.0 == 1003)) , "AFF"] <- -9

write.table(new_fam_file[,c("FID","IID","PID","MID","SEX","AFF")], file = paste0(scratch_dir,"GWAS_all_controls/",control_cond,'/',gout_cond,'/chrallimpv1.fam_',gout_cond), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = ' ')



#### fam gout hospital
#reset case/control as precaution
new_fam_file$AFF <- -9
new_fam_file$SEX <- -9
#in ukbio males are coded as female = 0, males = 1
new_fam_file[new_fam_file$f.31.0.0 == 1 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 1
new_fam_file[new_fam_file$f.31.0.0 == 0 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 2

# affstat for hosp gout
gout_cond <- 'hosp'
new_fam_file[new_fam_file$gouthosp == 1 & !is.na(new_fam_file$gouthosp), "AFF"] <- 2
new_fam_file[new_fam_file$control == 1 & !is.na(new_fam_file$control) & is.na(new_fam_file$diuretics),"AFF"] <- 1

table(new_fam_file$AFF, new_fam_file$f.21000.0.0, exclude=NULL)
# blank out non-white ethnicities
new_fam_file[!(!is.na(new_fam_file$f.21000.0.0) & (new_fam_file$f.21000.0.0 == 1001 | new_fam_file$f.21000.0.0 == 1002 | new_fam_file$f.21000.0.0 == 1003)) , "AFF"] <- -9

write.table(new_fam_file[,c("FID","IID","PID","MID","SEX","AFF")], file = paste0(scratch_dir,"GWAS_all_controls/",control_cond,'/',gout_cond,'/chrallimpv1.fam_',gout_cond), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = ' ')



#### fam gout hospital male
#reset case/control as precaution
new_fam_file$AFF <- -9
new_fam_file$SEX <- -9
#in ukbio males are coded as female = 0, males = 1
new_fam_file[new_fam_file$f.31.0.0 == 1 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 1
new_fam_file[new_fam_file$f.31.0.0 == 0 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- NA

# affstat for hosp gout
gout_cond <- 'hosp_male'
new_fam_file[new_fam_file$gouthosp == 1 & !is.na(new_fam_file$gouthosp), "AFF"] <- 2
new_fam_file[new_fam_file$control == 1 & !is.na(new_fam_file$control) & is.na(new_fam_file$diuretics),"AFF"] <- 1

table(new_fam_file$AFF, new_fam_file$f.21000.0.0, exclude=NULL)
# blank out non-white ethnicities
new_fam_file[!(!is.na(new_fam_file$f.21000.0.0) & (new_fam_file$f.21000.0.0 == 1001 | new_fam_file$f.21000.0.0 == 1002 | new_fam_file$f.21000.0.0 == 1003)) , "AFF"] <- -9

write.table(new_fam_file[,c("FID","IID","PID","MID","SEX","AFF")], file = paste0(scratch_dir,"GWAS_all_controls/",control_cond,'/',gout_cond,'/chrallimpv1.fam_',gout_cond), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = ' ')



#### fam gout winnard

#reset sex and affstat
#reset case/control as precaution
new_fam_file$AFF <- -9
new_fam_file$SEX <- -9
#in ukbio males are coded as female = 0, males = 1
new_fam_file[new_fam_file$f.31.0.0 == 1 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 1
new_fam_file[new_fam_file$f.31.0.0 == 0 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 2

# affstat for winnard gout
gout_cond <- 'winnard'
new_fam_file[new_fam_file$goutwinnard == 1 & !is.na(new_fam_file$goutwinnard), "AFF"] <- 2
new_fam_file[new_fam_file$control == 1 & !is.na(new_fam_file$control) & is.na(new_fam_file$diuretics),"AFF"] <- 1

table(new_fam_file$AFF, new_fam_file$f.21000.0.0, exclude=NULL)
# blank out non-white ethnicities
new_fam_file[!(!is.na(new_fam_file$f.21000.0.0) & (new_fam_file$f.21000.0.0 == 1001 | new_fam_file$f.21000.0.0 == 1002 | new_fam_file$f.21000.0.0 == 1003)) , "AFF"] <- -9

write.table(new_fam_file[,c("FID","IID","PID","MID","SEX","AFF")], file = paste0(scratch_dir,"GWAS_all_controls/",control_cond,'/',gout_cond,'/chrallimpv1.fam_',gout_cond), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = ' ')



#### fam gout self

#reset sex and affstat
#reset case/control as precaution
new_fam_file$AFF <- -9
new_fam_file$SEX <- -9
#in ukbio males are coded as female = 0, males = 1
new_fam_file[new_fam_file$f.31.0.0 == 1 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 1
new_fam_file[new_fam_file$f.31.0.0 == 0 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 2

# affstat for self reported gout
gout_cond <- 'self'
new_fam_file[new_fam_file$goutself == 1 & !is.na(new_fam_file$goutself), "AFF"] <- 2
new_fam_file[new_fam_file$control == 1 & !is.na(new_fam_file$control) & is.na(new_fam_file$diuretics),"AFF"] <- 1

table(new_fam_file$AFF, new_fam_file$f.21000.0.0, exclude=NULL)
# blank out non-white ethnicities
new_fam_file[!(!is.na(new_fam_file$f.21000.0.0) & (new_fam_file$f.21000.0.0 == 1001 | new_fam_file$f.21000.0.0 == 1002 | new_fam_file$f.21000.0.0 == 1003)) , "AFF"] <- -9

write.table(new_fam_file[,c("FID","IID","PID","MID","SEX","AFF")], file = paste0(scratch_dir,"GWAS_all_controls/",control_cond,'/',gout_cond,'/chrallimpv1.fam_',gout_cond), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = ' ')



### self defined or ULT
#reset case/control as precaution
new_fam_file$AFF <- -9
new_fam_file$SEX <- -9
#in ukbio males are coded as female = 0, males = 1
new_fam_file[new_fam_file$f.31.0.0 == 1 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 1
new_fam_file[new_fam_file$f.31.0.0 == 0 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 2

# affstat for self + ULT gout
gout_cond <- 'self_ult'
new_fam_file[new_fam_file$gout_self_ult == 1 & !is.na(new_fam_file$gout_self_ult), "AFF"] <- 2
new_fam_file[new_fam_file$control == 1 & !is.na(new_fam_file$control) & is.na(new_fam_file$diuretics),"AFF"] <- 1

table(new_fam_file$AFF, new_fam_file$f.21000.0.0, exclude=NULL)
# blank out non-white ethnicities
new_fam_file[!(!is.na(new_fam_file$f.21000.0.0) & (new_fam_file$f.21000.0.0 == 1001 | new_fam_file$f.21000.0.0 == 1002 | new_fam_file$f.21000.0.0 == 1003)) , "AFF"] <- -9

write.table(new_fam_file[,c("FID","IID","PID","MID","SEX","AFF")], file = paste0(scratch_dir,"GWAS_all_controls/",control_cond,'/',gout_cond,'/chrallimpv1.fam_',gout_cond), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = ' ')


### only self, not classified any other way
#reset sex and affstat
#reset case/control as precaution
new_fam_file$AFF <- -9
new_fam_file$SEX <- -9
#in ukbio males are coded as female = 0, males = 1
new_fam_file[new_fam_file$f.31.0.0 == 1 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 1
new_fam_file[new_fam_file$f.31.0.0 == 0 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 2

# affstat for self reported gout
gout_cond <- 'self_only'
new_fam_file[new_fam_file$goutself == 1 & !is.na(new_fam_file$goutself), "AFF"] <- 2
new_fam_file[new_fam_file$goutwinnard == 1 & !is.na(new_fam_file$goutwinnard), "AFF"] <- NA #includes ULT too
new_fam_file[new_fam_file$gouthosp == 1 & !is.na(new_fam_file$gouthosp), "AFF"] <- NA
new_fam_file[new_fam_file$control == 1 & !is.na(new_fam_file$control),"AFF"] <- 1

table(new_fam_file$AFF, new_fam_file$f.21000.0.0, exclude=NULL)
# blank out non-white ethnicities
new_fam_file[!(!is.na(new_fam_file$f.21000.0.0) & (new_fam_file$f.21000.0.0 == 1001 | new_fam_file$f.21000.0.0 == 1002 | new_fam_file$f.21000.0.0 == 1003)) , "AFF"] <- -9

write.table(new_fam_file[,c("FID","IID","PID","MID","SEX","AFF")], file = paste0(scratch_dir,"GWAS_all_controls/",control_cond,'/',gout_cond,'/chrallimpv1.fam_',gout_cond), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = ' ')


### ULT
#reset case/control as precaution
new_fam_file$AFF <- -9
new_fam_file$SEX <- -9
#in ukbio males are coded as female = 0, males = 1
new_fam_file[new_fam_file$f.31.0.0 == 1 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 1
new_fam_file[new_fam_file$f.31.0.0 == 0 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 2

# affstat for self + ULT gout
gout_cond <- 'ult'
new_fam_file[new_fam_file$goutult == 1 & !is.na(new_fam_file$goutult), "AFF"] <- 2
new_fam_file[new_fam_file$control == 1 & !is.na(new_fam_file$control) & new_fam_file$diuretics == 1 & !is.na(new_fam_file$diuretics),"AFF"] <- 1

table(new_fam_file$AFF, new_fam_file$f.21000.0.0, exclude=NULL)
# blank out non-white ethnicities
new_fam_file[!(!is.na(new_fam_file$f.21000.0.0) & (new_fam_file$f.21000.0.0 == 1001 | new_fam_file$f.21000.0.0 == 1002 | new_fam_file$f.21000.0.0 == 1003)) , "AFF"] <- -9

write.table(new_fam_file[,c("FID","IID","PID","MID","SEX","AFF")], file = paste0(scratch_dir,"GWAS_all_controls/",control_cond,'/',gout_cond,'/chrallimpv1.fam_',gout_cond), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = ' ')


# only need to create this once and needs to include everyone
co_var_file <- new_fam_file[,c("FID","IID","f.21003.0.0", "f.21001.0.0", 'waist_to_height_ratio', 'f.48.0.0', colnames(new_fam_file)[grep('22009', colnames(new_fam_file))])]
colnames(co_var_file) <- c("FID","IID","AGE","BMI", "WaistHeightRatio", 'Waist', paste0('PCA',1:15))
write.table(co_var_file, file = paste0(scratch_dir,"GWAS_all_controls/chrallimpv1.covar"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep = ' ')
