# Murray Cadzow
# University of Otago
# 4 March 2016

# script for creating diuretic controls gout GWAS

# run ukbiobank_gout.R first to create matches between gout and control - after pair_controls() function (about line 314)
## loop through 100 times to create different selected controls
# only change affection status based on gout criteria, everything else is -9 and 
# plink will handle subsetting during analysis



####
####
#### CONTROLS WITH ON DIURETICS
####

#requires gout_list, control_ids from when matching all controls above

## match people that currently have genotypes
#read in fam file again (just in case measure)
fam_file <- read.table("/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/UkBio/chr1impv1.fam", header=FALSE, stringsAsFactors = FALSE)
colnames(fam_file) <- c("FID","IID","PID","MID","SEX","AFF")
fam_file$sort <- as.numeric(rownames(fam_file))


## we keep diuretic controls inside the loop
ids_on_diuretics <- bd_refined[rows_w_diuretics,'f.eid']

#special matching for if can't pick 3, (loosen eth criteria)
pair_controls1 <- function(id){
  samples <- control_age_sex_eth_match_chr4[[id]]
  samples <- control_ids[which(control_ids$ids %in% samples), ]
  samples <- samples[samples$avail == TRUE, "ids"]
  if(length(samples) >= 3){
    s <- sample(samples, replace = FALSE, size = 3)
  } else {
    s <- samples
    case <- bd_refined[bd_refined$f.eid == id, c("f.eid", "f.31.0.0", "f.21000.0.0", "f.21003.0.0")]
    #select same age, and bmi similar eth
    v <- which( bd_refined$f.31.0.0 == case$f.31.0.0 & #same sex
                    bd_refined$f.21003.0.0 %in% c(case$f.21003.0.0-1,case$f.21003.0.0, case$f.21003.0.0+1) & #similar BMI
                    bd_refined$f.21000.0.0 %in% c(1001,1002,1003) & #similar eth
                    bd_refined$f.eid %in% control_ids[control_ids$avail == TRUE, 'ids' ] & # is an available id
                    bd_refined$goutaff == 0)  #make sure control 
    if(length(v) >= 3-length(s)){
    s <- c(s, sample(bd_refined[v,"f.eid"], replace=FALSE, size = 3 - length(s) ))
    }
  }
  control_ids[which(control_ids$ids %in% s),"avail"] <<- FALSE
  return(s)
  
}

## loop through 100 times to create different selected controls
# only change affection status based on gout criteria, everything else is -9 and 
# plink will handle subsetting during analysis

for (i in 0:100){
  print(i)
  ### set up start
  control_ids$avail <- FALSE
  control_ids[control_ids$ids %in% fam_file$IID, 'avail'] <-TRUE
  #remove non diuretic controls
  control_ids[!(control_ids$ids %in% ids_on_diuretics), 'avail'] <-FALSE 
  prior <-table(control_ids$avail)["FALSE"]
  
  #do pairing (randomly select 3 matched controls for each case)
  paired_list <- list()
  for(id in gout_list$ids){
    paired_list[[id]] <- pair_controls1(id)
  }
  
  # for each case include the number of controls selected to use as a check
  gout_list$pair_len <-NA
  for (case in genotyped[which(genotyped$goutaff == 1 & !is.na(genotyped$goutaff)),"f.eid"]){
    gout_list[gout_list$ids == case, "pair_len"] <- length(paired_list[[case]])
  }
  
  #sanity checks:
  test <- c(min(gout_list$pair_len) == 3,
            max(gout_list$pair_len) == 3,
            length(gout_list[,1]) * 3 == table(control_ids$avail)["FALSE"] - prior)
  names(test) <- c()
  print(test)
  length(unique(unlist(paired_list)))
  #merge fam file
  new_fam_file <- merge(fam_file, genotyped[,c("f.eid", "f.21000.0.0", "f.31.0.0", "f.21003.0.0","f.21001.0.0", "goutaff", "control","control_no_diuretics", "goutwinnard","goutself", "gouthosp","gout_winnard_self","goutall","gout_self_ult" )], by.x = "IID", by.y = "f.eid", all.x=TRUE)
  #resort
  new_fam_file <- new_fam_file[order(new_fam_file$sort),]
  #check lengths are equal
  length(fam_file[,1]) == length(new_fam_file[,1])
  #### set up end
  
  
  
  ###
  ### gout winnard self
  ###
  
  #reset case/control as precaution
  new_fam_file$AFF <- -9
  new_fam_file$SEX <- -9
  
  #in ukbio males are coded as female = 0, males = 1
  new_fam_file[new_fam_file$f.31.0.0 == 1 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 1
  new_fam_file[new_fam_file$f.31.0.0 == 0 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 2
  
  # affstat for all gout
  new_fam_file[new_fam_file$gout_winnard_self == 1 & !is.na(new_fam_file$gout_winnard_self), "AFF"] <- 2
  c_id <- unlist(lapply(gout_list[gout_list$ids %in% subset(new_fam_file, new_fam_file$gout_winnard_self == 1)$IID, "ids"] , function (x) {paired_list[[x]]}))
  new_fam_file[new_fam_file$IID %in% c_id,"AFF"] <- 1
  
  # make sure the number of controls is 3x gout
  # seems that some people are said to have been genoytped but do not exist in the fam file
  table(new_fam_file$AFF)
  length(c_id)
  table(new_fam_file$gout_winnard_self == 1 & !is.na(new_fam_file$gout_winnard_self))
  gout_list[!gout_list$ids %in%  new_fam_file[new_fam_file$goutaff == 1 & !is.na(new_fam_file$goutaff),"IID"], ]
  c_id[!c_id %in% new_fam_file$IID]
  
  table(new_fam_file$AFF, new_fam_file$f.21000.0.0, exclude=NULL)
  # blank out non-white ethnicities
  new_fam_file[!(!is.na(new_fam_file$f.21000.0.0) & (new_fam_file$f.21000.0.0 == 1001 | new_fam_file$f.21000.0.0 == 1002 | new_fam_file$f.21000.0.0 == 1003)) , "AFF"] <- -9
  write.table(new_fam_file[,c("FID","IID","PID","MID","SEX","AFF")], file = paste0("/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/UkBio/diuretics/all/chrallimpv1.fam_goutall_diuretic",i), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = ' ')
  
  
  ###
  ### fam gout hospital
  ###
  
  #reset sex and affstat
  new_fam_file$SEX <- -9
  new_fam_file$AFF <- -9
  
  #in ukbio males are coded as female = 0, males = 1
  new_fam_file[new_fam_file$f.31.0.0 == 1 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 1
  new_fam_file[new_fam_file$f.31.0.0 == 0 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 2
  
  # affstat for hosp criteria
  new_fam_file[new_fam_file$gouthosp == 1 & !is.na(new_fam_file$gouthosp), "AFF"] <- 2
  #find gout ids that match criteria then find the paired controls that go with them
  c_id <- unlist(lapply(gout_list[gout_list$ids %in% subset(new_fam_file, new_fam_file$gouthosp == 1)$IID, "ids"] , function (x) {paired_list[[x]]}))
  new_fam_file[new_fam_file$IID %in% c_id,"AFF"] <- 1
  
  table(new_fam_file$AFF, new_fam_file$f.21000.0.0, exclude=NULL)
  # blank out non-white ethnicities
  new_fam_file[!(!is.na(new_fam_file$f.21000.0.0) & (new_fam_file$f.21000.0.0 == 1001 | new_fam_file$f.21000.0.0 == 1002 | new_fam_file$f.21000.0.0 == 1003)) , "AFF"] <- -9
  write.table(new_fam_file[,c("FID","IID","PID","MID","SEX","AFF")], file = paste0("/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/UkBio/diuretics/hosp/chrallimpv1.fam_gouthosp_diuretic",i), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = ' ')
  rm(c_id)
  
  ###
  ### fam gout winnard
  ###
  
  #reset sex and affstat
  new_fam_file$SEX <- -9
  new_fam_file$AFF <- -9
  
  #in ukbio males are coded as female = 0, males = 1
  new_fam_file[new_fam_file$f.31.0.0 == 1 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 1
  new_fam_file[new_fam_file$f.31.0.0 == 0 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 2
  
  # affstat for winnard criteria
  new_fam_file[new_fam_file$goutwinnard == 1 & !is.na(new_fam_file$goutwinnard), "AFF"] <- 2
  #find gout ids that match criteria then find the paired controls that go with them
  c_id <- unlist(lapply(gout_list[gout_list$ids %in% subset(new_fam_file, new_fam_file$goutwinnard == 1)$IID, "ids"] , function (x) {paired_list[[x]]}))
  new_fam_file[new_fam_file$IID %in% c_id,"AFF"] <- 1
  
  table(new_fam_file$AFF, new_fam_file$f.21000.0.0, exclude=NULL)
  # blank out non-white ethnicities
  new_fam_file[!(!is.na(new_fam_file$f.21000.0.0) & (new_fam_file$f.21000.0.0 == 1001 | new_fam_file$f.21000.0.0 == 1002 | new_fam_file$f.21000.0.0 == 1003)) , "AFF"] <- -9
  write.table(new_fam_file[,c("FID","IID","PID","MID","SEX","AFF")], file = paste0("/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/UkBio/diuretics/winnard/chrallimpv1.fam_goutwinnard_diuretic",i), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = ' ')
  rm(c_id)
  
  ###
  ### fam gout self
  ###
  
  #reset sex and affstat
  new_fam_file$SEX <- -9
  new_fam_file$AFF <- -9
  
  #in ukbio males are coded as female = 0, males = 1
  new_fam_file[new_fam_file$f.31.0.0 == 1 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 1
  new_fam_file[new_fam_file$f.31.0.0 == 0 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 2
  
  # affstat for self criteria
  new_fam_file[new_fam_file$goutself == 1 & !is.na(new_fam_file$goutself), "AFF"] <- 2
  #find gout ids that match criteria then find the paired controls that go with them
  c_id <- unlist(lapply(gout_list[gout_list$ids %in% subset(new_fam_file, new_fam_file$goutself == 1)$IID, "ids"] , function (x) {paired_list[[x]]}))
  new_fam_file[new_fam_file$IID %in% c_id,"AFF"] <- 1
  
  table(new_fam_file$AFF, exclude=NULL)
  # blank out non-white ethnicities
  new_fam_file[!(!is.na(new_fam_file$f.21000.0.0) & (new_fam_file$f.21000.0.0 == 1001 | new_fam_file$f.21000.0.0 == 1002 | new_fam_file$f.21000.0.0 == 1003)) , "AFF"] <- -9
  table(new_fam_file$AFF, new_fam_file$f.21000.0.0, exclude=NULL)
  write.table(new_fam_file[,c("FID","IID","PID","MID","SEX","AFF")], file = paste0("/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/UkBio/diuretics/self/chrallimpv1.fam_goutself_diuretic",i), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = ' ')
  
  ###
  ### self defined or ULT
  ###
  
  #reset sex and affstat
  new_fam_file$SEX <- -9
  new_fam_file$AFF <- -9
  
  #in ukbio males are coded as female = 0, males = 1
  new_fam_file[new_fam_file$f.31.0.0 == 1 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 1
  new_fam_file[new_fam_file$f.31.0.0 == 0 & !is.na(new_fam_file$f.31.0.0), "SEX"] <- 2
  
  # affstat for self criteria
  new_fam_file[new_fam_file$gout_self_ult == 1 & !is.na(new_fam_file$gout_self_ult), "AFF"] <- 2
  #find gout ids that match criteria then find the paired controls that go with them
  c_id <- unlist(lapply(gout_list[gout_list$ids %in% subset(new_fam_file, new_fam_file$gout_self_ult == 1)$IID, "ids"] , function (x) {paired_list[[x]]}))
  new_fam_file[new_fam_file$IID %in% c_id,"AFF"] <- 1
  
  table(new_fam_file$AFF, exclude=NULL)
  # blank out non-white ethnicities
  new_fam_file[!(!is.na(new_fam_file$f.21000.0.0) & (new_fam_file$f.21000.0.0 == 1001 | new_fam_file$f.21000.0.0 == 1002 | new_fam_file$f.21000.0.0 == 1003)) , "AFF"] <- -9
  table(new_fam_file$AFF, new_fam_file$f.21000.0.0, exclude=NULL)
  write.table(new_fam_file[,c("FID","IID","PID","MID","SEX","AFF")], file = paste0("/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/UkBio/diuretics/self_ult/chrallimpv1.fam_goutself_ult_diuretic",i), col.names=FALSE, row.names=FALSE, quote=FALSE, sep = ' ')
  
}

# only need to create this once and needs to include everyone
co_var_file <- new_fam_file[,c("FID","IID","f.21003.0.0", "f.21001.0.0")]
colnames(co_var_file) <- c("FID","IID","AGE","BMI")
write.table(co_var_file, file = "/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/UkBio/diuretics/chrallimpv1.covar", col.names=TRUE, row.names=FALSE, quote=FALSE, sep = ' ')
