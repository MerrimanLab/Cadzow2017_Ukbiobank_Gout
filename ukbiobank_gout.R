# Murray Cadzow
# University of Otago
# 19 Jan 2016

# Script to select gout cases from the UKBioBank data
# criteria for classifying as gout as is follows (based on Winnard et al 2012)
# 
# Either Primary or secondary hospitalisation ICD-10 == M10
# or
# On allopurinol/probenecid/colchicine but exlcude if diagnosed with leukemia or lymphoma (ICD-10-AM C80-C96)
# 
# self-reported and on ULT?


# basic flow of script follows:
#   read in ukbiobank data
#   read in codings and drug info
#   work out rows ids for people
#     -based on drugs
#     -hospital records for various conditions eg KD or cancers
#     -gout affection based on various criteria
#   then
#     -work out people to exclude as controls
#   match cases with controls based on age/sex/eth
#   for each gout criteria
#     -randomly select 3 matched controls for each case
#     -write out fam file

# where the phenotype data lives
##Murray's Version
data_dir <- "/run/user/1000/gvfs/smb-share:server=biocldap,share=staff_groups/merrimanlab/Documents/Murray/ukbiobank_util/data/"
##Tanya's Version
#data_dir <- "/Volumes/staff_groups/merrimanlab/Merriman Documents/Murray/ukbiobank_util/data/"
# where the genotype data lives
##Murray's Version
scratch_dir <- "/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/UkBio/"
##Tanya's Version
#scratch_dir <- "/Volumes/scratch/merrimanlab/murray/working_dir/UkBio/"

# load reduced dataset (missing many food columns from right end of the original dataframe)
load(paste0(data_dir,"bd_refined.RData"))

# make sure withdrawn people are removed
withdrawn_ids <- read.csv(paste0(data_dir,'w1261_20160406.csv'), head=FALSE)
bd_refined <- bd_refined[!bd_refined$f.eid %in% withdrawn_ids[1,],]

bd_refined$gouthosp <- NA # meets hospital criteria for gout (1/NA)
bd_refined$goutwinnard <- NA # meets winnard criteria for gout (1/NA)
bd_refined$goutself <- NA # meets self report criteria for gout (1/NA)
bd_refined$gout_winnard_self <- NA # meets self report and winnard criteria for gout (1/NA)
bd_refined$control <- NA # meets conditions to be included as a control for gout (1/NA)
bd_refined$goutaff <- NA # gout affection, sum of all gout and control criteria (1 = case/0 = control/NA)
bd_refined$renaldisease <- NA # has renal disease (1/NA)
bd_refined$allopurinol <- NA # on allopurinol (1/NA)
bd_refined$colchicine <- NA # on colchicine (1/NA)
bd_refined$febuxostat <- NA # on febuxostat (1/NA)
bd_refined$sulphinpyrazone <- NA # on sulphinpyrazone (1/NA)
bd_refined$probenecid <- NA # on probenecid (1/NA)
bd_refined$cortico_prednisone <- NA # on cortico prenisone (1/NA)
bd_refined$nsaids <- NA # on nsaids (1/NA)
bd_refined$diuretics <- NA # on diuretics (1/NA)
bd_refined$gout_exclude <- NA # excluded from gout case control study for any gout reason (eg renal disease) (1/NA)
# gout_exclude_drug and gout_exclude_hosp are for establishing inclusion by ULT only
bd_refined$gout_exclude_drug <- NA # excluded because of gout_exclude_hosp and taking ULT (1/NA)
bd_refined$gout_exclude_hosp <- NA # exclude because have blood cancer (1/NA) 
bd_refined$control_exclude_drug <- NA # excluded from controls because on probenecid, cortico prenisone or nsaids (1/NA)
bd_refined$control_exclude <- NA # controls excluded because of renaldisease or control_exclude_drug (does not include people classified as gout) (1/NA)
bd_refined$qc_exclude <- NA # any reason to exclude a person from the analysis eg renal disease, or genotyping qc fail
#columns to do with non-cancer illnesses and medication
illness <- bd_refined[,grep('f.eid|f.2000[2389]', colnames(bd_refined))]

#columns for non-cancer illness code, self-reported: uses data-coding 6
ill_cols <- grep("f.20002", names(illness))
datacoding6 <- read.delim(paste0(data_dir,"codings/coding6.tsv"), stringsAsFactors = FALSE, sep="\t")


gout_code <- datacoding6[grep('gout', datacoding6[,"meaning"]), "coding"]

#find all the rows with gout from illness - currently finding 7337 which is less than ukbiobank reports (7387)
rows_w_self_gout <-unique(sort(unlist(sapply(ill_cols,function(x){ which(gout_code == illness[,x])}))))


#columns for medication/treatment: uses data-coding 4
med_cols <- grep("f.20003", names(illness))
datacoding4 <- read.delim(paste0(data_dir,"codings/coding4.tsv"), stringsAsFactors = FALSE, sep="\t")

# Nicola annotated datacoding4 file
nd_drugs <- read.csv(paste0(data_dir,"UKBio drugs ND formatted.csv"), header=TRUE, stringsAsFactors = FALSE)

# Allopurinol/probenecid/colchicine codes
apc_codes <- subset(nd_drugs, nd_drugs$Allopurinol == 'y' | nd_drugs$Colchicine == 'y' | nd_drugs$Febuxostat == 'y' | nd_drugs$Sulphinpyrazone == 'y')[,"coding"]

afs_only_codes <- subset(nd_drugs, nd_drugs$Allopurinol == 'y' | nd_drugs$Febuxostat == 'y' | nd_drugs$Sulphinpyrazone == 'y')[,"coding"]

colchicine_codes <- subset(nd_drugs,  nd_drugs$Colchicine == 'y')[,"coding"]
rows_w_colchicine <- unique(sort(unlist(sapply(med_cols,function(x){ which(illness[,x] %in% colchicine_codes)}))))
bd_refined[rows_w_colchicine,"colchicine"] <-1

allopurinol_codes <- subset(nd_drugs,  nd_drugs$Allopurinol == 'y')[,"coding"]
rows_w_allopurinol <- unique(sort(unlist(sapply(med_cols,function(x){ which(illness[,x] %in% allopurinol_codes)}))))
bd_refined[rows_w_allopurinol,"allopurinol"] <-1

febuxostat_codes <- subset(nd_drugs,  nd_drugs$Febuxostat == 'y')[,"coding"]
rows_w_febuxostat <- unique(sort(unlist(sapply(med_cols,function(x){ which(illness[,x] %in% febuxostat_codes)}))))
bd_refined[rows_w_febuxostat,"febuxostat"] <-1

sulphinpyrazone_codes <- subset(nd_drugs,  nd_drugs$Sulphinpyrazone == 'y')[,"coding"]
rows_w_sulphinpyrazone <- unique(sort(unlist(sapply(med_cols,function(x){ which(illness[,x] %in% sulphinpyrazone_codes)}))))
bd_refined[rows_w_sulphinpyrazone,"sulphinpyrazone"] <-1


rows_w_drug_crit <-unique(sort(unlist(sapply(med_cols,function(x){ which(illness[,x] %in% apc_codes)}))))
rows_w_afs <- unique(sort(unlist(sapply(med_cols,function(x){ which(illness[,x] %in% afs_only_codes)}))))

diuretics_codes <- subset(nd_drugs, nd_drugs$Diuretic == 'y' )[,"coding"]
rows_w_diuretics <- unique(sort(unlist(sapply(med_cols,function(x){ which(illness[,x] %in% diuretics_codes)}))))
bd_refined[rows_w_diuretics, "diuretics"] <- 1
#####
#### HOSPITAL BASED CRITERIA ####
#####

#columns to do with hospital admission
hospital <- bd_refined[,grep('f.eid|f.4120[24]', colnames(bd_refined))]

#ICD-10 uses data-coding 19
datacoding19 <-  read.delim(paste0(data_dir,"codings/coding19.tsv"), 
                            stringsAsFactors = FALSE, sep="\t")

icd_gout_codes <- datacoding19[grep('M10', datacoding19[,"coding"]), "coding"]

primary_cols <- grep("f.41202", names(hospital))
secondary_cols <- grep("f.41204", names(hospital))

#find all people with primary diagnosis of gout ICD-10 M10 (ukbiobank reports 283)
rows_w_hosp_prim_gout <-sort(unlist(sapply(primary_cols,function(x){ grep("M10", hospital[,x])})))
#find all people with secondary diagnosis of gout ICD-10 M10 (ukbiobank reports 1320)
rows_w_hosp_sec_gout <-sort(unlist(sapply(secondary_cols,function(x){ grep("M10", hospital[,x])})))

hosp_gout_rows <- unique(c(rows_w_hosp_prim_gout, rows_w_hosp_sec_gout))


# hospital criteria for exclusion
# have lymphoma or leukemia and on allopurinol/probenecid/colchicine
#find all people with primary diagnosis of ICD-10  C80 -> C96  (ukbiobank reports 280 + 3033)
prim_c81_c96 <- list()
for( code in  paste0("C", 81:96)){
  prim_c81_c96[[code]] <-sort(unlist(sapply(primary_cols,function(x){ grep(code, hospital[,x])})))
}
rows_w_hosp_prim_exclude <-unique(unlist(prim_c81_c96))
#find all people with secondary diagnosis of C81 -> C96 ICD-10  (ukbiobank reports 358 + 1109)
sec_c81_c96 <- list()
for( code in  paste0("C", 81:96)){
  sec_c81_c96[[code]] <-sort(unlist(sapply(secondary_cols,function(x){ grep(code, hospital[,x])})))
}
rows_w_hosp_sec_exclude <- unique(unlist(sec_c81_c96))
bd_refined[c(rows_w_hosp_prim_exclude,rows_w_hosp_sec_exclude), "gout_exclude_hosp"] <-1
# rows that have drug criteria but need to be excluded because the med might be for that disease not gout
exclude_drug <- intersect(c(rows_w_hosp_prim_exclude,rows_w_hosp_sec_exclude), rows_w_drug_crit)
bd_refined[exclude_drug, "gout_exclude_drug"] <- 1
# final rows for people to include based on drugs
rows_drug_crit_include <- rows_w_drug_crit[ !(rows_w_drug_crit %in% exclude_drug) ]
rows_w_afs_include <- rows_w_afs[ !(rows_w_afs %in% exclude_drug) ]
# remove renal disease as defined in: 
# Quan, Hude, et al. "Coding algorithms for defining comorbidities in ICD-9-CM and ICD-10 administrative data." Medical care (2005): 1130-1139.
# AND doi:10.1186/1471-2288-11-83
rd_icd_10_codes <- c( "I12", "I13", paste0("N0", 0:5), "N07", "N11", "N14", paste0("N",17:19), "Q61","N25.0", "Z49", "Z94.0", "Z99.2")
rd_list <- list()
for( code in  rd_icd_10_codes){
  rd_list[[code]] <-sort(unlist(sapply(primary_cols,function(x){ grep(code, hospital[,x])})))
}
rows_rd_prim <- unique(unlist(rd_list))

rd_list <- list()
for( code in  rd_icd_10_codes){
  rd_list[[code]] <-sort(unlist(sapply(secondary_cols,function(x){ grep(code, hospital[,x])})))
}
rows_rd_sec <- unique(unlist(rd_list))

# use this to remove renal disease from cases/controls prior to making subsets of bd_refined
rows_rd <- union(rows_rd_prim, rows_rd_sec)
bd_refined[rows_rd,"renaldisease"] <- 1

rows_has_chr4 <- which(!is.na(bd_refined[,"f.22104.0.0"]))

#exclude identified gout people with renal disease
hospital <- unique(c(rows_w_hosp_prim_gout, rows_w_hosp_sec_gout) )
hospital <- hospital[!(hospital %in% rows_rd)]
winnard  <- unique(c(rows_drug_crit_include, rows_w_hosp_prim_gout, rows_w_hosp_sec_gout))
winnard <- winnard[!(winnard %in% rows_rd)]
self <- unique(rows_w_self_gout)
self <- self[!(self %in% rows_rd)]
self_winnard <- unique(c(rows_drug_crit_include, rows_w_hosp_prim_gout, rows_w_hosp_sec_gout, rows_w_self_gout))
self_winnard <- self_winnard[!(self_winnard %in% rows_rd)]
# ult = allopurinol, febuxostat and sulfinpyrazone
self_ult <- unique(c(rows_w_afs_include, rows_w_self_gout))
self_ult <- self_ult[!(self_ult %in% rows_rd)]
ult <- unique(rows_w_afs_include)
ult <- ult[!(ult %in% rows_rd)]
drug <- unique(rows_drug_crit_include)
drug <- drug[!(drug %in% rows_rd)]

bd_refined[hospital, 'gouthosp'] <- 1
bd_refined[winnard, 'goutwinnard'] <- 1
bd_refined[self, 'goutself'] <- 1
bd_refined[self_winnard, 'gout_winnard_self'] <- 1
bd_refined[self_ult, 'gout_self_ult'] <- 1
bd_refined[ult, 'goutult'] <- 1
bd_refined[drug, 'goutdrug'] <-1
# make sure to add any further gout selecting conditions into this line:
bd_refined[unique(c(hospital, winnard,self, self_winnard, self_ult, ult, drug )),'goutall'] <- 1

#### NUMBERS ####
#total
length(bd_refined[,1])

#drug criteria
length(subset(bd_refined, bd_refined$goutdrug == 1)[,1])

# self reported gout
length(subset(bd_refined, bd_refined$goutself == 1)[,1])

# self reported or drug
length(subset(bd_refined, bd_refined$goutself == 1 | bd_refined$goutdrug)[,1])

# hospital total
length(subset(bd_refined, bd_refined$gouthosp == 1)[,1])

# self or hospital
length(subset(bd_refined, bd_refined$goutself == 1 | bd_refined$gouthosp == 1)[,1])

# hospital or drug (winnard)
length(subset(bd_refined, bd_refined$goutdrug == 1 | bd_refined$gouthosp == 1)[,1])

# hospital or drug or self
length(subset(bd_refined, bd_refined$goutself == 1 | bd_refined$gouthosp == 1 | bd_refined$goutdrug)[,1])

# ult
length(subset(bd_refined, bd_refined$goutult == 1)[,1])

#self or ult
length(subset(bd_refined, bd_refined$goutself == 1 | bd_refined$goutult == 1)[,1])

library(dplyr)
## NUMBERS with chr4 genotypes
#total
length(rows_has_chr4)
# drug
length(subset(bd_refined[rows_has_chr4,], bd_refined[rows_has_chr4,'goutdrug'] == 1)[,1])
# self
length(subset(bd_refined[rows_has_chr4,], bd_refined[rows_has_chr4,'goutself'] == 1)[,1])
#self drug
length(subset(bd_refined[rows_has_chr4,], bd_refined[rows_has_chr4,'goutdrug'] == 1 | bd_refined[rows_has_chr4,'goutself'] == 1)[,1])
#hospital
length(subset(bd_refined[rows_has_chr4,], bd_refined[rows_has_chr4,'gouthosp'] == 1)[,1])
#hospital and self
length(subset(bd_refined[rows_has_chr4,], bd_refined[rows_has_chr4,'gouthosp'] == 1 | bd_refined[rows_has_chr4,'goutself'] == 1)[,1])
#drug or hosp
length(subset(bd_refined[rows_has_chr4,], bd_refined[rows_has_chr4,'goutdrug'] == 1 | bd_refined[rows_has_chr4,'gouthosp'] == 1)[,1])
#drug or hosp or self
length(subset(bd_refined[rows_has_chr4,], bd_refined[rows_has_chr4,'goutdrug'] == 1 | bd_refined[rows_has_chr4,'goutself'] == 1 | bd_refined[rows_has_chr4,'gouthosp'] == 1)[,1])
#ult
length(subset(bd_refined[rows_has_chr4,], bd_refined[rows_has_chr4,'goutult'] == 1)[,1])
#ult or self
length(subset(bd_refined[rows_has_chr4,], bd_refined[rows_has_chr4,'goutult'] == 1 | bd_refined[rows_has_chr4,'goutself'] == 1)[,1])




### Controls ####

# exclude if have gout or on ULT, allopurinol, colchicine, prednisone, or NSAIDS
# 2 subsets, without diuretics and with diuretics

# self or hospital reported or based on drug use
rows_any_gout_inclusion <- which(bd_refined$goutall == 1)
# ULT and NSAIDS
ult_nsaids_codes <- subset(nd_drugs, nd_drugs$Allopurinol == 'y' | nd_drugs$Colchicine == 'y' | nd_drugs$Febuxostat == 'y' | nd_drugs$Sulphinpyrazone == 'y' | nd_drugs$Probenecid == 'y' | nd_drugs$Corticosteroids.prednisone == 'y' | nd_drugs$NSAIDs == 'y')[,"coding"]
probenecid_codes <- subset(nd_drugs,  nd_drugs$Probenecid == 'y')[,"coding"]
bd_refined[unique(sort(unlist(sapply(med_cols,function(x){ which(illness[,x] %in% probenecid_codes)})))), "probenecid"] <- 1
cortico_prednisone_codes <- subset(nd_drugs, nd_drugs$Corticosteroids.prednisone == 'y' )[,"coding"]
bd_refined[unique(sort(unlist(sapply(med_cols,function(x){ which(illness[,x] %in% cortico_prednisone_codes)})))), "cortico_prednisone"] <- 1
nsaids_codes <- subset(nd_drugs, nd_drugs$NSAIDs == 'y')[,"coding"]
bd_refined[unique(sort(unlist(sapply(med_cols,function(x){ which(illness[,x] %in% nsaids_codes)})))), "nsaids"] <- 1

diuretics_codes <- subset(nd_drugs, nd_drugs$Diuretic == 'y' )[,"coding"]
rows_w_diuretics <- unique(sort(unlist(sapply(med_cols,function(x){ which(illness[,x] %in% diuretics_codes)}))))

#remove all people that may match some sort of gout criteria or on ULT or NSAIDS or prednisone
rows_exclude_controls_bc_drugs <-unique(sort(unlist(sapply(med_cols,function(x){ which(illness[,x] %in% ult_nsaids_codes)}))))

bd_refined[rows_exclude_controls_bc_drugs,"control_exclude_drug"] <- 1
rows_exclude_controls_bc_gout <- rows_any_gout_inclusion


#exclude gouts, drugs and renal disease
rows_exclude_controls <- unique(sort(c(rows_exclude_controls_bc_gout, rows_exclude_controls_bc_drugs, rows_rd)))

#bd_refined[!(rownames(bd_refined) %in% rows_exclude_controls), "control"] <- 1
bd_refined[!( (bd_refined$goutall == 1 & !is.na(bd_refined$goutall)) |
                (bd_refined$control_exclude_drug == 1 &!is.na(bd_refined$control_exclude_drug)) |
                (bd_refined$renaldisease == 1 & ! is.na(bd_refined$renaldisease)) ) ,'control'] <- 1
#bd_refined[!(rownames(bd_refined) %in% c(rows_exclude_controls, rows_w_diuretics)), 'control_no_diuretics'] <- 1 # removed - should now use 'control' inconjuction with 'diuretics' column


# fill in final exclusion column
bd_refined[which(bd_refined$control_exclude_drug == 1 | bd_refined$renaldisease == 1), "control_exclude"] <- 1
bd_refined[which(bd_refined$renal_disease ==1),"gout_exclude"] <-1


#bd_refined2 = bd_refined #create 'restore point' before enacting exclusions

#remove related controls (gout cases that are related to other cases still left in)
#find all related pair ids
related_pairs  <- unlist(bd_refined[,grep('f.22011', colnames(bd_refined))])
# create a unique list
related_pairs  <- unique(related_pairs[!is.na(related_pairs)])
#for each control, see if they have a related pair id and if they do mark the related controls column '1'
bd_refined[bd_refined$control ==1 & !is.na(bd_refined$control) &is.na(bd_refined$control_exclude),"control_exclude"] <- apply(bd_refined[bd_refined$control == 1 & !is.na(bd_refined$control) & is.na(bd_refined$control_exclude),grep('f.22011', colnames(bd_refined))], 1, function(x){return(ifelse(sum(x %in% related_pairs) > 0,1,NA))})


# remove people who have reported and genetic gender disagreement
bd_refined[bd_refined$f.31.0.0 != bd_refined$f.22001.0.0 & !is.na(bd_refined$f.22001.0.0) & !is.na(bd_refined$f.31.0.0),"qc_exclude"] <- 1

# remove people who failed genotyping QC
bd_refined[bd_refined$f.22010.0.0 == 1 & !is.na(bd_refined$f.22010.0.0),'qc_exclude'] <- 1


# 'remove'/NA all controls with a 1 in the related_controls column
bd_refined[!is.na(bd_refined$control_exclude), "control" ] <- NA
bd_refined[!is.na(bd_refined$qc_exclude), "control" ] <- NA

# 'remove'/NA all gout cases that have a reason to exclude
bd_refined[!is.na(bd_refined$qc_exclude), "goutall" ] <- NA
bd_refined[!is.na(bd_refined$qc_exclude), "goutwinnard" ] <- NA
bd_refined[!is.na(bd_refined$qc_exclude), "goutself" ] <- NA
bd_refined[!is.na(bd_refined$qc_exclude), "gout_self_ult" ] <- NA
bd_refined[!is.na(bd_refined$qc_exclude), "gouthosp" ] <- NA
bd_refined[!is.na(bd_refined$qc_exclude), 'gout_winnard_self'] <- NA
bd_refined[!is.na(bd_refined$qc_exclude), 'goutult'] <- NA
bd_refined[!is.na(bd_refined$qc_exclude), 'goutdrug'] <-NA



# fill in 'most inclusive rows' for both gout and control
bd_refined[bd_refined$control == 1 & !is.na(bd_refined$control), "goutaff"] <- 0
bd_refined[bd_refined$goutall == 1 & !is.na(bd_refined$goutall), "goutaff"] <-1

#sanity check: Should all be NAs
table(bd_refined$control == 1 & bd_refined$goutall == 1 , exclude = NULL)
table(bd_refined$control == 1 & bd_refined$gouthosp == 1 , exclude = NULL)
table(bd_refined$control == 1 & bd_refined$goutwinnard == 1 , exclude = NULL)
table(bd_refined$control == 1 & bd_refined$goutself == 1 , exclude = NULL)
table(bd_refined$control == 1 & bd_refined$gout_winnard_self == 1 , exclude = NULL)



# test a few conditions that shouldn't be possible if correct filtering was applied:
table(bd_refined$goutaff == 1 & bd_refined$control == 1)
table(!is.na(bd_refined$goutaff) & bd_refined$renaldisease == 1 )
table(bd_refined$goutaff == 0 & bd_refined$gouthosp == 1)
table(bd_refined$goutaff == 0 & bd_refined$goutwinnard == 1)
table(bd_refined$goutaff == 0 & bd_refined$goutself == 1)
table(bd_refined$goutaff == 0 & bd_refined$gout_winnard_self == 1)
table(bd_refined$goutaff == 0 &(bd_refined$allopurinol == 1 | bd_refined$febuxostat ==1 | bd_refined$sulphinpyrazone == 1 | bd_refined$colchicine == 1))
table(bd_refined$goutaff == 0 & (bd_refined$nsaids == 1 | bd_refined$cortico_prednisone == 1 | bd_refined$probenecid == 1))
table(bd_refined$goutaff == 1 & bd_refined$gout_exclude == 1)
table(bd_refined$goutaff == 1 & bd_refined$gout_exclude_drug == 1) # the true results should hosp or self report classification
table(bd_refined$goutaff == 1 & bd_refined$gout_exclude_hosp == 1) # the true results should hosp or self report classification




#columns used that are relevant to gout
gout_cols <- unique(colnames(bd_refined)[which(colnames(bd_refined) %in% c(colnames(illness),colnames(hospital),"gouthosp", "goutwinnard", "goutself","gout_winnard_self","control","goutaff","f.20003.0.0","f.eid","f.31.0.0","f.21000.0.0","f.21001.0.0","f.21003.0.0",'f.49.0.0','f.48.0.0','f.50.0.0',"gout_self_ult", "goutult", "goutdrug","gout_self_ult","goutall","renaldisease","allopurinol","colchicine","febuxostat","sulphinpyrazone","diuretics","cortico_prednisone","nsaids","probenecid","gout_exclude","gout_exclude_drug","gout_exclude_hosp", "control_exclude", "control_exclude_drug", "qc_exclude"))])
bd_refined_gout <- bd_refined[,gout_cols]
save(bd_refined_gout, file=paste0('ukbiobank_gout_bd_refined',format(Sys.Date()),'.RData'))
# plot metrics
#boxplot( bd_refined$f.21003.0.0 ~ as.factor(bd_refined$goutaff), main = "age")
#boxplot( bd_refined$f.21001.0.0 ~ as.factor(bd_refined$goutaff), main = "bmi")
#plot( as.factor(bd_refined$f.31.0.0) ~ as.factor(bd_refined$goutaff), main = "sex")

# match all possible gout cases with controls
# method to match controls to cases

### USE FOR WHEN ALL GENOTYPES AVAILABLE
# control_age_sex_match <- list()
# control_age_sex_eth_match <- list()
# for (case_row in which(bd_refined$goutall == 1 )){
#   control_age_sex_match[[case_row]] <- which(bd_refined$f.21003.0.0 == bd_refined[case_row,"f.21003.0.0"] & bd_refined$f.31.0.0 == bd_refined[case_row,"f.31.0.0"])
#   control_age_sex_eth_match[[case_row]] <- which(bd_refined$f.21003.0.0 == bd_refined[case_row,"f.21003.0.0"] & bd_refined$f.31.0.0 == bd_refined[case_row,"f.31.0.0"] & bd_refined$f.21000.0.0 == bd_refined[case_row, "f.21000.0.0"])
#   
# } 
###


## match people that currently have genotypes
#read in fam file
fam_file <- read.table(paste0(scratch_dir,"chr1impv1.fam"), header=FALSE, stringsAsFactors = FALSE)
colnames(fam_file) <- c("FID","IID","PID","MID","SEX","AFF")
fam_file$sort <- as.numeric(rownames(fam_file))


#changed to use fam file 21 April 2016 - look at version history on git for old method
genotyped <- bd_refined[bd_refined$f.eid %in% fam_file$IID, ]
save(genotyped, file = paste0("ukbiobank_genotyped",format(Sys.Date()),".RData"))



####################################################################################################
################# LEGACY CODE FOR CASE-CONTROL MATCHING BEYOND THIS POINT ##########################
####################################################################################################
####################################################################################################
#control_age_sex_match_chr4 <- list()
#control_age_sex_eth_match_chr4 <- list()
#for (case in genotyped[which(genotyped$goutaff == 1),"f.eid"]){
#  control_age_sex_match_chr4[[case]] <- genotyped[which(genotyped$f.21003.0.0 == genotyped[which(case == genotyped$f.eid),"f.21003.0.0"] & genotyped$f.31.0.0 == genotyped[which(genotyped$f.eid == case) & genotyped$goutaff == 0,"f.31.0.0"]), "f.eid"]

#  control_age_sex_eth_match_chr4[[case]] <- genotyped[which(genotyped$f.21003.0.0 == genotyped[which(genotyped$f.eid == case),"f.21003.0.0"] & genotyped$f.31.0.0 == genotyped[which(genotyped$f.eid == case),"f.31.0.0"] & genotyped$f.21000.0.0 == genotyped[which(genotyped$f.eid == case) & genotyped$goutaff == 0, "f.21000.0.0"]), "f.eid"]
#} 

#make id lists for controls possibly exposed to diuretics
# control_ids <- data.frame(ids = genotyped[genotyped$goutaff == 0 & !is.na(genotyped$goutaff), "f.eid"])
# control_ids$avail <- FALSE
# control_ids[control_ids$ids %in% fam_file$IID, 'avail'] <-TRUE
# 
# gout_list <- data.frame(ids = genotyped[genotyped$goutaff == 1 & !is.na(genotyped$goutaff), "f.eid"])
# gout_list$len <- NA
# for (case in genotyped[which(genotyped$goutaff == 1 & !is.na(genotyped$goutaff)),"f.eid"]){
#   gout_list[gout_list$ids == case, "len"] <- length(control_age_sex_eth_match_chr4[[case]])
# }
# 
# #order by number of matches so least matches is paired up first
# gout_list <- gout_list[order(gout_list$len),]
# 
# 
# 
# 
# pair_controls <- function(id){
#   samples <- control_age_sex_eth_match_chr4[[id]]
#   samples <- control_ids[which(control_ids$ids %in% samples), ]
#   samples <- samples[samples$avail == TRUE, "ids"]
#   if(length(samples) >= 3){
#     s <- sample(samples, replace = FALSE, size = 3)
#   } else {
#     s <- samples
#   }
#   control_ids[which(control_ids$ids %in% s),"avail"] <<- FALSE
#   return(s)
#   
# }


### pca experimentation ###

# pca
# reduced <-bd_refined[rownames(bd_refined) %in% rows_has_chr4 & !is.na(bd_refined$goutaff),c("f.eid","f.21003.0.0", "f.31.0.0", "f.21000.0.0", "goutaff")]
# dim(reduced)
# 
# reduced$f.21000.0.0 <- as.integer(as.factor(as.character(reduced$f.21000.0.0)))
# reduced$agescale <- scale(reduced$f.21003.0.0)
# reduced$ethscale <- scale(reduced$f.31.0.0)
# pca <- prcomp(reduced[,c("agescale","ethscale", "f.31.0.0", "goutaff")])
# head(pca$x[,1:2])
# 
# library(ggplot2)
# ggplot(as.data.frame(cbind(pca$x[,1:2] )), aes(x=PC1, y=PC2)) + geom_point(colour = reduced$goutaff+1)
# summary(pca)
# 
# 
# 
# d <- as.data.frame(cbind(pca$x[,1:2], reduced))
# km <- kmeans(pca$x[,1:2], 44)
# d$group <- km$cluster %% 4
# d$cluster <- km$cluster
# ggplot(d , aes(x=PC1, y=PC2, colour = as.factor(cluster))) + geom_point() 
# 
# e <- list()
# for(id in d$f.eid[1:1000]){
#   e[[id]] <-d[which(d[d$f.eid == id,]$PC2 <= d$PC2 +0.1 & d[d$f.eid == id,]$PC2 >= d$PC2 -0.1 & d$f.31.0.0 == d[d$f.eid == id,]$f.31.0.0 ),"f.eid"]
# }
# 
# plot(NULL, xlim=c(min(d$f.eid), max(d$f.eid)), ylim=c(min(d$f.eid), max(d$f.eid)))
# for(id in d$f.eid[1:1000]) {
#   for(i in e[[id]]){
#     points(id, i, pch='.') }
# }
# distmat <- dist(d[, c('PC1', 'PC2')])
