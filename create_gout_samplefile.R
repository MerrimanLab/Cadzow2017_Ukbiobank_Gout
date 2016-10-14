fam <- read.table('/media/xsan/scratch/merrimanlab/murray/working_dir/UkBio/GWAS_all_controls/controls/all/chrallimpv1.fam_all', stringsAsFactors = FALSE)
sam <- read.table('/media/xsan/archive/merrimanlab/central_datasets/ukbiobank/bgen/impv1.sample', stringsAsFactors = FALSE, header= TRUE)

covar <- read.table('/media/xsan/scratch/merrimanlab/murray/working_dir/UkBio/GWAS_all_controls/chrallimpv1.covar', header=TRUE, stringsAsFactors = FALSE)

# check order of ids is the same between sample, fam, and covar file
table(sam$ID_1[2:nrow(sam)] == fam$V1)
table(sam$ID_2[2:nrow(sam)] == fam$V2)
table(sam$ID_2[2:nrow(sam)] == covar$IID)
table(sam$ID_2[2:nrow(sam)] == covar$FID)

# AGE BMI WaistToHeightRatio Waist PCA1:15
sam2 <- cbind(sam,rbind(c('D',rep('C',18), 'B'), cbind(covar[,3:ncol(covar)], fam$V6)))
names(sam2)[23] <-"Gout"
names(sam2)
sam2$Gout[2:length(sam2$Gout)] <- as.numeric(sam2$Gout[2:length(sam2$Gout)]) -1
sam2[sam2$Gout == -10,]$Gout <- NA

write.table(sam2, file = '/media/xsan/archive/merrimanlab/central_datasets/ukbiobank/bgen/impv1_gout_and_covar.sample', quote=FALSE, col.names=TRUE, row.names=FALSE, sep =' ')

#create list to remove people without gout case/control status
rem <- sam2[is.na(sam2$Gout), 'ID_1']
write.table(rem, file = '/media/xsan/archive/merrimanlab/central_datasets/ukbiobank/bgen/samples_remove_gout_gwas.txt', col.names=FALSE, row.names=FALSE, quote=FALSE)
