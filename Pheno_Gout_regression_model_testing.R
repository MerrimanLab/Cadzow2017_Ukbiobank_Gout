load('~/ukbiobank_gout_bd_refined.RData')
all_fam <- read.table("/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/UkBio/controls_d/all/chrallimpv1.fam_goutall1", header=FALSE)
colnames(all_fam) <- c("FID","IID","PID","MID", "SEX","AFF")
genotyped  <-  bd_refined_gout[bd_refined_gout$f.eid %in% all_fam$IID,]
genotyped$waist_height_ratio <- genotyped$f.48.0.0 / genotyped$f.50.0.0
genotyped$waist_hip_ratio <- genotyped$f.48.0.0 / genotyped$f.49.0.0

gout_glm  <- glm(as.factor(genotyped$goutaff) ~ genotyped$f.31.0.0, family='binomial')
summary(gout_glm)

gout_glm2  <- glm(as.factor(genotyped$goutaff) ~ as.factor(genotyped$f.31.0.0) + genotyped$f.21003.0.0 + genotyped$f.21000.0.0 , family='binomial')
step(gout_glm2)


gout_glm_wheight  <- glm(as.factor(genotyped$goutaff) ~ as.factor(genotyped$f.31.0.0) + genotyped$f.21003.0.0 + genotyped$f.48.0.0 +(genotyped$waist_height_ratio), family='binomial')
summary(gout_glm_wheight)
step(gout_glm_wheight)

gout_glm_whip  <- glm(as.factor(genotyped$goutaff) ~ as.factor(genotyped$f.31.0.0) + genotyped$f.21003.0.0  +(genotyped$waist_hip_ratio), family='binomial')
summary(gout_glm_whip)
step(gout_glm_whip)

gout_glm_weight  <- glm(as.factor(genotyped$goutaff) ~ as.factor(genotyped$f.31.0.0) + genotyped$f.21003.0.0  +genotyped$f.48.0.0, family='binomial')
summary(gout_glm_weight)
step(gout_glm_weight)



gout_glm_w  <- glm(as.factor(genotyped$goutaff) ~ as.factor(genotyped$f.31.0.0) * genotyped$f.21003.0.0, family='binomial')
summary(gout_glm_w)
step(gout_glm_w)
