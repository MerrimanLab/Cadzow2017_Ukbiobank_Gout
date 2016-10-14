library(data.table)
cond <- c('all','winnard','hosp','self','self_ult','ult')

for(gc in cond){
  tmp <- list()
  for(i in 1:22){
    tmp[[i]] <- fread(paste0('/media/xsan/scratch/merrimanlab/murray/working_dir/UkBio/GWAS_all_controls/controls/',gc,'/adjusted/controls',gc,'_age_sex_chr',i,'.assoc.logistic.tsv'))
  }
  gwas <- rbindlist(tmp)
  
  write.table(gwas[TEST == 'ADD' & P < 1e-5]$SNP, file = paste0('/media/xsan/scratch/merrimanlab/murray/working_dir/UkBio/GWAS_all_controls/controls/Heritability/',gc,'/',gc,'_age_sex_nominally_sig_snps.txt'),
              quote = FALSE, row.names = FALSE, col.names=FALSE, sep ='\t')
}

