# Kottgen loci table

#grep -w 'rs1471633\|rs1260326\|rs12498742\|rs2231142\|rs675209\|rs1165151\|rs1171614\|rs2078267\|rs478607\|rs3741414\|rs11264341\|rs17050272\|rs2307394\|rs6770152\|rs17632159\|rs729761\|rs1178977\|rs10480300\|rs17786744\|rs2941484\|rs10821905\|rs642803\|rs653178\|rs1394125\|rs6598541\|rs7193778\|rs7188445\|rs7224610\|rs2079742\|rs164009\|CHR' *logistic | grep -v 'BMI\|SEX' | sed 's/ \+/\t/g' > kottgen_results.txt


library(dplyr)
create_kottgen_table1 <- function(kottgen_loci = NULL){
  kot_t1 <- read.csv('kottgen_table1.csv', header=TRUE, stringsAsFactors = FALSE)
  kot_t1$sort <- as.numeric(row.names(kot_t1))
  a <- read.table(file = kottgen_loci, stringsAsFactors = FALSE)
  colnames(a) <- c('file',"CHR","SNP","BP","A1","TEST","NMISS","OR","SE","L95","U95","STAT","P")
  # remove 'CHR' from the chromosome column
  a <- a[grep("CHR", a$CHR,invert = TRUE),]
  # change the appropriate columns to be numeric
  a[,c("CHR","NMISS","OR","SE","L95","U95","STAT","P")] <- apply(a[,c("CHR","NMISS","OR","SE","L95","U95","STAT","P")], 2, function(x) {as.numeric(x)})
  
  #kott_table1 <-  as.data.frame(a  %>% group_by(SNP) %>% summarise(mean_OR = exp(mean(log(OR))), 
  #                                                                 max_SE = max(SE), 
  #                                                                 lower95 = mean_OR - 1.96 * max_SE,  
  #                                                                 upper95 = mean_OR + 1.96 * max_SE,  
  #                                                                 pcalc = 2* pnorm(-abs(log(mean_OR)/max_SE)), 
  #                                                                 n = length(OR)))
  
  # kott_tab1_snps <- data.frame(snps=c("rs1471633","rs1260326","rs12498742","rs2231142","rs675209","rs1165151","rs1171614","rs2078267","rs478607",
  #                     "rs3741414","rs11264341","rs17050272","rs2307394","rs6770152","rs17632159","rs729761","rs1178977","rs10480300",
  #                     "rs17786744","rs2941484","rs10821905","rs642803","rs653178","rs1394125","rs6598541","rs7193778","rs7188445",
  #                     "rs7224610","rs2079742","rs164009"))
  #kott_tab1_snps$sort <- as.numeric(row.names(kott_tab1_snps))
  kott_table1 <- a
  kott_table1 <- merge(kott_table1, kot_t1[,colnames(kot_t1) %in% c('SNP','CHR','BP', 'closest_Gene','grail_gene','sort')], by = "SNP")
  kott_table1 <- kott_table1[,!colnames(kott_table1) %in% c('CHR.y','BP.y')]
  colnames(kott_table1) <- c("SNP","file","CHR","BP","A1","TEST","NMISS","OR","SE", "L95","U95","STAT","P","grail_gene","sort")
  kott_table1 <- kott_table1[order(kott_table1$sort),]

  kott_table1$OR <- round(kott_table1$OR, digits = 3)
  kott_table1$SE <- round(kott_table1$SE, digits = 3)
  kott_table1$lower95 <- round(kott_table1$L95, digits = 3)
  kott_table1$upper95 <- round(kott_table1$U95, digits = 3)
  kott_table1$P <- signif(kott_table1$P,digits = 4)
  
  return(list(sum_table = kott_table1[,c("SNP","CHR", "BP", "grail_gene", "OR", "SE", "lower95", "upper95", "P" )], dat = a))
}


#kottgen_loci ="/run/user/1000/gvfs/smb-share:server=biocldap,share=scratch/merrimanlab/murray/working_dir/UkBio/GWAS_all_controls/controls/all/unadjusted/kottgen_results.txt"
