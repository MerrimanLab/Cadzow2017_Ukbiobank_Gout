library(data.table)
library(qqman)
library(dplyr)
l <-  list()
for(i in 1:22){
  print(i)
  l[[i]] <- fread(paste0("controlswinnard_male_age_sex_chr",i,".assoc.logistic.tsv"),header = TRUE, sep = "\t", col.names = c("CHR","SNP","BP","A1","TEST","NMISS","OR","SE","L95","U95","STAT","P"))
  
}

dt <- rbindlist(l = l )
setkey(dt, CHR, BP)
head(dt)

dt  %>% filter(P > 0.05)  %>% tally()
dt2 <- dt[P < 0.05 & TEST == 'ADD']
manhattan(dt2[P< 0.05], main = 'age adjusted')