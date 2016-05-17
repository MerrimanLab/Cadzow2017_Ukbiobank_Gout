library(dplyr)

read_fams <- function(fam_path = NULL){
  fam = data.frame()
  it_summary = list()
  drug_summary = list()
  case_ids = c()
  control_ids = c()
  
  fam <- read.table(paste0(fam_path), header=FALSE)
  colnames(fam) <- c("FID","IID","PID","MID", "SEX","AFF")
  genotyped$faminfo <- NA
  genotyped[genotyped$f.eid %in% fam$IID & 
              genotyped$f.eid %in% fam[which(fam$AFF==1 | fam$AFF == 2),"IID"], "faminfo"] <-1
  genotyped$fam <- NA
  genotyped[genotyped$f.eid %in% fam$IID,"fam"] <- 1
  
  case_ids <- genotyped  %>% select(faminfo, goutaff, f.eid)  %>% filter(goutaff==1 & !is.na(goutaff) & f.eid %in% fam$IID & faminfo == 1)  %>%  select(f.eid)
  control_ids <- genotyped  %>% select(faminfo, goutaff, f.eid)  %>% filter(goutaff==0 & !is.na(goutaff) & f.eid %in% fam$IID & faminfo == 1)  %>%  select(f.eid)
  
  drug_summary <- fam_drugs(cases = unlist(case_ids), controls = unlist(control_ids))
  
  age <- as.data.frame(genotyped  %>% 
                         select(faminfo, goutaff, f.21003.0.0 , f.31.0.0) %>% 
                         filter(faminfo == 1, !is.na(goutaff), !is.na(f.21003.0.0)) %>% 
                         mutate(sex = ifelse(f.31.0.0 == 0, "Female", ifelse(f.31.0.0 == 1, "Male",NA)), cc = ifelse(goutaff == 1, "Case", "Control")) %>% 
                         group_by(cc) %>% 
                         summarise(mean = mean(f.21003.0.0), sd = sd(f.21003.0.0)))
  
  bmi <-as.data.frame(genotyped  %>% 
                        select(faminfo, goutaff, f.21001.0.0 , f.31.0.0) %>% 
                        filter(faminfo == 1, !is.na(goutaff), !is.na(f.21001.0.0)) %>% 
                        mutate(sex = ifelse(f.31.0.0 == 0, "Female", ifelse(f.31.0.0 == 1, "Male",NA)), cc = ifelse(goutaff == 1, "Case", "Control")) %>% 
                        group_by(cc) %>% 
                        summarise(mean = mean(f.21001.0.0), sd = sd(f.21001.0.0) ))
  
  sex <- as.data.frame(genotyped  %>% 
                         select(goutaff, f.31.0.0, faminfo) %>% 
                         filter(!is.na(goutaff), faminfo == 1 ) %>% 
                         mutate(sex = ifelse(f.31.0.0 == 0, "Female", ifelse(f.31.0.0 == 1, "Male",NA)) , cc = ifelse(goutaff == 1, "Case", "Control")) %>% 
                         group_by(cc, sex)  %>% tally())
  per_male <- as.data.frame(t(data.frame(c("Cases",percent_male(genotyped  %>% 
                                                                  select(goutaff, f.31.0.0, faminfo) %>% 
                                                                  filter(!is.na(goutaff), faminfo == 1 ) %>% 
                                                                  mutate(sex = ifelse(f.31.0.0 == 0, "Female", ifelse(f.31.0.0 == 1, "Male",NA)) , cc = ifelse(goutaff == 1, "Case", "Control")) 
                                                                %>% filter(cc == 'Case')  %>% 
                                                                  group_by(sex)  %>% tally())),
                                         c("Controls", percent_male(genotyped  %>% 
                                                                      select(goutaff, f.31.0.0, faminfo) %>% 
                                                                      filter(!is.na(goutaff), faminfo == 1 ) %>% 
                                                                      mutate(sex = ifelse(f.31.0.0 == 0, "Female", ifelse(f.31.0.0 == 1, "Male",NA)) , cc = ifelse(goutaff == 1, "Case", "Control")) %>% filter(cc == 'Control')  %>% 
                                                                      group_by(sex)  %>% tally()))
  )))
  colnames(per_male) <- c('cc','percentMale')
  rownames(per_male) <- 1:2
  per_male$percentMale <-as.numeric(as.character(per_male$percentMale))
  
  eth_codes <- as.numeric(levels(as.factor(genotyped[!is.na(genotyped$faminfo) & !is.na(genotyped$goutaff), "f.21000.0.0" ])))
  eth_case <- sapply(eth_codes ,function (x) {nrow(genotyped[!is.na(genotyped$faminfo) & genotyped$goutaff == 1 & genotyped$f.21000.0.0 == x, ])})
  
  eth_control <- sapply(eth_codes ,function (x) {nrow(genotyped[!is.na(genotyped$faminfo) & genotyped$goutaff == 0 & genotyped$f.21000.0.0 == x,  ])})
  
  it_summary[['age']] <-  age
  it_summary[['bmi']] <-  bmi
  it_summary[['sex']] <-  sex
  it_summary[['per_male']] <-  per_male
  it_summary[['eth_case']] <-  eth_case
  it_summary[['eth_control']] <-  eth_control
  rm(eth_control,eth_case,eth_codes,sex,bmi,age)
  
  
  drugs <-  drug_summary
  
  summary_table <- data.frame(
    # case = [,1], control = [,2]
    n = c(nrow(case_ids),nrow(control_ids)) ,
    Sex = c(round(per_male[1,2], digits = 2), round(per_male[2,2], digits = 2) ),
    Age = c(paste(round(it_summary[['age']][1,2],digits =2), "+/-", round(it_summary[['age']][1,3], digits = 2) , sep = " ") , paste(round(it_summary[['age']][2,2],digits =2), "+/-", round(it_summary[['age']][2,3], digits = 2) , sep = " ")),
    Ethnicity = c(paste(it_summary[['eth_case']], collapse = " "), paste(it_summary[['eth_control']], collapse = " ")),
    Diuretics = c(drugs$diuretics_cases, drugs$diuretics_controls),
    Allopurinol = c(drugs$allopurinol_cases, drugs$allopurinol_controls),
    Febuxostat = c(drugs$febuxostat_cases, drugs$febuxostat_controls),
    Sulphinpyrazone = c(drugs$sulphinpyrazone_cases, drugs$sulphinpyrazone_controls),
    Colchicine = c(drugs$colchicine_cases, drugs$colchicine_controls),
    BMI = c(paste(round(it_summary[['bmi']][1,2],digits =2), "+/-", round(it_summary[['bmi']][1,3], digits = 2) , sep = " ") , paste(round(it_summary[['bmi']][2,2],digits =2), "+/-", round(it_summary[['bmi']][2,3], digits = 2) , sep = " "))
  )
  
  rownames(summary_table) <-c("Gout", "Control")
  summary_table <-  t(summary_table)
  
  return(list(summary_table, drugs))
}


percent_male  <- function(sex_n) {
  # expects df with 'sex' {'Female','Male'} and 'n' columns
	return(sex_n[sex_n$sex == "Male",]$n/ sum(sex_n$n))
} 


fam_drugs <- function(cases, controls){
drug_tally <- data.frame(
  allopurinol_cases = as.numeric(genotyped  %>% filter(allopurinol == 1 & f.eid %in% cases & goutaff == 1)  %>% select(f.eid) %>% tally()),
  allopurinol_controls = as.numeric(genotyped  %>% filter(allopurinol == 1 & f.eid %in% controls & goutaff == 0)  %>% select(f.eid) %>% tally()),
  febuxostat_cases = as.numeric(genotyped  %>% filter(febuxostat == 1 & f.eid %in% cases & goutaff ==1)  %>% select(f.eid) %>% tally()),
  febuxostat_controls = as.numeric(genotyped  %>% filter(febuxostat == 1 & f.eid %in% controls & goutaff == 0)  %>% select(f.eid) %>% tally()),
  colchicine_cases = as.numeric(genotyped  %>% filter(colchicine == 1 & f.eid %in% cases & goutaff == 1)  %>% select(f.eid) %>% tally()),
  colchicine_controls = as.numeric(genotyped  %>% filter(colchicine == 1 & f.eid %in% controls & goutaff == 0)  %>% select(f.eid) %>% tally()),
  sulphinpyrazone_cases = as.numeric(genotyped  %>% filter(sulphinpyrazone == 1 & f.eid %in% cases & goutaff == 1)  %>% select(f.eid) %>% tally()),
  sulphinpyrazone_controls = as.numeric(genotyped  %>% filter(sulphinpyrazone == 1 & f.eid %in% controls & goutaff == 0)  %>% select(f.eid) %>% tally()),
  colchicine_alone_cases = as.numeric(genotyped  %>% 
                                  filter(f.eid %in% cases & goutaff == 1 & colchicine ==1 & !(allopurinol == 1 | sulphinpyrazone == 1 | febuxostat == 1))  %>% 
                                  select(f.eid) %>% tally()),
  colchicine_alone_controls = as.numeric(genotyped  %>% 
                                           filter(goutaff == 0 & f.eid %in% controls & colchicine ==1 & !(allopurinol == 1 | sulphinpyrazone == 1 | febuxostat == 1))  %>% 
                                           select(f.eid) %>% tally()),
  diuretics_cases = as.numeric(genotyped  %>% filter(f.eid %in% cases & goutaff == 1 & diuretics == 1)  %>% select(f.eid) %>% tally()),
  diuretics_controls = as.numeric(genotyped  %>% filter(f.eid %in% controls & diuretics == 1 & goutaff == 0)  %>% select(f.eid) %>% tally())
)
}

#d <- read_fams(paste0(scratch_dir,"/GWAS_all_controls_old/controls_no_diuretics/all/chrallimpv1.fam_all"))
