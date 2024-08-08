# code for the main analysis 
library(readxl)
library(stringr)
library(ggplot2)
library(survminer)
library(gridExtra)
library(survival)
library(lubridate)
library(ThresholdROC) # for diagnostics 
library(tidyr)
library(dplyr)
library(DescTools)
source("Analysis/data_preparation_Jun2024.R")

# -------------------
# Data description (first paragraph of results section):
# (i) How old were they on average; 
mean(dat.incl$age)
median(dat.incl$age)

# (ii) how many had all 5 tests at baseline
dim(dat.incl[!is.na(dat.incl$baseline_cyt) & !is.na(dat.incl$hc2) & 
               !is.na(dat.incl$abb) & !is.na(dat.incl$aly) & !is.na(dat.incl$cob),])[1]

# (iii) distribution of CIN2/3 cases at baseline 
# with exclusion criteria (CIN2+ within first 6 months removed)
table(dat.full[dat.full$cin2plus.cens==1 & exclusion & dat.full$cin2plus.days.cens < 365.25/12*6 &
                 !is.na(dat.full$cin2plus.cens),]$worst.histo)

# without exclusion criteria (CIN2+ within first 6 months included in numbers)
table(dat.incl[dat.incl$cin2plus.cens==1 & !is.na(dat.incl$cin2plus.cens),]$worst.histo)

# (iv) # mean record of  visits (and range)
date.cols <-colnames(slovenia.dat)[str_detect(colnames(slovenia.dat), "Date")][-1] # remove birthdate column
visits <- (apply(slovenia.dat[inclusion & !exclusion, date.cols], 1, function(x) length(unique(x[!is.na(x)])) ))
range(visits); mean(visits)

# (v) how many women included in the follow up at 3, 6, 9 years; 
# check if they had a test result ANY date after X-years but doesn't mean they had a test result AT X-years
sum(apply(slovenia.dat[inclusion & !exclusion, date.cols], 1, function(x) any(ymd(x) >= ymd(x[1]) + years(3), na.rm = T) )) #3581 --> 0.86498
sum(apply(slovenia.dat[inclusion & !exclusion, date.cols], 1, function(x) any(ymd(x) >= ymd(x[1]) + years(6), na.rm = T) )) #2630 --> 0.63527
sum(apply(slovenia.dat[inclusion & !exclusion, date.cols], 1, function(x) any(ymd(x) >= ymd(x[1]) + years(9), na.rm = T) )) #1684 --> 0.40676

# -------------------
# (1) Are there statistically significant differences for prevalence assessment between assays? 

# total study population 
# using paired data structure 
dat.incl.prev <- dat.full[inclusion & !is.na(dat.full$hc2) & !is.na(dat.full$abb) & !is.na(dat.full$aly) &!is.na(dat.full$cob),]
prev.dat <- data.frame(tests = c('HC2', "Abbott", "Alinity", "Cobas"), N = rep(4107, 4), n.pos = c(sum(dat.incl.prev$hc2==1), sum(dat.incl.prev$abb==1), sum(dat.incl.prev$aly==1), sum(dat.incl.prev$cob==1)))

prev.dat$percent.pos <- round(prev.dat$n.pos/prev.dat$N*100, 2)
prev.dat$conf.int <- paste0("[", round(BinomCI( prev.dat$n.pos, prev.dat$N)[,2]*100, 2), "-", round(BinomCI(prev.dat$n.pos, prev.dat$N)[,3]*100,2), "]")

prev.dat$conf.int <- paste0("[", round(BinomCI( prev.dat$n.pos, prev.dat$N)[,2]*100, 2), "-", round(BinomCI(prev.dat$n.pos, prev.dat$N)[,3]*100,2), "]")


prev.dat$p.val <- c("Reference", c(mcnemar.test(table(dat.incl.prev$hc2, dat.incl.prev$abb))$p.val, 
                                   mcnemar.test(table(dat.incl.prev$hc2, dat.incl.prev$aly))$p.val,
                                   mcnemar.test(table(dat.incl.prev$hc2, dat.incl.prev$cob))$p.val))
prev.dat
# women 30 or older
dat.incl.prev.30 <- dat.full[inclusion & dat.full$age >=30 & !is.na(dat.full$hc2) & !is.na(dat.full$abb) & !is.na(dat.full$aly) &!is.na(dat.full$cob),]
prev.dat30 <- data.frame(tests = c('HC2', "Abbott", "Alinity", "Cobas"), N = rep(2900, 4), n.pos = c(sum(dat.incl.prev.30$hc2==1), sum(dat.incl.prev.30$abb==1), sum(dat.incl.prev.30$aly==1), sum(dat.incl.prev.30$cob==1)))

prev.dat30$percent.pos <- round(prev.dat30$n.pos/prev.dat30$N*100, 2)
prev.dat30$conf.int <- paste0("[", round(BinomCI( prev.dat30$n.pos, prev.dat30$N)[,2]*100, 2), "-", round(BinomCI(prev.dat30$n.pos, prev.dat30$N)[,3]*100,2), "]")
prev.dat30$p.val <- c("Reference", c(mcnemar.test(table(dat.incl.prev.30$hc2, dat.incl.prev.30$abb))$p.val, 
                                     mcnemar.test(table(dat.incl.prev.30$hc2, dat.incl.prev.30$aly))$p.val,
                                     mcnemar.test(table(dat.incl.prev.30$hc2, dat.incl.prev.30$cob))$p.val))
prev.dat30

# Table 1. Baseline screening test positivity rates by four different clinically validated HPV assays used 
# for HPV testing in women ≥30 years old and in total study population (women aged 20-64 years).
data.frame(rbind(prev.dat30, prev.dat)) %>% write.csv(., "Results/Jun2024/Table_1.csv")

# -------------------
# (2a) Show that the difference between women with normal baseline cytology and negative HPV result
# is maintained throughout all 9 years (p<0,05) 
cin2plus.fit <- survfit(Surv(cin2plus.days.cens/365.25, cin2plus.cens) ~test, data=dat.long)
cin3plus.fit <- survfit(Surv(cin3plus.days.cens/365.25, cin3plus.cens) ~test, data=dat.long)
cin2plus.fit30 <- survfit(Surv(cin2plus.days.cens/365.25, cin2plus.cens) ~test, data=dat.long.30)
cin3plus.fit30 <- survfit(Surv(cin3plus.days.cens/365.25, cin3plus.cens) ~test, data=dat.long.30)

# (2b) determine when 3-year safety of normal cytology cross HPV negative curve to determine/suggest 
# length of interval screening
safety.hpv_cross_cytoloy <- function(fit, df){
  surv_summary(fit, df) %>%  filter(., test != 'Cytology' & round(surv, 3) == round(summary(fit, time=3)$surv[1], 3)) %>%
    group_by(test) %>% summarise(min.time = min(time)) %>% as.data.frame() %>% summarise(avg = mean(min.time)) %>% as.numeric()
}

cin2plus.fit.safety <- safety.hpv_cross_cytoloy(cin2plus.fit, dat.long)
cin3plus.fit.safety <- safety.hpv_cross_cytoloy(cin3plus.fit, dat.long)
cin2plus.fit30.safety <- safety.hpv_cross_cytoloy(cin2plus.fit30, dat.long.30)
cin3plus.fit30.safety <- safety.hpv_cross_cytoloy(cin3plus.fit30, dat.long.30)

# Table 2: The cumulative risk of cervical intraepithelial neoplasia grade 2 or worse (CIN2+) and grade 3 or worse (CIN3+) in women with 
# normal baseline cytology result or negative HPV baseline result at 3-years, 6-years and 9-years of follow-up, stratified by age groups.
risk.table <- function(x){
  return(data.frame(test=x$strata, time=x$time, risk=round((1-x$surv)*100, 2), lower=round((1-x$upper)*100,2), upper=round((1-x$lower)*100,2)))
}

table.function <- function(df, age, days, events){
  df %>% survfit(Surv(df[[days]]/365.25, df[[events]]) ~test, data=.) %>% summary(., times=c(3, 6, 9)) %>% risk.table() %>% 
    mutate(time= factor(time), risk = paste(risk, " (", lower, " - ", upper, ")", sep='')) %>% 
    select(test, time, risk) %>% tidyr::spread(time, risk) %>% mutate(age = age)
}

year.risks <- function(data, year){
  return(data %>% mutate(age = factor(age, levels =c("All", "20-29", "30-59", "60+"))) %>% select(test, "age", year) %>% group_by(test) %>%
           mutate(Id= row_number()) %>% ungroup %>% tidyr::spread(test, year) %>% select(!Id))
}

risk.results.cin2plus <- rbind(table.function(dat.long, "All", "cin2plus.days.cens", "cin2plus.cens"),
                               table.function(dat.long.20to30, "20-29", "cin2plus.days.cens", "cin2plus.cens"),
                               table.function(dat.long.30to59, "30-59", "cin2plus.days.cens", "cin2plus.cens"),
                               table.function(dat.long.60plus, "60+", "cin2plus.days.cens", "cin2plus.cens"))

data.frame(cbind(year.risks(risk.results.cin2plus, "3"), 
                 year.risks(risk.results.cin2plus, "6"), 
                 year.risks(risk.results.cin2plus, "9"))) %>% write.csv(., "Results/Jun2024/risksCIN2+_Jun2021.csv") #CIN2+ 3/6/9 year risks

risk.results.cin3plus <- rbind(table.function(dat.long, "All", "cin3plus.days.cens", "cin3plus.cens"),
                               table.function(dat.long.20to30, "20-29", "cin3plus.days.cens", "cin3plus.cens"),
                               table.function(dat.long.30to59, "30-59", "cin3plus.days.cens", "cin3plus.cens"),
                               table.function(dat.long.60plus, "60+", "cin3plus.days.cens", "cin3plus.cens"))

data.frame(cbind(year.risks(risk.results.cin3plus, "3"), 
                 year.risks(risk.results.cin3plus, "6"),
                 year.risks(risk.results.cin3plus, "9"))) %>%  write.csv(., "Results/Jun2024/risksCIN3+_Jun2024.csv") #CIN3+ 3/6/9 year risks



# -------------------
# (3)	NPV of normal baseline cytology and negative HPV result at 9 years; regardless of the screening test used, 
# they all offer comparable safety if they are clinically validated

# --> include women with baseline CIN2+ (previous exclusion criteria)
dat.incl.npv <- dat.full[inclusion, ]
dat.incl.npv.30 <- dat.full[inclusion & dat.full$age >= 30,]

diagnostic.function <- function(df, col, event){
  dat <- format(round(diagnostic(table(df[[col]], df[[event]]), method="exact")*100, 2), nsmall=2)
  return(noquote(apply(trimws(dat), 1, function(x) paste0(x[1], " (", x[2], " - ", x[3], ")")))[1:4])
}

cin2plus.whole <- data.frame(t(sapply(list("baseline_cyt", "hc2", "abb", "cob", "aly"), function(y) diagnostic.function(dat.incl.npv, y, "cin2plus.cens.factor"))))
cin3plus.whole <- data.frame(t(sapply(list("baseline_cyt", "hc2", "abb", "cob", "aly"), function(y) diagnostic.function(dat.incl.npv, y, "cin3plus.cens.factor"))))
cin2plus.30 <-  data.frame(t(sapply(list("baseline_cyt", "hc2", "abb", "cob", "aly"), function(y) diagnostic.function(dat.incl.npv.30, y, "cin2plus.cens.factor"))))
cin3plus.30 <- data.frame(t(sapply(list("baseline_cyt", "hc2", "abb", "cob", "aly"), function(y) diagnostic.function(dat.incl.npv.30, y, "cin3plus.cens.factor"))))

# Supplementary Table S1: Longitudinal test characteristics using cervical intraepithelial neoplasia grade 2 or worse (CIN2+) and grade 3 or worse (CIN3+) as outcome 
# in women by their baseline cytology and HPV result at 9 years of follow-up in women ≥ 30 years old and in total study population.
diagnostics <- data.frame(endpoint = rep(rep(c("CIN2+", "CIN3+"), each=5), 2),
                          rbind(cin2plus.30, cin3plus.30, cin2plus.whole, cin3plus.whole)) # results table

data.frame(t(apply(diagnostics[,2:5], 1, function(x) {
  temp <- strsplit(x, "\\(")
  return(unlist(lapply(temp, function(y) c(y[1], substr(y[2], 1, nchar(y[2])-1)))))}))) %>% 
  mutate(population=rep(c("30+", "total pop"), each=10), .before = Sensitivity1) %>% 
  mutate(endpoint=rep(rep(c("CIN2+", "CIN3+"), each=5), 2), .before = Sensitivity1) %>% 
  mutate(test = rep(c("cyt", "hc2", "abb", "cob", "aly"), 4), .before = Sensitivity1) %>% 
  write.csv(., "Results/Jun2024/Table_S1.csv")

# p-values:

# The patients with a (+, +) result and the patients with a ( - , - ) result do not distinguish between the two 
# diagnostic tests. The only information for comparing the sensitivities of the two diagnostic tests comes from 
# those patients with a (+, - ) or ( - , +) result.
pairs <- t(combn(c("baseline_cyt", "hc2", "abb", "aly", "cob"), 2))

# Endpoint CIN2+
# sensitivities of each test compared to each other - only consider CIN2+ women and then compare pairs of tests:
cin2plus.1st.round <- dat.incl.npv[dat.incl.npv$cin2plus.1st.round.cens ==1, ]
cin2plus.overall <- dat.incl.npv[dat.incl.npv$cin2plus.cens ==1 , ]

# specificities of each test compared to each other - only consider women without CIN2+ and then compare pairs of tests:
no.cin2plus.1st.round <- dat.incl.npv[dat.incl.npv$cin2plus.1st.round.cens ==0 , ]
no.cin2plus.overall <- dat.incl.npv[dat.incl.npv$cin2plus.cens ==0 , ]

mcnemar.pval <- data.frame( 
  sens.first.round = sapply(1:10, function(i) mcnemar.test(table(cin2plus.1st.round[[pairs[i,1]]], cin2plus.1st.round[[pairs[i,2]]], useNA='no'))$p.value), 
  spec.first.round = sapply(1:10, function(i) mcnemar.test(table(no.cin2plus.1st.round[[pairs[i,1]]], no.cin2plus.1st.round[[pairs[i,2]]], useNA='no'))$p.value),
  sens.other = sapply(1:10, function(i) mcnemar.test(table(cin2plus.overall[[pairs[i,1]]], cin2plus.overall[[pairs[i,2]]], useNA='no'))$p.value),
  spec.other = sapply(1:10, function(i) mcnemar.test(table(no.cin2plus.overall[[pairs[i,1]]], no.cin2plus.overall[[pairs[i,2]]], useNA='no'))$p.value))

cbind(pairs, round(mcnemar.pval, 3))

# Endpoint CIN3+ 
# sensitivities of each test compared to each other - only consider CIN3+ women and then compare pairs of tests:
cin3plus.1st.round <- dat.incl.npv[dat.incl.npv$cin3plus.1st.round.cens ==1, ]
cin3plus.overall <- dat.incl.npv[dat.incl.npv$cin3plus.cens ==1 , ]

# specificities of each test compared to each other - only consider women without CIN3+ and then compare pairs of tests:
no.cin3plus.1st.round <- dat.incl.npv[dat.incl.npv$cin3plus.1st.round.cens ==0 , ]
no.cin3plus.overall <- dat.incl.npv[dat.incl.npv$cin3plus.cens ==0 , ]

mcnemar.pvalCIN3 <- data.frame( 
  sens.first.round = sapply(1:10, function(i) mcnemar.test(table(cin3plus.1st.round[[pairs[i,1]]], cin3plus.1st.round[[pairs[i,2]]], useNA='no'))$p.value), 
  spec.first.round = sapply(1:10, function(i) mcnemar.test(table(no.cin3plus.1st.round[[pairs[i,1]]], no.cin3plus.1st.round[[pairs[i,2]]], useNA='no'))$p.value),
  sens.other = sapply(1:10, function(i) mcnemar.test(table(cin3plus.overall[[pairs[i,1]]], cin3plus.overall[[pairs[i,2]]], useNA='no'))$p.value),
  spec.other = sapply(1:10, function(i) mcnemar.test(table(no.cin3plus.overall[[pairs[i,1]]], no.cin3plus.overall[[pairs[i,2]]], useNA='no'))$p.value))

cbind(pairs, round(mcnemar.pvalCIN3, 3))

# In the column for sensitivity, all the pairs with NA are because they completely agreed with each other (see for loop below), 
# this was because we excluded cin2+ cases within 6 months after baseline --> no difference in sens/spec (perfect agreement), so p-value=1
# for (i in 8:10){
#   print(table(cin3plus.1st.round[[pairs[i,1]]], cin3plus.1st.round[[pairs[i,2]]], useNA='no'))
# }

# p-values for difference in NPV 
npv.pval <- function(dat, test, outcome){
  return(prop.test(x=c(sum(dat[[test]]==0 & dat[[outcome]]==0, na.rm=T),
                       sum(dat$baseline_cyt==0 & dat[[outcome]]==0, na.rm=T)),
                   n=c(sum(dat[[test]]==0 & dat[[outcome]]==0, na.rm=T) + sum(dat[[test]]==0 & dat[[outcome]]==1, na.rm=T),
                       sum(dat$baseline_cyt==0 & dat[[outcome]]==0, na.rm=T) + sum(dat$baseline_cyt==0 & dat[[outcome]]==1, na.rm=T)))$p.value
  )
}

npv.cin2 <- data.frame(npv = diagnostics$Neg.Pred.Val.[diagnostics$endpoint == 'CIN2+'],
                       p.val = c("ref", round(c(npv.pval(dat.incl.npv.30, "hc2", "cin2plus.cens"),
                                                npv.pval(dat.incl.npv.30, "abb", "cin2plus.cens"),
                                                npv.pval(dat.incl.npv.30, "cob", "cin2plus.cens"),
                                                npv.pval(dat.incl.npv.30, "aly", "cin2plus.cens")), 5),
                                 "ref", round(c(npv.pval(dat.incl.npv, "hc2", "cin2plus.cens"),
                                                npv.pval(dat.incl.npv, "abb", "cin2plus.cens"),
                                                npv.pval(dat.incl.npv, "cob", "cin2plus.cens"),
                                                npv.pval(dat.incl.npv, "aly", "cin2plus.cens")), 5)))

npv.cin3 <- data.frame(npv = diagnostics$Neg.Pred.Val.[diagnostics$endpoint == 'CIN3+'],
                       p.val = c("ref", round(c(npv.pval(dat.incl.npv.30, "hc2", "cin3plus.cens"),
                                                npv.pval(dat.incl.npv.30, "abb", "cin3plus.cens"),
                                                npv.pval(dat.incl.npv.30, "cob", "cin3plus.cens"),
                                                npv.pval(dat.incl.npv.30, "aly", "cin3plus.cens")), 5),
                                 "ref", round(c(npv.pval(dat.incl.npv, "hc2", "cin3plus.cens"),
                                                npv.pval(dat.incl.npv, "abb", "cin3plus.cens"),
                                                npv.pval(dat.incl.npv, "cob", "cin3plus.cens"),
                                                npv.pval(dat.incl.npv, "aly", "cin3plus.cens")), 5)))
cbind(npv.cin2, npv.cin3) 


# -------------------
# (4) co-testing results 

surv.9year <- function(survmodel){
  df <- do.call(data.frame, lapply(c(10, 6, 15, 16) , function(x) summary(survmodel, times=9)[x])) %>% 
    mutate_if(is.numeric, function(x) round((1-x)*100, 2))
  colnames(df) <- c("strata", "surv", "lower", "upper")
  return(df[,c(1,2,4,3)])
}
# endpoint CIN2+
hc2.cotest.fit <- survfit(Surv(cin2plus.days.cens/365.25, cin2plus.cens) ~test, data=hc2.cotest)
abb.cotest.fit <- survfit(Surv(cin2plus.days.cens/365.25, cin2plus.cens) ~test, data=abb.cotest)
aly.cotest.fit <- survfit(Surv(cin2plus.days.cens/365.25, cin2plus.cens) ~test, data=aly.cotest)
cob.cotest.fit <- survfit(Surv(cin2plus.days.cens/365.25, cin2plus.cens) ~test, data=cob.cotest)

# 9-year risk of CIN2+ after co-test negative
rbind(surv.9year(hc2.cotest.fit)[3,],
      surv.9year(abb.cotest.fit)[3,],
      surv.9year(aly.cotest.fit)[3,],
      surv.9year(cob.cotest.fit)[3,])

# 9-year risk of CIN2+ after HPV test negative
rbind(surv.9year(hc2.cotest.fit)[2,],
      surv.9year(abb.cotest.fit)[2,],
      surv.9year(aly.cotest.fit)[2,],
      surv.9year(cob.cotest.fit)[2,])

# 9-year risk of CIN2+ after cytology negative (same for all assays because its based on cytology so just take 1 value)
surv.9year(cob.cotest.fit)[1,] 

# endpoint CIN3+
# 9-year risk of CIN3+ after co-test negative
rbind(surv.9year(survfit(Surv(cin3plus.days.cens/365.25, cin3plus.cens) ~test, data=hc2.cotest))[3,],
      surv.9year(survfit(Surv(cin3plus.days.cens/365.25, cin3plus.cens) ~test, data=abb.cotest))[3,],
      surv.9year(survfit(Surv(cin3plus.days.cens/365.25, cin3plus.cens) ~test, data=aly.cotest))[3,],
      surv.9year(survfit(Surv(cin3plus.days.cens/365.25, cin3plus.cens) ~test, data=cob.cotest))[3,])

# 9-year risk of CIN2+ after HPV test negative
rbind(surv.9year(survfit(Surv(cin3plus.days.cens/365.25, cin3plus.cens) ~test, data=hc2.cotest))[2,],
      surv.9year(survfit(Surv(cin3plus.days.cens/365.25, cin3plus.cens) ~test, data=abb.cotest))[2,],
      surv.9year(survfit(Surv(cin3plus.days.cens/365.25, cin3plus.cens) ~test, data=aly.cotest))[2,],
      surv.9year(survfit(Surv(cin3plus.days.cens/365.25, cin3plus.cens) ~test, data=cob.cotest))[2,])

# 9-year risk of CIN3+ after cytology negative (same for all assays because its based on cytology so just take 1 value)
surv.9year(survfit(Surv(cin3plus.days.cens/365.25, cin3plus.cens) ~test, data=cob.cotest))[1,]

# -------------------
# (5) 9-year risks after HPV+ test 
#cin2plus.pos.fit <- survfit(Surv(cin2plus.days.cens/365.25, cin2plus.cens) ~test, data=dat.long.pos)
#cin3plus.pos.fit <- survfit(Surv(cin3plus.days.cens/365.25, cin3plus.cens) ~test, data=dat.long.pos)
cin2plus.pos.fit30 <- survfit(Surv(cin2plus.days.cens/365.25, cin2plus.cens) ~test, data=dat.long.pos.30)
cin3plus.pos.fit30 <- survfit(Surv(cin3plus.days.cens/365.25, cin3plus.cens) ~test, data=dat.long.pos.30)

surv.9year(cin2plus.pos.fit30)[2:5,]
surv.9year(cin3plus.pos.fit30)[2:5,]

# -------------------
# (6) Genotype analysis - Extended genotyping can offer better risk stratification; for all genotypes if data are sufficient,
#  otherwise grouping HPV genotypes according to their risk 
#  
cin2plus.genotype.fit <- survfit(Surv(cin2plus.days.cens/365.25, cin2plus.cens) ~genotype, data=dat.incl)
cin3plus.genotype.fit <- survfit(Surv(cin3plus.days.cens/365.25, cin3plus.cens) ~genotype, data=dat.incl)
cin2plus.genotype.fit30 <- survfit(Surv(cin2plus.days.cens/365.25, cin2plus.cens) ~genotype, data=dat.incl.30)
cin3plus.genotype.fit30 <- survfit(Surv(cin3plus.days.cens/365.25, cin3plus.cens) ~genotype, data=dat.incl.30)

surv.9year(cin2plus.genotype.fit30)
surv.9year(cin3plus.genotype.fit30)

# -------------------
# (7) 9-year risk after Npos test results
ids <- rownames(dat.incl[!is.na(dat.incl$hc2) & !is.na(dat.incl$abb) & !is.na(dat.incl$aly) &!is.na(dat.incl$cob), ])
hpv.tests<- apply(select(dat.incl[!is.na(dat.incl$hc2) & !is.na(dat.incl$abb) & !is.na(dat.incl$aly) &!is.na(dat.incl$cob), ],
                         "hc2", "abb", "aly", "cob", "cin2plus.days.cens", "cin2plus.cens", "cin2plus",  "cin3plus.days.cens", "cin3plus.cens", "age"), 2, as.numeric) %>% 
  as.data.frame() %>% mutate(n.pos = rowSums(.[1:4], na.rm=T))

hpv.tests$n.pos <- factor(case_when(
  hpv.tests$n.pos == 0 ~ "0 tests positive",
  hpv.tests$n.pos == 1 ~ "1 test positive",
  hpv.tests$n.pos == 2 ~ "2 tests positive",
  hpv.tests$n.pos == 3 ~ "3 tests positive",
  hpv.tests$n.pos == 4 ~ "4 tests positive",
), levels = c("0 tests positive", "1 test positive", "2 tests positive", "3 tests positive", "4 tests positive"))
hpv.tests$id <- ids

n.pos.CIN2 <- survfit(Surv(cin2plus.days.cens/365.25, cin2plus.cens) ~n.pos, data=hpv.tests)
n.pos.CIN3 <- survfit(Surv(cin2plus.days.cens/365.25, cin3plus.cens) ~n.pos, data=hpv.tests)
n.pos.CIN2.30 <- survfit(Surv(cin2plus.days.cens/365.25, cin2plus.cens) ~n.pos, data=hpv.tests[hpv.tests$age >=30,])
n.pos.CIN3.30 <- survfit(Surv(cin2plus.days.cens/365.25, cin3plus.cens) ~n.pos, data=hpv.tests[hpv.tests$age >=30,])

# 9year survival
surv.9year(n.pos.CIN2)[,2:4]
surv.9year(n.pos.CIN3)
surv.9year(n.pos.CIN2.30)[,2:4]
surv.9year(n.pos.CIN3.30)

# concordance 
hpv.pos.30 <- hpv.tests[hpv.tests$n.pos != '0 tests positive' & hpv.tests$age >=30,] #hpv.pos[hpv.pos$age >= 30,]
concordance.30plus <- hpv.pos.30 %>% group_by(n.pos) %>% summarise(n = n(), N= dim(hpv.pos.30)[1],
                                             perc = round(n/N*100, 1), 
                                             lower = round(BinomCI(n, N, method='wilson')[2]*100,1), 
                                             upper = round(BinomCI(n, N, method='wilson')[3]*100, 1), 
                                             n.cin2pl = sum(cin2plus.cens, na.rm=T),
                                             N.cin2pl = sum(hpv.pos.30$cin2plus.cens==1, na.rm=T),
                                             perc.cin2pl = round(n.cin2pl/N.cin2pl*100, 1),
                                             lower.cin2pl = round(BinomCI(n.cin2pl, N.cin2pl, method='wilson')[2]*100,1), 
                                             upper.cin2pl = round(BinomCI(n.cin2pl, N.cin2pl, method='wilson')[3]*100, 1)) 

hpv.pos <- hpv.tests[hpv.tests$n.pos != '0 tests positive',]
concordance.total.pop <- hpv.pos %>% group_by(n.pos) %>% summarise(n = n(), N= dim(hpv.pos)[1],
                                          perc = round(n/N*100, 1), 
                                          lower = round(BinomCI(n, N, method='wilson')[2]*100,1), 
                                          upper = round(BinomCI(n, N, method='wilson')[3]*100, 1),
                                          n.cin2pl = sum(cin2plus.cens, na.rm=T),
                                          N.cin2pl = sum(hpv.pos$cin2plus.cens==1, na.rm=T),
                                          perc.cin2pl = round(n.cin2pl/N.cin2pl*100, 1),
                                          lower.cin2pl = round(BinomCI(n.cin2pl, N.cin2pl, method='wilson')[2]*100,1), 
                                          upper.cin2pl = round(BinomCI(n.cin2pl, N.cin2pl, method='wilson')[3]*100, 1))


# PCR tests only 
hpv.tests.PCR <- apply(select(dat.incl[!is.na(dat.incl$hc2) & !is.na(dat.incl$abb) & !is.na(dat.incl$aly) &!is.na(dat.incl$cob), ],
                         "abb", "aly", "cob", "cin2plus.days.cens", "cin2plus.cens", "cin2plus",  "cin3plus.days.cens", "cin3plus.cens", "age"), 2, as.numeric) %>% 
  as.data.frame() %>% mutate(n.pos = rowSums(.[1:3], na.rm=T))

hpv.tests.PCR$n.pos <- factor(case_when(
  hpv.tests.PCR$n.pos == 0 ~ "0 tests positive",
  hpv.tests.PCR$n.pos == 1 ~ "1 test positive",
  hpv.tests.PCR$n.pos == 2 ~ "2 tests positive",
  hpv.tests.PCR$n.pos == 3 ~ "3 tests positive"
), levels = c("0 tests positive", "1 test positive", "2 tests positive", "3 tests positive"))
hpv.tests.PCR$id <- ids

hpv.pos.PCR.30 <- hpv.tests.PCR[hpv.tests.PCR$age >= 30 & hpv.tests.PCR$n.pos != '0 tests positive',]
concordance.30plus.PCR <- hpv.pos.PCR.30 %>% group_by(n.pos) %>% summarise(n = n(), N= dim(hpv.pos.PCR.30)[1],
                                                                   perc = round(n/N*100, 1), 
                                                                   lower = round(BinomCI(n, N, method='wilson')[2]*100,1), 
                                                                   upper = round(BinomCI(n, N, method='wilson')[3]*100, 1), 
                                                                   n.cin2pl = sum(cin2plus.cens, na.rm=T),
                                                                   N.cin2pl = sum(hpv.pos.PCR.30$cin2plus.cens==1, na.rm=T),
                                                                   perc.cin2pl = round(n.cin2pl/N.cin2pl*100, 1),
                                                                   lower.cin2pl = round(BinomCI(n.cin2pl, N.cin2pl, method='wilson')[2]*100,1), 
                                                                   upper.cin2pl = round(BinomCI(n.cin2pl, N.cin2pl, method='wilson')[3]*100, 1)) 

hpv.pos.PCR <- hpv.tests.PCR[hpv.tests.PCR$n.pos != '0 tests positive',]
concordance.total.pop.PCR <- hpv.pos.PCR %>% group_by(n.pos) %>% summarise(n = n(), N= dim(hpv.pos.PCR)[1],
                                                                           perc = round(n/N*100, 1), 
                                                                           lower = round(BinomCI(n, N, method='wilson')[2]*100,1), 
                                                                           upper = round(BinomCI(n, N, method='wilson')[3]*100, 1), 
                                                                           n.cin2pl = sum(cin2plus.cens, na.rm=T),
                                                                           N.cin2pl = sum(hpv.pos.PCR$cin2plus.cens==1, na.rm=T),
                                                                           perc.cin2pl = round(n.cin2pl/N.cin2pl*100, 1),
                                                                           lower.cin2pl = round(BinomCI(n.cin2pl, N.cin2pl, method='wilson')[2]*100,1), 
                                                                           upper.cin2pl = round(BinomCI(n.cin2pl, N.cin2pl, method='wilson')[3]*100, 1)) 

# Table 3:  Concordance between four clinically validated HPV assays by number of concurrent baseline HPV-positive results and cervical 
# intraepithelial neoplasia grade 2 or worse (CIN2+) identified over 9-year follow-up in women ≥ 30 years old and in total study population.
rbind(rbind(concordance.30plus, concordance.total.pop) %>% mutate(assays = rep("All 4 assays", 8), .before = n.pos) %>% mutate(endpoint = rep(c("30+", "total pop"), each=4), .before = n.pos ), 
      rbind(concordance.30plus.PCR, concordance.total.pop.PCR) %>% mutate(assays = rep("PCR assays", 6), .before = n.pos) %>%  mutate(endpoint = rep(c("30+", "total pop"), each=3), .before = n.pos )) %>% 
  write.csv(., "Results/Jun2024/Table_3.csv")

# -------------------
# (8) Response to reviewers 
# 
# Table 2, and Figures 1 - 3 refer risk of CIN2+ at various points beyond baseline. Within the text it is noted that there are 32 CIN2, 61 CIN3, 4 SCC
# and an adenocarcinoma identified through the follow up period but it is unclear at what timepoint. Information on how many of these abnormal findings 
# were identified, e.g. post baseline  to 3 years, post 3 year to 6 year, and so on.  
dat.incl$baseline_to_3years <- ifelse(ymd(dat.incl$worst.histo.date) < ymd(dat.incl$baseline) + years(3), 1, 0)
dat.incl$three_to_6years <- ifelse((ymd(dat.incl$worst.histo.date) > ymd(dat.incl$baseline) + years(3)) & ymd(dat.incl$worst.histo.date) < ymd(dat.incl$baseline) + years(6), 1, 0)
dat.incl$six_to_9years <- ifelse((ymd(dat.incl$worst.histo.date) > ymd(dat.incl$baseline) + years(6)) & ymd(dat.incl$worst.histo.date) < ymd(dat.incl$baseline) + years(9), 1, 0)
dat.incl$nine_plus <- ifelse((ymd(dat.incl$worst.histo.date) > ymd(dat.incl$baseline) + years(9)), 1, 0)

dat.incl$histo.time <- case_when(
  !is.na(dat.incl$baseline_to_3years) & dat.incl$baseline_to_3years==1 ~ "post baseline  to 3 years",
  !is.na(dat.incl$three_to_6years) & dat.incl$three_to_6years==1 ~ "post 3 years to 6 years",
  !is.na(dat.incl$six_to_9years) & dat.incl$six_to_9years==1 ~ "post 6 years to 9 years",
  !is.na(dat.incl$nine_plus) & dat.incl$nine_plus==1 ~ "post 9 years",
  T ~ NA
)

dat.incl$histo.time <- factor(dat.incl$histo.time, levels = c("post baseline  to 3 years", "post 3 years to 6 years", "post 6 years to 9 years", "post 9 years", NA))

table(dat.incl$histo.time[dat.incl$worst.histo >=7], dat.incl$worst.histo[dat.incl$worst.histo >=7])

# Table 3 is very interesting, can the authors do the same for CIN3+ :  
# concordance 
concordance.30plus.cin3pl <- hpv.pos.30 %>% group_by(n.pos) %>% summarise(n = n(), N= dim(hpv.pos.30)[1],
                                                                   perc = round(n/N*100, 1), 
                                                                   lower = round(BinomCI(n, N, method='wilson')[2]*100,1), 
                                                                   upper = round(BinomCI(n, N, method='wilson')[3]*100, 1), 
                                                                   n.cin3pl = sum(cin3plus.cens, na.rm=T),
                                                                   N.cin3pl = sum(hpv.pos.30$cin3plus.cens==1, na.rm=T),
                                                                   perc.cin3pl = round(n.cin3pl/N.cin3pl*100, 1),
                                                                   lower.cin3pl = round(BinomCI(n.cin3pl, N.cin3pl, method='wilson')[2]*100,1), 
                                                                   upper.cin3pl = round(BinomCI(n.cin3pl, N.cin3pl, method='wilson')[3]*100, 1)) 

concordance.total.pop.cin3pl <- hpv.pos %>% group_by(n.pos) %>% summarise(n = n(), N= dim(hpv.pos)[1],
                                                                   perc = round(n/N*100, 1), 
                                                                   lower = round(BinomCI(n, N, method='wilson')[2]*100,1), 
                                                                   upper = round(BinomCI(n, N, method='wilson')[3]*100, 1),
                                                                   n.cin3pl = sum(cin3plus.cens, na.rm=T),
                                                                   N.cin3pl = sum(hpv.pos$cin3plus.cens==1, na.rm=T),
                                                                   perc.cin3pl = round(n.cin3pl/N.cin3pl*100, 1),
                                                                   lower.cin3pl = round(BinomCI(n.cin3pl, N.cin3pl, method='wilson')[2]*100,1), 
                                                                   upper.cin3pl = round(BinomCI(n.cin3pl, N.cin3pl, method='wilson')[3]*100, 1))

concordance.30plus.PCR.cin3pl <- hpv.pos.PCR.30 %>% group_by(n.pos) %>% summarise(n = n(), N= dim(hpv.pos.PCR.30)[1],
                                                                           perc = round(n/N*100, 1), 
                                                                           lower = round(BinomCI(n, N, method='wilson')[2]*100,1), 
                                                                           upper = round(BinomCI(n, N, method='wilson')[3]*100, 1), 
                                                                           n.cin3pl = sum(cin3plus.cens, na.rm=T),
                                                                           N.cin3pl = sum(hpv.pos.PCR.30$cin3plus.cens==1, na.rm=T),
                                                                           perc.cin3pl = round(n.cin3pl/N.cin3pl*100, 1),
                                                                           lower.cin3pl = round(BinomCI(n.cin3pl, N.cin3pl, method='wilson')[2]*100,1), 
                                                                           upper.cin3pl = round(BinomCI(n.cin3pl, N.cin3pl, method='wilson')[3]*100, 1)) 

concordance.total.pop.PCR.cin3pl <- hpv.pos.PCR %>% group_by(n.pos) %>% summarise(n = n(), N= dim(hpv.pos.PCR)[1],
                                                                           perc = round(n/N*100, 1), 
                                                                           lower = round(BinomCI(n, N, method='wilson')[2]*100,1), 
                                                                           upper = round(BinomCI(n, N, method='wilson')[3]*100, 1), 
                                                                           n.cin3pl = sum(cin3plus.cens, na.rm=T),
                                                                           N.cin3pl = sum(hpv.pos.PCR$cin3plus.cens==1, na.rm=T),
                                                                           perc.cin3pl = round(n.cin3pl/N.cin3pl*100, 1),
                                                                           lower.cin3pl = round(BinomCI(n.cin3pl, N.cin3pl, method='wilson')[2]*100,1), 
                                                                           upper.cin3pl = round(BinomCI(n.cin3pl, N.cin3pl, method='wilson')[3]*100, 1)) 


# Table 3 CIN3+
rbind(rbind(concordance.30plus.cin3pl, concordance.total.pop.cin3pl) %>% mutate(assays = rep("All 4 assays", 8), .before = n.pos) %>% mutate(endpoint = rep(c("30+", "total pop"), each=4), .before = n.pos ), 
      rbind(concordance.30plus.PCR.cin3pl, concordance.total.pop.PCR.cin3pl) %>% mutate(assays = rep("PCR assays", 6), .before = n.pos) %>%  mutate(endpoint = rep(c("30+", "total pop"), each=3), .before = n.pos )) %>% 
  mutate(endpoint = 'cin3plus', .before = assays) %>% 
  write.csv(., "Results/Jun2024/Table_3_CIN3+.csv")
