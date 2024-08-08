# code for Data preparation 
slovenia.dat <- read.csv("Data/2022-10-22-MASTER_SEND.csv")[1:4510,]
#slovenia.dat <- as.data.frame(readxl::read_excel("Data/2022-10-22-MASTER_SEND.xlsx"))
rownames(slovenia.dat) <- slovenia.dat$ID
# create an indicator variable for normal cytology at baseline (1 = abnormal, 0 = normal cytology)
slovenia.dat$baseline_cyt <- ifelse(slovenia.dat$Cytology_result_1 == 1, 0, 1)

# there are multiple formats of dates in some columns so we convert them so they are all matching formats 
fix.date.col <- function(column){
  d <- parse_date_time(column, orders = "dmy") # c(c('%d-%m-%Y', '%Y/%m/%d', "%d/%m/%Y", "%d.%m.%Y"))
  d[is.na(d)] <- as.POSIXct(as.Date(as.numeric(column[is.na(d)]), origin ="1899-12-30"))
  return(format(as.Date(d), '%Y/%m/%d'))
} 

# apply the fix.date.col() function to all columns which contain dates 
# ignore warning message: that is because some dates are NA
slovenia.dat[,colnames(slovenia.dat)[str_detect(colnames(slovenia.dat), "Date")]] <- apply(slovenia.dat[,colnames(slovenia.dat)[str_detect(colnames(slovenia.dat), "Date")]], 2, fix.date.col)

# store the hpv test data and dates
hpvs <-  slovenia.dat[, colnames(slovenia.dat)[str_detect(colnames(slovenia.dat),"_res_" )]]
hpvs <- hpvs[, colnames(hpvs)[str_detect(colnames(hpvs), "Gen", negate = T)]] # remove genotype columns
hpvs.date <- slovenia.dat[,colnames(slovenia.dat)[str_detect(colnames(slovenia.dat), "DateHPV")]]

cyto <- slovenia.dat[, colnames(slovenia.dat)[str_detect(colnames(slovenia.dat),"Cytology_" )]]
cyto.date <-  slovenia.dat[, colnames(slovenia.dat)[str_detect(colnames(slovenia.dat),"Date_visit_" )]]

# First we use create_data() which formats the full data into the information we need like cin2+/cin3+ indicators and times since baseline 
# --> makes it easy to do the analysis 
# -------------------
create_data <- function(){
  # find date of last Pap 1 which can be used as substitute of CIN0/1 date
  # find the index of the first repeat cytology by finding the first non-empty entry in each row 
  last.cyt.results <- trimws(apply(cyto, 1, function(i) ifelse(sum(is.na(i)) == 13, NA, tail(na.omit(i), 1))))
  last.cyt.date <- apply(cyto.date, 1, function(i) ifelse(sum(is.na(i)) == 13, NA, tail(na.omit(i), 1)))
  
  last.pap.1 <- cyto[cbind(seq_along(apply(cyto, 1, function(x) max(which(x == 1)))),  apply(cyto, 1, function(x) max(which(x == 1))))]
  last.pap.1.date <- cyto.date[cbind(seq_along(apply(cyto, 1, function(x) max(which(x == 1)))),  apply(cyto, 1, function(x) max(which(x == 1))))]
  
  # create a data frame which includes date and time of worst histology result (to determine if CIN2/CIN3+)
  histo <- slovenia.dat[, colnames(slovenia.dat)[str_detect(colnames(slovenia.dat),"Histo_" )]]
  histo.date <-  slovenia.dat[, colnames(slovenia.dat)[str_detect(colnames(slovenia.dat),"Date_histo_" )]]
  which.worst <- max.col(replace(histo, is.na(histo), -Inf), ties.method="first") * NA^!rowSums(!is.na(histo)) 
  worst.histo <- do.call(pmax, c(histo, list(na.rm=TRUE))) # find worst histology result
  worst.histo.date <- histo.date[cbind(seq_along(which.worst), which.worst)] # find worst histology date
  
  # replace unknown histology values with last Pap 1 result and date (if available) - since women were only referred for colposcopy after abnormal cytology results 
  # (if you don't do this then the risk estimated for CIN2+ are round 40% at 9 years rather than 3%)
  # -> basically assuming that CIN2/3+ was not present because they were not referred 
  worst.histo.cyt <- ifelse((is.na(worst.histo) & is.na(worst.histo.date))  & !is.na(last.pap.1), 1 ,worst.histo)
  worst.histo.cyt.date <- ifelse((is.na(worst.histo) & is.na(worst.histo.date)) & !is.na(last.pap.1), last.pap.1.date, worst.histo.date)
  
  # the only ones with missing histology results now are those who had abnormal cytology and never returned for histology 
  
  # make an indicator variable for whether the worst cytology result was CIN2+ or not (using definitions from the data dictionary)
  cin2plus <- ifelse(worst.histo.cyt >=7 & worst.histo.cyt!=13 & worst.histo!=11, 1, 0)
  # store the time difference between the 
  cin2plus.days <- difftime(worst.histo.cyt.date, slovenia.dat$Date_visit_1, units="days")
  cin2plus.days[cin2plus.days<=0] <- 0
  # add right-censoring at 9 years ---> if their CIN2+ days is higher than 9 years we set CIN2+ indicator to 0 and set days to 9 years  
  cin2plus.cens <- ifelse(cin2plus.days > 365.25*9, 0, cin2plus)
  cin2plus.days.cens <- ifelse(cin2plus.days > 365.25*9, 365.25*9, cin2plus.days)
  
  cin2plus.1st.round.cens <- ifelse(cin2plus.days > 365.25*3, 0, cin2plus)
  cin2plus.1st.round.cens.days <- ifelse(cin2plus.days > 365.25*3, 365.25*3, cin2plus.days)
  
  # repeat the process above for CIN3+ 
  cin3plus <- ifelse(worst.histo.cyt >=8 & worst.histo.cyt!=13 & worst.histo!=11, 1, 0)
  cin3plus.days <- difftime(worst.histo.cyt.date, slovenia.dat$Date_visit_1, units="days")
  cin3plus.days[cin3plus.days<=0] <- 0
  # add censoring
  cin3plus.cens <- ifelse(cin3plus.days > 365.25*9, 0, cin3plus)
  cin3plus.days.cens <- ifelse(cin3plus.days > 365.25*9, 365.25*9, cin3plus.days)
  
  cin3plus.1st.round.cens <- ifelse(cin3plus.days > 365.25*3, 0, cin3plus)
  cin3plus.1st.round.cens.days <- ifelse(cin3plus.days > 365.25*3, 365.25*3, cin3plus.days)
  
  # create indicator variables for specific baseline genotypes 
  hpv.16 <- str_detect(slovenia.dat$Genotype_1, "16")
  hpv.18 <- str_detect(slovenia.dat$Genotype_1, "18")
  hpv.31 <- str_detect(slovenia.dat$Genotype_1, "31")
  
  genotype <- case_when(
    str_detect(slovenia.dat$Genotype_1, "16") ~ "HPV16+",
    str_detect(slovenia.dat$Genotype_1, "18|45")~ "HPV18/45+",
    str_detect(slovenia.dat$Genotype_1, "31|33|35|52|58") ~ "HPV31/33/35/52/58+",
    str_detect(slovenia.dat$Genotype_1, "39|51|56|59") ~ "HPV39/51/56/59+",
    TRUE ~ "HPV-"
  )
  
  hpv.2nd.dates <- unlist(apply(hpvs.date, 1, function(x) { x <- x[!is.na(x)] ; diff <- difftime(x, x[1], units='days') ; min(which(diff >= 3*365.25 & diff < 6*365.25))}))
  hpv.2nd.dates[hpv.2nd.dates > 3] <- Inf # there are only 3 Abbott visits so we're only interested in those
  abb.2nd.round <- hpvs[,c("Abb_RT_res_1", "Abb_RT_res_2", "Abb_RT_res_3")][cbind(seq_along(hpv.2nd.dates), hpv.2nd.dates)]
  
  baseline_cyt_pap <- case_when(
    str_detect(slovenia.dat$Cytology_result_1, "1") ~ "Normal",
    str_detect(slovenia.dat$Cytology_result_1, "2") ~ "ASCUS",
    str_detect(slovenia.dat$Cytology_result_1, "3") ~ "LSIL",
    str_detect(slovenia.dat$Cytology_result_1, "4") ~ "HSIL",
    TRUE ~ NA
  )
  
  # format Abbott baseline genotype because some 18 and neg were written as 18*/neg*
  abb_genotype1 = case_when(
    str_detect(slovenia.dat$Abb_RT_1, "18") ~ "18",
    str_detect(slovenia.dat$Abb_RT_1, "neg") ~ "neg",
    TRUE ~ slovenia.dat$Abb_RT_1,
  )
  # combine all new variables in a dataframe which is easier to use for analysis 
  dat.full <- data.frame(baseline= slovenia.dat$Date_visit_1, 
                         worst.histo, worst.histo.date,
                         worst.histo.cyt, worst.histo.cyt.date,
                         cin2plus, cin2plus.days, 
                         cin2plus.cens = cin2plus.cens, 
                         cin2plus.cens.factor = factor(cin2plus.cens, levels=c(1, 0)), 
                         cin2plus.days.cens,
                         cin2plus.1st.round.cens, 
                         cin2plus.1st.round.cens.days, 
                         cin3plus, cin3plus.days, 
                         cin3plus.cens = cin3plus.cens, 
                         cin3plus.cens.factor = factor(cin3plus.cens, levels=c(1, 0)), cin3plus.days.cens, 
                         cin3plus.1st.round.cens, cin3plus.1st.round.cens.days,
                         hpv.16, hpv.18, hpv.31, hrhpv = slovenia.dat$Gen_hr_res_1, 
                         baseline_cyt= factor(slovenia.dat$baseline_cyt, levels=c(1,0)),
                         hc2 = factor(slovenia.dat$HC2_res_1, levels=c(1,0)),
                         abb = factor(slovenia.dat$Abb_RT_res_1, levels=c(1,0)),
                         aly = factor(slovenia.dat$Aly_RT_res_1, levels=c(1,0)),
                         cob = factor(slovenia.dat$Cob_RT_res_1, levels=c(1,0)),
                         genotype = factor(genotype, levels=c("HPV16+", "HPV18/45+", "HPV31/33/35/52/58+", "HPV39/51/56/59+", "HPV-")),
                         abb.2nd.round = abb.2nd.round,
                         baseline_cyt_pap = factor(baseline_cyt_pap, levels=c("Normal", "ASCUS", "LSIL", "HSIL")),
                         abb_genotype1,
                         aly_genotype1 = slovenia.dat$Aly_RT_1,
                         cob_genotype1 = slovenia.dat$Cob_RT_1,
                         hc2_genotype1 = slovenia.dat$HC2_res_1,
                         age = slovenia.dat$Age) #as.numeric(difftime(slovenia.dat$Date_visit_1, slovenia.dat$Date_birth)/365.25))
  rownames(dat.full) <- slovenia.dat$ID
  return(dat.full)
}
dat.full <- create_data()

#-------------------
# Inclusion criteria:
#   -	valid cytology or HPV result at baseline
#   -	at least one repeat cytology or HPV result in follow-up
# Exclusion criteria:
#   -	CIN2+ diagnosed within 6 months of baseline enrollment
valid.baseline.hpvs <- apply(hpvs[,1:4],1, function(x) sum(is.na(x)) ==0 ) # all baseline hpv tests --> need matching data to compare sens/spec with p-values 
valid.baseline.cyt <- slovenia.dat$Cytology_result_1 !=99 & !is.na(dat.full$baseline_cyt) # from data dictionary: cytology=99 means not available (lost/destroyed)
cin2plus.diagnosis <- dat.full$cin2plus.cens ==1 & ! is.na( dat.full$cin2plus.cens)

# check that there is any repeat cytology or repeat HPV test after baseline by looking at the
# data frame of HPV test dates and checking if there is a non-NA date value recorded (which means there was a test performed)
repeat.hpv.test <- apply(hpvs.date, 1, function(x){ x <- x[-1] ; any(!is.na(x))})
repeat.cyt.test <- apply(cyto.date, 1, function(x){ x <- x[-1] ; any(!is.na(x))})

inclusion <- ifelse((valid.baseline.cyt | valid.baseline.hpvs) & (cin2plus.diagnosis | (repeat.cyt.test | repeat.hpv.test)), TRUE, FALSE)
exclusion <- ifelse(cin2plus.diagnosis & dat.full$cin2plus.days.cens < 365.25/12*6, 1, 0) # excluded if CIN2+ within 6 months of baseline enrolment
#-------------------
# subset data
dat.incl <- dat.full[inclusion & !exclusion,]
dat.incl.30 <- dat.full[inclusion & !exclusion & dat.full$age >= 30,]

#-------------------
# convert data to long format of different HPV tests for the survival analysis functions 
# --> so that test can be used as a strata variable 
create.long.data <- function(data, outcome=0){
  cols <- c("cin2plus.cens", "cin2plus.days.cens", "cin3plus.cens", "cin3plus.days.cens")
  cyt.neg <- data[data$baseline_cyt == outcome & !is.na(data$baseline_cyt), cols]
  cyt.neg$test <- rep("Cytology", dim(cyt.neg)[1])
  cyt.neg$id <- rownames(cyt.neg)
  
  hc2.neg <- data[data$hc2 == outcome & !is.na(data$hc2), cols]
  hc2.neg$test <- rep("HC2", dim(hc2.neg)[1])
  hc2.neg$id <- rownames(hc2.neg)
  
  abb.neg <- data[data$abb == outcome & !is.na(data$abb), cols]
  abb.neg$test <- rep("Abbott", dim(abb.neg)[1])
  abb.neg$id <- rownames(abb.neg)
  
  aly.neg <- data[data$aly == outcome & !is.na(data$aly), cols]
  aly.neg$test <- rep("Alinity", dim(aly.neg)[1])
  aly.neg$id <- rownames(aly.neg)
  
  cob.neg <- data[data$cob == outcome & !is.na(data$cob), cols]
  cob.neg$test <- rep("Cobas", dim(cob.neg)[1])
  cob.neg$id <- rownames(cob.neg)
  
  new_dat <- rbind(cyt.neg, hc2.neg, abb.neg, aly.neg, cob.neg)
  new_dat$test <- factor(new_dat$test, levels = c("Cytology", "HC2", "Abbott", "Cobas", "Alinity"))
  return(new_dat)
}
# women who were negative at baseline !!
dat.long <- create.long.data(dat.incl, 0)
dat.long.30 <- create.long.data(dat.incl.30, 0)
dat.long.20to30 <- create.long.data(dat.incl[dat.incl$age >=20 & dat.incl$age<30,], 0)
dat.long.30to59 <- create.long.data(dat.incl[dat.incl$age >=30 & dat.incl$age<60,], 0)
dat.long.60plus <- create.long.data(dat.incl[dat.incl$age >=60,], 0)

dat.long.pos <- create.long.data(dat.incl, 1)
dat.long.pos.30 <- create.long.data(dat.incl.30, 1)

rm(hpvs, hpvs.date, cyto, cyto.date, repeat.hpv.test, repeat.cyt.test, valid.baseline.cyt, valid.baseline.hpvs, cin2plus.diagnosis)


#-------------------
# convert data to long format of different HPV tests for the survival analysis functions 
# --> so that test can be used as a strata variable 
create.cotest.data <- function(data, test){
  cols <- c("cin2plus.cens", "cin2plus.days.cens", "cin3plus.cens", "cin3plus.days.cens")
  cyt.neg <- data[data$baseline_cyt == 0 & !is.na(data$baseline_cyt), cols]
  cyt.neg$test <- rep("Cytology", dim(cyt.neg)[1])
  cyt.neg$id <- rownames(cyt.neg)
  cyt.neg$label <- "Cytology normal"
  
  cotest.neg <- data[data[[test]] == 0 & !is.na(data[[test]]) & 
                       data$baseline_cyt == 0 & !is.na(data$baseline_cyt), cols]
  cotest.neg$test <- rep("Co-test", dim(cotest.neg)[1])
  cotest.neg$id <- rownames(cotest.neg)
  cotest.neg$label <- "Co\255 test negative"
  
  hpv.neg <- data[data[[test]] == 0 & !is.na(data[[test]]), cols]
  hpv.neg$test <- rep(test, dim(hpv.neg)[1])
  hpv.neg$id <- rownames(hpv.neg)
  hpv.neg$label <- case_when(
    test == "hc2" ~ "hc2 negative",
    test == "abb" ~ "RealTi*m*e negative",
    test == "aly" ~ "Alinity negative",
    test == "cob" ~ "cobas 4800 negative"
  )
  
  new_dat <- rbind(cyt.neg, hpv.neg, cotest.neg)
  new_dat$test <- factor(new_dat$test, levels = c("Cytology", test, "Co-test"))
  new_dat$label <- factor(new_dat$label, levels = c("Cytology normal", hpv.neg$label[1], "Co-test negative"))
  return(new_dat)
}

hc2.cotest <- create.cotest.data(dat.incl.30, "hc2")
abb.cotest <- create.cotest.data(dat.incl.30, "abb")
aly.cotest <- create.cotest.data(dat.incl.30, "aly")
cob.cotest <- create.cotest.data(dat.incl.30, "cob")
