# code for producing plots for the manuscript - first have to run "analysis.R"
library(ggplot2)
library(survminer)
library(gdata)
library(grid)
library(ggtext)# remotes::install_github("clauswilke/ggtext")
library(dplyr)

null_ggsurvplot <- ggsurvplot(survfit(Surv(0, 0) ~ 1), data=data.frame(time=0, surv=0), risk.table = FALSE, legend = "none", conf.int = FALSE, color = 'white')
null_ggsurvplot$plot <- null_ggsurvplot$plot + theme_void() +theme(legend.position='none')

# -------------------
# (1) HPV test comparison - safety after negative HPV/ cytology test
# Figure 1: Cumulative risks of cervical intraepithelial neoplasia grade 2 or worse (CIN2+) and grade 3 or worse (CIN3+) by normal baseline cytology and negative baseline 
# HPV result according to HPV assay used among women ≥ 30 years old (A, B) and in total study population (C, D). Risk bands are 95% confidence intervals. 
ggsurvplot.risk <- function(fit, ylab, safety, tag, agegroup){
  newplot <- ggsurvplot(fit,  fun="event", break.x.by=3, surv.scale = "percent", risk.table ="nrisk_cumevents",
                         legend.labs=c("Cytology", "hc2", "RealTi*m*e", "cobas 4800", "Alinity"),xlim=c(0,9.1),
                         legend.title="", ylab=ylab, xlab='Time (years)', tables.theme = theme_cleantable(),
                         conf.int=T, ylim=c(0,0.032), font.x = 20, font.y=20, font.tickslab=20, risk.table.fontsize=5)
  
  newplot$plot <- newplot$plot +  geom_vline(xintercept = c(3, safety), lty=2, col='grey') +
    geom_hline(yintercept = 1-summary(fit, times=3)$surv[1], lty=2, col='grey') + 
    scale_x_continuous(breaks=c(0, 3, 6, round(safety, 3), 9), labels = c(0, 3, 6, round(safety, 3), 9)) +
    theme(legend.position = c(0.1, 0.9), text=element_text(size=25), plot.tag = element_text(face='bold'), 
          legend.text = element_markdown(size = 14), plot.margin = margin(5.5, 10, 5.5, 5.5, "points") ) +
    labs(title=agegroup, tag=tag)
  #newplot$table <- newplot$table + theme( plot.margin = margin(5.5, 50, 5.5, 5.5, "points") ) +
  return(newplot)
}

cin2plot <- ggsurvplot.risk(cin2plus.fit, "Cumulative risk of CIN2+", cin2plus.fit.safety, "C", "CIN2+ total study population")
cin3plot <- ggsurvplot.risk(cin3plus.fit, "Cumulative risk of CIN3+", cin3plus.fit.safety, "D", "CIN3+ total study population")
cin2plot.30 <- ggsurvplot.risk(cin2plus.fit30, "Cumulative risk of CIN2+", cin2plus.fit30.safety, "A", expression("CIN2+" >="30 years old"))
cin3plot.30 <-   ggsurvplot.risk(cin3plus.fit30, "Cumulative risk of CIN3+", cin2plus.fit30.safety, "B", expression("CIN3+" >="30 years old"))

plot.list <- list(cin2plot.30, cin2plot, cin3plot.30, cin3plot)
res <- arrange_ggsurvplots(plot.list, print = F,ncol = 2, nrow = 2, risk.table.height = 0.25) 
ggsave("Plots/Figure_1.pdf", res, width=20, height=20, device=cairo_pdf)

# -------------------
# (2) co-test results
# Figure 2: Cumulative risks of cervical intraepithelial neoplasia grade 2 or worse (CIN2+) by normal baseline cytology (red), negative baseline HPV result (green), 
# and co-testing (both normal baseline cytology and negative baseline HPV result; blue) according to HPV assay in women ≥ 30 years old.
ggsurvplot.cotest <- function(fit, ylab, tag, legendlabel){
  rightmargin <- 5.5
  if (tag =="B" | tag =='D') {rightmargin <- 100}
  newplot <- ggsurvplot(fit, ylab=ylab,break.x.by=3, legend.labs = levels(legendlabel),legend.title="",
                        fun = 'event', conf.int=F, xlab="Time (years)", xlim=c(0,9.2),ylim=c(0, 0.02), risk.table ="nrisk_cumevents", 
                        tables.theme = theme_cleantable(), surv.scale='percent',  font.x = 20, font.y=20, font.tickslab=20, 
                        risk.table.fontsize = 5, size=1, censor.size=2)
  newplot$plot <- newplot$plot + labs(title="", tag=tag) + 
    theme(legend.position = c(0.15, 0.9), text=element_text(size=25), plot.tag = element_text(face='bold'), 
          legend.text = element_markdown(size = 14), plot.margin = margin(5.5, 10, 5.5, 5.5, "points"))
  return(newplot)
}

hc2.cotest.plot <- ggsurvplot.cotest(hc2.cotest.fit, "Cumulative risk of CIN2+", "A", hc2.cotest$label)
abb.cotest.plot <- ggsurvplot.cotest(abb.cotest.fit, "Cumulative risk of CIN2+", "B", abb.cotest$label)
aly.cotest.plot <- ggsurvplot.cotest(aly.cotest.fit, "Cumulative risk of CIN2+", "D", aly.cotest$label)
cob.cotest.plot <- ggsurvplot.cotest(cob.cotest.fit, "Cumulative risk of CIN2+", "C", cob.cotest$label)

cotest.plots.list <- list(hc2.cotest.plot, cob.cotest.plot, abb.cotest.plot,  aly.cotest.plot)
res2 <- arrange_ggsurvplots(cotest.plots.list, print = F,ncol = 2, nrow = 2, risk.table.height = 0.15) 
ggsave("Plots/Figure_2 (No CI).pdf", res2, width=25, height=25, dpi=500, device=cairo_pdf)

# -------------------
# (3) Genotype analysis - Extended genotyping can offer better risk stratification
# Figure 3: Cumulative risks of cervical intraepithelial neoplasia grade 2 or worse (CIN2+) and grade 3 or worse (CIN3+) by HPV genotype during 9 years of 
# follow-up among women ≥ 30 years old (A, B) and in total study population (C, D).
geontype.ggsurvplot <- function(fit, ytitle, tag, agegroup, ci=T){
  newplot <- ggsurvplot(fit, ylab=ytitle,break.x.by=3, legend.labs = levels(dat.incl$genotype),legend.title="",
                        fun = 'event', conf.int=ci, xlab="Time (years)", xlim=c(0,9.2),ylim=c(0, 0.5), risk.table ="nrisk_cumevents", tables.theme = theme_cleantable(),
                        surv.scale='percent',  font.x = 15, font.y=15, font.tickslab=15, risk.table.fontsize=5 )
  newplot$plot <- newplot$plot + labs(title=agegroup, tag=tag) + 
    theme(legend.position = c(0.16, 0.9), text=element_text(size=25),plot.tag = element_text(face='bold'),
          plot.margin = unit(c(5.5, 10, 5.5, 5.5), "points"), legend.text = element_markdown(size = 14))
  #newplot$table <- newplot$table + theme( plot.margin = unit(c(5.5, 50, 5.5, 5.5), "points") ) +
  return(newplot)
}
cin2plus.genotype.plot30 <- geontype.ggsurvplot(cin2plus.genotype.fit30, "Cumulative risk of CIN2+", "A", expression("CIN2+" >="30 years old"), ci=F)
cin2plus.genotype.plot <- geontype.ggsurvplot(cin2plus.genotype.fit, "Cumulative risk of CIN2+", "C", "CIN2+ total study population", ci=F)
cin3plus.genotype.plot30 <- geontype.ggsurvplot(cin3plus.genotype.fit30, "Cumulative risk of CIN3+", "B", expression("CIN3+" >="30 years old"), ci=F)
cin3plus.genotype.plot <- geontype.ggsurvplot(cin3plus.genotype.fit, "Cumulative risk of CIN3+", "D", "CIN3+ total study population", ci=F)

genotype.plots.list <- list(cin2plus.genotype.plot30, cin2plus.genotype.plot, cin3plus.genotype.plot30, cin3plus.genotype.plot)
res3 <- arrange_ggsurvplots(genotype.plots.list, print = F,ncol = 2, nrow = 2, risk.table.height = 0.2) 
ggsave("Plots/Figure_3 (no CI).pdf", res3, width=20, height=20,  device=cairo_pdf)


# -------------------
# (4) Supplementary Figures: Venn diagram of number of HPV pos test results
source("ggvennXY.R")

venn.func <- function(data){
  return(list("hc2"=data$id[data$hc2==1 & data$n.pos != '0 tests positive'], 
              "Alinity"=data$id[data$aly==1 & data$n.pos != '0 tests positive'], 
              "cobas 4800"=data$id[data$cob==1 & data$n.pos != '0 tests positive'], 
              "Realti*m*e"= data$id[data$abb==1 & data$n.pos != '0 tests positive']))
}
hpv.tests.list <- venn.func(hpv.tests)
cin2.tests.list <- venn.func(hpv.tests[hpv.tests$cin2plus.cens == 1,])
hpv.tests.30list <- venn.func(hpv.tests[hpv.tests$age >= 30,])
cin2.tests.30list <- venn.func(hpv.tests[hpv.tests$age >= 30 & hpv.tests$cin2plus.cens == 1,])

fills <- c("#999999", "#E69F00", "#56B4E9", "#009E73")

total.pop.venn <- my.ggvenn(hpv.tests.list, cin2.tests.list, show_percentage_X = T, fill_color = fills, stroke_size = 0.25, text_size = 4)
older30.venn <- my.ggvenn(hpv.tests.30list, cin2.tests.30list, show_percentage_X = T, fill_color = fills, stroke_size = 0.25, text_size = 4)
ggsave(total.pop.venn, file="Plots/Figure_S1.png", width=7, height=5)
ggsave(older30.venn, file="Plots/Figure_S3.png", width=7, height=5)

PCR.total.pop.venn <- my.ggvenn(hpv.tests.list[2:4], cin2.tests.list[2:4], show_percentage_X = T, text_size = 3, fill_color = fills[2:4], stroke_size = 0.25)
PCR.older30.venn <- my.ggvenn(hpv.tests.30list[2:4], cin2.tests.30list[2:4], show_percentage_X = T, text_size = 3, fill_color = fills[2:4], stroke_size = 0.25)
ggsave(PCR.total.pop.venn, file="Plots/Figure_S2.png", width=7, height=5)
ggsave(PCR.older30.venn, file="Plots/Figure_S4.png", width=7, height=5)
