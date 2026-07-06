# Libraries
## Data Manip
suppressMessages(library(tidyverse))
## Modeling
suppressMessages(library(brms))
suppressMessages(library(cmdstanr))
set_cmdstan_path("C:/Users/colem/.cmdstan/cmdstan-2.39.0")
## Model Diagnostics
suppressMessages(library(DHARMa))
## Data Vis
suppressMessages(library(gt))
suppressMessages(library(gridGraphics))
suppressMessages(library(cowplot))
suppressMessages(library(bayesplot))
suppressMessages(library(ggbreak))

Sys.setenv(CHROMOTE_CHROME = "C://Program Files (x86)/Microsoft/Edge/Application/msedge.exe")

"%nin%"<- Negate("%in%")

# Data
data<- read.csv("Data/FinalData_12-11-25_Miller/TotalSystem_4-24-2026_Miller.csv") %>%
  filter(Tank %nin% c(11,14)) %>% 
  ungroup() %>% 
  group_by(Tank, Position, Treatment)

## Data with t-1 time lag
dataLag <- data %>%
  ungroup() %>% 
  arrange(Tank,Position, time) %>% 
  group_by(Position) %>%
  mutate(TempLag = lag(Temp, n =1)) %>%
  mutate(ChlLag = lag(ChlAConcentration, n = 1)) %>%
  mutate(logChl = log(ChlAConcentration)) %>%
  mutate(logChlLag = lag(logChl, n =1)) %>%
  mutate(TadLag = lag(TotalTadMass, n = 1)) %>%
  mutate(AnuraLag = lag(TotalAnuranMass, n = 1)) %>% 
  mutate(CopeLag = lag(CopepodPerLiter, n = 1)) %>%
  mutate(DiploLag = lag(DiploPerLiter, n = 1)) %>%
  mutate(OstraLag = lag(OstracodPerLiter, n = 1)) %>%
  mutate(CopeLag10 = lag(CopepodPerLiter/100, n = 1)) %>%
  mutate(DiploLag10 = lag(DiploPerLiter/100, n = 1)) %>%
  mutate(OstraLag10 = lag(OstracodPerLiter/100, n = 1)) %>%
  mutate(ZoopLag = lag(TotalCrust, n = 1)) %>%
  mutate(tadNumLag = lag(Density, n = 1)) %>%
  mutate(DiploPres = if_else(DiploPerLiter == 0, 0, 1)) %>% 
  mutate(Time = as.factor(time))  %>%
  mutate(Density.at.t.1 = lag(Density, n = 1)) %>%
  mutate(Treatment = as.numeric(Treatment)) %>%
  filter(time != 1) %>%
  ungroup() %>%
  group_by(Time)%>%
  mutate(TempLag2 = scale(TempLag)) %>% ungroup()
#mutate(Treatment = factor(Treatment, levels=c('0.4', '3', '5'))) %>%


## Global Prior Specification
priors<- c(
  prior(normal(0,1), class = "b"),
  prior(normal(0,10), class = "Intercept"),
  prior(exponential(1), class = "sd")
)

#### Data
zoopData<- dataLag %>% 
  ungroup() %>%  
  group_by(Tank, Time, Treatment) %>% 
  summarise(across(c(Diplostraca, Ostracoda, Copepoda, Collected.Volume..mL.),
                   sum)) %>% mutate(Volume = Collected.Volume..mL./1000)
predictorData<- dataLag %>%
  ungroup() %>%
  group_by(Tank, Time, Treatment) %>%
  summarise_at(vars(DiploPres, TadLag, logChlLag, ChlLag, CopeLag10, DiploLag10, TempLag, 
                    pH), mean, .groups="keep")

zoopData<- left_join(zoopData, predictorData, join_by(Tank, Time ,Treatment))
rm(predictorData)



# Modeling
#---------------------------------------------------------------------------
### Algae/Chlorophyll Modeling
### Formulas
chlModList<- list(
  bf(ChlAConcentration ~ (1|Tank) + (1|Position) +
       Treatment + CopeLag10 + DiploLag10 + TadLag + ChlLag),
  bf(ChlAConcentration ~ (1|Tank) + (1|Position) +
       Treatment* CopeLag10 + DiploLag10 + TadLag + ChlLag),
  bf(ChlAConcentration ~ (1|Tank) + (1|Position) +
       Treatment* DiploLag10 + CopeLag10 + TadLag + ChlLag), 
  bf(ChlAConcentration ~ (1|Tank) + (1|Position) +
       Treatment* TadLag + DiploLag10 + CopeLag10 + ChlLag)
)

#Processing
chlFitList<- list()
chlSumList<- list()
chlLooList<- list()
chlParams <- list()


for(i in seq_along(chlModList)){
  chlFitList[[i]] <- brm(chlModList[[i]],
                         data    = dataLag,
                         family  = lognormal(),
                         chains  = 2,
                         cores   = 2,
                         iter    = 6000,
                         prior   = priors,
                         file    = paste("Code/Final Code/Saved Models/savedChl/m_", i, sep = ""),
                         backend = 'cmdstanr'
  )
  chlSumList[[i]]<- summary(chlFitList[[i]])
  chlLooList[[i]]<- loo(chlFitList[[i]], reloo = T)
  chlParams[[i]] <- nrow(fixef(chlFitList[[i]]))
  names(chlFitList)[i] <- paste("g", i, sep="")
  names(chlLooList)[i] <- paste("g", i, sep="")
  names(chlParams)[i]  <- paste("g", i, sep="")
  names(chlModList)[i]  <- paste("g", i, sep="")
}

chlParams<- as.data.frame(t(as.data.frame(chlParams))) %>% 
  rename_at(1, ~"Params")
chlLoo<- loo_compare(chlLooList) %>%
  as.data.frame() %>%
  round(.,2)
chlLoo<- merge(chlLoo, chlParams, by=0) %>% arrange(desc(elpd_diff))

bestChlMod<- chlLoo %>% filter(elpd_diff + 1.96*se_diff >= 0) %>%
  arrange(Params) %>% .[1,1] # Chl Mod 2 was used due to poor fit metrics 

chlModSum<- summary(chlFitList[["g2"]], robust = T)
chlModSum<- chlModSum$fixed %>% round(., 2) %>%
  mutate(Converge = if_else(Rhat <= 1.05, TRUE, FALSE)) %>%
  select(!c(Bulk_ESS, Tail_ESS, Rhat)) %>%
  mutate(Meaningful = if_else(`l-95% CI` <= 0 & `u-95% CI` >= 0, FALSE, TRUE)) %>%
  rownames_to_column("Term")

chlPPFit<- pp_check(chlFitList[["g2"]], type = 'scatter_avg', stat="mean") 
chlPPDist<- pp_check(chlFitList[["g2"]], ndraws = 100) + scale_x_continuous(limits = c(0,45)) +
  xlab("Density") 
plot_grid(chlPPFit, chlPPDist, labels = "AUTO", align = "h", scale = 0.9)

chlPPFit2<- pp_check(chlFitList[["g3"]], type = 'scatter_avg', stat="mean") 
chlPPDist2<- pp_check(chlFitList[["g3"]], ndraws = 100) + scale_x_continuous(limits = c(0,45)) +
  xlab("Density") 
plot_grid(chlPPFit2, chlPPDist2, labels = "AUTO", align = "h", scale = 0.9)
#---------------------------------------------------------------------------
### Tadpole Mass Modeling
### Formulas
tadModList<- list(
  bf(TotalTadMass ~ (1|Tank) + (1|Position) + TadLag + 
       Treatment + CopeLag10 + DiploLag10 + ChlLag,
     hu ~ tadNumLag + Treatment + CopeLag10 + DiploLag10 + ChlLag + (1|Tank) + (1|Position)),
  
  bf(TotalTadMass ~ (1|Tank) + (1|Position) + TadLag +
       Treatment*CopeLag10 + DiploLag10 + ChlLag,
     hu ~ tadNumLag + Treatment + CopeLag10 + DiploLag10 + ChlLag + (1|Tank) + (1|Position)),
  
  bf(TotalTadMass ~ (1|Tank) + (1|Position) + TadLag +
       Treatment*DiploLag10 + CopeLag10+ ChlLag,
     hu ~ tadNumLag + Treatment + CopeLag10 + DiploLag10 + ChlLag + (1|Tank) + (1|Position)),
  
  bf(TotalTadMass  ~ (1|Tank) + (1|Position) + TadLag +
       Treatment*logChlLag + DiploLag10 + CopeLag10,
     hu ~ tadNumLag + Treatment + CopeLag10 + DiploLag10 + ChlLag +(1|Tank) + (1|Position))
)


#Processeing
tadFitList<- list()
tadLooList<- list()
tadParams<- list()

for(i in seq_along(tadModList)){
  tadFitList[[i]] <- brm(tadModList[[i]],
                         data    = dataLag ,
                         family  = hurdle_lognormal(),
                         chains  = 2,
                         cores   = 2,
                         iter    = 6000,
                         prior    = priors,
                         save_pars = save_pars(all = TRUE),
                         file = paste("Code/Final Code/Saved Models/savedTad/m_", i, sep = ""),
                         backend = 'cmdstanr'
  )
  tadLooList[[i]]<- loo(tadFitList[[i]], reloo=T)
  tadParams[[i]] <- nrow(fixef(tadFitList[[i]]))
  names(tadFitList)[i] <- paste("g", i, sep="")
  names(tadLooList)[i] <- paste("g", i, sep="")
  names(tadParams)[i]  <- paste("g", i, sep="")
  names(tadModList)[i] <- paste("g", i, sep="")
}

tadParams<- as.data.frame(t(as.data.frame(tadParams))) %>% rename_at(1, ~"Params")
tadLoo<- loo_compare(tadLooList) %>%
  as.data.frame() %>%
  round(.,2)
tadLoo<- merge(tadLoo, tadParams, by=0) %>% arrange(desc(elpd_diff))

bestTadMod<- tadLoo %>% filter(elpd_diff + 1.96*se_diff >= 0) %>%
  arrange(Params) %>% .[1,1]

tadModSum<- summary(tadFitList[[bestTadMod]], robust = T)
tadModSum<- tadModSum$fixed %>% round(., 2) %>%
  mutate(Converge = if_else(Rhat <= 1.05, TRUE, FALSE)) %>%
  select(!c(Bulk_ESS, Tail_ESS, Rhat)) %>%
  mutate(Meaningful = if_else(`l-95% CI` <= 0 & `u-95% CI` >= 0, FALSE, TRUE)) %>%
  rownames_to_column("Term")


tadPPFit<- pp_check(tadFitList[[bestTadMod]], type = 'scatter_avg')
tadPPDist<-pp_check(tadFitList[[bestTadMod]], ndraws = 100) + scale_x_continuous(limits = c(0,20))  +
  xlab("Density") 
plot_grid(tadPPFit, tadPPDist, labels = "AUTO", align = "h", scale = 0.9)

#---------------------------------------------------------------------------
### Copepod Abundance Modeling
### Formulas
copeModList<- list(
  bf(Copepoda ~  (1|Tank) + offset(log(Volume)) +
       Treatment + TadLag + DiploLag10 + ChlLag+ CopeLag10),
  bf(Copepoda ~  (1|Tank) + offset(log(Volume)) +
       Treatment*TadLag + DiploLag10 + ChlLag + CopeLag10),
  bf(Copepoda ~  (1|Tank) + offset(log(Volume)) +
       Treatment*DiploLag10 + TadLag + ChlLag + CopeLag10),
  bf(Copepoda ~  (1|Tank) + offset(log(Volume)) +
       Treatment*ChlLag + DiploLag10 + TadLag + CopeLag10)
)
#Processeing
copeFitList<- list()
copeLooList<- list()
copeParams<- list()

for(i in seq_along(copeModList)){
  copeFitList[[i]] <- brm(copeModList[[i]],
                          data    = zoopData,
                          family  = negbinomial(),
                          chains  = 2,
                          cores   = 2,
                          iter    = 6000,
                          prior   = priors,
                          file    = paste("Code/Final Code/Saved Models/savedCope/m_", i, sep = ""),
                          save_pars = save_pars(all = TRUE),
                          backend = 'cmdstanr'
  )
  copeLooList[[i]]<- loo(copeFitList[[i]], reloo = T)
  copeParams[[i]] <- nrow(fixef(copeFitList[[i]]))
  names(copeFitList)[i] <- paste("g", i, sep="")
  names(copeLooList)[i] <- paste("g", i, sep="")
  names(copeParams)[i]  <- paste("g", i, sep="")
  names(copeModList)[i] <- paste("g", i, sep="")
}

copeParams<- as.data.frame(t(as.data.frame(copeParams))) %>% rename_at(1, ~"Params")
copeLoo<- loo_compare(copeLooList) %>%
  as.data.frame() %>%
  round(.,2)
copeLoo<- merge(copeLoo, copeParams, by=0) %>% arrange(desc(elpd_diff))
bestCopeMod<- copeLoo %>% filter(elpd_diff + 1.96*se_diff >= 0) %>%
  arrange(Params) %>% .[1,1]


copeModSum<- summary(copeFitList[[bestCopeMod]], robust = T)
copeModSum<- copeModSum$fixed %>% round(., 2) %>%
  mutate(Converge = if_else(Rhat <= 1.05, TRUE, FALSE)) %>%
  select(!c(Bulk_ESS, Tail_ESS, Rhat)) %>%
  mutate(Meaningful = if_else(`l-95% CI` <= 0 & `u-95% CI` >= 0, FALSE, TRUE)) %>%
  rownames_to_column("Term")

copePPFit<- pp_check(copeFitList[[bestCopeMod]], type = 'scatter_avg')
copePPDist<-pp_check(copeFitList[[bestCopeMod]], ndraws = 100) + scale_x_continuous(limits = c(0,3000))  +
  xlab("Density") 
plot_grid(copePPFit, copePPDist, labels = "AUTO", align = "h", scale = 0.9)


#---------------------------------------------------------------------------
### Diplostracan Abundance Modeling
diploModList<- list(
  bf(Diplostraca ~  (1|Tank) + offset(log(Volume)) +
       Treatment + TadLag + CopeLag10 + ChlLag+ DiploLag10),
  bf(Diplostraca ~  (1|Tank) + offset(log(Volume)) +
       Treatment*TadLag + CopeLag10 + ChlLag + DiploLag10),
  bf(Diplostraca ~  (1|Tank) + offset(log(Volume)) +
       Treatment*CopeLag10 + TadLag + ChlLag + DiploLag10),
  bf(Diplostraca ~  (1|Tank) + offset(log(Volume)) +
       Treatment*ChlLag + CopeLag10 + TadLag + DiploLag10)
)

#Processeing
diploFitList<- list()
diploLooList<- list()
diploParams<- list()

for(i in seq_along(diploModList)){
  diploFitList[[i]] <- brm(diploModList[[i]],
                           data    = zoopData,
                           family  = negbinomial(),
                           chains  = 2,
                           cores   = 2,
                           iter    = 6000,
                           prior   = priors,
                           file    = paste("Code/Final Code/Saved Models/savedDiploAbun/m_", i, sep = ""),
                           save_pars = save_pars(all = TRUE),
                           backend = 'cmdstanr'
  )
  diploLooList[[i]]<- loo(diploFitList[[i]], reloo=T)
  diploParams[[i]] <- nrow(fixef(diploFitList[[i]]))
  names(diploFitList)[i] <- paste("g", i, sep="")
  names(diploLooList)[i] <- paste("g", i, sep="")
  names(diploParams)[i]  <- paste("g", i, sep="")
  names(diploModList)[i] <- paste("g", i, sep="")
}

diploParams<- as.data.frame(t(as.data.frame(diploParams))) %>% rename_at(1, ~"Params")
diploLoo<- loo_compare(diploLooList) %>%
  as.data.frame() %>%
  round(.,2)
diploLoo<- merge(diploLoo, diploParams, by=0) %>% arrange(desc(elpd_diff))

bestDiploMod<- diploLoo %>% filter(elpd_diff + 1.96*se_diff >= 0) %>%
  arrange(Params) %>% .[1,1]

diploModSum<- summary(diploFitList[[bestDiploMod]], robust = T)
diploModSum<- diploModSum$fixed %>% round(., 2) %>%
  mutate(Converge = if_else(Rhat <= 1.05, TRUE, FALSE)) %>%
  select(!c(Bulk_ESS, Tail_ESS, Rhat)) %>%
  mutate(Meaningful = if_else(`l-95% CI` <= 0 & `u-95% CI` >= 0, FALSE, TRUE)) %>%
  rownames_to_column("Term")

diploPPFit <-pp_check(diploFitList[[bestDiploMod]], type = 'scatter_avg')
diploPPDist<-pp_check(diploFitList[[bestDiploMod]], ndraws = 100) + scale_x_continuous(limits = c(0,750))  +
  xlab("Density") 
plot_grid(diploPPFit, diploPPDist, labels = "AUTO", align = "h", scale = 0.9)


#-----Figs and Tables------------------------------------------------------------------
# Tables
## chl model summary
gt(chlModSum, rowname_col = "Term") %>% 
  tab_stubhead(label = "Term") %>%
  cols_move(columns  = Estimate,  after = `l-95% CI`) %>%
  cols_label(Estimate = "Median", `l-95% CI` = "2.5%", `u-95% CI`  = "97.5%") %>%
  tab_style(style = list(cell_text(weight = "bold")), locations = cells_body(rows = Meaningful == TRUE)) %>%
  tab_style(style = list(cell_text(weight = "bold")), locations = cells_stub(Meaningful == TRUE))  %>%
  cols_hide(columns  = c(Est.Error, Converge, Meaningful)) %>%
  gtsave("Figures/ThesisFigures/ChlModSum.docx")

## tad model summary
gt(tadModSum, rowname_col = "Term") %>% 
  tab_stubhead(label = "Term") %>%
  cols_move(columns  = Estimate,  after = `l-95% CI`) %>%
  cols_label(Estimate = "Median", `l-95% CI` = "2.5%", `u-95% CI`  = "97.5%") %>%
  tab_style(style = list(cell_text(weight = "bold")), locations = cells_body(rows = Meaningful == TRUE)) %>%
  tab_style(style = list(cell_text(weight = "bold")), locations = cells_stub(Meaningful == TRUE))  %>%
  cols_hide(columns  = c(Est.Error, Converge, Meaningful)) %>%
  gtsave("Figures/ThesisFigures/TadModSum.docx")

# cope model summary
gt(copeModSum, rowname_col = "Term") %>% 
  tab_stubhead(label = "Term") %>%
  cols_move(columns  = Estimate,  after = `l-95% CI`) %>%
  cols_label(Estimate = "Median", `l-95% CI` = "2.5%", `u-95% CI`  = "97.5%") %>%
  tab_style(style = list(cell_text(weight = "bold")), locations = cells_body(rows = Meaningful == TRUE)) %>%
  tab_style(style = list(cell_text(weight = "bold")), locations = cells_stub(Meaningful == TRUE))  %>%
  cols_hide(columns  = c(Est.Error, Converge, Meaningful)) %>%
  gtsave("Figures/ThesisFigures/CopeModSum.docx")

# diplo model summary
gt(diploModSum, rowname_col = "Term") %>% 
  tab_stubhead(label = "Term") %>%
  cols_move(columns  = Estimate,  after = `l-95% CI`) %>%
  cols_label(Estimate = "Median", `l-95% CI` = "2.5%", `u-95% CI`  = "97.5%") %>%
  tab_style(style = list(cell_text(weight = "bold")), locations = cells_body(rows = Meaningful == TRUE)) %>%
  tab_style(style = list(cell_text(weight = "bold")), locations = cells_stub(Meaningful == TRUE))  %>%
  cols_hide(columns  = c(Est.Error, Converge, Meaningful)) %>%
  gtsave("Figures/ThesisFigures/diploModSum.docx")

#Figures
## Chlorophyll 
bestChlMod<- "g2"
chlPlotList<- list()
chlPlotList[[1]]<- plot(conditional_effects(chlFitList[[bestChlMod]], effects = "Treatment", resolution =1000, prob = 0.8, plot=F))[[1]] +
  theme_classic() + xlab("Salinity (ppt)") + ylab("Chlorophyll a Concentration (\u03bcg/L)") +
  geom_point(data = dataLag, mapping = aes(Treatment, ChlAConcentration), inherit.aes = FALSE,  size = 2, alpha = 0.3) + 
  geom_ribbon(alpha = .4) + geom_smooth(color = "black") +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI")) +
  coord_trans(y = "log")
chlPlotList[[2]]<- plot(conditional_effects(chlFitList[[bestChlMod]], effects = "TadLag", resolution =1000, prob = 0.8, plot=F))[[1]] +
  theme_classic() + xlab("Total Tadpole Biomass (g)") + ylab("Chlorophyll a Concentration (\u03bcg/L)") +
  geom_point(data = dataLag, mapping = aes(TadLag, ChlAConcentration), inherit.aes = FALSE,  size = 2, alpha = 0.3) + 
  geom_ribbon(alpha = .4) + geom_smooth(color = "black") +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI"))+
  coord_trans(y = "log")
chlPlotList[[3]]<- plot(conditional_effects(chlFitList[[bestChlMod]], effects = "DiploLag10", resolution =1000,
                                            prob = 0.8, int_conditions = list(Treatment = c(0.4,3,5)), plot=F))[[1]] +
  theme_classic() + xlab(expression(bold(paste("Diplostracan Concentration (No. hL"^"-1",")")))) + 
  ylab("Chlorophyll a Concentration (\u03bcg/L)") +
  geom_point(data = dataLag, mapping = aes(DiploLag10, ChlAConcentration),
             inherit.aes = FALSE,  size = 2, alpha = 0.3) + 
  geom_ribbon(alpha = .4) + geom_smooth(color = "black") +
  theme(axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI")) +
  coord_trans(y = "log")



## Interaction 
drawMatrix<- as_draws_df(chlFit) %>% select(b_CopeLag10, `b_Treatment:CopeLag10`) %>%
  mutate(`0.4` =(b_CopeLag10+(`b_Treatment:CopeLag10`*0.4))) %>%
  mutate(`3` =(b_CopeLag10+(`b_Treatment:CopeLag10`*3))) %>%
  mutate(`5` =(b_CopeLag10+(`b_Treatment:CopeLag10`*5))) %>%
  .[,-c(1,2)] %>% 
  pivot_longer(1:3, names_to = "Treatment", values_to = "value") %>% 
  as.data.frame() %>%
  group_by(Treatment) %>%
  mutate(q5 = quantile(value, probs= 0.025)) %>%
  mutate(q95 = quantile(value, probs= 0.975))

chlSlopesPlot<- ggplot(data= drawMatrix, aes(value, fill = Treatment)) + geom_histogram(color = "black") + facet_wrap(~Treatment, ncol=1) +
  geom_vline(xintercept = 0, lwd = 2, color = "firebrick") + 
  geom_vline(data = drawMatrix, aes(group = Treatment, xintercept = q5+.01), lwd = 2, linetype = "dashed") +
  geom_vline(data = drawMatrix, aes(group = Treatment, xintercept = q95), lwd = 2, linetype = "dashed") +
  ylab("Draws") + 
  xlab("Copepod Interaction Slope Estimate") + 
  labs(fill = "Salinity Treatment") +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12)) +
  scale_fill_hue(direction= -1)

chlPlotList[[4]]<- plot(conditional_effects(chlFitList[[bestChlMod]], effects = "CopeLag10:Treatment", resolution =1000,
                         prob = 0.8, int_conditions = list(Treatment = c(0.4,3,5)), plot=F))[[1]] +
  theme_classic() + xlab(expression(bold(paste("Copepod Concentration (No. hL"^"-1",")")))) + 
  ylab("Chlorophyll a Concentration (\u03bcg/L)") +
  geom_point(data = dataLag, mapping = aes(CopeLag10, ChlAConcentration, color = as.factor(Treatment)),
             inherit.aes = FALSE,  size = 2, alpha = 0.3) + 
  geom_ribbon(alpha = .4) +
  theme(axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI")) +
  coord_trans(y = "log")
plot_grid(
  plot_grid(chlPlotList[[1]], chlPlotList[[2]], chlPlotList[[3]], scale = 0.9, nrow =1, labels = "AUTO"),
  plot_grid(chlPlotList[[4]]+theme(legend.position = "none"), chlSlopesPlot, scale = 0.9, nrow = 1, labels = c("D", "E"),
            rel_widths = c(1, 1.2)),
  ncol = 1
  )


## Tadpoles
tadPlotListh<- list()
dataLag2<- dataLag %>% mutate(TadPres = if_else(TotalTadMass == 0, 1, 0))
tadPlotListh[[1]]<- print(plot(conditional_effects(tadFitList[[bestTadMod]], resolution =1000, dpar = "hu", effects = "tadNumLag",
                                             prob = 0.8, plot=F))[[1]] +
  theme_classic() + xlab("Tadpole Density (Num)") + ylab("Non-Zero Total Tadpole Biomass (g)")+
  geom_jitter(data = dataLag2, mapping = aes(tadNumLag, TadPres), inherit.aes = FALSE,  size = 2, alpha = 0.3,
              height = 0, width = 0.4) + 
  geom_ribbon(alpha = .4) + geom_line(color = "black", lwd = 1) +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI")))
tadPlotListh[[2]]<- print(plot(conditional_effects(tadFitList[[bestTadMod]], resolution =1000, dpar = "hu", effects = "Treatment",
                                             prob = 0.8, plot=F))[[1]] +
  theme_classic() + xlab("Salinity Treatment (ppt)") + ylab("Non-Zero Total Tadpole Biomass (g)")+
  geom_jitter(data = dataLag2, mapping = aes(Treatment, TadPres), inherit.aes = FALSE,  size = 2, alpha = 0.3,
              height = 0, width = 0.4) + 
  geom_ribbon(alpha = .4) + geom_line(color = "black", lwd = 1) +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI")) +
  scale_y_break(c(0.004, .999), ticklabels = 1, symbol = "slash", expand= T, space = 0.1))
tadPlotListh[[3]]<- print(plot(conditional_effects(tadFitList[[bestTadMod]], resolution =1000, dpar = "hu", effects = "CopeLag10",
                                             prob = 0.8, plot=F))[[1]] +
  theme_classic() + xlab(expression(bold(paste("Copepod Concentration (No. hL"^"-1",")")))) + ylab("Non-Zero Total Tadpole Biomass (g)")+
  geom_jitter(data = dataLag2, mapping = aes(CopeLag10, TadPres), inherit.aes = FALSE,  size = 2, alpha = 0.3,
              height = 0, width = 0.4) + 
  geom_ribbon(alpha = .4) + geom_line(color = "black", lwd = 1) +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI")) +
  scale_y_break(c(0.03, .999), ticklabels = 1, symbol = "slash", expand= T, space = 0.1))
tadPlotListh[[4]]<- print(plot(conditional_effects(tadFitList[[bestTadMod]], resolution =1000, dpar = "hu", effects = "DiploLag10",
                                             prob = 0.8, plot=F))[[1]] +
  theme_classic() + xlab(expression(bold(paste("Diplostracan Concentration (No. hL"^"-1",")")))) + ylab("Non-Zero Total Tadpole Biomass (g)") +
  geom_jitter(data = dataLag2, mapping = aes(DiploLag10, TadPres), inherit.aes = FALSE,  size = 2, alpha = 0.3,
              height = 0, width = 0) + 
  geom_ribbon(alpha = .4) + geom_line(color = "black", lwd = 1) +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI")) +
  scale_y_break(c(0.0075, .999), ticklabels = 1, symbol = "slash", expand= T, space = 0.1))
tadPlotListh[[5]]<- print(plot(conditional_effects(tadFitList[[bestTadMod]], resolution =1000, dpar = "hu", effects = "ChlLag",
                                             prob = 0.8, plot=F))[[1]] +
  theme_classic() + xlab("Chlorophyll a Concentration (\u03bcg/L)") + ylab("Non-Zero Total Tadpole Biomass (g)") +
  geom_jitter(data = dataLag2, mapping = aes(ChlLag, TadPres), inherit.aes = FALSE,  size = 2, alpha = 0.3,
              height = 0, width = 0) + 
  geom_ribbon(alpha = .4) + geom_line(color = "black", lwd = 1) +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI")) +
  scale_y_break(c(0.075, .999), ticklabels = 1, symbol = "slash", expand= T, space = 0.1))
hurdlePlots<- plot_grid(plotlist = tadPlotListh, labels = "AUTO", scale = 0.9)
hurdlePlots

tadPlotList<- list()
tadPlotList[[1]]<- plot(conditional_effects(tadFitList[[bestTadMod]], resolution =1000,  effects = "Treatment",
                                             prob = 0.8, plot=F))[[1]] +
  theme_classic() + xlab("Salinity Treatment (ppt)") + ylab("Total Tadpole Biomass (g)")+
  geom_point(data = dataLag2, mapping = aes(Treatment, TotalTadMass), inherit.aes = FALSE,  size = 2, alpha = 0.3) + 
  geom_ribbon(alpha = .4) + geom_line(color = "black", lwd = 1) +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI")) 
tadPlotList[[2]]<- plot(conditional_effects(tadFitList[[bestTadMod]], resolution =1000, effects = "CopeLag10",
                                             prob = 0.8, plot=F))[[1]] +
  theme_classic() + xlab(expression(bold(paste("Copepod Concentration (No. hL"^"-1",")")))) + ylab("Total Tadpole Biomass (g)")+
  geom_point(data = dataLag2, mapping = aes(CopeLag10, TotalTadMass), inherit.aes = FALSE,  size = 2, alpha = 0.3) + 
  geom_ribbon(alpha = .4) + geom_line(color = "black", lwd = 1) +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI")) 
tadPlotList[[3]]<- plot(conditional_effects(tadFitList[[bestTadMod]], effects = "DiploLag10", resolution =1000, prob = 0.8, 
                                            plot=F))[[1]] +
  theme_classic() + xlab(expression(bold(paste("Diplostracan Concentration (No. hL"^"-1",")")))) + ylab("Total Tadpole Biomass (g)")+
  geom_point(data = dataLag2, mapping = aes(DiploLag10, TotalTadMass), inherit.aes = FALSE,  size = 2, alpha = 0.3) + 
  geom_ribbon(alpha = .4) + geom_line(color = "black", lwd = 1) +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI"))
tadPlotList[[4]]<- plot(conditional_effects(tadFitList[[bestTadMod]], effects = "ChlLag", resolution =1000, prob = 0.8, 
                                            plot=F))[[1]] +
  theme_classic() + xlab("Chlorophyll a Concentration (\u03bcg/L)") + ylab("Total Tadpole Biomass (g)")+
  geom_point(data = dataLag2, mapping = aes(ChlLag, TotalTadMass), inherit.aes = FALSE,  size = 2, alpha = 0.3) + 
  geom_ribbon(alpha = .4) + geom_line(color = "black", lwd = 1) +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI"))
nonZeroPlots<- plot_grid(plotlist = tadPlotList, labels = c("AUTO"), scale = 0.9)
nonZeroPlots


## Copepods
copePlotList<- list()
copePlotList[[1]]<- plot(conditional_effects(copeFitList[[bestCopeMod]], effects = "Treatment", resolution =1000, prob = 0.8, plot=F))[[1]] +
  theme_classic() + xlab("Salinity (ppt)") + ylab("Copepod Abundance (No.)") +
  geom_point(data = zoopData, mapping = aes(Treatment, Copepoda), inherit.aes = FALSE,  size = 2, alpha = 0.3) + 
  geom_ribbon(alpha = .4) + geom_smooth(color = "black") +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI"))
copePlotList[[2]]<- plot(conditional_effects(copeFitList[[bestCopeMod]], effects = "TadLag", resolution =1000, prob = 0.8, plot=F))[[1]] +
  theme_classic() + xlab("Total Tadpole Biomass (g)") + ylab("Copepod Abundance (No.)") +
  geom_point(data = zoopData, mapping = aes(TadLag, Copepoda), inherit.aes = FALSE,  size = 2, alpha = 0.3) + 
  geom_ribbon(alpha = .4) + geom_smooth(color = "black") +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI"))
copePlotList[[3]]<- plot(conditional_effects(copeFitList[[bestCopeMod]], effects = "DiploLag10", resolution =1000, prob = 0.8, plot=F))[[1]] +
  theme_classic() + xlab(expression(bold(paste("Diplostracan Concentration (No. hL"^"-1",")")))) + ylab("Copepod Abundance (No.)") +
  geom_point(data = zoopData, mapping = aes(DiploLag10, Copepoda), inherit.aes = FALSE,  size = 2, alpha = 0.3) + 
  geom_ribbon(alpha = .4) + geom_smooth(color = "black") +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI"))
copePlotList[[4]]<- plot(conditional_effects(copeFitList[[bestCopeMod]], effects = "ChlLag", resolution =1000, prob = 0.8, plot=F))[[1]] +
  theme_classic() + xlab("Chlorophyll a Concentration (\u03bcg/L)") + ylab("Copepod Abundance (No.)") +
  geom_point(data = zoopData, mapping = aes(ChlLag, Copepoda), inherit.aes = FALSE,  size = 2, alpha = 0.3) + 
  geom_ribbon(alpha = .4) + geom_smooth(color = "black") +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI"))
plot_grid(plotlist = copePlotList, labels = "AUTO", scale = 0.9)

## DiploStracans
diploPlotList<- list()
diploPlotList[[1]]<- plot(conditional_effects(diploFitList[[bestDiploMod]], effects = "Treatment", resolution =1000, prob = 0.8, 
                                              plot=F))[[1]] +
  theme_classic() + xlab("Salinity (ppt)") + ylab("Diplostraca Abundance (No.)") +
  geom_jitter(data = zoopData, mapping = aes(Treatment, Diplostraca), inherit.aes = FALSE,  size = 2, alpha = 0.3, width = .2) + 
  geom_ribbon(alpha = .4) + geom_smooth(color = "black", lwd = 1) +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI"))
diploPlotList[[2]]<- plot(conditional_effects(diploFitList[[bestDiploMod]], effects = "TadLag", resolution =1000, prob = 0.8, 
                                              plot=F))[[1]] +
  theme_classic() + xlab("Total Tadpole Mass (g)") + ylab("Diplostraca Abundance (No.)") +
  geom_point(data = zoopData, mapping = aes(TadLag, Diplostraca), inherit.aes = FALSE,  size = 2, alpha = 0.3) + 
  geom_ribbon(alpha = .4) + geom_smooth(color = "black", lwd = 1) +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI"))
diploPlotList[[3]]<- plot(conditional_effects(diploFitList[[bestDiploMod]], effects = "CopeLag10", resolution =1000, prob = 0.8, plot=F))[[1]] +
  theme_classic() + xlab(expression(bold(paste("Copepod Concentration (No. hL"^"-1",")")))) + ylab("Diplostraca Abundance (No.)") +
  geom_point(data = zoopData, mapping = aes(DiploLag10, Diplostraca), inherit.aes = FALSE,  size = 2, alpha = 0.3) + 
  geom_ribbon(alpha = .4) + geom_smooth(color = "black") +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI"))
diploPlotList[[4]]<- plot(conditional_effects(diploFitList[[bestDiploMod]], effects = "ChlLag", resolution =1000, prob = 0.8, plot=F))[[1]] +
  theme_classic() + xlab("Chlorophyll a Concentration (\u03bcg/L)") + ylab("Diplostraca Abundance (No.)") +
  geom_point(data = zoopData, mapping = aes(ChlLag, Diplostraca), inherit.aes = FALSE,  size = 2, alpha = 0.3) + 
  geom_ribbon(alpha = .4) + geom_smooth(color = "black") +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12), text = element_text(family = "Segoe UI"))
plot_grid(plotlist = diploPlotList, labels = "AUTO", scale = 0.9)


# ---- Interaction Validation ---
suppressMessages(library(interactions))
chlFit<- chlFitList[[bestChlMod]]


plot_grid(
  plot_grid(chlPlotList[[1]],chlPlotList[[2]],chlPlotList[[3]], nrow = 1, labels = LETTERS[1:3], scale = .9),
  plot_grid(chlPlotList[[4]]+theme(legend.position = "none"), chlSlopesPlot, nrow = 1, labels = LETTERS[4:5], scale = .9), ncol = 1
)

