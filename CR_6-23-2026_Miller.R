## Restart R
.rs.restartR()

## Clear environment
rm(list=ls(all=TRUE)) 

# Seed for reproductibilbity
set.seed("1234")

# Libraries
## Data Manip and visualization 
suppressMessages(library(tidyverse))
suppressMessages(library(tidyverse))
suppressMessages(library(ggtext))
suppressMessages(library(cowplot))
suppressMessages(library(bayesplot))
suppressMessages(library(gt))

## Modeling
suppressMessages(library(rstan))
suppressMessages(library(cmdstanr))
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Data 
anura_person_period<- read.csv("Data/FinalData_12-11-25_Miller/anuraTTEFormatedData_6-12-2026_Miller.csv")

## Removing T1 
anura_person_period<- anura_person_period %>% 
  group_by(TankEn, ID) %>%
  mutate(tadNumLag = lag(Density, n = 1)) %>%
  ungroup() %>%
  filter(time != 1) %>%
  mutate(Event_num = ifelse(Event == "Live", 1, 
                            ifelse(Event == "Metamorph", 2, 3))) %>%
  group_by(Tank) %>%
  mutate(Tank = cur_group_id()) %>%
  ungroup() %>% 
  group_by(EnID) %>%
  mutate(EnID = cur_group_id()) %>%
  ungroup() %>% 
  group_by(time) %>%
  mutate(time = cur_group_id()) %>%
  ungroup() %>%
  mutate(CopeLag10 = CopeLag/100) %>%
  mutate(DiploLag10 = DiploLag/100) %>% 
  arrange(TankEn, time)

range(anura_person_period$CopeLag10)

#Stan Model
#crModel<- stan_model("Code/Final Code/Stan Code/Event Analysis/CompetingRisks_RandomEffectsandTime.stan")
crModelcmd<- cmdstan_model("Code/Final Code/Stan Code/Event Analysis/CompetingRisks_RandomEffectsandTime_cmd.stan")

#-----Competing Risks Integration-------------------------------------------------
# Matrix
generalMatrix<-model.matrix( ~ 0 + log(tadNumLag) + Treatment + logChlLag + 
                               CopeLag10 + DiploLag10, anura_person_period)
crData<- 
  list(
    T         = length(unique(anura_person_period$time)),
    N         = nrow(anura_person_period),
    P         = ncol(generalMatrix),
    Time      = anura_person_period$time,
    K         = length(unique(anura_person_period$Event_num)),
    Event     = anura_person_period$Event_num,
    X         = generalMatrix,
    G         = length(unique(anura_person_period$Tank)),
    tank      = as.numeric(anura_person_period$Tank),
    E         = length((unique(anura_person_period$EnID))),
    enclosure = as.numeric(anura_person_period$EnID)
  )

# Model Fit
crModelGeneral <- sampling(crModel, data = crData, chains = 1,
                           iter = 10000,
                           init = "0")
crModelGeneral@stanmodel@dso <- new("cxxdso")
saveRDS(crModelGeneral, file = "Code/Final Code/Saved Models/savedCR/GTF_General.rds")
crModelGeneral<- readRDS("Code/Final Code/Saved Models/savedCR/GTF_General.rds")

crModelGeneralcmd<- crModelcmd$sample(
  data = crData, chains = 1, iter_warmup=5000, iter_sampling = 5000, init = 0
)
crModelGeneralcmd$save_object("Code/Final Code/Saved Models/savedCR/GTF_Generalcmd.rds")

#crModelGeneral<- readRDS("Code/Final Code/Saved Models/savedCR/GTF_General.rds")
#crGeneralModSum<- summary(crModelGeneral, pars = c("alpha", "beta", "theta", "gamma"),
 #                         probs = c(0.025,.1,.9, 0.975))$summary

crGeneralModSum<- crModelGeneralcmd$summary(variables = c("alpha", "beta", "theta", "gamma"))
crGeneralModSumTab<- crGeneralModSum %>% as.data.frame() %>%  
  filter(variable == gsub("alpha","", variable)) %>%
  filter(variable == gsub("theta","", variable)) %>%
  mutate(coef = rep(colnames(generalMatrix), 2)) %>%
  mutate(id = rep(c("M","D"), each =5)) %>%
  mutate(Converge = if_else(rhat < 1.05, TRUE, FALSE)) %>% 
  mutate(Meaningful = if_else(q5 <= 0 & q95 >= 0, FALSE, TRUE)) %>% 
  select(!c(variable, mean, sd, mad, rhat, ess_bulk, ess_tail)) %>%
  mutate_at(c("median", "q5", "q95"), round, 2)
gt(crGeneralModSumTab, rowname_col = "coef") %>%
  tab_stubhead(label = "Term") %>%
  tab_row_group(label = "Metamorphosis", rows = id == "M") %>%
  tab_options(row_group.default_label = "Death") %>%
  cols_move(columns  = median,  after = q5) %>%
  cols_label(median = "Median", q5 = "2.5%", q95  = "97.5%") %>%
  tab_style(style = list(cell_text(weight = "bold")), locations = cells_body(rows = Meaningful == 1)) %>%
  tab_style(style = list(cell_text(weight = "bold")), locations = cells_stub(Meaningful == 1))  %>%
  cols_hide(columns =c(Meaningful, Converge, id)) %>%
  gtsave("Figures/ThesisFigures/crModSum.docx")

#-----Posterior Predictive Checking8-------------------------------------------------
# General PP checks
#M_rep<-extract(crModelGeneral, pars = "M_rep")$M_rep
#M_rep_mat<- as.matrix(M_rep)
M_rep_mat<- crModelGeneral$draws("M_rep", format = "matrix")
M_obs <- as.vector(crData$Event)
crPPDist<- bayesplot::ppc_dens_overlay(y = M_obs, yrep = M_rep_mat[1:100,]) +
  xlab("Density") + theme(axis.title.x = element_text(face="bold"))

M_rep_mat_meta<- M_rep_mat -1
M_rep_mat_meta[M_rep_mat_meta==2]<-0
prop_rep_meta<- rowMeans(M_rep_mat_meta)
prob_obs_meta<- length(which(crData$Event==2))/length(crData$Event)
crPPMeta<- ggplot() + 
  geom_histogram(aes(prop_rep_meta), color = "black") + 
  geom_vline(xintercept = prob_obs_meta, linetype = "dashed", linewidth = 2) +
  xlab("Posterior Predictive Rate of Metamorphosis")
M_rep_mat_dead<-  M_rep_mat -1
M_rep_mat_dead[M_rep_mat_dead==1]<-0
M_rep_mat_dead[M_rep_mat_dead==2]<-1
prop_rep_dead<- rowMeans(M_rep_mat_dead)
prob_obs_dead<- length(which(crData$Event==3))/length(crData$Event)
crPPDeath<- ggplot() + 
  geom_histogram(aes(prop_rep_dead), color = "black") + 
  geom_vline(xintercept = prob_obs_dead, linetype = "dashed", linewidth = 2) +
  xlab("Posterior Predictive Rate of Death")
crProps<- cowplot::plot_grid(crPPMeta, crPPDeath, labels = c("B", "C"))
cowplot::plot_grid(crPPDist, crProps, labels = c("A",""), ncol=1, scale = .9)

# Marginal PP checks
M_rep<- M_rep_mat
E <- nrow(M_rep)
obs_by_time <- anura_person_period %>% ungroup() %>%
  group_by(time) %>%
  summarise(p_met = mean(Event_num == 2),
            p_die = mean(Event_num == 3),
            n = n())
time_index <- anura_person_period$time 

## This function extracts the posterior predictive checks
ppc_by_time <- function(event_code) {
  ## for each draw, get mean occurrence of event_code by day -> matrix (days x D)
  sapply(1:E, function(d) tapply(M_rep[d, ] == event_code, time_index, mean)) %>%
    ## transpose: rows = draws, columns = days
    t() %>% 
    ## convert matrix to data frame
    as.data.frame() %>% 
    ## reshape wide (one col/day) to long (one row per draw-day)
    pivot_longer(everything(), names_to = "days", values_to = "p") %>%
    ## column names came in as strings; convert back to numeric
    mutate(days = as.numeric(days)) %>%
    ## group draws by day for summarizing
    group_by(days) %>% 
    ## 95% posterior predictive interval + median per day
    summarise(lo = quantile(p, 0.025), med = quantile(p, 0.5), hi = quantile(p, 0.975))  
}

## Run marginal PPC for each event
ppc_met_time <- ppc_by_time(2)
ppc_die_time <- ppc_by_time(3)

## Compare observed proportions of each event occurrance to predicted values
## from simulated draws
observed<- "y"
predicted<- "y predicted"
Met_plot <- ggplot() +
  geom_ribbon(data = ppc_met_time, aes(days, ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.3) +
  geom_line(data = ppc_met_time, aes(days, med, col = predicted), lwd = 1) +
  geom_point(data = obs_by_time, aes(time, p_met, fill = observed), size = 5, alpha = 0.7) +
  labs(y = "P(Metamorphosis)", x = "Time (fortnights)") +
  scale_color_manual(values = "steelblue") + theme_classic() +
  theme(legend.title = element_blank(), axis.title = element_text(face = "bold", size = 12),
        legend.text = element_text(face = "bold"), axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.spacing.y = unit(.05,"cm")) 

Die_plot <- ggplot() +
  geom_ribbon(data = ppc_die_time, aes(days, ymin = lo, ymax = hi),  fill = "steelblue",alpha = 0.3) +
  geom_line(data = ppc_die_time, aes(days, med, col = predicted), lwd = 1) +
  geom_point(data = obs_by_time, aes(time, p_die, fill = observed), size = 5, alpha = 0.7) +
  labs(y = "P(Death)", x = "Time (fortnights)") +
  scale_color_manual(values = "steelblue") + theme_classic() +
  theme(legend.title = element_blank(), axis.title = element_text(face = "bold", size = 12),
        legend.text = element_text(face = "bold")) 
## Combine plots
legend<- get_legend(Met_plot)
plot_grid(
  plot_grid(Met_plot+theme(legend.position="none"), Die_plot+theme(legend.position="none"), 
          ncol = 1, align = "v", labels = "AUTO"),
  legend, rel_widths = c(3,.4)
  )

#-----Counterfactual Sim-------------------------------------------------
Counterfactual <- expand.grid(Time = 2:9,
                              Density = log(c(10, 50, 90)),
                              Treatment = c(0.4, 3, 5),
                              Chl = c(min(anura_person_period$logChl),
                                      mean(anura_person_period$logChl),
                                      max(anura_person_period$logChl)),
                              Cope = c(min(anura_person_period$CopepodPerLiter),
                                       mean(anura_person_period$CopepodPerLiter),
                                       max(anura_person_period$CopepodPerLiter))/100,
                              Diplo = c(min(anura_person_period$DiploPerLiter),
                                        mean(anura_person_period$DiploPerLiter),
                                        max(anura_person_period$DiploPerLiter))/100
)

## Posterior distributions from Stan
#posterior <- rstan::extract(crModelGeneralMix)      ## Stanfit object
## --- Metamorphosis --- ##
#alpha_post <- posterior$alpha    ## matrix: iterations x T
#beta_post  <- posterior$beta     ## matrix: iterations x P
alpha_post <- crModelGeneral$draws(variables = "alpha", format = "matrix") ## matrix: iterations x T
beta_post  <- crModelGeneral$draws(variables = "beta", format = "matrix")    ## matrix: iterations x P

## --- Death --- ##
#theta_post <- posterior$theta    ## matrix: iterations x T
#gamma_post <- posterior$gamma    ## matrix: iterations x K
theta_post <- crModelGeneral$draws(variables = "theta", format = "matrix")    ## matrix: iterations x T
gamma_post <- crModelGeneral$draws(variables = "gamma", format = "matrix")    ## matrix: iterations x K

D<- 50 # Number of posterior draws

## Create a matrix for each event type
Pi_Non <- matrix(NA, nrow = nrow(Counterfactual), ncol = D)
Pi_Met <- matrix(NA, nrow = nrow(Counterfactual), ncol = D)
Pi_Die <- matrix(NA, nrow = nrow(Counterfactual), ncol = D)


for (d in 1:D) {
  for (i in 1:nrow(Counterfactual)) {
    t = Counterfactual$Time[i]-1
    lambda <- exp(alpha_post[d, t] + sum(Counterfactual[i,-1] * beta_post[d, ]))
    eta <- exp(theta_post[d, t] + sum(Counterfactual[i,-1] * gamma_post[d, ]))
    Pi_Non[i,d] <- exp(-lambda - eta)
    Pi_Met[i,d] <- lambda/(lambda + eta)*(1-exp(-lambda - eta))
    Pi_Die[i,d] <- eta/(lambda + eta)*(1-exp(-lambda - eta))
  }
}

## Take summary statistics
Pi_Non_Sum <- apply(Pi_Non, 1, function(x){quantile(x, c(0.1, 0.5, 0.9))})%>%t()
Pi_Met_Sum <- apply(Pi_Met, 1, function(x){quantile(x, c(0.1, 0.5, 0.9))})%>%t()
Pi_Die_Sum <- apply(Pi_Die, 1, function(x){quantile(x, c(0.1, 0.5, 0.9))})%>%t()


## Take full data frame
CF_Data_Non <- data.frame(Time = Counterfactual$Time,
                          Salinity = Counterfactual$Treatment,
                          Density = exp(Counterfactual$Density),
                          Chl = exp(Counterfactual$Chl),
                          Cope = Counterfactual$Cope,
                          Diplo = Counterfactual$Diplo,
                          Min = Pi_Non_Sum[,1],
                          Med = Pi_Non_Sum[,2],
                          Max = Pi_Non_Sum[,3],
                          Category = "None")

CF_Data_Met <- data.frame(Time = Counterfactual$Time,
                          Salinity = Counterfactual$Treatment,
                          Density = exp(Counterfactual$Density),
                          Chl = exp(Counterfactual$Chl),
                          Cope = Counterfactual$Cope,
                          Diplo = Counterfactual$Diplo,
                          Min = Pi_Met_Sum[,1],
                          Med = Pi_Met_Sum[,2],
                          Max = Pi_Met_Sum[,3],
                          Category = "Metamorphosis")

CF_Data_Die <- data.frame(Time = Counterfactual$Time,
                          Salinity = Counterfactual$Treatment,
                          Density = exp(Counterfactual$Density),
                          Chl = exp(Counterfactual$Chl),
                          Cope = Counterfactual$Cope,
                          Diplo = Counterfactual$Diplo,
                          Min = Pi_Die_Sum[,1],
                          Med = Pi_Die_Sum[,2],
                          Max = Pi_Die_Sum[,3],
                          Category = "Death")
CF_Data_Full <- rbind(CF_Data_Non,
                      CF_Data_Met,
                      CF_Data_Die)
CF_Data_Full$Density = as.factor(CF_Data_Full$Density)

# Cumulative Probabilities
Q <- as.factor(paste(Counterfactual$Treatment, Counterfactual$Density, 
                     Counterfactual$Chl, Counterfactual$Cope, 
                     Counterfactual$Diplo, 
                     sep = "_"))

S_all <- NULL
M_all <- NULL
D_all <- NULL

for (q in 1:length(Q)){ ## For each unique predictor combination q
  ## Extract out specific CF matrix slice
  Pi_Non_q <- Pi_Non[which(Q == unique(Q)[q]),] ## Pr(Survival)
  Pi_Met_q <- Pi_Met[which(Q == unique(Q)[q]),] ## Pr(Metamorphosis)
  Pi_Die_q <- Pi_Die[which(Q == unique(Q)[q]),] ## Pr(Death)
  
  S_q <- matrix(NA, nrow = nrow(Pi_Non_q), ncol = D) ## Cumulative Pr(Survival to t)
  M_q <- matrix(NA, nrow = nrow(Pi_Non_q), ncol = D) ## Cumulative Pr(Metamorphosis to t)
  D_q <- matrix(NA, nrow = nrow(Pi_Non_q), ncol = D) ## Cumulative Pr(Death to t)
  
  for (d in 1:D){ ## For each posterior draw d
    ## Step 1: survival through t
    S_full <- cumprod(Pi_Non_q[, d])
    
    ## Step 2: survival just before t (S_{s-1}) 
    ## prepend 1 (S0 = 1) and drop last element to align dimensions
    S_prev <- c(1, head(S_full, -1))
    
    ## Step 3: cumulative sums
    M_q[, d] <- cumsum(S_prev * Pi_Met_q[, d])
    D_q[, d] <- cumsum(S_prev * Pi_Die_q[, d])
    S_q[, d] <- S_full  ## optional, keep as survival through t
  }
  
  S_all <- rbind(S_all, S_q)
  M_all <- rbind(M_all, M_q)
  D_all <- rbind(D_all, D_q)
}

## Take summary statistics
S_Sum <- apply(S_all, 1, function(x){quantile(x, c(0.1, 0.5, 0.9))})%>%t()
M_Sum <- apply(M_all, 1, function(x){quantile(x, c(0.1, 0.5, 0.9))})%>%t()
D_Sum <- apply(D_all, 1, function(x){quantile(x, c(0.1, 0.5, 0.9))})%>%t()




## Take full data frame
CF_Data_S <- data.frame(Time = Counterfactual$Time,
                        Salinity = Counterfactual$Treatment,
                        Density = exp(Counterfactual$Density),
                        Chl = exp(Counterfactual$Chl),
                        Cope = Counterfactual$Cope,
                        Diplo = Counterfactual$Diplo,
                        Min = S_Sum[,1],
                        Med = S_Sum[,2],
                        Max = S_Sum[,3],
                        Category = "None")

CF_Data_M <- data.frame(Time = Counterfactual$Time,
                        Salinity = Counterfactual$Treatment,
                        Density = exp(Counterfactual$Density),
                        Chl = exp(Counterfactual$Chl),
                        Cope = Counterfactual$Cope,
                        Diplo = Counterfactual$Diplo,
                        Min = M_Sum[,1],
                        Med = M_Sum[,2],
                        Max = M_Sum[,3],
                        Category = "Metamorphosis")

CF_Data_D <- data.frame(Time = Counterfactual$Time,
                        Salinity = Counterfactual$Treatment,
                        Density = exp(Counterfactual$Density),
                        Chl = exp(Counterfactual$Chl),
                        Cope = Counterfactual$Cope,
                        Diplo = Counterfactual$Diplo,
                        Min = D_Sum[,1],
                        Med = D_Sum[,2],
                        Max = D_Sum[,3],
                        Category = "Death")
CF_Data_Cum <- rbind(CF_Data_S,
                     CF_Data_M,
                     CF_Data_D)

CF_Data_Cum$Density = as.factor(CF_Data_Cum$Density)

#-----Counterfactual Vis-------------------------------------------------
mean_chl = exp(mean((anura_person_period$logChl))) %>% round(.,2)
mean_cope = (mean(anura_person_period$CopepodPerLiter)/100)%>% round(.,2)
mean_diplo = (mean(anura_person_period$DiploPerLiter)/100)  %>% round(.,2)

# Denisty Plots
Density <- CF_Data_Full %>% filter(round(Chl,2) == mean_chl, round(Cope,2) == mean_cope,
                         round(Diplo,2) == mean_diplo, Salinity == 3) %>%
    ggplot()+
    theme_classic() +
    geom_line(aes(Time, Med, col = Density), lwd=1, position = position_dodge(width=.5))+
    geom_point(aes(Time, Med, col = Density), size=3, position = position_dodge(width=.5))+
    geom_linerange(aes(Time, ymin = Min, ymax = Max, color = Density), lwd=4, alpha = 0.3, position = position_dodge(width=.5))+
    scale_color_manual(values = c("goldenrod","#b73779","#000004"), name = "Tadpole Density") +
    theme(legend.position = "top")+labs(y = "Hazard", x = "Time (fortnights)") + 
    facet_wrap(~as.factor(Category), labeller = as_labeller(c(`Death` = "Death*", `Metamorphosis` = "Metamorphosis*", `None`="None"))) +
    theme(axis.title = element_text(face = "bold", size = 12), axis.title.x = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank()) #+
    #scale_colour_manual(values=my_orange)
CumDensity <-CF_Data_Cum %>% filter(round(Chl,2) == mean_chl, round(Cope,2) == mean_cope,
                        round(Diplo,2) == mean_diplo, Salinity == 3) %>%
  ggplot()+
  theme_classic() +
  geom_line(aes(Time, Med, col = Density), lwd=1, position = position_dodge(width=.5))+
  geom_point(aes(Time, Med, col = Density), size=3, position = position_dodge(width=.5))+
  geom_linerange(aes(Time, ymin = Min, ymax = Max, color = Density), lwd=4, alpha = 0.3, position = position_dodge(width=.5))+
  #scale_colour_manual(values=my_orange) +
  scale_color_manual(values = c("goldenrod","#b73779","#000004")) +
  theme(legend.position = "top")+labs(y = "Cumulative Probability", x = "Time (fortnights)") + 
  facet_wrap(~as.factor(Category)) +
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12),
        strip.text.x = element_blank())
plot_grid(Density, CumDensity, ncol = 1,  align = "v", labels = "AUTO")  

#Salinity Plots
Salinity <- CF_Data_Full %>% filter(round(Chl,2) == mean_chl, round(Cope,2) == mean_cope, round(Diplo,2) == mean_diplo, 
                         Density == 50) %>% mutate(Salinity = as.factor(Salinity)) %>%
  ggplot()+
  theme_classic() +
  geom_line(aes(Time, Med, col = Salinity), lwd=1, position = position_dodge(width=.5))+
  geom_point(aes(Time, Med, col = Salinity), size=3, position = position_dodge(width=.5))+
  geom_linerange(aes(x =Time, Med, ymin = Min, ymax = Max, color = Salinity), lwd=4, alpha = 0.3, position = position_dodge(width=.5))+
  theme(legend.position = "top")+labs(y = "Hazard", x = "Time (fortnights)") + 
  facet_wrap(~as.factor(Category), labeller = as_labeller(c(`Death` = "Death*", `Metamorphosis` = "Metamorphosis", `None`="None"))) +
  theme(axis.title = element_text(face = "bold", size = 12), axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
CumSalinity <- CF_Data_Cum %>% filter(round(Chl,2) == mean_chl, round(Cope,2) == mean_cope, round(Diplo,2) == mean_diplo, 
                                      Density == 50) %>% mutate(Salinity = as.factor(Salinity)) %>%
  ggplot()+
  theme_classic() +
  geom_line(aes(Time, Med, col = Salinity), lwd=1, position = position_dodge(width=.5))+
  geom_point(aes(Time, Med, col = Salinity), size=3, position = position_dodge(width=.5))+
  geom_linerange(aes(Time, ymin = Min, ymax = Max, color = Salinity), lwd=4, alpha = 0.3, position = position_dodge(width=.5))+
  labs(y = "Cumulative Probability", x = "Time (fortnights)") + 
  facet_wrap(~as.factor(Category)) + 
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12),
        strip.text.x = element_blank())
plot_grid(Salinity, CumSalinity, ncol = 1,  align = "v", labels = "AUTO")

#Food Plots
Chl <- CF_Data_Full %>% filter(Salinity == 3, round(Cope,2) == mean_cope, round(Diplo,2) == mean_diplo, 
                                    Density == 50) %>% mutate(Chl = round(Chl,2)) %>%
                                    mutate(Chl = as.factor(Chl)) %>%
  ggplot()+
  theme_classic() +
  geom_line(aes(Time, Med, col = Chl), lwd=1, position = position_dodge(width=.5))+
  geom_point(aes(Time, Med, col = Chl), size=3, position = position_dodge(width=.5))+
  geom_linerange(aes(Time, ymin = Min, ymax = Max, color = Chl), lwd=4, alpha = 0.3,position = position_dodge(width=.5))+
  #scale_colour_manual(values=my_green) +
  scale_color_manual(values = c("#5ec962","#3b528b","#440154"), 
                     name = expression(paste("Chlorophyll-A Concentration(\u03bcg L"^"-1",")"))) +
  theme(legend.position = "top")+labs(y = "Hazard", x = "Time (fortnights)") + 
  facet_wrap(~as.factor(Category), labeller = as_labeller(c(`Death` = "Death*", `Metamorphosis` = "Metamorphosis*", `None`="None"))) +
  theme(axis.title = element_text(face = "bold", size = 12), axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
CumChl <- CF_Data_Cum %>%  filter(Salinity == 3, round(Cope,2) == mean_cope, round(Diplo,2) == mean_diplo, 
                                  Density == 50) %>% mutate(Chl = round(Chl,2)) %>%
                                  mutate(Chl = as.factor(Chl)) %>%
  ggplot()+
  theme_classic() +
  geom_line(aes(Time, Med, col = Chl), lwd=1, position = position_dodge(width=.5))+
  geom_point(aes(Time, Med, col = Chl), size=3, position = position_dodge(width=.5))+
  geom_linerange(aes(Time, ymin = Min, ymax = Max, color = Chl), lwd=4, alpha = 0.3, position = position_dodge(width=.5))+
  #scale_colour_manual(values=my_green) +
  scale_color_manual(values = c("#5ec962","#3b528b","#440154")) +
  labs(y = "Cumulative Probability", x = "Time (fortnights)") + 
  facet_wrap(~as.factor(Category)) + 
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12),
        strip.text.x = element_blank()) 
plot_grid(Chl, CumChl, ncol = 1,  align = "v", labels = "AUTO")

#Cope Plots
Cope <- CF_Data_Full %>% filter(Salinity == 3, round(Chl,2) == mean_chl, round(Diplo,2) == mean_diplo, 
                               Density == 50) %>% mutate(Cope = round(Cope,2)) %>%
  mutate(Cope = as.factor(Cope)) %>%
  ggplot()+
  theme_classic() +
  geom_line(aes(Time, Med, col = Cope), lwd=1, position = position_dodge(width=.5))+
  geom_point(aes(Time, Med, col = Cope), size=3, position = position_dodge(width=.5))+
  geom_linerange(aes(Time, ymin = Min, ymax = Max, color = Cope), lwd=4, alpha = 0.3, position = position_dodge(width=.5))+
  #scale_colour_manual(values=my_red) +
  scale_color_manual(values = c("#fc8961","#b73779","#51127c"), 
                     name = expression(paste("Copepod Concentration (No. hL"^"-1",")"))) +
  theme(legend.position = "top")+labs(y = "Hazard", x = "Time (fortnights)") + 
  facet_wrap(~as.factor(Category), labeller = as_labeller(c(`Death` = "Death*", `Metamorphosis` = "Metamorphosis", `None`="None"))) +
  theme(axis.title = element_text(face = "bold", size = 12), axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
CumCope<- CF_Data_Cum %>% filter(Salinity == 3, round(Chl,2) == mean_chl, round(Diplo,2) == mean_diplo, 
                                 Density == 50) %>% mutate(Cope = round(Cope,2)) %>%
  mutate(Cope = as.factor(Cope)) %>%
  ggplot()+
  theme_classic() +
  geom_line(aes(Time, Med, col = Cope), lwd=1, position = position_dodge(width=.5))+
  geom_point(aes(Time, Med, col = Cope), size=3, position = position_dodge(width=.5))+
  geom_linerange(aes(Time, ymin = Min, ymax = Max, color = Cope), lwd=4, alpha = 0.3, position = position_dodge(width=.5))+
  #scale_colour_manual(values=my_red) +
  scale_color_manual(values = c("#fc8961","#b73779","#51127c"), 
                     name = expression(paste("Copepod Concentration (No. hL"^"-1",")"))) +
  labs(y = "Cumulative Probability", x = "Time (fortnights)") + 
  facet_wrap(~as.factor(Category)) + 
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12),
        strip.text.x = element_blank()) 
plot_grid(Cope, CumCope, ncol = 1,  align = "v", labels = "AUTO")

#Diplo Plots
Diplo <- CF_Data_Full %>% filter(Salinity == 3, round(Chl,2) == mean_chl, round(Cope,2) == mean_cope, 
                                Density == 50) %>% mutate(Diplo = round(Diplo,2)) %>%
  mutate(Diplo = as.factor(Diplo)) %>%
  ggplot()+
  theme_classic() +
  geom_line(aes(Time, Med, col = Diplo), lwd=1, position = position_dodge(width=.5))+
  geom_point(aes(Time, Med, col = Diplo), size=3, position = position_dodge(width=.5))+
  geom_linerange(aes(Time, ymin = Min, ymax = Max, color = Diplo), lwd=4, alpha = 0.3, position = position_dodge(width=.5))+
  #scale_colour_manual(values=my_purple) +
  scale_color_manual(values = c("#A6BDD7","#C10020","#CEA262"), 
                     name = expression(paste("Diplostracan Concentration (No. hL"^"-1",")"))) +
  theme(legend.position = "top")+labs(y = "Hazard", x = "Time (fortnights)") + 
  facet_wrap(~as.factor(Category), labeller = as_labeller(c(`Death` = "Death*", `Metamorphosis` = "Metamorphosis", `None`="None"))) +
  theme(axis.title = element_text(face = "bold", size = 12), axis.title.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
CumDiplo<- CF_Data_Cum %>% filter(Salinity == 3, round(Chl,2) == mean_chl, round(Cope,2) == mean_cope, 
                                  Density == 50) %>% mutate(Diplo = round(Diplo,2)) %>%
  mutate(Diplo = as.factor(Diplo)) %>%
  ggplot()+
  theme_classic() +
  geom_line(aes(Time, Med, col = Diplo), lwd=1, position = position_dodge(width=.5))+
  geom_point(aes(Time, Med, col = Diplo), size=3, position = position_dodge(width=.5))+
  geom_linerange(aes(Time, ymin = Min, ymax = Max, color = Diplo), lwd=4, alpha = 0.3, position = position_dodge(width=.5))+
  #scale_colour_manual(values=my_purple) +
  scale_color_manual(values = c("#A6BDD7","#C10020","#CEA262"), 
                     name = expression(paste("Diplostracan Concentration (No. hL"^"-1",")"))) +
  labs(y = "Cumulative Probability", x = "Time (fortnights)") + 
  facet_wrap(~as.factor(Category)) + 
  theme(legend.position = "none", axis.title = element_text(face = "bold", size = 12),
        strip.text.x = element_blank()) 
plot_grid(Diplo, CumDiplo, ncol = 1,  align = "v", labels = "AUTO")


#-----Alternate Model-------------------------------------------------
generalMatrix2<-model.matrix( ~ 0 + log(tadNumLag) + Treatment + ChlLag + 
                               CopeLag10 + DiploLag10, anura_person_period)
crData2<- 
  list(
    T         = length(unique(anura_person_period$time)),
    N         = nrow(anura_person_period),
    P         = ncol(generalMatrix2),
    Time      = anura_person_period$time,
    K         = length(unique(anura_person_period$Event_num)),
    Event     = anura_person_period$Event_num,
    X         = generalMatrix2,
    G         = length(unique(anura_person_period$Tank)),
    tank      = as.numeric(anura_person_period$Tank),
    E         = length((unique(anura_person_period$EnID))),
    enclosure = as.numeric(anura_person_period$EnID)
  )


crModelGeneralcmd2<- crModelcmd$sample(
  data = crData2, chains = 1, iter_warmup=5000, iter_sampling = 5000, init = 0
)
crModelGeneralcmd2$save_object("Code/Final Code/Saved Models/savedCR/GTF_Generalcmd2.rds")
crModelGeneralcmd2$summary(variables = c("beta", "gamma")) %>% 
  as.data.frame()  %>%  
  mutate(coef = rep(colnames(generalMatrix2), 2)) %>%
  mutate(id = rep(c("M","D"), each =5)) %>%
  mutate(Converge = if_else(rhat < 1.05, TRUE, FALSE)) %>% 
  mutate(Meaningful = if_else(q5 <= 0 & q95 >= 0, FALSE, TRUE)) %>% 
  select(!c(variable, mean, sd, mad, rhat, ess_bulk, ess_tail)) %>%
  mutate_at(c("median", "q5", "q95"), round, 2)

