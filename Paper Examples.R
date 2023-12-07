
library("nimble");
library("dplyr");
library("survival");
library("runjags");
library("ggplot2")
pathway <- "~/Change-Point Simulation Studies/Simulation Study 2023/"
source(paste0(pathway, "Nimble Functions.R"))
source(paste0(pathway, "Jags_codes2.R"))
Ccp_LL<- compileNimble(cp_LL)


n.thin <- 5
n.burnin <- 5000
n.samp <- 10000


E1690.dat <- read.table(paste0(pathway,"e1690.missing.dat"),
                        header=TRUE)

# Ibhrahim data
#Drop PFS events with time equal zero
E1690.dat <- E1690.dat[-which(E1690.dat$FAILTIME ==0),]

#Convert to the correct notation for survival objects
E1690.dat[which(E1690.dat$FAILCENS == 1),"FAILCENS"] <-0
E1690.dat[which(E1690.dat$FAILCENS == 2),"FAILCENS"] <-1
E1690.dat[which(E1690.dat$SURVCENS == 1),"SURVCENS"] <-0
E1690.dat[which(E1690.dat$SURVCENS == 2),"SURVCENS"] <-1
E1690.dat$AGE <- as.numeric(E1690.dat$AGE)
E1690.dat$AGE[is.na(E1690.dat$AGE)] <- mean(E1690.dat$AGE, na.rm = T)
E1690.dat$AGE_scale <- as.numeric(scale(E1690.dat$AGE))
E1690.dat$AGE_scale_abs <- abs(E1690.dat$AGE_scale)
E1690.dat$TRT <- E1690.dat$TRT -1


E1690.dat <- E1690.dat %>% mutate(Treatment = ifelse(TRT == 1, "INF", "OBS"))
E1690.dat$Treatment <- as.factor(E1690.dat$Treatment)
E1690.dat$Treatment <- relevel(E1690.dat$Treatment, ref ="OBS")

E1690.dat_final <- E1690.dat[,c("SURVTIME","SURVCENS","TRT", "AGE_scale")]
E1690.dat_final <- cbind(tstart = 0,E1690.dat_final)


fit.OS <- survfit(Surv(SURVTIME, SURVCENS)~Treatment,
                  data = E1690.dat)

#Define Model matrix

ff <- Surv(SURVTIME,SURVCENS)~as.factor(Treatment)+AGE_scale
utils::str(m <- model.frame(ff, E1690.dat))
mat <- model.matrix(ff, m)
mat[,3] <- 0
mat <- mat[!duplicated(mat),]


#-- Sceanrio 1 Change-point model assuming a change-point in HR and baseline hazards
#Important defines the model that you are going to fit.
beta_1 <- beta_2 <- matrix(1, byrow = FALSE, ncol = 2, nrow = 3)
beta_2[2:3,1:2] <-0


data_jags <- list()
data_jags$N <- nrow(E1690.dat)
data_jags$time <- E1690.dat$SURVTIME
data_jags$status <- E1690.dat$SURVCENS
data_jags$MAXX <- max(E1690.dat$SURVTIME)
data_jags$N_CP <- 1
data_jags$X_mat<- as.matrix(cbind(1 ,E1690.dat$TRT, E1690.dat$AGE_scale))
data_jags$beta_1_ind <- beta_1
data_jags$beta_2_ind <- beta_2

#If you wanted a fixed change-point... Set to zero if you want the change-point to be estimated
data_jags$cp_fix <- c(0, 1, 5)
data_jags$cp_fixed <- 0
mod <- "wei"
#Change-inits


mod.cp_scenario1 <- runjags::run.jags(
  model = get(paste0("jags.piecewise_",mod,"_chng")),
  #model = jags.piecewise.expo_chng,
  data = data_jags,
  n.chains = 2,
  monitor = c("cp", "beta_cp","beta_cp_anc", "total_loglik","loglik"),
  sample=n.samp,
  thin = n.thin,
  burnin = n.burnin,
  inits = inits,
  method ='rjparallel',
  summarise = FALSE)


#runjags::failed.jags(c('model','data','inits'))[["model]]
#flexsurv::flexsurvreg(Surv(SURVTIME,SURVCENS)~as.factor(TRT)+AGE_scale,data =E1690.dat_final, dist = "weibullPH" )



df_origin <- E1690.dat_final %>% rename(time = SURVTIME, status = SURVCENS, arm = TRT)

scenario1_rmst <- export_outcomes(fit.km = fit.OS,mod.cp_scenario1, data_jags = data_jags,
                pred_mod_mat = mat,
                pred_time= seq(0,15, by = 0.1),
                plt_title = "Weibull 1 change-point model (Proportional Hazards)",
                xlab_title = "Time (Years)",
                pathway_export= gen_pathway(paste0(pathway,"Scenario 1/")),
                df.origin =  df_origin)


scenario1_rmst[1,1] <- "Change-point for HR"


cp_all <- data.frame(Scenario = "Change-point for HR",
                     CP = as.matrix(mod.cp_scenario1$mcmc)[,"cp[2]"])

mod.cp_scenario1_mat <- as.matrix(mod.cp_scenario1$mcmc)
df_params_weib <- rbind(quantile(mod.cp_scenario1_mat[, "cp[2]"], probs = c(0.025,0.5,0.975)),
                       quantile(exp(mod.cp_scenario1_mat[, "beta_cp[2,1]"]), probs = c(0.025,0.5,0.975)),
                       quantile(exp(mod.cp_scenario1_mat[, "beta_cp[2,2]"]), probs = c(0.025,0.5,0.975)),
                       quantile(exp(mod.cp_scenario1_mat[, "beta_cp[3,1]"]), probs = c(0.025,0.5,0.975)))

df_params_final <- cbind(Param = c("Change-point", "HR 1", "HR 2", "HR Age"),
                         as.data.frame(df_params_weib))

write.csv(df_params_final,paste0(pathway,"summary_output_df_weib_cp.csv"))
print(xtable::xtable(df_params_final %>% mutate_if(is.numeric, round, digits = 2)),   include.rownames = FALSE)


#Scenario 2 -- Common Hazards for final interval

data_jags$beta_1_ind <- beta_1
data_jags$beta_2_ind <- beta_2

data_jags$beta_1_ind[2,2] <- 0


mod.cp_scenario2 <- runjags::run.jags(
  model = get(paste0("jags.piecewise_",mod,"_chng")),
  #model = jags.piecewise.expo_chng,
  data = data_jags,
  n.chains = 2,
  monitor = c("cp", "beta_cp","beta_cp_anc", "total_loglik","loglik"),
  sample=n.samp,
  thin = n.thin,
  burnin = n.burnin,
  inits = inits,
  method ='rjparallel',
  summarise = FALSE)

scenario2_rmst <- export_outcomes(fit.km = fit.OS,mod.cp_scenario2, data_jags = data_jags,
                                  pred_mod_mat = mat,
                                  pred_time= seq(0,15, by = 0.1),
                                  plt_title = "Weibull 1 change-point model (Equal Final Hazards)",
                                  xlab_title = "Time (Years)",
                                  pathway_export= gen_pathway(paste0(pathway,"Scenario 2/")),
                                  df.origin =  df_origin)


scenario2_rmst[1,1] <- "Change-point: Equal Final Hazards"

cp_all <- rbind(cp_all, data.frame(Scenario = "Change-point for HR (Equal Final Hazards)",
                     CP = as.matrix(mod.cp_scenario2$mcmc)[,"cp[2]"]))

rm(mod.cp_scenario2)
## Scenario 3 -- Converging Hazards


data_jags$beta_1_ind <- beta_1
data_jags$beta_2_ind <- beta_2
data_jags$beta_1_ind[2,2] <- 0

trt_ind <- 2
data_jags$trt_ind <- trt_ind

mod.cp_scenario3 <- runjags::run.jags(
  model = get(paste0("jags.piecewise_",mod,"_wane")),
  #model = jags.piecewise.expo_chng,
  data = data_jags,
  n.chains = 2,
  monitor = c("cp", "beta_cp","beta_cp_anc", "total_loglik", "loglik", "lambda_wane"),
  sample=n.samp,
  thin = n.thin,
  burnin = n.burnin,
  inits = inits,
  method ='rjparallel',
  summarise = FALSE)


scenario3_rmst <- export_outcomes(fit.km = fit.OS,mod.cp_scenario3, data_jags = data_jags,
                                  pred_mod_mat = mat,
                                  pred_time= seq(0,15, by = 0.1),
                                  plt_title = "Weibull 1 change-point model (Converging hazards)",
                                  xlab_title = "Time (Years)",
                                  pathway_export= gen_pathway(paste0(pathway,"Scenario 3/")),
                                  df.origin =  df_origin)

scenario3_rmst[1,1] <- "Change-point: Converging Hazards"


scenario_all <- rbind(scenario1_rmst, scenario2_rmst[1,],scenario3_rmst[1,])
write.csv(scenario_all, file = paste0(pathway,"RMST_scenario_all.csv"))


cp_all <- rbind(cp_all, data.frame(Scenario = "Change-point for HR (Converging hazards)",
                                   CP = as.matrix(mod.cp_scenario3$mcmc)[,"cp[2]"]))



mu <- cp_all %>% group_by(Scenario) %>% summarise(grp.mean=mean(CP))

gg <- ggplot(cp_all, aes(x = CP,fill = Scenario))+
  geom_density(alpha=0.4)+
  xlab("Time (Years)")+
  ggtitle("Posterior Distribution of Change-point")+

  #geom_vline(data=mu, aes(xintercept=grp.mean, color=Model),
  #           linetype="dashed")+
  theme_bw()+
  theme(legend.position="bottom")


ggsave(filename = paste0(pathway, "Changepoint-Posterior-Distribution.png"),
       plot = gg,
       width = 10, height = 5, units = 'in')


rm(mod.cp_scenario2)

# Cox Model Example ---

#Alternative More flexible code not currently implemented

if(FALSE){
  library("survival")
  library("nimble")
  #undebug(LL_cox_ph)
  #debug(LL_cox_ph)
  # LL_cox_ph(params_cp = c(0,2,10000),
  #           params_beta = c(-0.3,0.3,0.2), 
  #           data = data_nimble$data,
  #           penalty = 0)
  
  LL_cox_ph <- function(params_cp,params_beta, data, penalty = 0){
    params_cp <- head(params_cp, n = -1)
    params_cp <- tail(params_cp, n = -1)    
    params_cp <- sort(params_cp)
    
    ncp <- length(params_cp)
    data <- as.data.frame(data)
    ncol <- dim(data)[2]
    colnames(data)[1:2] <- c("time", "status")
    colnames(data)[3:ncol] = LETTERS[1:(ncol(data)-2)]
    data_split <- try({survival::survSplit(data = data,
                                      cut = params_cp,
                                      end = "time",
                                      start = "start",
                                      event = "status", 
                                      episode = "period")
      })
    
    if(class(data_split) == "try-error"){
      return(-999999)
    }
    
    if(ncol(data)>3){
      formula.eval <- paste0(paste0("I((A==1)&(period==",1:(ncp+1),"))",collapse = "+"), "+",
                             paste0(LETTERS[2:(ncol(data)-3)],collapse = "+"))
    }else{
      formula.eval <- paste0("I((A==1)&(period==",1:(ncp+1),"))",collapse = ",")
    }
    
    formula_all <- as.formula(paste0("Surv(start,time,status) ~", formula.eval))
    
    m <- model.frame(formula_all, data_split)
    mat <- model.matrix(formula_all,m)
    mat_final <- mat[,-1]#remove intercept
    index_trt <- grep("^A$",colnames(mat_final))
    mat_final <- mat_final[,-index_trt]#remove redundant treatment; we don't need it
    colnames(mat_final) <- LETTERS[1:(ncol(mat_final))]
    
    
    mod.matrix <- cbind(start = data_split$start, time = data_split$time, status = data_split$status, mat_final)
   
    if(penalty == 1){
      formula_post_tidle <- paste0("ridge(", paste0(colnames(mat_final), collapse = ","),", theta = 3")
    }else{
      formula_post_tidle <- paste0(colnames(mat_final), collapse = "+")
    }
    
    formula_final <- as.formula(paste0("Surv(start, time, status) ~", formula_post_tidle))
    
    #browser()
    cox_fit <- try({ coxph(formula_final,
                           data = data.frame(mod.matrix),
                           init = params_beta,
                           control = coxph.control(iter.max = 0))
      
    })
    
    if(class(cox_fit)[1] == "try-error"){
      return(-999999)
    }else{
        return(cox_fit$loglik[2])
    }
    
  }
  
  R_LL_cox_ph <- nimbleRcall(function(params_cp = double(1),
                                      params_beta = double(1),
                                      data = double(2),
                                      penalty = double(0)){},
                             Rfun = 'LL_cox_ph',
                             returnType = double(0)
                             )
  
  sort_vec <- nimbleFunction(
    #X matrix tstart; tstop,status, X var all
    run = function(X = double(1)){
        lx      <- dim(X)[1]
        lxMin1  <- lx-1
        madeChange <- TRUE        
        while(madeChange) {
          madeChange <- FALSE
          for (ii in 1:lxMin1) {
            ii1 <- ii+1
            if (X[ii] > X[ii1]) {
              tmp    <- X[ii]
              X[ii]  <- X[ii1]
              X[ii1] <- tmp
              madeChange <- TRUE
            }
          }
        }
        return(X)
        returnType(double(1))
      }
    )
    
#Csort_vec <-   nimble::compileNimble(sort_vec)
 
 
  mod <- nimbleCode({
    C <- 10000
    cp_unsort[1] <- 0
    for(i in 1:N_CP){
      params_cp[i] ~ dunif(0, max_cp)
      cp_unsort[i+1] <- params_cp[i]
    }
    cp_unsort[N_CP+2] <- 1000
    cp[1:(N_CP+2)] <- sort_vec(cp_unsort[1:(N_CP+2)])
    
    for(i in 2:(N_CP+1)){
      diff[i-1] <- cp[i]-cp[i-1]
      #log_diff[i-1] <- log(diff[i])
    }
    
    diff[N_CP+1] <- 1 #need to padd but needs to evaluate to 0
    
    #sum_log_diff <- sum(log_diff[1:N_CP])
    
    log_prior <- loggam(2*N_CP +2) + sum(log(diff[1:(N_CP+1)])) - ((2*N_CP +1)*log(MAXX))
    zero_prior ~ dpois(C - log_prior)
    
    
    for(i in 1:N_covar){
      params_beta[i] ~ dnorm(0, 5)
     
    } 
    
    LL <- R_LL_cox_ph(cp[1:(N_CP+2)],
                      params_beta[1:N_covar],
                      data[1:df_nrow,1:df_ncol],
                      0)
    
    zero ~ dpois(-LL + C)
    
  })

    
  
data_mod <- as.matrix(E1690.dat[, c("SURVTIME","SURVCENS", "TRT", "AGE_scale")])
#params_beta = c(1,1) #coefficient for hazard(s) and changepoint
#params_cp = c(10, -999) #change-point padding the last val

data_nimble = list(zero = 0,
                   zero_prior = 0,
                   data =data_mod, max_cp = max(data_mod[,1])*0.7)
N_CP <- 1

constants_nimble = list(N_CP = N_CP,
                        df_nrow = nrow(data_mod),
                        df_ncol = ncol(data_mod),
                        N_covar = 1+N_CP+1,
                        MAXX = max(data_mod[,1])) #Age plus the TRT


#https://oliviergimenez.github.io/nimble-workshop/#25
mod_Cox <- nimbleMCMC(mod,
                    data = data_nimble,
                    constants = constants_nimble,
                    monitors = c("params_beta", "cp", "LL", "log_prior"),
                    niter = n.samp*10,
                    inits = inits,
                    nchains =2,
                    nburnin = n.burnin*10)

#coxph_fit <- coxph(Surv(SURVTIME,SURVCENS)~ TRT+AGE_scale,data = data.frame(data_mod))
#coxph_fit$loglik

mod_Cox_mat <- do.call("rbind", mod_Cox)

index_save <- grep("cp|beta",colnames(mod_Cox_mat))

colnames_save <- index_save[apply(mod_Cox_mat[,index_save],2, function(x){ifelse(length(unique(x)) ==1, FALSE,TRUE)})]

mcmc_list_eval <- lapply(mod_Cox,function(x){x[,colnames_save]})
mcmc_list_eval <- lapply(mcmc_list_eval,as.mcmc)

ggmcmc(ggs(coda::as.mcmc.list(mcmc_list_eval)), file = gen_pathway(paste0(pathway,"Cox-CP/","MCMC_diag.pdf")))



df_params_cox <- rbind(quantile(mod_Cox_mat[, "cp[2]"], probs = c(0.025,0.5,0.975)),
                   quantile(exp(mod_Cox_mat[, "params_beta[1]"]), probs = c(0.025,0.5,0.975)),
                   quantile(exp(mod_Cox_mat[, "params_beta[2]"]), probs = c(0.025,0.5,0.975)),
                   quantile(exp(mod_Cox_mat[, "params_beta[3]"]), probs = c(0.025,0.5,0.975)))

df_params_final <- cbind(Param = c("Change-point", "HR 1", "HR 2", "HR Age"),
                         as.data.frame(df_params_cox))

write.csv(df_params_final,paste0(pathway, "summary_output_df_cox_partial.csv"))

# plot(density(model_res_final[, "LL"]))
# plot(density(model_res_final[, "params_cp[1]"]))
# hist(model_res_final[, "params_cp[1]"])


cp_all <- rbind(cp_all, data.frame(Scenario = "Change-point for HR (Cox Model)",
                                   CP = mod_Cox_mat[,"cp[2]"]))


# res.df <- data.frame(changepoint = c(model_output_compare[,"cp[2]"],model_output_compare2[, "params_cp[1]"]),
#                      Model = rep(c("Weibull", "Cox"), c(nrow(model_output_compare),nrow(model_output_compare2))))


#mu <- res.df %>% group_by(Model) %>% summarise(grp.mean=mean(changepoint))

cp_all <- cp_all %>% mutate(Scenario = ifelse(Scenario == "Change-point for HR", "Change-point for HR (Weibull)",
                                               Scenario))
gg <- ggplot(cp_all %>% filter(Scenario %in% c("Change-point for HR (Cox Model)","Change-point for HR (Weibull)") ), aes(x = CP,fill = Scenario))+
  geom_density(alpha=0.4)+
  xlab("Time (Years)")+
  ggtitle("Posterior Distribution of Change-point")+
  
  #geom_vline(data=mu, aes(xintercept=grp.mean, color=Model),
  #           linetype="dashed")+
  theme_bw()+
  theme(legend.position="bottom")


ggsave(filename = paste0(pathway, "Changepoint-Posterior-Distribution (Cox vs Weibull).png"),
       plot = gg,
       width = 10, height = 5, units = 'in')


}

#HTA Examples



# Treatment Delay
#TA347

path_mod_eval <- pathway

TA347_OS_N_df <- read.table(paste0(path_mod_eval,"TA347_N_OS_Initial","/IPDdata.txt"), header = T)
TA347_OS_N_df$arm <- 1#"Nintedanib"
TA347_OS_Doce_df <- read.table(paste0(path_mod_eval,"TA347_Doce_OS_Initial","/IPDdata.txt"), header = T)
TA347_OS_Doce_df$arm <- 0#"Docetaxel"
TA347_df_all <- rbind(TA347_OS_N_df,TA347_OS_Doce_df)

TA347_df_all <- TA347_df_all %>% mutate(Treatment = ifelse(arm == 1,
                                                           "Nintedanib + Docetaxel",
                                                           "Docetaxel"))
TA347_df_all$Treatment <- as.factor(TA347_df_all$Treatment)
TA347_df_all$Treatment <- relevel(TA347_df_all$Treatment, ref ="Docetaxel")


TA347.km <- survfit(Surv(time,status)~Treatment, data = TA347_df_all)


ff <- Surv(time,status)~Treatment
utils::str(m <- model.frame(ff, TA347_df_all))
mat <- model.matrix(ff, m)
mat <- mat[!duplicated(mat),]


if(FALSE){
  survp <- ggsurvplot(TA347.km, legend.labs = c("Docetaxel", "Nintedanib + Docetaxel"), xlab = "Time (Months)", conf.int = F)
  ggsave(file = paste0(pathway, "TA347-Kaplan-Meier.png"), width = 10)
  ggsurvplot(TA347.km, fun = "cumhaz")

}


data_jags <- list()
data_jags$N <- nrow(TA347_df_all)
data_jags$time <- TA347_df_all$time
data_jags$status <- TA347_df_all$status
data_jags$MAXX <- max(TA347_df_all$time)
data_jags$N_CP <- 1
data_jags$X_mat<- matrix(c(rep(1, data_jags$N), TA347_df_all$arm),  ncol = 2)

data_jags$cp_fix <- c(0, 1, 5)
data_jags$cp_fixed <- 0


beta_1 <- beta_2 <- matrix(1, byrow = FALSE, ncol = 2, nrow = 2)
beta_1[2,1] <- 0
beta_2[2,1:2] <-0

data_jags$beta_1_ind <- beta_1
data_jags$beta_2_ind <- beta_2


mod.cp_TA347 <- runjags::run.jags(
  model = get("jags.piecewise_wei_chng"),
  data = data_jags,
  n.chains = 2,
  monitor = c("cp", "beta_cp","beta_cp_anc", "total_loglik","loglik"),
  sample=n.samp,
  thin = n.thin,
  burnin = n.burnin,
  inits = inits,
  method ='rjparallel',
  summarise = FALSE)




TA347_rmst <- export_outcomes(fit.km = TA347.km,mod.cp_TA347, data_jags = data_jags,
                                  pred_mod_mat = mat,
                                  pred_time= seq(0,50, by = 1),
                                  plt_title = "Weibull 1 change-point model (Common hazards Before Change-point)",
                                  xlab_title = "Time (Months)",
                                  pathway_export= gen_pathway(paste0(pathway,"TA347/")),
                                  df.origin =  TA347_df_all)

rm(mod.cp_TA347)

# BRIM - Convergence of hazards


#Bagust and Beale hazard convergence.


path_BRIM <- paste0(pathway,"BRIM-3/BRIM-3/")

TA269_OS_D_df <- xlsx::read.xlsx(paste0(path_BRIM,"Dacarbazine","/Pseudo_IPD.xlsx"), header = T, sheetIndex = 1)[,-1]
TA269_OS_D_df$arm <- 0

TA269_OS_V_df <- xlsx::read.xlsx(paste0(path_BRIM,"Vemurafenib","/Pseudo_IPD.xlsx"), header = T, sheetIndex = 1)[,-1]
TA269_OS_V_df$arm <- 1

TA269_df_all <- rbind(TA269_OS_V_df,TA269_OS_D_df)

TA269_df_all <- TA269_df_all %>% mutate(Treatment = ifelse(arm == 1,
                                                           "Vemurafenib",
                                                           "Dacarbazine"))
TA269_df_all$Treatment <- as.factor(TA269_df_all$Treatment)
TA269_df_all$Treatment <- relevel(TA269_df_all$Treatment, ref ="Dacarbazine")

TA269.km <- survfit(Surv(time,status)~Treatment, data = TA269_df_all)


ff <- Surv(time,status)~Treatment
utils::str(m <- model.frame(ff, TA269_df_all))
mat <- model.matrix(ff, m)
mat <- mat[!duplicated(mat),]


data_jags <- list()
data_jags$N <- nrow(TA269_df_all)
data_jags$time <- TA269_df_all$time
data_jags$status <- TA269_df_all$status
data_jags$MAXX <- max(TA269_df_all$time)
data_jags$N_CP <- 1
data_jags$X_mat<- matrix(c(rep(1, data_jags$N), TA269_df_all$arm),  ncol = 2)

data_jags$cp_fix <- c(0, 1, 5)
data_jags$cp_fixed <- 0


beta_1 <- beta_2 <- matrix(1, byrow = FALSE, ncol = 2, nrow = 2)

beta_2[1,1:2] <- 2 #same intercept for the shape across both intervals
beta_1[1,1:2] <- 2 #same intercept for the scale across both intervals
beta_1[2,2] <- 0 # allow a treatment effect in the first interval but not the second
beta_2[2,2] <- 0 # allow a treatment effect in the first interval but not the second


data_jags$beta_1_ind <- beta_1
data_jags$beta_2_ind <- beta_2

mod.cp_TA269 <- runjags::run.jags(
  model = get("jags.piecewise_wei_chng"),
  data = data_jags,
  n.chains = 2,
  monitor = c("cp", "beta_cp","beta_cp_anc", "total_loglik","loglik"),
  sample=n.samp,
  thin = n.thin,
  burnin = n.burnin,
  inits = inits,
  method ='rjparallel',
  summarise = FALSE)



TA269_rmst_weibull <- export_outcomes(fit.km = TA269.km,mod.cp_TA269, data_jags = data_jags,
                              pred_mod_mat = mat,
                              pred_time= seq(0.01,50, by = 1),
                              inc_con_km = TRUE,
                              plt_title = "Weibull 1 change-point model (Common Hazards after Change-point)",
                              xlab_title = "Time (Months)",
                              pathway_export= gen_pathway(paste0(pathway,"TA269/Weibull/")),
                              df.origin =  TA269_df_all)
rm(mod.cp_TA269)

#Now assume constant hazard model

data_jags$beta_2_ind <- matrix(0, byrow = FALSE, ncol = 2, nrow = 2)

mod.cp_TA269_expo <- runjags::run.jags(
  model = get("jags.piecewise_wei_chng"),
  data = data_jags,
  n.chains = 2,
  monitor = c("cp", "beta_cp","beta_cp_anc", "total_loglik","loglik"),
  sample=n.samp,
  thin = n.thin,
  burnin = n.burnin,
  inits = inits,
  method ='rjparallel',
  summarise = FALSE)

TA269_rmst_expo <- export_outcomes(fit.km = TA269.km,mod.cp_TA269_expo, data_jags = data_jags,
                                      pred_mod_mat = mat,
                                      pred_time= seq(0.01,1, by = .1),
                                      inc_con_km = TRUE,
                                      plt_title = "Exponential 1 change-point model (Common Hazards after Change-point)",
                                      xlab_title = "Time (Months)",
                                      pathway_export= gen_pathway(paste0(pathway,"TA269/Exponential/")),
                                      df.origin =  TA269_df_all)
rm(mod.cp_TA269_expo)



#Output results

rmst_TA347 <- read.csv(paste0(pathway,"TA347/rmst.csv"))[,-1]

print(xtable::xtable(rmst_TA347 %>% mutate_if(is.numeric,round, digits = 2) %>% arrange(WAIC)),
             include.rownames=FALSE)



rmst_TA269_expo <- read.csv(paste0(pathway,"TA269/Exponential/rmst.csv"))[,-1]
rmst_TA269_expo[1,1] <- "Change-point Exponential"

rmst_TA269_weibull <- read.csv(paste0(pathway,"TA269/Weibull/rmst.csv"))[,-1]
rmst_TA269_weibull[1,1] <- "Change-point Weibull"

rmst_TA269 <- rbind(rmst_TA269_expo,rmst_TA269_weibull[1,] )

print(xtable::xtable(rmst_TA269 %>% mutate_if(is.numeric,round, digits = 2) %>% arrange(WAIC)),
      include.rownames=FALSE)


rmst_E190 <- read.csv(file = paste0(pathway,"RMST_scenario_all.csv"))[,-1]

print(xtable::xtable(rmst_E190 %>% mutate_if(is.numeric,round, digits = 2) %>% arrange(WAIC)),
      include.rownames=FALSE)
