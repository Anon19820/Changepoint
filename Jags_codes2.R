
library("loo")
library("flexsurv")
library("ggmcmc")
library("coda")

inits <- function(){
  list(unif = runif(data_jags$N_CP,max = data_jags$MAXX/data_jags$N_CP)
  )
}
gen_scaled_dens <- function(cp, scale = NULL){

  dens <- density(cp)
  if(!is.null(scale)){
    scaled <- max(dens$y)/scale
  }else{
    scaled <- 1
  }

  return(data.frame(x = dens$x, y = dens$y/scaled))
}

# WAIC gives almost exact same results as AIC.
# library("expertsurv")
# require("dplyr")
# param_expert_example1 <- list()
# 
# param_expert_example1[[1]] <- data.frame(dist = c("beta"),
#                                          wi = c(1), # Ensure Weights sum to 1
#                                          param1 = c(1),
#                                          param2 = c(1),
#                                          param3 = c(NA))
# 
# 
# timepoint_expert <- 14
# 
# data2 <- data %>% rename(status = censored) %>% mutate(time2 = ifelse(time > 10, 10, time),
#                                                        status2 = ifelse(time> 10, 0, status))
# 
# example1  <- fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
#                                distr=c("wph", "gomp"),
#                                method="hmc",
#                                pool_type = "log pool", 
#                                opinion_type = "survival",
#                                times_expert = timepoint_expert, 
#                                param_expert = param_expert_example1)
# 
# 
# example1.mle  <- fit.models.expert(formula=Surv(time2,status2)~1,data=data2,
#                                    distr=c("wph", "gomp"),
#                                    method="mle",
#                                    pool_type = "log pool", 
#                                    opinion_type = "survival",
#                                    times_expert = timepoint_expert, 
#                                    param_expert = param_expert_example1)
# 
# model.fit.plot(example1.mle, type = "aic")
# 
# model.fit.plot(example1, type = "waic")

get_rmst_all <- function(df_all, rmst_diff_cp, time_vec, pathway_export = NULL ){


  t_max <- max(time_vec)

  df_res <- data.frame(Model = c( "Changepoint"),
                       Est= c(rmst_diff_cp), WAIC = NA)
  param_mods <- c("exp","weibullPH",  "lnorm", "gamma",
                  "gompertz", "llogis","gengamma", "rps", "rps-anc")
  param_mods_full <- c("Exponential", "Weibull", "Log-Normal", "Gamma", "Gompertz", "Log-Logistic", "Generalized-Gamma", "Royston Parmar (PH)",
                       "Royston Parmar (non-PH)")


  for(p in 1:length(param_mods)){

    if(param_mods[p] == "rps"){

      mle.mod <- try(flexsurvspline(Surv(time, status)~as.factor(arm),
                                    data = df_all,  k = 1), silent = TRUE)

    }else if(param_mods[p] == "rps-anc"){

      mle.mod <- try(flexsurvspline(Surv(time, status)~as.factor(arm),
                                    data = df_all,  k = 1, anc =list(gamma1 = ~ as.factor(arm))), silent = TRUE)

    }else{
      mle.mod <- flexsurvreg(Surv(time, status)~as.factor(arm),
                             data = df_all, dist = param_mods[p],)
    }

    if(class(mle.mod) != "try-error"){
      rmst_diff_mod <- summary(mle.mod, t = t_max, type = "rmst")
      df_res <- rbind(df_res,c(param_mods_full[p], rmst_diff_mod[[1]]$est- rmst_diff_mod[[2]]$est,
                               AIC(mle.mod)))
      
      if(!is.null(pathway_export)){
        png(paste0(pathway_export, param_mods_full[p], ".png"))
        
        plot(mle.mod, t = seq(0:t_max), xlim = c(0, t_max), main = paste0("Model : ", param_mods_full[p]))
        
        dev.off()
      }

    }else{

      df_res <- rbind(df_res,c(param_mods_full[p], NA,NA))

    }

  }

  df_res$Est <- round(as.numeric(df_res$Est), digits = 2)
  return(df_res)

}

export_outcomes <- function(fit.km, mod.cp, data_jags, pred_mod_mat,
                            pred_time= seq(0,10, by = 0.25),
                            inc_con_km = FALSE,
                            inc_cred = FALSE,
                            interval_width =0.95,
                            plt_title = NULL,
                            xlab_title = NULL,
                            pathway_export= getwd(),
                            df.origin = NULL){

  n_last <- 1                                # Specify number of characters to extract

  last_letter<-   substr(pathway_export, nchar(pathway_export) - n_last + 1, nchar(pathway_export))

  if(last_letter != "/"){
    pathway_export <- paste0(pathway_export,"/")
  }

  n_pred <- length(pred_time)
  pred_mod_mat_rep<- do.call("rbind", replicate(n_pred, pred_mod_mat, simplify = FALSE))
  t_pred <- rep(pred_time, by = nrow(pred_mod_mat))

  pred_mod_mat_rep <-  pred_mod_mat_rep[order(pred_mod_mat_rep[,2]),]

  pred_mat_final <- cbind(tstart = 0, t_pred,status= 1,pred_mod_mat_rep)

  s <- summary(fit.km, times = fit.km$time, extend = T)

  plot_data <- tibble(
    'time' = s$time,
    'n.risk' = s$n.risk,
    'n.event' = s$n.event,
    'n.censor' = s$n.censor,
    'Survival' = s$surv,
    'std.error' = s$std.err,
    'Treatment' = s$strata,
    'lower' = s$lower,
    'upper' = s$upper,
  )

  plot_data$Treatment<-  gsub("Treatment=", "" ,plot_data$Treatment)

  #browser()
  # Extract Model Diagnostics.

  index_save <- grep("cp|lambda_wane",colnames(mod.cp$mcmc[[1]]))

  colnames_save <- index_save[apply(mod.cp$mcmc[[1]][,index_save],2, function(x){ifelse(length(unique(x)) ==1, FALSE,TRUE)})]

  mcmc_list_eval <- lapply(mod.cp$mcmc,function(x){x[,colnames_save]})

  lapply(mcmc_list_eval, as.numeric)
  #For some reason the output of runjags is a charachter need to run as.mcmc.list
  ggmcmc(ggs(coda::as.mcmc.list(mcmc_list_eval)), file = paste0(pathway_export,"MCMC_diag.pdf"))

    mat_res <- as.matrix(mod.cp$mcmc)
  param_names <- colnames(mat_res)
  loglik_mat <- mat_res[,grep("loglik",param_names)]
  loglik_mat <- loglik_mat[,-grep("total_loglik",colnames(loglik_mat))]
  mean_cp <- mean(mat_res[,grep("cp[2]",param_names,fixed = T)])
  waic <- loo::waic(loglik_mat)[["estimates"]][3, 1]
  pml <- -2 * sum(log(nrow(loglik_mat)/colSums(1/exp(loglik_mat))))

  nrow_beta <- nrow(data_jags$beta_1_ind)
  ncol_beta <- ncol(data_jags$beta_1_ind)
  nsims <- nrow(mat_res)

  beta_1_array <- array(as.numeric(t(mat_res[,grep( "beta_cp[",colnames(mat_res), fixed = TRUE)])), c(nrow_beta,ncol_beta,nsims))
  beta_2_array <- array(as.numeric(t(mat_res[,grep( "beta_cp_anc[",colnames(mat_res), fixed = TRUE)])), c(nrow_beta,ncol_beta,nsims))
  cp_vec <- mat_res[,"cp[2]"]

  wane_flag <- any(colnames(mat_res) %in% "lambda_wane")

  if(wane_flag){
    lambda_wane <- mat_res[,"lambda_wane"]
  }



  haz_mat <- St_mat<- matrix(nrow = nrow(pred_mat_final),ncol =  nsims)
  #browser()

  for(i in 1:nsims){
    if(wane_flag){
      St_mat[,i] <- Ccp_LL(as.matrix(pred_mat_final),c( cp_vec[i], Inf), beta_1_array[,,i], beta_2_array[,,i], type = 2,trt_wane = 1, lambda_wane= lambda_wane[i]) #survival
      haz_mat[,i] <- Ccp_LL(as.matrix(pred_mat_final),c( cp_vec[i], Inf), beta_1_array[,,i], beta_2_array[,,i], type = 3,trt_wane = 1, lambda_wane= lambda_wane[i]) #survival

    }else{
      St_mat[,i] <- Ccp_LL(as.matrix(pred_mat_final),c( cp_vec[i], Inf), beta_1_array[,,i], beta_2_array[,,i], type = 2,trt_wane = 0, lambda_wane= 0) #survival
      haz_mat[,i] <- Ccp_LL(as.matrix(pred_mat_final),c( cp_vec[i], Inf), beta_1_array[,,i], beta_2_array[,,i], type = 3,trt_wane = 0, lambda_wane= 0) #survival

    }

  }
  lower_int <- (1-interval_width)/2
  upper_int <- interval_width + (1-interval_width)/2

  St_median <- apply(St_mat,1,quantile , probs =  c(lower_int, 0.5, upper_int), na.rm = TRUE) %>% t()
  haz_median <- apply(haz_mat,1,quantile , probs =  c(lower_int, 0.5, upper_int), na.rm = TRUE) %>% t()
  hr_mat <- haz_mat[pred_mat_final[,5]==1,]/haz_mat[pred_mat_final[,5]==0,]
  hr_median <-  apply(hr_mat,1,quantile , probs =  c(lower_int, 0.5, upper_int), na.rm = TRUE) %>% t()
  colnames(hr_median) <- c("lower", "Hazard", "upper")

  mat_predict<- cbind(pred_mat_final[,c(2,5)],St_median,haz_median) %>% data.frame()

  # mat_predict[,2]<- as.numeric(as.character(mat_predict[,2]))
  #
  # mat_predict <- mat_predict %>%
  #   arrange(desc(.[[2]]))

  mat_predict[,2] <- factor(mat_predict[,2], labels  = unique(plot_data$Treatment))
  colnames(mat_predict) <- c("time", "Treatment", "St_lower","Survival","St_upper",
                             "Haz_lower", "Hazard", "Haz_upper")

#browser()
  p_St <- ggplot(plot_data, aes(x =time, y = Survival, color = Treatment))+
    geom_step( size = 1)+
    ylim(c(0,1))+
    xlab(xlab_title)+
    ggtitle(plt_title)+
    geom_line(data = mat_predict,
              aes(x = time, y = Survival,color = Treatment))+
    # geom_point(data = data.frame(time = mean.cp,
    #                              Survival = mean(pull(Surv_extrap,Survival)[which(dev == min(dev))])),
    #            aes(x = time, Survival), shape = 23, fill = "green",
    #            color = "darkred", size = 4, inherit.aes = F, stroke  = 1.5)+
    geom_ribbon(data = gen_scaled_dens( mat_res[,"cp[2]"], scale = 0.4),
                aes(x = x,ymin = 0,ymax = y),
                inherit.aes = FALSE, color="darkred", fill="lightgreen", alpha = 0.4)+

    theme_bw()


 if(inc_con_km){
   p_St <- p_St+
     geom_line(data = plot_data,
               aes(x = time, y = lower,color = Treatment), linetype = "dashed")+
     geom_line(data = plot_data,
               aes(x = time, y = upper,color = Treatment), linetype = "dashed")
 }


  p_haz <- ggplot(data = mat_predict,
                  aes(x = time, y = Hazard,color = Treatment))+
    geom_line()+
    xlab(xlab_title)+
    ggtitle(plt_title)+
    theme_bw()


  if(inc_cred){
    p_St <- p_St+
      geom_line(data = mat_predict,
                aes(x = time, y = St_lower,color = Treatment), linetype = 2)+
      geom_line(data = mat_predict,
                aes(x = time, y = St_upper,color = Treatment), linetype = 2)

    p_haz <- p_haz+
      geom_line(data = mat_predict,
                aes(x = time, y = Haz_lower,color = Treatment), linetype = 2)+
      geom_line(data = mat_predict,
                aes(x = time, y = Haz_upper,color = Treatment), linetype = 2)
  }

  p_HR <-  ggplot(data = cbind(time = pred_time,hr_median) %>% data.frame(),
                  aes(x = time, y = Hazard))+
    geom_line()+
    geom_line(aes(x = time, y = upper), linetype = 2)+
    geom_line(aes(x = time, y = lower), linetype = 2)+
    ylim(0,NA)+
    xlab(xlab_title)+
    ggtitle(plt_title)+
    theme_bw()



  ggsave(filename = gen_pathway(paste0(pathway_export,"hr_plot.png")), plot = print(p_HR), width = 10, height = 5, units = 'in')
  ggsave(filename = gen_pathway(paste0(pathway_export,"St_plot.png")), plot = p_St, width = 10, height = 5, units = 'in')
  ggsave(filename = gen_pathway(paste0(pathway_export,"haz_plot.png")), plot = p_haz, width = 10, height = 5, units = 'in')

  mat_predict_surv <- mat_predict[,c("time","Treatment","Survival")]

  rmst_cp_df <- tidyr::spread(mat_predict_surv, Treatment, Survival) %>% mutate(diff = .[[3]]-.[[2]])

  rmst_diff_cp <- sfsmisc::integrate.xy(x = rmst_cp_df$time,rmst_cp_df$diff)


  rmst_all <- get_rmst_all(df_all =df.origin,
               rmst_diff_cp = rmst_diff_cp,
               time_vec =rmst_cp_df$time,
               pathway_export = pathway_export )

  rmst_all[1,3] <- waic

  rmst_all[,3] <- as.numeric(rmst_all[,3] )

  write.csv(x = rmst_all,file= gen_pathway(paste0(pathway_export,"rmst.csv")))

  return(rmst_all)


}


jags.piecewise_wei_chng <-"

data{
  for(i in 1:N){
    zero[i] <- 0
  }
  zero_prior <- 0
  D_cp <- dim(beta_1_ind)
  ncovar_cp <- D_cp[1]
}

model {

  #Constant for zerors trick
  C <- 10000

  #Prior for change-point locations - See Chapell
  cp[1] = 0  # mcp helper value. Should be zero
  cp[N_CP+2] = 9999  # mcp helper value. very large number

  for(k in 1:N_CP){
    unif[k] ~ dunif(0.001, min(max(time),MAXX))
  }
  cp_x[2:(N_CP+1)] <- sort(unif)

  for(i in 2:(N_CP+1)){
    cp[i] <- cp_x[i]*equals(cp_fixed,0) + cp_fix[i]*equals(cp_fixed,1)
  }

  for(i in 2:(N_CP+1)){
    diff[i-1] <- cp[i]-cp[i-1]
  }

  log_prior <- loggam(2*N_CP +2) + sum(log(diff)) - ((2*N_CP +1)*log(MAXX))
  zero_prior ~ dpois(C - log_prior)

  #Prior for the model parameters

  for(i in 1:2){
    sd[i]  <- 5 # ~ dunif(0.2,5)
    prec[i] <- pow(sd[i], -2)
  }

  for(k in 1:(N_CP+1)){
    for(j in 1:ncovar_cp){
      beta_cp_x[j,k] ~ dnorm(0,prec[1])
      beta_cp_anc_x[j,k] ~ dnorm(0,prec[1])
      #if 2 we have a common value of the covariate across the time-points;
      #otherwise we have different change-point models.
      beta_cp[j,k] <- beta_cp_x[j,1]*equals(beta_1_ind[j,k],2) + beta_cp_x[j,k]*(1-equals(beta_1_ind[j,k],2))
      beta_cp_anc[j,k] <- beta_cp_anc_x[j,1]*equals(beta_2_ind[j,k],2) + beta_cp_anc_x[j,k]*(1-equals(beta_2_ind[j,k],2))
    }

    param_1_ln[1:N,k] <-  X_mat %*% beta_cp[1:ncovar_cp,k]
    param_2_ln[1:N,k] <-  X_mat %*% beta_cp_anc[1:ncovar_cp,k]

  }


  # Model and likelihood
  for (i in 1:N) {

    for(k in 1:(N_CP+1)){
      #variable which gives the difference between the two intervals if time[i]>cp[k+1]
      #(i.e. cp[k+1] -cp[k]) or time between time[i] and cp[k]
      X[i,k] = max(min(time[i], cp[k+1]) - cp[k],0)

      #Indicator variable which highlights which interval time is in
      X_ind[i,k] = step(time[i]-cp[k])*step(cp[k+1]-time[i])

      param_1[i,k] <- exp(param_1_ln[i,k])
      param_2[i,k] <- exp(param_2_ln[i,k])

      log_haz_seg[i,k] <-  log(param_2[i,k]*param_1[i,k]*pow(time[i],param_2[i,k]-1))*X_ind[i,k]
      cum_haz_seg[i,k] <- param_1[i,k]*pow(X[i,k]+cp[k],param_2[i,k]) -  param_1[i,k]*pow(cp[k],param_2[i,k])

    }

    log_haz[i] <- sum(log_haz_seg[i,])
    cum_haz[i] <- sum(cum_haz_seg[i,])

    loglik[i] = status[i]*log_haz[i] - cum_haz[i]
    zero[i] ~ dpois(C - loglik[i])

  }

  total_loglik <- sum(loglik)
}

"


jags.piecewise_wei_wane <-"

  data{
    for(i in 1:N){
    zero[i] <- 0
    }
    zero_prior <- 0
    D_cp <- dim(beta_1_ind)
    ncovar_cp <- D_cp[1]
  }

  model {
  #Constant for zerors trick
  C <- 10000

 #Prior for change-point locations - See Chapell
  cp[1] = 0  # mcp helper value. Should be zero
  cp[N_CP+2] = 9999  # mcp helper value. very large number

  for(k in 1:N_CP){
   unif[k] ~ dunif(0.001, max(time))
  }
  cp[2:(N_CP+1)] <- sort(unif)

  for(i in 2:(N_CP+1)){
  diff[i-1] <- cp[i]-cp[i-1]
  }

  log_prior <- loggam(2*N_CP +2) + sum(log(diff)) - ((2*N_CP +1)*log(MAXX))
  zero_prior ~ dpois(C - log_prior)


  #Prior for the model parameters

  for(i in 1:2){
    sd[i]  <- 5 # ~ dunif(0.2,5)
    prec[i] <- pow(sd[i], -2)
  }


  for(k in 1:(N_CP+1)){
    for(j in 1:ncovar_cp){
      beta_cp_x[j,k] ~ dnorm(0,prec[1])
      beta_cp_anc_x[j,k] ~ dnorm(0,prec[1])
      #if 2 we have a common value of the covariate across the time-points;
      #otherwise we have different change-point models.
      beta_cp[j,k] <- beta_cp_x[j,1]*equals(beta_1_ind[j,k],2) + beta_cp_x[j,k]*(1-equals(beta_1_ind[j,k],2))
      beta_cp_anc[j,k] <- beta_cp_anc_x[j,1]*equals(beta_2_ind[j,k],2) + beta_cp_anc_x[j,k]*(1-equals(beta_2_ind[j,k],2))
    }

    param_1_ln[1:N,k] <-  X_mat %*% beta_cp[1:ncovar_cp,k]  #beta_cp[trt_flg,ncp] must be zero
    param_2_ln[1:N,k] <-  X_mat %*% beta_cp_anc[1:ncovar_cp,k]

  }
  ln_lambda_wane ~ dunif(-3,1.1)
  lambda_wane <- exp(ln_lambda_wane)

  initial_HR <-   exp(beta_cp[trt_ind,N_CP])

  # Model and likelihood
  for (i in 1:N) {

    for(k in 1:N_CP){
      #variable which gives the difference between the two intervals if time[i]>cp[k+1]
      #(i.e. cp[k+1] -cp[k]) or time between time[i] and cp[k]
      X[i,k] = max(min(time[i], cp[k+1]) - cp[k],0)

      #Indicator variable which highlights which interval time is in
      X_ind[i,k] = step(time[i]-cp[k])*step(cp[k+1]-time[i])

      param_1[i,k] <- exp(param_1_ln[i,k])
      param_2[i,k] <- exp(param_2_ln[i,k])

      log_haz_seg[i,k] <-  log(param_2[i,k]*param_1[i,k]*pow(time[i],param_2[i,k]-1))*X_ind[i,k]
      cum_haz_seg[i,k] <- param_1[i,k]*pow(X[i,k]+cp[k],param_2[i,k]) -  param_1[i,k]*pow(cp[k],param_2[i,k])

    }


  for(k in (N_CP+1):(N_CP+1)){ #Treatment waning interval - Last interval

    X[i,k] = max(min(time[i], cp[k+1]) - cp[k],0)
    X_ind[i,k] = step(time[i]-cp[k])*step(cp[k+1]-time[i])

    HR_wane[i] <- 1-(1-initial_HR)*exp(-lambda_wane*X[i,k])

    param_1[i,k] <- exp(param_1_ln[i,k])
    param_2[i,k] <- exp(param_2_ln[i,k])


    #log_haz_stn[i] <-  log(param_2[i,k]*param_1[i,k]*pow(time[i],param_2[i,k]-1))*X_ind[i,k]
    #(log(haz)+log(haz_wane)) == log(haz*haz_waner)
    log_haz_seg[i,k] <-  (log(param_2[i,k]*param_1[i,k]*pow(time[i],param_2[i,k]-1)) + log(HR_wane[i])*equals(X_mat[i,trt_ind],1))*X_ind[i,k]

    #Need to integrate this
    #Treatment waning... a is shape; m is scale
    #(amt^{a-1})*(1-(1-r)exp(-wt)
    #drop am -- need to add it back in at end
    #(t^{a-1})*(1-(1-r)exp(-w(t-q))
    #(1/lambda2-exp(-t*lambda2)/lambda2)*r*lambda

    #incomplete_gamma_upper<- function(s,x){
    #  exp(lgamma(s))*(1-pgamma(x,s,1))
    #}

    #int_1 <- (t_0^a)/a-(incomplete_gamma_upper(a,hr_wane_lambda*t_0)*(r-1)*exp(q*hr_wane_lambda)*(t_0^a))/((hr_wane_lambda*t_0)^a)
    #int_2 <- (t_1^a)/a-(incomplete_gamma_upper(a,hr_wane_lambda*t_1)*(r-1)*exp(q*hr_wane_lambda)*(t_1^a))/((hr_wane_lambda*t_1)^a)

    upper_inc_gamma1[i] <- exp(loggam(param_2[i,k]))*(1-pgamma(lambda_wane*time[i],param_2[i,k],1))
    upper_inc_gamma2[i] <- exp(loggam(param_2[i,k]))*(1-pgamma(lambda_wane*cp[k],param_2[i,k],1))
    upper_int[i] <- (time[i]^param_2[i,k])/param_2[i,k]-(upper_inc_gamma1[i]*(initial_HR-1)*exp(cp[k]*lambda_wane)*(time[i]^param_2[i,k]))/((lambda_wane*time[i])^param_2[i,k])
    lower_int[i] <- (cp[k]^param_2[i,k])/param_2[i,k]-(upper_inc_gamma2[i]*(initial_HR-1)*exp(cp[k]*lambda_wane)*(cp[k]^param_2[i,k]))/((lambda_wane*cp[k])^param_2[i,k])

    cum_haz_stn[i] <- param_1[i,k]*pow(X[i,k]+cp[k],param_2[i,k]) -  param_1[i,k]*pow(cp[k],param_2[i,k])
    cum_haz_wane[i] <- ((upper_int[i] - lower_int[i])*param_2[i,k]*param_1[i,k])*step(time[i]-cp[k]) #Has to be treatment and cp < time
    cum_haz_seg[i,k] <- cum_haz_wane[i]*equals(X_mat[i,trt_ind],1) +  cum_haz_stn[i]*equals(X_mat[i,trt_ind],0)


 }

  log_haz[i] <- sum(log_haz_seg[i,])
  cum_haz[i] <- sum(cum_haz_seg[i,])

  loglik[i] = status[i]*log_haz[i] - cum_haz[i]
  zero[i] ~ dpois(C - loglik[i])

  }


  total_loglik <- sum(loglik)
  }

"
# Nimble Cox Models

LL_cox_ph <-  function(params_beta,params_cp,  data){
  library("survival")
  params_cp <- head(params_cp, n= -1)
  data <- as.data.frame(data)
  colnames(data) <- c("time","status","trt", "age")
  data$trt <- as.factor(data$trt)

  data_split <- survival::survSplit(data,cut=params_cp,
                          end="time",start="start",
                          event="status", episode="period")

  cox_fit <- coxph(Surv(start, time, status) ~
                     I((trt==2)&(period==1)) + I((trt==2)&(period==2)) + age,
                   data=data_split,
                   init =params_beta, #
                   control = coxph.control(iter.max = 0))
  #model.matrix(cox_fit)

  return(cox_fit$loglik[2])

}

R_LL_cox_ph <- nimbleRcall(function(params_beta = double(1),
                                    params_cp = double(1),
                                    data = double(2)){},
                           Rfun = 'LL_cox_ph',
                           returnType = double(0))


mod <- nimbleCode({

  for(i in 1:n_param_covar){
    params_beta[i] ~ dunif(-10,10)
    # params_beta[i] <- 1
    # params_cp[i] <- 10
  }


  for(i in 1:n_param){
    params_cp[i] ~ dunif(0,max_cp)
    # params_beta[i] <- 1
    # params_cp[i] <- 10
  }

  C <- 1000

  LL <- R_LL_cox_ph(params_beta[1:n_param_covar],
                    params_cp[1:n_param],
                    data[1:df_nrow,1:df_ncol])

  zero ~ dpois(-LL + C)


})
