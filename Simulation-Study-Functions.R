gen_cp_df <- function(n.samp = 100, plot = FALSE, t_cens, St_trt, St_comp){
  
  sim_unif1 <- runif(n.samp)#seq(0,1, length.out = n.samp_size)#runif(n.samp)
  
  sim_unif2 <- runif(n.samp) #sim_unif1 #runif(n.samp)
  
  time_trt_samp <- pred_time_grid[sapply(sim_unif1, function(x){which.min(abs(St_trt-x))})]
  time_comp_samp <- pred_time_grid[sapply(sim_unif2, function(x){which.min(abs(St_comp-x))})]
  
  df_all <- data.frame(time = c(time_trt_samp,time_comp_samp),
                       arm = c(rep(1, length(time_trt_samp)),
                               rep(0, length(time_comp_samp)))) 
  df_all <- dplyr::mutate(df_all,
           time_event = time,
           status_true = 1,
           status = ifelse(time > t_cens, 0,1),
           time = ifelse(time > t_cens, t_cens, time))
  
  df_all$time <- df_all$time +0.0001
  if(plot){
    km.all <- survfit(survival::Surv(time, status)~arm, data = df_all)
    plot(km.all)
  }
  
  return(df_all)
  
  
}
#To parrallize

make_data_jags <- function(df, beta_1_ind, beta_2_ind){
  data_jags <- list()
  data_jags$N <- nrow(df)
  data_jags$time <- df$time
  data_jags$status <- df$status
  data_jags$MAXX <- max(df$time)
  data_jags$N_CP <- 1
  data_jags$X_mat<- as.matrix(cbind(1,df$arm))
  data_jags$beta_1_ind <- beta_1_ind
  data_jags$beta_2_ind <- beta_2_ind
  
  #If you wanted a fixed change-point... Set to zero if you want the change-point to be estimated
  data_jags$cp_fix <- c(0, 1, 5)
  data_jags$cp_fixed <- 0
  
  data_jags
}
# inits <- function(data_jags){
#   list(unif = runif(data_jags$N_CP,max = data_jags$MAXX/data_jags$N_CP)
#   )
# }

inits <- function(true_cp, true_beta_cp, true_beta_cp_anc, true_lambda_wane){
  init_list <- list(unif = true_cp,
                    beta_cp_x = true_beta_cp,
                    beta_cp_anc_x = true_beta_cp_anc)
  
  if(true_lambda_wane != 0){
    init_list$ln_lambda_wane = log(true_lambda_wane)
  }
  
  list(init_list)
}




cp_LL_nim <- function (X, cut, beta_1, beta_2, type, trt_wane, lambda_wane){
  X_res <- survSplit(X, cut)
  n <- nimble::nimDim(X)[1]
  nrow_X = nimble::nimDim(X_res)[1]
  ncol_X = nimble::nimDim(X_res)[2]
  n_int = nimble::nimDim(beta_1)[2]
  param_1 <- nimble::nimMatrix(nrow = nrow_X, ncol = n_int)
  param_2 <- nimble::nimMatrix(nrow = nrow_X, ncol = n_int)
  param_1_final <- nimble::nimNumeric(nrow_X)
  param_2_final <- nimble::nimNumeric(nrow_X)
  res <- nimble::nimNumeric(nrow_X)
  for (i in 1:n_int) {
    param_1[, i] <- X_res[, 6:ncol_X] %*% beta_1[, i]
    param_2[, i] <- X_res[, 6:ncol_X] %*% beta_2[, i]
  }
  for (i in 1:nrow_X) {
    param_1_final[i] <- exp(param_1[i, X_res[i, 4]])
    param_2_final[i] <- exp(param_2[i, X_res[i, 4]])
    if (type == 1) {
      res[i] <- hweibullPH(X_res[i, 2], param_1_final[i], 
                           param_2_final[i], 1) * X_res[i, 3] - (HweibullPH(X_res[i, 
                                                                                  2], param_1_final[i], param_2_final[i], 0) - 
                                                                   HweibullPH(X_res[i, 1], param_1_final[i], param_2_final[i], 
                                                                              0))
    }
    else if (type == 2) {
      if (trt_wane == 1 & X_res[i, 7] == 1 & X_res[i, 4] == 
          n_int) {
        initial_HR <- exp(beta_1[2, n_int - 1])
        t_end <- X_res[i, 2]
        t_start <- X_res[i, 1]
        upper_inc_gamma1 <- exp(loggam(param_2_final[i])) * 
          (1 - pgamma(lambda_wane * t_end, param_2_final[i], 
                      1))
        upper_inc_gamma2 <- exp(loggam(param_2_final[i])) * 
          (1 - pgamma(lambda_wane * t_start, param_2_final[i], 
                      1))
        upper_int <- (t_end^param_2_final[i])/param_2_final[i] - 
          (upper_inc_gamma1 * (initial_HR - 1) * exp(t_start * 
                                                       lambda_wane) * (t_end^param_2_final[i]))/((lambda_wane * 
                                                                                                    t_end)^param_2_final[i])
        lower_int <- (t_start^param_2_final[i])/param_2_final[i] - 
          (upper_inc_gamma2 * (initial_HR - 1) * exp(t_start * 
                                                       lambda_wane) * (t_start^param_2_final[i]))/((lambda_wane * 
                                                                                                      t_start)^param_2_final[i])
        res[i] <- ((upper_int - lower_int) * param_2_final[i] * 
                     param_1_final[i])
      }
      else {
        res[i] <- (HweibullPH(X_res[i, 2], param_1_final[i], 
                              param_2_final[i], 0) - HweibullPH(X_res[i, 
                                                                      1], param_1_final[i], param_2_final[i], 0))
      }
    }
    else {
      res[i] <- hweibullPH(X_res[i, 2], param_1_final[i], 
                           param_2_final[i], 0) * X_res[i, 3]
      if (trt_wane == 1 & X_res[i, 7] == 1 & X_res[i, 4] == 
          n_int) {
        t_end <- X_res[i, 2]
        t_start <- X_res[i, 1]
        initial_HR <- exp(beta_1[2, n_int - 1])
        HR_wane <- 1 - (1 - initial_HR) * exp(-lambda_wane * 
                                                (t_end - t_start))
        res[i] <- res[i] * HR_wane
      }
    }
  }
  res_final <- nimble::nimNumeric(n)
  for (i in 1:n) {
    res_final[i] <- sum(res[X_res[, 5] == i])
  }
  if (type == 2) {
    res_final <- exp(-res_final)
  }
  return(res_final)
}

export_outcomes_sim_study <- function(mod.cp, 
                                      data_jags, 
                                      pred_mod_mat,
                                      pred_time= seq(0,15, by = 0.1),
                                      interval_width =0.95,
                                      df.origin = NULL,
                                      true_trt_diff_vec,
                                      path_export){
  
  library("nimble")
  n_pred <- length(pred_time)
  pred_mod_mat_rep<- do.call("rbind", replicate(n_pred, pred_mod_mat, simplify = FALSE))
  t_pred <- rep(pred_time, by = nrow(pred_mod_mat))
  
  pred_mod_mat_rep <-  pred_mod_mat_rep[order(pred_mod_mat_rep[,2]),]
  
  pred_mat_final <- cbind(tstart = 0, t_pred,status= 1,pred_mod_mat_rep)
  
  
  index_save <- grep("cp|lambda_wane",colnames(mod.cp$mcmc[[1]]))

  colnames_save <- index_save[apply(mod.cp$mcmc[[1]][,index_save],2, function(x){ifelse(length(unique(round(x, digits = 5))) ==1, FALSE,TRUE)})]
  
  mcmc_list_eval <- lapply(mod.cp$mcmc,function(x){x[,colnames_save]})
  mcmc_med_parm   <-  apply(do.call(rbind,mcmc_list_eval),2, median)
  
  lapply(mcmc_list_eval, as.numeric)
  #For some reason the output of runjags is a charachter need to run as.mcmc.list
  #ggmcmc(ggs(coda::as.mcmc.list(mcmc_list_eval)), file = paste0(pathway_export,"MCMC_diag.pdf"))
  
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
  
  
  
  St_mat<- matrix(nrow = nrow(pred_mat_final),ncol =  nsims)
  #browser()
  
  for(i in 1:nsims){
    if(wane_flag){
      St_mat[,i] <- cp_LL_nim(as.matrix(pred_mat_final),c( cp_vec[i], Inf), beta_1_array[,,i], beta_2_array[,,i], type = 2,trt_wane = 1, lambda_wane= lambda_wane[i]) #survival
    }else{
      St_mat[,i] <- cp_LL_nim(as.matrix(pred_mat_final),c( cp_vec[i], Inf), beta_1_array[,,i], beta_2_array[,,i], type = 2,trt_wane = 0, lambda_wane= 0) #survival
    }
    
  }
  lower_int <- (1-interval_width)/2
  upper_int <- interval_width + (1-interval_width)/2
  
  St_median <- t(apply(St_mat,1,quantile , probs =  c(lower_int, 0.5, upper_int), na.rm = TRUE))
  
  mat_predict<- data.frame(cbind(pred_mat_final[,c(2,5)],St_median))
  
  # mat_predict[,2]<- as.numeric(as.character(mat_predict[,2]))
  #
  # mat_predict <- mat_predict %>%
  #   arrange(desc(.[[2]]))
  
  mat_predict[,2] <-  factor(mat_predict[,2], labels  = c("Comparator","Treatment"))
  colnames(mat_predict) <- c("time", "arm", "St_lower","Survival","St_upper")
  
  mat_predict_surv <- mat_predict[,c("time","arm","Survival")]
  
  rmst_cp_df <- dplyr::mutate(tidyr::spread(mat_predict_surv, arm, Survival), diff = Treatment-Comparator)
  
  survKM<- survival::survfit(survival::Surv(time_event,status_true) ~ arm, data = df.origin)
  
  #if(FALSE){
    png(paste0(path_export,"Change-point.png"))
    plot(survKM, xlab = c(0,max(pred_time)))
    lines(x = rmst_cp_df$time,rmst_cp_df$Treatment ,col = "blue")
    lines(x = rmst_cp_df$time,rmst_cp_df$Comparator ,col = "red")
    dev.off()
  #}
  
  rmst_diff_cp <- sfsmisc::integrate.xy(x = rmst_cp_df$time,rmst_cp_df$diff)
  rmst_trt_cp <- sfsmisc::integrate.xy(x = rmst_cp_df$time,rmst_cp_df$Treatment)
  rmst_comp_cp <- sfsmisc::integrate.xy(x = rmst_cp_df$time,rmst_cp_df$Comparator)
  
  
  rmst_all <- get_rmst_all_sim(df_all =df.origin,
                           time_vec =rmst_cp_df$time,
                           pathway_export = NULL,
                           true_trt_diff = true_trt_diff_vec,
                           path_export = path_export)
  
  rmst_all_1 <- rmst_all[[1]]
  
  rmst_all_1[1,2:5] <- c(rmst_diff_cp,rmst_trt_cp,rmst_comp_cp, waic)
  
  #rmst_all[,3] <- as.numeric(rmst_all[,3] )
  
  rmst_all_1[,2:5] <- apply(rmst_all_1[,2:5], 2, as.numeric)
  
  int_Surv <- rmst_all[[2]]
  
  int_Surv[1,2] <- sfsmisc::integrate.xy(fx =  abs(true_trt_diff_vec - rmst_cp_df$diff), x = rmst_cp_df$time)
  
  int_Surv[,2] <- as.numeric(int_Surv[,2])
  
  #This is the difference between the KM and the predicted RMST
  rmst_diff <- sweep(as.matrix(rmst_all_1[-2,2:4]), 2, STATS = as.numeric(rmst_all_1[2,2:4]))
  rmst_diff <- data.frame(Model = rmst_all_1[-2,1],rmst_diff)
  
  return(list(rmst_all_1,rmst_diff,mcmc_med_parm,int_Surv))
  
  
}


km_AUC <- function(surv_boot_km.df, time_vec){
  surv_km_time <- rep(NA, length(time_vec))
  for (i in 1:length(time_vec)) {
    if (time_vec[i] < surv_boot_km.df$time[1]) {
      surv_km_time[i] <- 1
    }
    else if (time_vec[i] > surv_boot_km.df$time[nrow(surv_boot_km.df)]) {
      surv_km_time[i] <- NA
    }
    else {
      surv_km_time[i] <- surv_boot_km.df$survival[max(which(surv_boot_km.df$time <=
                                                              time_vec[i]))]
    }
  }
  t_eval <- !is.na(surv_km_time)
  AUC_true <- sfsmisc::integrate.xy(time_vec[t_eval], surv_km_time[t_eval])
  
  return(list(AUC_true, df_km = data.frame(time = time_vec[t_eval],surv_km_time[t_eval])))
}



get_rmst_all_sim <- function(df_all,  time_vec, pathway_export = NULL, true_trt_diff, path_export){
  
  n.boots <- 100
  #Get non-parametric estimator
  
  df_trt <-   dplyr::filter(df_all,arm == 1)
  df_comp <-   dplyr::filter(df_all,arm == 0)
  
  boot_ids <- replicate(n.boots,sample(nrow(df_trt), size = nrow(df_trt), replace = TRUE))
  
  AUC_true_comp <- AUC_true_trt <- rep(NA, n.boots)
  
  surv_km <- survival::survfit(survival::Surv(time_event,status_true) ~ arm, data = df_all)
  for (b in 1:n.boots) {
    
    surv_boot_km <- survival::survfit(survival::Surv(time_event, status_true) ~ arm,
                            data = df_all[c(boot_ids[,b],boot_ids[,b]+nrow(df_trt) ),])
    
    surv_boot_km.df <- data.frame(cbind(surv_boot_km[[c("time")]],
                                        surv_boot_km[[c("surv")]],
                                        surv_boot_km[[c("upper")]],
                                        surv_boot_km[[c("lower")]]))
    
    surv_boot_km.df$arm <- NA
    surv_boot_km.df$arm[1:surv_boot_km$strata[1]] <-0
    index_end <- surv_boot_km$strata[1]+1
    surv_boot_km.df$arm[index_end:(index_end+surv_boot_km$strata[2]-1)] <-1
    
    colnames(surv_boot_km.df) <- c("time", "survival", "upper",
                                   "lower", "arm")
    surv_boot_km_trt <-  dplyr::filter(surv_boot_km.df,arm == 1)
    surv_boot_km_comp <- dplyr::filter(surv_boot_km.df,arm == 0)
    
    res_trt <- km_AUC(surv_boot_km.df = surv_boot_km_trt, time_vec = time_vec)
    res_comp <- km_AUC(surv_boot_km.df = surv_boot_km_comp, time_vec = time_vec)
    AUC_true_trt[b]<-  res_trt[[1]]
    AUC_true_comp[b]<-  res_comp[[1]]
    
  }
  
  # parametric results
  
  t_max <- max(time_vec)
  
  df_res <- data.frame(Model = c("Changepoint", "KM-Data"),
                       Diff_Est= c(NA, mean(AUC_true_trt-AUC_true_comp)),
                       AUC_trt = c(NA, mean(AUC_true_trt)),
                       AUC_comp = c(NA, mean(AUC_true_comp)),
                       WAIC = c(NA,NA))
  
  df_res2 <- data.frame(Model = c("Changepoint"),
                       Int_Diff= c(NA))
  
  
  
  
  param_mods <- c("exp","weibullPH",  "lnorm", "gamma",
                  "gompertz", "llogis","gengamma", "rps", "rps-anc")
  param_mods_full <- c("Exponential", "Weibull", "Log-Normal", "Gamma", "Gompertz", "Log-Logistic", "Generalized-Gamma", "Royston Parmar (PH)",
                       "Royston Parmar (non-PH)")
  
  
  for(p in 1:length(param_mods)){
    
    if(param_mods[p] == "rps"){
      
      mle.mod <- try(flexsurv::flexsurvspline(survival::Surv(time, status)~as.factor(arm),
                                    data = df_all,  k = 1), silent = TRUE)
      
    }else if(param_mods[p] == "rps-anc"){
      
      mle.mod <- try(flexsurv::flexsurvspline(survival::Surv(time, status)~as.factor(arm),
                                    data = df_all,  k = 1, anc =list(gamma1 = ~ as.factor(arm))), silent = TRUE)
      
    }else{
      mle.mod <- flexsurv::flexsurvreg(survival::Surv(time, status)~as.factor(arm),
                             data = df_all, dist = param_mods[p],)
    }
    
    if(class(mle.mod) != "try-error"){
      rmst_diff_mod <- summary(mle.mod, t = t_max, type = "rmst")
      Surv_mod <- summary(mle.mod, t = time_vec)
      surv_diff_mod <- Surv_mod[[1]]$est -Surv_mod[[2]]$est
      
      Int_error<-  sfsmisc::integrate.xy( fx =  abs(true_trt_diff -surv_diff_mod), x = time_vec)
      
      
      
      df_res <- rbind(df_res,c(param_mods_full[p], 
                               rmst_diff_mod[[1]]$est- rmst_diff_mod[[2]]$est,
                               rmst_diff_mod[[1]]$est,
                               rmst_diff_mod[[2]]$est,
                               AIC(mle.mod)))
      
      df_res2 <- rbind(df_res2,c(param_mods_full[p], 
                                Int_error))
      
      
      if(!is.null(path_export)){
        png(paste0(path_export, param_mods_full[p], ".png"))
        
        plot(mle.mod, t = seq(0:t_max), xlim = c(0, t_max), main = paste0("Model : ", param_mods_full[p]))
        
        dev.off()
      }
      
    }else{
      
      df_res <- rbind(df_res,c(param_mods_full[p], NA,NA, NA,NA))
      df_res2 <- rbind(df_res2,c(param_mods_full[p], NA))
      
    }
    
  }
  
  #df_res$Est <- round(as.numeric(df_res$Est), digits = 2)
  return(list(df_res, df_res2))
  
}



FUN <- function(x,df, beta_1_ind, beta_2_ind, mod_name, monitor_vec, n.samp = 500,
                n.burnin = 100,n.thin = 1 ,adapt = 0, true_trt_diff_vec, pred_time = seq(0, 15, by = 0.1),
                path_export, true_cp, true_beta_cp, true_beta_cp_anc, true_lambda_wane, export_jags = FALSE){
  #Change-inits
  data_jags <- make_data_jags(df = df[[x]],
                              beta_1_ind,
                              beta_2_ind)
  
  if(mod_name == "jags.piecewise_wei_wane"){
    data_jags$trt_ind <- 2
  }
  
  jags_mod <- #try({
    runjags::run.jags(
      model = get(mod_name),
      data = data_jags,
      n.chains = 2,
      monitor = monitor_vec,
      sample=n.samp,
      thin = n.thin,
      burnin = n.burnin,
      adapt = 0,
      inits = inits(true_cp, true_beta_cp, true_beta_cp_anc, true_lambda_wane),
      method ='rjags',
      summarise = FALSE)
  #})

  if(class(jags_mod)=="try-error"){
    return(NA)
  }else{
    mat <- cbind(1,c(1,0))
     sim_res <-export_outcomes_sim_study(mod.cp =jags_mod, 
                                         data_jags =data_jags, 
                                         pred_mod_mat=mat,
                                         df.origin = df[[x]],
                                         pred_time = pred_time,
                                         true_trt_diff_vec = true_trt_diff_vec,
                                         path_export = gen_pathway(paste0(path_export,"Simulation ",x,"/")))
    if(export_jags){
      return(list(jags_mod, sim_res))
    }else{
      return(sim_res)
    }
  }
  # we keep the posterior survival and compare that to the true sample survival...
}
