library("nimble")
library("dplyr")
library("survival")
library("parallel")
pathway <- "~/Change-Point Simulation Studies/Simulation Study 2023/"
source(paste0(pathway, "Jags_codes2.R"))
source(paste0(pathway, "Nimble Functions.R"))
source(paste0(pathway, "Simulation-Study-Functions.R"))


param_mods_full <- c("Change-point","Exponential", "Weibull", "Log-Normal", "Gamma", "Gompertz", "Log-Logistic", "Generalized-Gamma", "Royston Parmar (PH)",
                     "Royston Parmar (non-PH)")

#scenario_names <- c("Treatment Delay", "Loss of Treatment Effect", "Converging Treatment Effect")
#model_names <- c("jags.piecewise_wei_chng","jags.piecewise_wei_chng","jags.piecewise_wei_wane")

#scenario_names <- c("Converging Treatment Effect Redo")
#model_names <- c("jags.piecewise_wei_wane")


monitor_vec <- c("cp", "beta_cp","beta_cp_anc", "total_loglik","loglik")

monitor_list <- list()
monitor_list[[1]] <- monitor_vec
monitor_list[[2]] <- monitor_vec
monitor_list[[3]] <- c(monitor_vec, "lambda_wane")

list_beta1 <- list()
list_beta2 <- list()

pred_time_grid <- seq(0, 30, by = 0.01)
pred_time_eval <- seq(0, 15, by = 0.1 )
n_sims <- 50

run_flag <- TRUE # Run the analysis or just plot survival and hazard plots
#pred_mod_mat <- cbind(1,c(1,0)) # Intercept and treatment
pred_mat_final_trt <- cbind(tstart = 0, pred_time_grid,status= 1,Intercept = 1,trt = 1)
pred_mat_final_comp <- cbind(tstart = 0, pred_time_grid,status= 1,Intercept =1,trt = 0)
pred_mat_final_trt2 <- cbind(tstart = 0, pred_time_eval,status= 1,Intercept = 1,trt = 1)
pred_mat_final_comp2 <- cbind(tstart = 0, pred_time_eval,status= 1,Intercept =1,trt = 0)

#Set up Parallel Code

ncores <- min(12, n_sims) # 6 seems to correct;; 48 is much too large.

export.cluster = c('make_data_jags','inits','jags.piecewise_wei_chng','jags.piecewise_wei_wane','make_data_jags','export_outcomes_sim_study','get_rmst_all_sim','km_AUC', "cp_LL",
                   "survSplit","pweibullPH","HweibullPH","mat_adjust", "hweibullPH","sort_vec","cp_LL_nim", "gen_pathway")

parallel.options <- list()
parallel.options$X <- 1:n_sims
parallel.options$fun <- FUN


parallel.options$pred_time = seq(0,15, 0.1)
#parallel.options$mod_name <- "jags.piecewise_wei_chng"
# parallel.options$beta_1_ind <-beta_1_ind 
# parallel.options$beta_2_ind <-beta_2_ind 
#parallel.options$monitor_vec <- c("cp", "beta_cp","beta_cp_anc", "total_loglik","loglik")
parallel.method <- parLapply
parallel.options$cl <- parallel::makeCluster(ncores)

clusterExport(parallel.options$cl, varlist = export.cluster, 
              envir = .GlobalEnv)

#helper function
`%!in%` <- Negate(`%in%`)
abs_mean <- function(x){mean(abs(x), na.rm= TRUE)}

#generate the datasets

shape_1 <- c(1.3,0.75)
scale_1 <- 0.3
HR <- c(0.25,0.5, 0.75)
cp <- c(1)
t_cens = 5
n.samp_size <- c(100, 300, 500)
lambda_wane <- c(1,0.75,0.5, .25)

max_eval_time <- 15

param_vals_overall <- expand.grid(n.samp = n.samp_size, cp = cp,t_cens = t_cens,
                                  scale_1=scale_1, HR =HR,shape_1=shape_1, lambda_wane =lambda_wane)


param_vals_overall<- param_vals_overall %>% mutate(t_cens = ifelse(shape_1 == 1.3, 3, 5)) 

#TD Scenario
beta_1_ind <- beta_2_ind <- matrix(1, byrow = FALSE, ncol = 2, nrow = 2)
beta_1_ind[1,1:2] <- 2
beta_1_ind[2,1] <- 0
beta_2_ind[1,1:2] <- 2
beta_2_ind[2,1:2] <-0

list_beta1[[1]] <- beta_1_ind
list_beta2[[1]] <- beta_2_ind

#LTE
beta_1_ind <- matrix(1, byrow = FALSE, ncol = 2, nrow = 2)
beta_1_ind[1,1:2] <- 2
beta_1_ind[2,2] <- 0
list_beta1[[2]] <- beta_1_ind
list_beta2[[2]] <- beta_2_ind

#Converge Treatment effect - Same as the LTE
list_beta1[[3]] <- beta_1_ind
list_beta2[[3]] <- beta_2_ind


for(scen in 1:length(scenario_names)){
  
  parallel.options$n.samp = 1000
  parallel.options$n.burnin = 500
  parallel.options$n.thin = 1
  parallel.options$adapt = 100
  
  if(scenario_names[scen] == "Converging Treatment Effect"){
    parallel.options$n.samp <- parallel.options$n.samp*2
    parallel.options$n.burnin <- parallel.options$n.burnin*2
    parallel.options$n.thin <- parallel.options$n.thin*2
    parallel.options$adapt <- parallel.options$adapt*2
  }

  parallel.options$mod_name <- model_names[scen]
  parallel.options$monitor_vec <-  monitor_list[[scen]]
  parallel.options$beta_1_ind <-list_beta1[[scen]] 
  parallel.options$beta_2_ind <-list_beta2[[scen]] 
  
  
  if(scenario_names[scen] != "Converging Treatment Effect"){
    param_vals <- distinct(param_vals_overall[,-ncol(param_vals_overall)])
  }else{
    param_vals <- param_vals_overall
  }
  
  write.csv(param_vals, gen_pathway(paste0(pathway,scenario_names[scen],"/param_vals.csv")))
  param_vals_res2 <-param_vals_res <- matrix(nrow = nrow(param_vals), ncol = ncol(param_vals)) #Evaluate the Mean Param
  param_vals_res2[,1:2] <-param_vals_res[,1:2] <- as.matrix(param_vals[,1:2]) #Evaluate the SE
  
  param_surv_mods <- matrix(nrow = nrow(param_vals), ncol = ncol(param_vals)+10 )
  colnames(param_surv_mods) <- c(names(param_vals),param_mods_full )
  param_surv_mods[1:nrow(param_vals), 1:ncol(param_vals)] <- as.matrix(param_vals)
  
  param_surv_mods2 <-param_surv_mods
  
  
 for(i in 1:nrow(param_vals_overall)){
 # i = 1
   
 scenario_pathway <- paste0(pathway,scenario_names[scen],"/Scenario ", i,"/" )
  
  beta_1_sim <- beta_2_sim <- matrix(0, byrow = FALSE, ncol = 2, nrow = 2)
  beta_1_sim[1,] <- log(param_vals$scale_1[i])
  if(scenario_names[scen] == "Treatment Delay"){
    beta_1_sim[2,] <- c(0,log(param_vals$HR[i]))# Treatment Delay
  }else{
    beta_1_sim[2,] <- c(log(param_vals$HR[i]),0) #Other scenarios
  }
  beta_2_sim[1,] <- log(param_vals$shape_1[i])  
  
  if(is.null(param_vals[i,"lambda_wane"])){ #IF lambda_wane doesn't exist
    lambda_wane_sim = 0
    trt_wane_sim = 0
  }else{
    lambda_wane_sim = param_vals$lambda_wane[i]
    trt_wane_sim = 1
  }
  
  

  parallel.options$true_cp <- param_vals$cp[i]
  parallel.options$true_beta_cp <- beta_1_sim
  parallel.options$true_beta_cp_anc <- beta_2_sim
  parallel.options$true_lambda_wane <- lambda_wane_sim
  

  St_trt <- Ccp_LL(pred_mat_final_trt, c(param_vals$cp[i],Inf),beta_1 = beta_1_sim,beta_2 = beta_2_sim,type = 2, lambda_wane = lambda_wane_sim, trt_wane = trt_wane_sim)
  St_comp <- Ccp_LL(pred_mat_final_comp, c(param_vals$cp[i],Inf),beta_1 = beta_1_sim,beta_2 = beta_2_sim,type = 2, lambda_wane = lambda_wane_sim, trt_wane = trt_wane_sim)
  
  St_trt2 <- Ccp_LL(pred_mat_final_trt2, c(param_vals$cp[i],Inf),beta_1 = beta_1_sim,beta_2 = beta_2_sim,type = 2, lambda_wane = lambda_wane_sim, trt_wane = trt_wane_sim)
  St_comp2 <- Ccp_LL(pred_mat_final_comp2, c(param_vals$cp[i],Inf),beta_1 = beta_1_sim,beta_2 = beta_2_sim,type = 2, lambda_wane = lambda_wane_sim, trt_wane = trt_wane_sim)
  
  
  parallel.options$true_trt_diff_vec <- St_trt2-St_comp2
  parallel.options$path_export <- scenario_pathway
  
  png(gen_pathway(paste0(scenario_pathway,"St_true.png")))
  plot(y = St_trt, x = pred_time_grid, ylim  = c(0,1), type = "l")
  lines(y = St_comp, x = pred_time_grid, col = "red")
  abline(v =param_vals$t_cens[i], lty  = 2 )
  dev.off()
  
  haz_trt <- Ccp_LL(pred_mat_final_trt, c(param_vals$cp[i],Inf),beta_1 = beta_1_sim,beta_2 = beta_2_sim,type = 3, lambda_wane = lambda_wane_sim, trt_wane = trt_wane_sim)
  haz_comp <- Ccp_LL(pred_mat_final_comp, c(param_vals$cp[i],Inf),beta_1 = beta_1_sim,beta_2 = beta_2_sim,type = 3, lambda_wane = lambda_wane_sim, trt_wane = trt_wane_sim)
  
  png(gen_pathway(paste0(scenario_pathway,"haz_true.png")))
  plot(y = haz_trt, x = pred_time_grid, xlim  = c(0,10),ylim = c(0, 2), type = "l")
  lines(y = haz_comp, x = pred_time_grid, col = "red")
  abline(v =param_vals$t_cens[i], lty  = 2 )
  dev.off()
  
  if(run_flag){
  
  #True Values -- Important
  true_AUC_trt <- sfsmisc::integrate.xy(x = pred_time_grid[pred_time_grid<max_eval_time], fx = St_trt[pred_time_grid<max_eval_time])
  true_AUC_comp <- sfsmisc::integrate.xy(x = pred_time_grid[pred_time_grid<max_eval_time], fx = St_comp[pred_time_grid<max_eval_time])
  
  
  
  data_list <- list()
  for(j in 1:n_sims){
      data_list[[j]] <- gen_cp_df(param_vals$n.samp[i],t_cens =param_vals$t_cens[i], plot = T, St_trt =St_trt, St_comp = St_comp )
  }
  
  parallel.options$df <- data_list
  

  success <- try({
      output <- capture.output({
        results <- do.call("parallel.method", args = parallel.options)
      })
  })
  
  # Check to see if there are any errors....
  index_keep <- sapply(results, function(x){!(is.nan(x[[2]][1,2]) |is.na(x[[2]][1,2])) })
  if(any(index_keep==FALSE)){
    warning(paste0(n_sims - sum(index_keep), " Simulations produced errors"))
  }
   
  results <- results[index_keep]

  diff_mat <-  do.call(rbind,args = lapply(results, function(x){x[[2]]})) %>% group_by(Model) %>%
    summarise_if(is.numeric, abs_mean)
  
  int_mat <- do.call(rbind,args = lapply(results, function(x){x[[4]]})) %>% group_by(Model) %>%
    summarise_if(is.numeric, mean)
  
  diff_mat <-  diff_mat %>%  arrange(match(Model, c("Changepoint", "Exponential", "Weibull", "Log-Normal", "Gamma",
                                        "Gompertz", "Log-Logistic", "Generalized-Gamma", "Royston Parmar (PH)",
                                        "Royston Parmar (non-PH)")))
  int_mat <-  int_mat %>%  arrange(match(Model, c("Changepoint", "Exponential", "Weibull", "Log-Normal", "Gamma",
                                                    "Gompertz", "Log-Logistic", "Generalized-Gamma", "Royston Parmar (PH)",
                                                    "Royston Parmar (non-PH)")))
  
  
  write.csv(diff_mat, paste0(scenario_pathway, "RMST_abs.csv")) # This is the RMST_diff of the model
  
  diff_mat_prop_error <- diff_mat
  diff_mat_prop_error[,2:4] <-  sweep(diff_mat_prop_error[,2:4],MARGIN = 2, STATS = c(true_AUC_trt-true_AUC_comp,true_AUC_trt,true_AUC_comp), FUN = "/")
  
  
  write.csv(diff_mat_prop_error, paste0(scenario_pathway, "RMST_prop.csv")) # This is the RMST_diff/true estimate of the model
  
  beta_mat <- do.call(rbind,args = lapply(results, function(x){x[[3]]}))
  if(scenario_names[scen] != "Converging Treatment Effect"){
    beta_mat[,2:ncol(beta_mat)] <- exp(beta_mat[,2:ncol(beta_mat)])
  }else{
    beta_mat[,2:(ncol(beta_mat)-1)] <- exp(beta_mat[,2:(ncol(beta_mat)-1)])
  }

  beta_mat_final <- beta_mat[,!duplicated(beta_mat[1,])]
  
  avg.post_median <- colMeans(beta_mat_final)
  se.post_median <- apply(beta_mat_final,2,sd)
  
 true_vals <- t(param_vals[i,colnames(param_vals) %!in% c("n.samp","t_cens")])
 df_sims <-  data.frame(Mean = avg.post_median, SE =se.post_median, True_val =as.numeric(true_vals))
 
 write.csv(df_sims,
           paste0(scenario_pathway, "Param_error.csv")) 
 
  
 if(i == 1){
   colnames(param_vals_res2)  <- colnames(param_vals_res) <- c("n.samp","t_cens", names(avg.post_median))
 }
 param_vals_res[i,3:ncol(param_vals_res)] <- avg.post_median
 param_vals_res2[i,3:ncol(param_vals_res)] <- se.post_median

 param_surv_mods[i,(ncol(param_vals)+1):ncol(param_surv_mods)] <-  int_mat$Int_Diff
 param_surv_mods2[i,(ncol(param_vals)+1):ncol(param_surv_mods2)] <-  diff_mat_prop_error$Diff_Est
 
 write.csv(int_mat,
           paste0(scenario_pathway, "int_mat.csv")) 
 
 
  check_km <- do.call(rbind,args = lapply(results, function(x){x[[1]][2,1:4]})) %>% mutate_if(is.numeric, mean) %>% distinct()
  check_km <- rbind(check_km, data.frame(Model = "True", 
                                         Diff_Est = true_AUC_trt-true_AUC_comp,
                                         AUC_trt = true_AUC_trt,
                                         AUC_comp = true_AUC_comp)) %>%
              mutate_if(is.numeric, round, digits = 3)
  
  
  write.csv(check_km,
            paste0(scenario_pathway, "qc_compare_km.csv")) 
  
  
  # crashed <- sapply(results, class) == "try-error"
  # error <- sapply(results[crashed], as.character)
  # if (length(error) != 0) {
  #   names(error) <- paste("simulation.", which(crashed),
  #                         sep = "")
  # 
  # }
  
   }

 }
  
  if(run_flag){
    write.csv(param_vals_res ,paste0(pathway,scenario_names[scen], "/Mean_param.csv"))
    write.csv(param_vals_res2 ,paste0(pathway,scenario_names[scen], "/SE_param.csv"))
    write.csv(param_surv_mods ,paste0(pathway,scenario_names[scen], "/Int_error_param.csv"))
    write.csv(param_surv_mods2 ,paste0(pathway,scenario_names[scen], "/Diff_RMST.csv"))
  }
  
}  
#model




#test code for parallel

md <- FUN(1, df =parallel.options$df,
          beta_1_ind = parallel.options$beta_1_ind,
          beta_2_ind = parallel.options$beta_2_ind,
          mod_name = parallel.options$mod_name,
          monitor_vec = parallel.options$monitor_vec,
          true_trt_diff_vec = parallel.options$true_trt_diff_vec,
          path_export = parallel.options$path_export,
          true_cp = parallel.options$true_cp,
          true_beta_cp = parallel.options$true_beta_cp,
          true_beta_cp_anc = parallel.options$true_beta_cp_anc,
          true_lambda_wane = parallel.options$true_lambda_wane,
          export_jags = TRUE)

str(inits(true_cp, true_beta_cp, true_beta_cp_anc, true_lambda_wane))

inits(true_cp, true_beta_cp, true_beta_cp_anc, true_lambda_wane)
apply(md[[1]][,2:5], 2, as.numeric)

fit.km <- survfit(Surv(time,status)~arm, data = parallel.options$df[[1]])
fit.km_final <- survfit(Surv(time_event ,status_true )~arm, data = parallel.options$df[[1]])
scenario1_rmst <- export_outcomes(fit.km = fit.km_final,md[[1]], 
                                  data_jags = make_data_jags(parallel.options$df[[1]],
                                                             beta_1_ind = parallel.options$beta_1_ind, 
                                                             beta_2_ind = parallel.options$beta_2_ind),
                                  pred_mod_mat = cbind(1,c(1,0)),
                                  pred_time= seq(0,15, by = 0.1),
                                  plt_title = "Test Sim",
                                  xlab_title = "Time (Years)",
                                  pathway_export= gen_pathway(paste0(pathway,"test/")),
                                  df.origin =  parallel.options$df[[1]])



