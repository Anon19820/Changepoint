pathway_R <- "~/Change-Point Simulation Studies/Simulation Study 2023/Treatment Delay/"
#<- read.csv(paste0(pathway_R, "param_vals.csv"))

Mean_param <- read.csv(paste0(pathway_R, "Mean_param.csv"))[, - 1]
SE_param <- read.csv(paste0(pathway_R, "SE_param.csv"))[,-1]

Int_error_param <- read.csv(paste0(pathway_R, "/Int_error_param.csv"))[,-1]
Diff_RMST <- read.csv(paste0(pathway_R, "/Diff_RMST.csv"))[,-1]
nsims <- nrow(Mean_param)

for(i in 1:nsims){
 df_temp <- read.csv(paste0(pathway_R,"Scenario ", i,"/Param_error.csv"))
 Mean_param[i,3:ncol(Mean_param)] <- df_temp$Mean
 SE_param[i,3:ncol(SE_param)] <- df_temp$SE
   
 df_temp <- read.csv(paste0(pathway_R,"Scenario ", i,"/int_mat.csv")) %>% 
   arrange(match(Model, c("Changepoint", "Exponential", "Weibull", "Log-Normal", "Gamma",
                          "Gompertz", "Log-Logistic", "Generalized-Gamma", "Royston Parmar (PH)",
                          "Royston Parmar (non-PH)")))
 index_start <- which(colnames(Int_error_param) == "Change.point")
 
 Int_error_param[i,index_start:ncol(Int_error_param)] <- df_temp$Int_Diff
 
 df_temp <- read.csv(paste0(pathway_R,"Scenario ", i,"/RMST_prop.csv")) %>% 
   arrange(match(Model, c("Changepoint", "Exponential", "Weibull", "Log-Normal", "Gamma",
                          "Gompertz", "Log-Logistic", "Generalized-Gamma", "Royston Parmar (PH)",
                          "Royston Parmar (non-PH)")))
 Diff_RMST[i,index_start:ncol(Diff_RMST)] <- df_temp$Diff_Est
 
}

write.csv(Mean_param,paste0(pathway_R, "/Mean_param.csv"))
write.csv(SE_param,paste0(pathway_R, "/SE_param.csv"))

write.csv(Int_error_param,paste0(pathway_R, "/Int_error_param_correct.csv"))
write.csv(Diff_RMST,paste0(pathway_R, "/Diff_RMST_correct.csv"))

# 
# Int_error_param_compare <- read.csv(paste0(pathway_R, "/Int_error_param.csv"))[,-1]
# Diff_RMST_compare <- read.csv(paste0(pathway_R, "/Diff_RMST.csv"))[,-1]

