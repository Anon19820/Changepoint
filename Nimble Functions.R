survSplit <- nimbleFunction(
  buildDerivs = list(run = list(ignore = c('i','j','p','k'))),
  #X matrix tstart; tstop,status, X var all
  run = function(X = double(2),
                 cut = double(1)){ # type declarations

    n = dim(X)[1];
    ncut = length(cut);
    ncol = dim(X)[2]
    tstart = X[,1]
    tstop = X[,2]
    status_orig = X[,3]
    extra =0;
    extra_cov = ncol- 3 #must have tstart; tstop,status,
    for (i in 1:n) {
      for (j in 1:ncut){
        if (cut[j] > tstart[i] & cut[j] < tstop[i]){
          extra = extra + 1;
        }
      }
    }
    n2 <- n + extra

    start <- numeric(n2);
    interval <- numeric(n2);
    end <- numeric(n2);
    status <- numeric(n2);
    id <- numeric(n2);
    X_res <- matrix(nrow = n2, ncol = ncol+2)

    k =1;

    for (i in 1:n) {

      j_eval <- 1

      for(j in 1:ncut){
        if(cut[j] <= tstart[i]){
          j_eval <- j_eval + 1
        }
      }

      start[k] = tstart[i];
      interval[k] = j_eval;
      id[k] = i;
      for(p in 4:ncol){
        X_res[k,p+2] <-X[i,p]
      }

      for (j in j_eval:ncut){

        if(cut[j] < tstop[i]){

          if (cut[j] > tstart[i]){
            end[k] = cut[j];
            status[k] = 0;
            k =k+1;
            start[k] = cut[j];
            interval[k] = j+1;
            id[k] = i;
            for(p in 4:ncol){
              X_res[k,p+2] <-X[i,p]
            }
          }
        }
      }
      end[k] = tstop[i];
      status[k] =status_orig[i];
      for(p in 4:ncol){
        X_res[k,p+2] <-X[i,p]
      }
      k <- k +1 ;
    }

    X_res[,1] <- start
    X_res[,2] <- end
    X_res[,3] <- status
    X_res[,4] <- interval
    X_res[,5] <- id
    #X_res[,6] <- 1 #Intercept

    return(X_res)
    returnType(double(2)) # return type declaration
  })

pweibullPH <- nimbleFunction(
  #X matrix tstart; tstop,status, X var all
  run = function(q = double(0),

                 scale = double(0),
                 shape = double(0),
                 lower.tail = double(0), #logical; if TRUE (default), probabilities are P(X \le x)P(X≤x), otherwise, P(X > x)P(X>x).
                 log = double(0)){
    #shape <- exp(ln_shape)
    #scale <- exp(ln_scale)

    if(lower.tail==1){
      res <- 1-exp(-scale*(q^shape))

    }else{
      res <-  exp(-scale*(q^shape)) #ST
    }

    if(log ==1){
      return(log(res))
      #returnType(double(0))
    }else{
      return(res)

    }
    returnType(double(0))
  }
)

HweibullPH <- nimbleFunction(
  #X matrix tstart; tstop,status, X var all
  run = function(q = double(0),

                 scale = double(0),
                 shape = double(0),
                 log = double(0)){



      res <-  scale*(q^shape)

    if(log ==1){
      return(log(res))
      #returnType(double(0))
    }else{
      return(res)

    }
    returnType(double(0))
  }
)


mat_adjust <- nimbleFunction(
  buildDerivs = list(run = list(ignore = c('i','j'))),
  #X matrix tstart; tstop,status, X var all
  run = function(mat = double(2),
                 mat2 = double(2)){
    n_col = dim(mat)[2];
    n_row = dim(mat)[1];
    mat_final = matrix(nrow = n_row, ncol = n_col)
   for(i in 1:n_row){

     for(j in 1:n_col){

       if(mat2[i,j] == 2){

         mat_final[i,1:n_col] <- mat[i,j]

       }else{
         mat_final[i,j] <- mat[i,j]*mat2[i,j]
       }
     }

   }


     return(mat_final)
     returnType(double(2))
  }
)



 #beta_2 <- matrix(0, byrow = FALSE, ncol = 2, nrow = 3)


 #beta_2[1,1] <- 1

 #beta_2[2,2] <- 2
 #beta_1 <- matrix(rnorm(6), byrow = FALSE, ncol = 2, nrow = 3)

#Cmat_adjust<-  compileNimble(mat_adjust)
 #undebug(mat_adjust)
 #Cmat_adjust(beta_1,beta_2)

# log(flexsurv::pweibullPH(0,exp(0),exp(0), lower.tail = T))
# pweibullPH <- function(q, ln_shape, ln_scale,
#                        lower.tail, log) {
#
#   shape <- exp(ln_shape)
#   scale <- exp(ln_scale)
#
#   if(lower.tail==1){
#     res <-  1-exp(-shape*(q^scale))
#   }else{
#     res <- exp(-shape*(q^scale))
#   }
#
#   if(log ==1){
#     log(res)
#   }else{
#     res
#   }
#
# }


# hweibullPH <- function(x, ln_shape, ln_scale){
#
#   shape <- exp(ln_shape)
#   scale <- exp(ln_scale)
#
#   scale*shape*(x^(shape-1))
# }

hweibullPH <- nimbleFunction(
  #X matrix tstart; tstop,status, X var all
  run = function(x = double(0),
                 scale = double(0),
                 shape = double(0),
                 log = double(0)){

    res <- scale*shape*(x^(shape-1))

    if(log ==1){
      return(log(res))
      # returnType(double(0))
    }else{
      return(res)
    }
    returnType(double(0))

  }
)

#Old Version
# cp_LL <- nimbleFunction(
#   run = function(X = double(2),
#                  cut = double(1),
#                  beta_1 = double(2),
#                  beta_2 = double(2),
#                  type = double(0)){
#
#     X_res <- survSplit(X, cut)
#
#     n <- dim(X)[1]
#     nrow_X = dim(X_res)[1];
#     ncol_X = dim(X_res)[2];
#     n_int = dim(beta_1)[2];
#
#     param_1 <- matrix(nrow = nrow_X, ncol = n_int)
#     param_2 <- matrix(nrow = nrow_X, ncol = n_int)
#
#     param_1_final <- numeric(nrow_X)
#     param_2_final <- numeric(nrow_X)
#     LL <- numeric(nrow_X)
#
#
#     for(i in 1:n_int){
#       param_1[,i] <-   X_res[,6:ncol_X] %*% beta_1[,i];
#       param_2[,i] <-   X_res[,6:ncol_X] %*% beta_2[,i];
#
#     }
#
#     for(i in 1:nrow_X){
#       param_1_final[i] <-  param_1[i,X_res[i,4]]
#       param_2_final[i] <-  param_2[i,X_res[i,4]]
#
#       if(type == 1){
#         LL[i] <- hweibullPH(X_res[i,2],param_1_final[i],param_2_final[i],1)*X_res[i,3]-
#           (HweibullPH(X_res[i,2],param_1_final[i],param_2_final[i],0)-HweibullPH(X_res[i,1],param_1_final[i],param_2_final[i],0))
#
#       }else{ #cumulative hazard
#         LL[i] <- (HweibullPH(X_res[i,2],param_1_final[i],param_2_final[i],0)-HweibullPH(X_res[i,1],param_1_final[i],param_2_final[i],0))
#
#       }
#
#     }
#
#     LL_final <- numeric(n);
#     for(i in 1:n){
#       LL_final[i] <- sum(LL[X_res[,5]==i])
#     }
#
#     if(type != 1){
#       LL_final <- exp(-LL_final)
#     }
#
#     return(LL_final)
#     returnType(double(1)) # return type declaration
#
#
#   }
# )


cp_LL <- nimbleFunction(
  run = function(X = double(2),
                 cut = double(1),
                 beta_1 = double(2),
                 beta_2 = double(2),
                 type = double(0),
                 trt_wane = double(0),
                 lambda_wane = double(0)){

    X_res <- survSplit(X, cut)

    n <- dim(X)[1]
    nrow_X = dim(X_res)[1];
    ncol_X = dim(X_res)[2];
    n_int = dim(beta_1)[2];

    param_1 <- matrix(nrow = nrow_X, ncol = n_int)
    param_2 <- matrix(nrow = nrow_X, ncol = n_int)

    param_1_final <- numeric(nrow_X)
    param_2_final <- numeric(nrow_X)
    res <- numeric(nrow_X)


    for(i in 1:n_int){
      param_1[,i] <-   X_res[,6:ncol_X] %*% beta_1[,i];
      param_2[,i] <-   X_res[,6:ncol_X] %*% beta_2[,i];

    }

    for(i in 1:nrow_X){
      param_1_final[i] <-  exp(param_1[i,X_res[i,4]])
      param_2_final[i] <-  exp(param_2[i,X_res[i,4]])

      if(type == 1){
        res[i] <- hweibullPH(X_res[i,2],param_1_final[i],param_2_final[i],1)*X_res[i,3]-
          (HweibullPH(X_res[i,2],param_1_final[i],param_2_final[i],0)-HweibullPH(X_res[i,1],param_1_final[i],param_2_final[i],0))

      }else if(type == 2){ #Survival
        if(trt_wane == 1 & X_res[i, 7] == 1 & X_res[i, 4] == n_int){ #if we have a waning of treatment, treatment arm, and final interval



          initial_HR <- exp(beta_1[2,n_int-1])
          t_end <- X_res[i,2]
          t_start <- X_res[i,1] #change-point
          #HR_wane <- 1-(1-initial_HR)*exp(-lambda_wane*(t_end-t_start))

          upper_inc_gamma1 <- exp(loggam(param_2_final[i]))*(1-pgamma(lambda_wane*t_end,param_2_final[i],1))
          upper_inc_gamma2 <- exp(loggam(param_2_final[i]))*(1-pgamma(lambda_wane*t_start,param_2_final[i],1))
          upper_int <- (t_end^param_2_final[i])/param_2_final[i]-(upper_inc_gamma1*(initial_HR-1)*exp(t_start*lambda_wane)*(t_end^param_2_final[i]))/((lambda_wane*t_end)^param_2_final[i])
          lower_int <- (t_start^param_2_final[i])/param_2_final[i]-(upper_inc_gamma2*(initial_HR-1)*exp(t_start*lambda_wane)*(t_start^param_2_final[i]))/((lambda_wane*t_start)^param_2_final[i])


          res[i] <- ((upper_int - lower_int)*param_2_final[i]*param_1_final[i])#*step(time[i]-cp[k]) #Has to be treatment and cp < time




        }else{
          res[i] <- (HweibullPH(X_res[i,2],param_1_final[i],param_2_final[i],0)-HweibullPH(X_res[i,1],param_1_final[i],param_2_final[i],0))

        }

      }else{ #hazard

        res[i] <- hweibullPH(X_res[i,2],param_1_final[i],param_2_final[i],0)*X_res[i,3] #must be final interval for hazard

        if(trt_wane == 1 & X_res[i, 7] == 1 & X_res[i, 4] == n_int){
          t_end <- X_res[i,2]
          t_start <- X_res[i,1] #change-point
          initial_HR <- exp(beta_1[2,n_int-1])
          HR_wane <- 1-(1-initial_HR)*exp(-lambda_wane*(t_end-t_start))
          res[i] <- res[i]*HR_wane
        }


      }

    }

    res_final <- numeric(n);
    for(i in 1:n){
      res_final[i] <- sum(res[X_res[,5]==i])
    }

    if(type == 2){#survival
      res_final <- exp(-res_final)
    }

    return(res_final)
    returnType(double(1)) # return type declaration


  }
)

Ccp_LL<- compileNimble(cp_LL)


calcLoglike <- nimbleFunction(
  setup = function(model, wrt, nodes){},
  run = function(parVals = double(1)){
    values(model, wrt) <<- parVals
    ll <- model$calculate(nodes)
    return(ll)
    returnType(double())
  }
)


gen_pathway<- function(string){
  folder_bits <- strsplit(string, "\\/")[[1]]

  file_extensions <- c("R", "png", "xlsx", "csv", "jpg", "svg", "txt", "pptx","pdf", "Rmd", "bug")
  file_extensions_grep <- paste0(paste0("\\.",file_extensions), collapse ="|")


  exclude.last <- grepl(file_extensions_grep,folder_bits[length(folder_bits)], ignore.case = T)


  if(exclude.last){
    folder_num <- length(folder_bits) -1

  }else{
    folder_num <- length(folder_bits)
  }

  for(i in 1:folder_num){
    folder_path_temp <- paste0(folder_bits[1:i], collapse = "/")

    if(!file.exists(folder_path_temp)){
      dir.create(folder_path_temp)
    }
  }

  return(string) #Need to return string for actual pasting

}

sort_vec <- nimbleFunction(
  buildDerivs = list(run = list(ignore = c('ii','ii1'))),
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
