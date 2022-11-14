#####This is the complementary functions file to our simulation file#####


################################
######Estimate Beta#############

est.beta <- function(X,Y){ ###standard OLS estimator
  beta_est <- c()
  beta_est = solve(t(X)%*%X)%*%t(X)%*%Y
  return(beta_est)
}

###########################################################
######Create dummy to compare Beta_hat to true values######

comparison_dummy <- function(truetreatmenteffect, trueTimeeffects, trueGroupeffects, G, T){
  comparison <- c()
  for(t in 1:(T-1)){
    comparison[t] = trueTimeeffects[t+1] - trueTimeeffects[1] 
  }
  for(g in 1:(G-1)){
    comparison[(T-1)+g] = trueGroupeffects[g+1] - trueGroupeffects[1]
  }
  comparison[(T+G-1)] = truetreatmenteffect
  return(comparison)
}

################################
######Estimate r2###############
est.r2 <- function(X,Y,beta_hat){
  Y_hat <- c()
  Y_hat = X%*%beta_hat
  
  s_pred = 0
  s_tot = 0
  for(i in 1:length(Y)){
    s_pred = s_pred + (Y[i]-Y_hat[i])^2
    s_tot = s_tot + (Y[i]-mean(Y))^2
  }
  
  r2 = 1 - s_pred/s_tot
  return(r2)
}

###############################################
###### Standard deviation and T_statistic######

standard_deviation <- function(X,Y,beta_hat){ ###estimates the variance matrix for beta
  n = dim(X)[1]
  k = dim(X)[2]
  
  Y_hat = X%*%beta_hat
  sigma = 0
  for(i in 1:length(Y)){
    sigma = sigma + (Y[i] - Y_hat[i])**2
  }
  sigma = sqrt(sigma/(n-k-1))
  return((solve(t(X)%*%X)) * sigma)
}

test_t_stat <- function(X,Y,beta_hat,hypothesis,to_test){  #####Assumes N > 1000
  sd_beta = standard_deviation(X,Y,beta_hat)
  test_stat = (beta_hat[to_test] - hypothesis)/sd_beta[to_test,to_test]
  significance = 0
  t_stat_0.05 = 1.96
  t_stat_0.01 = 2.326
  t_stat_0.001 = 3.291
  
  if(abs(test_stat) >= t_stat_0.05){ 
    significance = 1
    if(abs(test_stat >= t_stat_0.01)){
      significance = 2
      if(abs(test_stat >= t_stat_0.001)){
        significance = 3
      }
    }
  }
  results <- list(significance,sd_beta[to_test,to_test])
  return(results)
}


##############################################################
######Estimate N1 = total number of treated observations######

est.N1 <- function(X,d){ ### d refers to the collum containing the Treatment information
  N1 = 0
  for(i in 1:dim(X)[1]){
    if(X[i,d] > 0){
      N1 = N1 + 1
    }
  }
  return(N1)
}


#########################
######Estimate N_gt######

est.N_gt <- function(X,N,G,T){
N_gt <- matrix(0, nrow = G, ncol = T)

for(i in 1:N){
  for(t in 1:T){
    for(g in 2:G){
      
      if(X[((i-1)*T + t),(T-1+(g-1))] == 1){  
        N_gt[g,2:T] = N_gt[g,2:T] + X[((i-1)*T +t),1:(T-1)]
      }
      
      
      if(X[((i-1)*T +t),(T-1+(g-1))] == 1 && sum(X[((i-1)*T +t),1:(T-1)]) == 0){ ###border case t = 1
        N_gt[g,1] = N_gt[g,1] + 1 
      }
    }
    if(sum(X[((i-1)*T +t),1:(T+G-2)]) == 0){ ###border case t = 1, g = 1
      N_gt[1,1] = N_gt[1,1] + 1
    }
    for(t2 in 1:(T-1)){
      if(X[((i-1)*T + t),t2] == 1 && sum(X[((i-1)*T + t),T:(T+G-2)])==0){ ###border case g = 1
        N_gt[1,(t2+1)] = N_gt[1,(t2+1)] + 1
      }
    }
  }
}
return(N_gt)
}

#####################################
#####Estimate LDelta_TR (ATT)########

est.LDelta_TR <- function(X,Y,d,beta_hat){
Y_1 <- c()
Y_0 <- c()
j=1
for (i in 1:dim(X)[1]){
  if (X[i,d] > 0){
    Y_1[j]=Y[i]          
    dummyvector <- c()
    for(h in 1:dim(X)[2]){
      dummyvector[h]=X[i,h]
    }
    dummyvector[d]=0
    Y_0[j]=dummyvector%*%beta_hat
    j = j +1
  }
}

LDelta_TR = 1/est.N1(X,d) * sum(Y_1 - Y_0)

return(LDelta_TR)
}

###################################
#####Estimate LDelta_gt (ATE)######

est.LDelta_gt <- function (X,Y,N,G,T,d){
  
LDelta_gt <- matrix(0, nrow = G, ncol = T)

N_gt <- est.N_gt(X,N,G,T)
beta_mod <- est.beta(X,Y)
beta_mod[d] = 0 ###used to estimate Y_igt(0)

for(i in 1:N){
  for(t in 1:T){
    for(g in 2:G){
      if(X[(i-1)*T+t,d] > 0 && X[(i-1)*T+t,(T-1+(g-1))] == 1){
        LDelta_gt[g,t] = LDelta_gt[g,t] + (1/N_gt[g,t]) * (Y[((i-1)*T+t)] - (X[((i-1)*T+t),]%*%beta_mod))
      }
    }
    if(sum(X[(i-1)*T+t,T:(T+G-2)]) == 0 && X[(i-1)*T+t,d] > 0){
      LDelta_gt[1,t] = LDelta_gt[1,t] + (1/N_gt[1,t]) * (Y[((i-1)*T+t)] - (X[((i-1)*T+t),]%*%beta_mod))
    }
  }
}
return(LDelta_gt)
}


###############################
######Estimate delta_TR########

est.delta_TR <- function(X,Y,N,G,T,d){ 
  delta_TR = 0
  N_gt <- est.N_gt(X,N,G,T)
  N1 <- est.N1(X,d)
  LDelta_gt <- est.LDelta_gt(X,Y,N,G,T,d)
  
  for (t in 1:T){
    for (g in 1:G){
      delta_TR = delta_TR + N_gt[g,t]/N1 * LDelta_gt[g,t] 
    }
  }
  return(delta_TR)
}

################################################################################################
#####create D_gt matrix containing information on average treatment in group g at period t######

est.D_gt <- function(X,Y,N,G,T,d){
D_gt <- matrix(0, nrow = G, ncol = T)
N_gt = est.N_gt(X,N,G,T)
for (i in 1:N){
  for(t in 1:T){
    for(g in 2:G){
      if(X[(i-1)*T + t,d] > 0 && X[(i-1)*T + t,(T-1+(g-1))] == 1){
        D_gt[g,t] = D_gt[g,t] + 1/N_gt[g,t] * X[(i-1)*T + t,d]
      }
    }
    if(X[(i-1)*T + t,d] > 0 && sum(X[(i-1)*T + t,T:(T+G-2)]) == 0){ #####border case group 1
      D_gt[1,t] = D_gt[1,t] + 1/N_gt[1,t] * X[(i-1)*T + t,d]
    }
  }
}
return(D_gt)
}

##############################################
#####estimate Two-Way weights version 1#######

est.weights_v1 <- function(X,Y,N,G,T,d){
  
  beta_hat = est.beta(X,Y)
  D_gt = est.D_gt(X,Y,N,G,T,d)
  N_gt = est.N_gt(X,N,G,T)
  N1 = est.N1(X,d)
  lambda <- c(0, beta_hat[1:(T-1)])          #####vector containing the time fixed effects estimates
  gamma <- c(0, beta_hat[(T-1):(T+G-2)])     #####vector containing the group fixed effects estimates

  alpha = mean(D_gt) - mean(lambda) - mean(gamma) 

  ## compute epsilon_v1[g,t] = D[g,t] - gamma[g] - lambda[t]
  epsilon <- matrix(NA, ncol = T, nrow = G)
  
  for(t in 1:T){
    for(g in 1:G){
      epsilon[g,t] = D_gt[g,t] - lambda [t] - gamma[g] - alpha
    }
  }

  w <- matrix(0, ncol = T, nrow = G)
  sum_eps = 0
  
  for (t in 1:T){
    for (g in 1:G){
      sum_eps = sum_eps + D_gt[g,t] *  N_gt[g,t]/N1 * epsilon[g,t]
    }
  }
  
  for(t in 1:T){
    for(g in 1:G){
      w[g,t] = epsilon[g,t]/sum_eps
    }
  }
  return(w)
}


###############################################
#####estimate Two-Way weights version 2########

est.weights_v2 <- function(X,Y,N,G,T,d){  ####### Note that for this method to be accurate N_gt[g,t]/N_gt[g,t-1] may not vary across g over time
  epsilon <- matrix(0, nrow = G, ncol = T)
  D_gt <- est.D_gt(X,Y,N,G,T,d)
  N_gt <- est.N_gt(X,N,G,T)
  N1 = est.N1(X,d)

 
  D_g_sum <- c()
  for(g in 1:G){
    D_g_sum[g] = 0
  }
  D_t_sum <- c()
  for(t in 1:T){
    D_t_sum[t] = 0
  }
  D_tot_sum = 0
  
  for(t in 1:T){
    for(g in 1:G){
      D_g_sum[g] = D_g_sum[g] + (N_gt[g,t]/sum(N_gt[g,]))*D_gt[g,t]
      D_t_sum[t] = D_t_sum[t] + (N_gt[g,t]/sum(N_gt[,t]))*D_gt[g,t]
      D_tot_sum = D_tot_sum + (N_gt[g,t]/sum(N_gt[,]))*D_gt[g,t]
    }
  }


  
  for(t in 1:T){
    for(g in 1:G){
      epsilon[g,t] = D_gt[g,t] - D_g_sum[g] - D_t_sum[t] + D_tot_sum  
    }
  }


  sum_eps = 0

  for(t in 1:T){
    for(g in 1:G){
      if(D_gt[g,t] > 0){
      sum_eps = sum_eps + N_gt[g,t]/N1 * epsilon[g,t]
      }
    }
  }
  
  
  w <- matrix(NA, ncol = T, nrow = G)
 
  for(t in 1:T){
    for(g in 1:G){
      w[g,t] = epsilon[g,t] / sum_eps 
    }
  }
  return(w)
}


###################################################
####### Control for sum(w_gt | D_gt = 1) = 1 ######

check.sum_weights <- function (X,Y,N,G,T,d,weights) {  ###### The sole purpose of this function is to be able to check, whether the computed w_gt[g,t] sum up to 1 once scaled properly.
  D_gt <- est.D_gt(X,Y,N,G,T,d)
  N1 <- est.N1(X,d)
  N_gt <- est.N_gt(X,N,G,T)
  sum_w = 0
  for(g in 1:G){
    for(t in 1:T){
      if(D_gt[g,t] > 0){
        sum_w = sum_w + ((N_gt[g,t]/N1) * w_gt[g,t])
      }
    }
  }
}


##########################################################
###### Robustness to heterogenous treatment effects ######

est.sigma_delta <- function(X,Y,N,G,T,d){ #####Estimates the standard deviation of the ATEs from the ATT
  sigma_delta = 0
  D_gt <- est.D_gt(X,Y,N,G,T,d)
  N_gt <- est.N_gt(X,N,G,T)
  N1 <- est.N1(X,d)
  delta_gt <- est.LDelta_gt(X,Y,N,G,T,d)
  LDelta_TR <- est.LDelta_TR(X,Y,d,est.beta(X,Y))
  
  for(g in 1:G){
    for(t in 1:T){
      if(D_gt[g,t] > 0){
        sigma_delta = sigma_delta + ((N_gt[g,t]/N1) * (delta_gt[g,t] - LDelta_TR)**2)
      }
    }
  }
  sigma_delta = sqrt(sigma_delta)
  return(sigma_delta)
}


est.sigma_weights <- function(X,Y,N,G,T,d,v_weights = 1){ ######Estimates the standard deviation of the weights
  sigma_weights = 0
  D_gt <- est.D_gt(X,Y,N,G,T,d)
  N_gt <- est.N_gt(X,N,G,T)
  N1 <- est.N1(X,d)
  if(v_weights == 1){
    w <- est.weights_v1(X,Y,N,G,T,d)
  }
  if(v_weights == 2){
    w <- est.weights_v2(X,Y,N,G,T,d)
  }
  
  
  for(g in 1:G){
    for(t in 1:T){
      if(D_gt[g,t] > 0){
        sigma_weights = sigma_weights + ((N_gt[g,t]/N1) * ((w[g,t] - 1)**2))
      }
    }
  }
  sigma_weights = sqrt(sigma_weights)
  return(sigma_weights)
}

inference <- function(X,Y,N,G,T,d,v_weights = 1){
  D_gt <- est.D_gt(X,Y,N,G,T,d)
  N_gt <- est.N_gt(X,N,G,T)
  N1 <- est.N1(X,d)
  LDelta_gt <- est.LDelta_gt(X,Y,N,G,T,d)
  beta_hat <- est.beta(X,Y)
  LDelta_TR <- est.LDelta_TR(X,Y,d,beta_hat)
  sigma_delta <- est.sigma_delta(X,Y,N,G,T,d)
  sigma_weights <- est.sigma_weights(X,Y,N,G,T,d,v_weights)
  beta_fe <- est.beta_fe(X,Y,N,G,T,d,v_weights)
  
  if(v_weights == 1){
    w <- est.weights_v1(X,Y,N,G,T,d)
  }
  if(v_weights == 2){
    w <- est.weights_v2(X,Y,N,G,T,d)
  }
  
  
  count_treated_gt_cells = 0
  for(g in 1:G){
    for(t in 1:T){
      if(D_gt[g,t] > 0){
        count_treated_gt_cells = count_treated_gt_cells + 1
      }
    }
  }
  
  sorted_by_weights <- matrix(NA, nrow = count_treated_gt_cells, ncol = 3)
  j = 1
  n = 0
  
  for(g in 1:G){
    for(t in 1:T){
      if(D_gt[g,t]  > 0){
         sorted_by_weights[j,1] = w[g,t]
         sorted_by_weights[j,2] = N_gt[g,t]
         sorted_by_weights[j,3] = LDelta_gt[g,t]
         n = n + 1
         j = j + 1
      }
    }
  }
  sorted_by_weights[order(sorted_by_weights[,1], decreasing = TRUE),]
  
  P_k <- c()
  S_k <- c()
  T_k <- c()
  for(k in 1:count_treated_gt_cells){
    P_k[k] = 0
    S_k[k] = 0
    T_k[k] = 0
  }
  for(k in 1:count_treated_gt_cells){
    i = count_treated_gt_cells
    while(i >= k){
      P_k[k] = P_k[k] + (sorted_by_weights[k,2]/N1)
      S_k[k] = S_k[k] + ((sorted_by_weights[k,2]/N1)*sorted_by_weights[k,1])
      T_k[k] = T_k[k] + ((sorted_by_weights[k,2]/N1)*(sorted_by_weights[k,1]**2))
      i = i - 1
    }
  }
  min_comp_sigma_weights = 0
  
  
  
  ######Test for Assumption 7
  assumption_7 = 0
  for(g in 1:G){
    for(t in 1:T){
      if(D_gt[g,t] > 0){
        assumption_7 = assumption_7 + ((N_gt[g,t]/N1)*(w[g,t]-1)*(LDelta_gt[g,t]-LDelta_TR))
      }
    }
  }
  
  
  
  
  if(sigma_weights > 0){
    mim_comp_sigma_weights = abs(beta_fe)/sigma_weights #####First lower bound for sigma_fe
    
    if(beta_fe != 0 && sorted_by_weights[count_treated_gt_cells,1] < 0){ #####Second lower bound for sigma_fe
      s = 1
      min2_comp_sigma_weights = 0
      for(i in 1:count_treated_gt_cells){
        if((sorted_by_weights[i,1] < (-S_k[i]/(1-P_k[i]))) && (sorted_by_weights[i,1] < sorted_by_weights[s,1])){
          s = i
        }
      }
      min2_comp_sigma_weights = abs(beta_fe) / sqrt(T_k[s] + (S_k[s]**2)/(1-T_k[s]))
      results <- data.frame(min_comp_sigma_weights, min2_comp_sigma_weights,assumption_7)
      return(results)
    } else {
      results <- data.frame(min_comp_sigma_weights, assumption_7)
      return(results)
    }
  } else {
    return(assumption_7)
  }
}


####################################################
######Estimate Beta_fe for both variants of w ######

est.beta_fe <- function(X,Y,N,G,T,d,v_weights = 1){ ##### estimates the decomposition result
  D_gt <- est.D_gt(X,Y,N,G,T,d)
  N_gt <- est.N_gt(X,N,G,T)
  N1 = est.N1(X,d)
  
  if(v_weights == 1){
    w <- est.weights_v1(X,Y,N,G,T,d) 
  }
  
  if(v_weights == 2){
    w <- est.weights_v2(X,Y,N,G,T,d)
  }
  LDelta_gt <- est.LDelta_gt(X,Y,N,G,T,d)
  
  
  beta_fe = 0
  for (t in 1:T){
    for (g in 1:G){
      beta_fe = beta_fe + (D_gt[g,t] * N_gt[g,t]/N1 * w[g,t] * LDelta_gt[g,t]) 
    }
  }
  return(beta_fe)
}



###################################################
###### Alternative Estimator ######################
###################################################

est.DID_M <- function(X,Y,N,G,T,d){
  DID_plus <- c()
  DID_minus <- c()
  N_gt = est.N_gt(X,N,G,T)
  D_gt <- est.D_gt(X,Y,N,G,T,d)
  Y_gt <- matrix(0, nrow = G, ncol = T)
  
  for(i in 1:N){
    for(t in 1:T){
      for(g in 2:G){
        if(X[(i-1)*T+t,(T-1)+(g-1)]==1){
          Y_gt[g,t] = Y_gt[g,t] + (1/N_gt[g,t])*Y[(i-1)*T+t]
        }
      }
      if(sum(X[(i-1)*T+t,T:(T+G-2)]) == 0){
        Y_gt[1,t] = Y_gt[1,t] + (1/N_gt[g,t])*Y[(i-1)*T+t]
      }
      if(sum(X[(i-1)*T+t,1:(T+G-2)]) == 0){
        Y_gt[1,1] = Y_gt[1,1] + (1/N_gt[g,t])*Y[(i-1)*T+t]
      }
    }
  }
  


  
  N_ddt <- matrix(0, nrow = 4, ncol = (T-1))
  
  #  N_ddt[case,t] = sum(N_gt[g,t]) conditional on g:D_gt = d, D_gt-1 = d' 
  #  counts switchers and nonswitchers in period t
  
  N_S = 0 ###Total number of switchers 
  
  for(g in 1:G){
    for(t in 1:(T-1)){
      if(D_gt[g,(t+1)] == 0 && D_gt[g,t] == 0){ #### case 1: d = d' = 0
        N_ddt[1,t] = N_ddt[1,t] + N_gt[g,(t+1)]
      }
      if(D_gt[g,(t+1)] == 0 && D_gt[g,t] > 0){ #### case 2: d = 0, d' = 1, leaver
        N_ddt[2,t] = N_ddt[2,t] + N_gt[g,(t+1)]
        N_S = N_S + N_gt[g,(t+1)]
      }
      if(D_gt[g,(t+1)] > 0 && D_gt[g,t] == 0){ #### case 3: d = 1, d' = 0, joiner
        N_ddt[3,t] = N_ddt[3,t] + N_gt[g,(t+1)]
        N_S = N_S + N_gt[g,(t+1)]
      }
      if(D_gt[g,(t+1)] > 0 && D_gt[g,t] > 0){ #### case 4: d = d' = 1
        N_ddt[4,t] = N_ddt[4,t] + N_gt[g,(t+1)]
      }
    }
  }
  
  
  Y_gt_sum <- matrix(0, nrow = 4, ncol = (T-1)) ##### collects N[g,t]/N[d,d',t] * (Y[g,t] - Y[g,t-1]) terms conditional on the appearance of each case.
  for(g in 1:G){
    for(t in 1:(T-1)){
      
      if(N_ddt[1,t] != 0 && D_gt[g,(t+1)] == 0 && D_gt[g,t] == 0){
          Y_gt_sum[1,t] = Y_gt_sum[1,t] + (N_gt[g,t]/N_ddt[1,t])*(Y_gt[g,(t+1)] - Y_gt[g,t])
      }
      
      if(N_ddt[2,t] != 0 && D_gt[g,(t+1)] == 0 && D_gt[g,t] > 0){
          Y_gt_sum[2,t] = Y_gt_sum[2,t] + (N_gt[g,t]/N_ddt[2,t])*(Y_gt[g,(t+1)] - Y_gt[g,t])
      }
      
      if(N_ddt[3,t] != 0 && D_gt[g,(t+1)] > 0 && D_gt[g,t] == 0){
          Y_gt_sum[3,t] = Y_gt_sum[3,t] + (N_gt[g,t]/N_ddt[3,t])*(Y_gt[g,(t+1)] - Y_gt[g,t])
      }
      
      if(N_ddt[4,t] != 0 && D_gt[g,(t+1)] > 0 && D_gt[g,t] > 0){
          Y_gt_sum[4,t] = Y_gt_sum[4,t] + (N_gt[g,t]/N_ddt[4,t])*(Y_gt[g,(t+1)] - Y_gt[g,t])
      }
    }
  }
  
  for(t in 1:(T-1)){ 
    if(N_ddt[3,t] != 0 && N_ddt[1,t] != 0){ #### estimates DID_[Plus,t]
      DID_plus[t] = Y_gt_sum[3,t] - Y_gt_sum[1,t]
    } else {
      DID_plus[t] = 0
    }
    
    if(N_ddt[4,t] != 0 && N_ddt[2,t] != 0){ #### estimates DID_[minus,t]
      DID_minus[t] = Y_gt_sum[4,t] - Y_gt_sum[2,t]
    } else {
      DID_minus[t] = 0
    }
  }
  
  DID_M = 0
  for(t in 1:(T-1)){
    if(DID_plus[t] != 0){
      DID_M = DID_M + (N_ddt[3,t]/N_S)*DID_plus[t]
    }
    if(DID_minus[t] != 0){
      DID_M = DID_M + (N_ddt[2,t]/N_S)*DID_minus[t]
    }
  }
  return(DID_M)
}


##################################
###### Placebo estimator #########
##################################


est.DID_M_PL <- function(X,Y,N,G,T,d){ ###funtions similar to DID_M, changes are commented
  D_gt <- est.D_gt(X,Y,N,G,T,d)
  N_gt <- est.N_gt(X,N,G,T)
  Y_gt <- matrix(0, nrow = G, ncol = T)
  
  for(i in 1:N){
    for(t in 1:T){
      for(g in 2:G){
        if(X[(i-1)*T+t,(T-1)+(g-1)]==1){
          Y_gt[g,t] = Y_gt[g,t] + (1/N_gt[g,t])*Y[(i-1)*T+t]
        }
      }
      if(sum(X[(i-1)*T+t,T:(T+G-2)]) == 0){
        Y_gt[1,t] = Y_gt[1,t] + (1/N_gt[g,t])*Y[(i-1)*T+t]
      }
      if(sum(X[(i-1)*T+t,1:(T+G-2)]) == 0){
        Y_gt[1,1] = Y_gt[1,1] + (1/N_gt[g,t])*Y[(i-1)*T+t]
      }
    }
  }
  
  N_dddt <- matrix(0,nrow = 6, ncol = T-2) ##### 6 cases,  d = d'' = 0, d' = 1 and d = d'' = 1, d' = 0 are left out
  
  N_S = 0 ###Total number of switchers 
  
  for(g in 1:G){
    for(t in 1:(T-2)){
      if(D_gt[g,(t+2)] == 0 && D_gt[g,(t+1)] == 0 && D_gt[g,t] == 0){ #### case 1: d = d' = d'' = 0
        N_dddt[1,t] = N_dddt[1,t] + N_gt[g,(t+1)]
      }
      if(D_gt[g,t+2] == 0 && D_gt[g,(t+1)] == 0 && D_gt[g,t] > 0){ #### case 2: d = d' = 0, d'' = 1 
        N_dddt[2,t] = N_dddt[2,t] + N_gt[g,(t+1)]
      }
      if(D_gt[g,t+2] == 0 && D_gt[g,(t+1)] > 0 && D_gt[g,t] > 0){ #### case 3: d = 0, d' = d'' = 1 (period d leaver)
        N_dddt[3,t] = N_dddt[3,t] + N_gt[g,(t+1)]
        N_S = N_S + N_gt[g,(t+1)]
      }
      if(D_gt[g,t+2] > 0 && D_gt[g,(t+1)] > 0 && D_gt[g,t] > 0){ #### case 4: d = d' = d'' = 1 
        N_dddt[4,t] = N_dddt[4,t] + N_gt[g,(t+1)]
      }
      if(D_gt[g,t+2] > 0 && D_gt[g,(t+1)] > 0 && D_gt[g,t] == 0){ #### case 5: d = d' = 1, d'' = 0
        N_dddt[5,t] = N_dddt[5,t] + N_gt[g,(t+1)]
      }
      if(D_gt[g,t+2] > 0 && D_gt[g,(t+1)] == 0 && D_gt[g,t] == 0){ #### case 6: d = 1, d' = d'' = 0 (period d joiner)
        N_dddt[6,t] = N_dddt[6,t] + N_gt[g,(t+1)]
        N_S = N_S + N_gt[g,(t+1)]
      }
    }
  }
  
  DID_plus <- c()
  DID_minus <- c()
  
  for(t in 1:(T-2)){
    DID_plus[t] = 0
    DID_minus[t] = 0
  }
  for(t in 3:T){
    for(g in 1:G){
      if(N_dddt[6,t-2] != 0 && N_dddt[1,t-2] != 0){ #### Estimates DID_[Plus,t-1] for period t joiners 
        if(D_gt[g,t] > 0 && D_gt[g,(t-1)] == 0 && D_gt[g,t-2] == 0){
          DID_plus[t-2] = DID_plus[t-2] + (N_gt[g,t]/N_dddt[6,t-2])*(Y_gt[g,t-1]-Y_gt[g,t-2])
        }
        if(D_gt[g,(t)] == 0 && D_gt[g,(t-1)] == 0 && D_gt[g,t-2] == 0){
          DID_plus[t-2] = DID_plus[t-2] - (N_gt[g,t]/N_dddt[1,t-2])*(Y_gt[g,t-1]-Y_gt[g,t-2])
        }
      } else {
        DID_plus[t-2] = 0
      }
      
      if(N_dddt[4,t-2] != 0 && N_dddt[3,t-2] != 0){ #### Estimates DID_[minus,t-1] for period t leavers
        if(D_gt[g,t] > 0 && D_gt[g,(t-1)] > 0 && D_gt[g,t-2] > 0){
          DID_minus[t-2] = DID_minus[t-2] - (N_gt[g,t]/N_dddt[4,t-2])*(Y_gt[g,t-1]-Y_gt[g,t-2])
        }
        if(D_gt[D_gt[g,t] == 0 && g,(t-1)] > 0 && D_gt[g,t-2] > 0){
          DID_minus[t-2] = DID_minus[t-2] - (N_gt[g,t]/N_dddt[3,t-2])*(Y_gt[g,t-1]-Y_gt[g,t-2])
        }
      } else {
        DID_minus[t-2] = 0
      }
    }
  }
  DID_M_PL = 0
  for(t in 3:T){
    DID_M_PL = DID_M_PL + ((N_dddt[6,t-2]/N_S)*DID_plus[t-2] + (N_dddt[3,t-2]/N_S)*DID_minus[t-2])
  }
  
  
  return(DID_M_PL)
}





