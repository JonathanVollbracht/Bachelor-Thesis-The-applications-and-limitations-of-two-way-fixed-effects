install.packages("stargazer")
library(stargazer)

##Let number of groups be g = 6, let number of periods be t = 5, let number of individuals be 300, individuals are evenly distributed throughout the 6 groups
set.seed(42)

########## Setup ###################
G = 6 
T = 10 
N = 300 
d = T+G-1 ##This refers to the collum which will contain our Treatment information within our modified datamatrix
beta_2 = 0.03 ##height effect
beta_3 = 1.1 ##experience effect


groupeffects <- c(2.5, 4, 3.5, 4, 2, 3.5)  ### True group effects vector
timeeffects <- c(8, 5, 8.7, 10, 12, 9.4, 7.8, 5.5, 6.5, 6) ### True time effects vector
treatmenteffect = 8 ### True treatment effect



###Construction of outcomevector
Randomdummy <- matrix(NA, nrow = N, ncol = T)
for(t in 1:T){ 
  Randomdummy[,] = c(rnorm(N*T, mean = 0, sd = 1.5))
}

Timedummy <- matrix(NA, nrow = N, ncol = T)
for(t in 1:T){
  Timedummy[,t]=timeeffects[t]
}
Groupdummy <- matrix(NA, nrow = N, ncol = T)
for(g in 1:G){
  j = (g-1)*(N/G)
for(i in 1:(N/G)){
  Groupdummy[j+i,]=groupeffects[g]
}  
}


controls_1 <- c(rnorm(N, mean = 175, sd = 10)) ##height
Controlsdummy_1 = controls_1*beta_2

controls_2 <- c(rnorm(N, mean = 8, sd = 1.5)) ##years of previous experience
Controlsdummy_2 = controls_2*beta_3


###### Construct homogeneous treatment effect data samples 1,2 and 3 ######

Treatmentdummy_S1 <- matrix(NA, nrow = N, ncol = T)
Treatmentdummy_S2 <- matrix(NA, nrow = N, ncol = T)
Treatmentdummy_S3 <- matrix(NA, nrow = N, ncol = T)

for(i in 1:(N/G)){ ##in all three samples group 1 serves as a control group
  
  Treatmentdummy_S1[i,] = 0
  
  Treatmentdummy_S2[i,] = 0
  
  Treatmentdummy_S3[i,] = 0
}


for(i in ((N/G)+1):(2*(N/G))){ ##in samples 1 and 2 group 2 also serves as a control group, in sample 3 group 2 enters treatment at period 9 
  
  Treatmentdummy_S1[i,] = 0
  
  Treatmentdummy_S2[i,] = 0
  
  Treatmentdummy_S3[i,1:8] = 0
  Treatmentdummy_S3[i,9:10] = 1
}


for(i in (2*(N/G)+1):(3*(N/G))){
  
  Treatmentdummy_S1[i,1:5] = 0
  Treatmentdummy_S1[i,6:10] = 1 ##for sample 1: group 3 joins the treatment in period 6
  
  Treatmentdummy_S2[i,1:5] = 0
  Treatmentdummy_S2[i,6:10] = 1 ##for sample 2: group 3 joins the treatment in period 6
  
  Treatmentdummy_S3[i,1:6] = 0
  Treatmentdummy_S3[i,7:10] = 1 ##for sample 3: group 3 joins the treatment in period 7
}


for(i in (3*(N/G)+1):(4*(N/G))){
  
  Treatmentdummy_S1[i,1:5] = 0
  Treatmentdummy_S1[i,6:10] = 1 ##for sample 1: group 4 joins the treatment in period 6
  
  Treatmentdummy_S2[i,1:5] = 0
  Treatmentdummy_S2[i,6:10] = 1 ##for sample 2: group 4 joins the treatment in period 6
  
  Treatmentdummy_S3[i,1:4] = 0
  Treatmentdummy_S3[i,5:10] = 1 ##for sample 3: group 4 joins the treatment in period 5
}


for(i in (4*(N/G)+1):(5*(N/G))){
  
  Treatmentdummy_S1[i,] = 1 ##for sample 1: group 5 is always treated
  
  Treatmentdummy_S2[i,1] = 0
  Treatmentdummy_S2[i,2:10] = 1 ##for sample 2: group 5 is treated in every period except for period 1
  
  Treatmentdummy_S3[i,1:2] = 0
  Treatmentdummy_S3[i,3:10] = 1 ##for sample 3: group 5 joins the treatment in period 3
}


for(i in (5*(N/G)+1):N){
  Treatmentdummy_S1[i,] = 1 ##for sample 1: group 6 is always treated
  
  Treatmentdummy_S2[i,1] = 0
  Treatmentdummy_S2[i,2:10] = 1 ##for sample 2: group 6 is treated in every period except for period 1
  
  Treatmentdummy_S3[i,1] = 0
  Treatmentdummy_S3[i,2:10] = 0 ##for sample 3: group 6 is treated in every period except for period 1
}

Treatmentdummy_S1 = Treatmentdummy_S1*treatmenteffect
Treatmentdummy_S2 = Treatmentdummy_S2*treatmenteffect
Treatmentdummy_S3 = Treatmentdummy_S3*treatmenteffect


Outcome_matrix_S1 <- matrix(Randomdummy, nrow = N, ncol = T)
Outcome_matrix_S1 = Outcome_matrix_S1 + Timedummy + Groupdummy + Treatmentdummy_S1
for (t in 1:T){
  Outcome_matrix_S1[,t] = Outcome_matrix_S1[,t] + Controlsdummy_1 + Controlsdummy_2
}

Outcome_matrix_S2 <- matrix(Randomdummy, nrow = N, ncol = T)
Outcome_matrix_S2 = Outcome_matrix_S2 + Timedummy + Groupdummy + Treatmentdummy_S2
for (t in 1:T){
  Outcome_matrix_S2[,t] = Outcome_matrix_S2[,t] + Controlsdummy_1 + Controlsdummy_2
}

Outcome_matrix_S3 <- matrix(Randomdummy, nrow = N, ncol = T)
Outcome_matrix_S3 = Outcome_matrix_S3 + Timedummy + Groupdummy + Treatmentdummy_S3
for (t in 1:T){
  Outcome_matrix_S3[,t] = Outcome_matrix_S3[,t] + Controlsdummy_1 + Controlsdummy_2
}

####### Restructuring the data ##########
Y_S1 <- c() ###Y = Output per worker, our dependent variable
Y_S2 <- c()
Y_S3 <- c()
for (i in 1:N){
  for (t in 1:T){
  Y_S1[(i-1)*T+t]=Outcome_matrix_S1[i,t]
  Y_S2[(i-1)*T+t]=Outcome_matrix_S2[i,t]
  Y_S3[(i-1)*T+t]=Outcome_matrix_S3[i,t]
  }
}




######Time and group information dummies for X######

Timedata <- matrix(0, nrow = N*T, ncol = (T-1))
for (i in 1:N){
  for (t in 2:T){
    Timedata[((i-1)*T+t),(t-1)]=1
  }
}

Groupdata <- matrix(0, nrow = N*T, ncol = (G-1))
for (g in 2:G){
  for (i in ((g-1)*(N/G)+1):(g*N/G)){
  Groupdata[((i-1)*T+1):(i*T),(g-1)] = 1
  }
}

######Treatment information collum for X######

Treatmentdata_S1 <- c()
for (i in 1:(4*(N/G)*T)){
  Treatmentdata_S1[i]=0
}
for (i in (2*(N/G)+1):(4*(N/G))){
  for (j in 6:10){
  Treatmentdata_S1[((i-1)*T+j)]=1
  }
}
Treatmentdata_S1[(4*(N/G)*T+1):(N*T)]=1




Treatmentdata_S2 <- c()
for (i in 1:(2*(N/G)*T)){
  Treatmentdata_S2[i]=0
}
for (i in (2*(N/G)+1):(4*(N/G))){
  for (j in 1:5){
    Treatmentdata_S2[((i-1)*T+j)]=0
  }
  for (j in 6:10){
    Treatmentdata_S2[((i-1)*T+j)]=1
  }
}
for (i in (4*(N/G)+1):N){
  Treatmentdata_S2[((i-1)*T+1)] = 0
  for (j in 2:10){
    Treatmentdata_S2[((i-1)*T+j)]=1
  }
}




Treatmentdata_S3 <- c()
for (i in 1:((N/G)*T)){
  Treatmentdata_S3[i]=0
}
for (i in ((N/G)+1):(2*(N/G))){
  for (j in 1:8){
    Treatmentdata_S3[((i-1)*T+j)]=0
  }
  for (j in 9:10){
    Treatmentdata_S3[((i-1)*T+j)]=1
  }
}
for (i in (2*(N/G)+1):(3*(N/G))){
  for (j in 1:6){
    Treatmentdata_S3[((i-1)*T+j)]=0
  }
  for (j in 7:10){
    Treatmentdata_S3[((i-1)*T+j)]=1
  }
}
for (i in (3*(N/G)+1):(4*(N/G))){
  for (j in 1:4){
    Treatmentdata_S3[((i-1)*T+j)]=0
  }
  for (j in 5:10){
    Treatmentdata_S3[((i-1)*T+j)]=1
  }
}
for (i in (4*(N/G)+1):(5*(N/G))){
  for (j in 1:2){
    Treatmentdata_S3[((i-1)*T+j)]=0
  }
  for (j in 3:10){
    Treatmentdata_S3[((i-1)*T+j)]=1
  }
}
for (i in (5*(N/G)+1):N){
  for (j in 1){
    Treatmentdata_S3[((i-1)*T+j)]=0
  }
  for (j in 2:10){
    Treatmentdata_S3[((i-1)*T+j)]=1
  }
}





#######Construction of X######
## our datamatrix containing all information concerning groups, periods, treatment and controls that belong to a certain observation


X_S1 <- matrix(NA, nrow = N*T, ncol = T + G + 1) ##Sample 1
X_S1[1:(N*T),1:(T-1)]=Timedata[,]
X_S1[1:(N*T),(T):(T+G-2)]=Groupdata[,]
X_S1[1:(N*T),d]=Treatmentdata_S1[]
for (i in 1:N){
  for (j in ((i-1)*T+1):(i*T)){
  X_S1[j,(T+G)] = controls_1[i]
  }
}
for (i in 1:N){
  for (j in ((i-1)*T+1):(i*T)){
    X_S1[j,(T+G+1)] = controls_2[i]
  }
}


X_S2 <- matrix(NA, nrow = N*T, ncol = T + G + 1) ##Sample 2
X_S2[1:(N*T),1:(T-1)]=Timedata[,]
X_S2[1:(N*T),(T):(T+G-2)]=Groupdata[,]
X_S2[1:(N*T),d]=Treatmentdata_S2[]
for (i in 1:N){
  for (j in ((i-1)*T+1):(i*T)){
    X_S2[j,(T+G)] = controls_1[i]
  }
}
for (i in 1:N){
  for (j in ((i-1)*T+1):(i*T)){
    X_S2[j,(T+G+1)] = controls_2[i]
  }
}


X_S3 <- matrix(NA, nrow = N*T, ncol = T + G + 1) ##Sample 3
X_S3[1:(N*T),1:(T-1)]=Timedata[,]
X_S3[1:(N*T),(T):(T+G-2)]=Groupdata[,]
X_S3[1:(N*T),d]=Treatmentdata_S3[]
for (i in 1:N){
  for (j in ((i-1)*T+1):(i*T)){
    X_S3[j,(T+G)] = controls_1[i]
  }
}
for (i in 1:N){
  for (j in ((i-1)*T+1):(i*T)){
    X_S3[j,(T+G+1)] = controls_2[i]
  }
}

####### TWFE Estimations #########

beta_hat_S1 <- est.beta(X_S1,Y_S1)
beta_hat_S2 <- est.beta(X_S2,Y_S2)
beta_hat_S3 <- est.beta(X_S3,Y_S3)

######Precision measurements for fixed effects estimates######

TWFE_beta_S1_deviation <- c()
for(t in 1:(T-1)){
  TWFE_beta_S1_deviation[t] = (timeeffects[1] + beta_hat_S1[t] - timeeffects[t+1])/timeeffects[t+1]
}
for(g in 1:(G-1)){
  TWFE_beta_S1_deviation[(T-1)+g] = (groupeffects[1] + beta_hat_S1[(T-1)+g] - groupeffects[g+1])/groupeffects[g+1]
}

avg_TE_deviation_S1 = mean(TWFE_beta_S1_deviation[1:(T-1)])
avg_GE_deviation_S1 = mean(TWFE_beta_S1_deviation[T:length(TWFE_beta_S1_deviation)])



TWFE_beta_S2_deviation <- c()
for(t in 1:(T-1)){
  TWFE_beta_S2_deviation[t] = (timeeffects[1] + beta_hat_S2[t] - timeeffects[t+1])/timeeffects[t+1]
}
for(g in 1:(G-1)){
  TWFE_beta_S2_deviation[(T-1)+g] = (groupeffects[1] + beta_hat_S2[(T-1)+g] - groupeffects[g+1])/groupeffects[g+1]
}

avg_TE_deviation_S2 = mean(TWFE_beta_S2_deviation[1:(T-1)])
avg_GE_deviation_S2 = mean(TWFE_beta_S2_deviation[T:length(TWFE_beta_S2_deviation)])



TWFE_beta_S3_deviation <- c()
for(t in 1:(T-1)){
  TWFE_beta_S3_deviation[t] = (timeeffects[1] + beta_hat_S3[t] - timeeffects[t+1])/timeeffects[t+1]
}
for(g in 1:(G-1)){
  TWFE_beta_S3_deviation[(T-1)+g] = (groupeffects[1] + beta_hat_S3[(T-1)+g] - groupeffects[g+1])/groupeffects[g+1]
}

avg_TE_deviation_S3 = mean(TWFE_beta_S3_deviation[1:(T-1)])
avg_GE_deviation_S3 = mean(TWFE_beta_S3_deviation[T:length(TWFE_beta_S3_deviation)])


###### Precision assessment of beta_hat ######
precision_beta_hat <- c()
precision_beta_hat[1] = (beta_hat_S1[d] - treatmenteffect)/treatmenteffect
precision_beta_hat[2] = (beta_hat_S2[d] - treatmenteffect)/treatmenteffect
precision_beta_hat[3] = (beta_hat_S3[d] - treatmenteffect)/treatmenteffect



###### Average treatment effect across all treated units (ATT) ##############

LDelta_TR_S1 <- est.LDelta_TR(X_S1,Y_S1,d,beta_hat_S1)
LDelta_TR_S2 <- est.LDelta_TR(X_S2,Y_S2,d,beta_hat_S2)
LDelta_TR_S3 <- est.LDelta_TR(X_S3,Y_S3,d,beta_hat_S3)

####### Average treatment effect in cell (g,t) ##########

LDelta_gt_S1 <- est.LDelta_gt(X_S1,Y_S1,N,G,T,d)
LDelta_gt_S2 <- est.LDelta_gt(X_S2,Y_S2,N,G,T,d)
LDelta_gt_S3 <- est.LDelta_gt(X_S3,Y_S3,N,G,T,d)

####### Small delta ######

delta_TR_S1 <- est.delta_TR(X_S1,Y_S1,N,G,T,d)
delta_TR_S2 <- est.delta_TR(X_S2,Y_S2,N,G,T,d)
delta_TR_S3 <- est.delta_TR(X_S3,Y_S3,N,G,T,d)

######Estimate beta_fe for both versions of w[g,t]#######

beta_fe_v1_S1 <- est.beta_fe(X_S1,Y_S1,N,G,T,d,1)
beta_fe_v1_S2 <- est.beta_fe(X_S2,Y_S2,N,G,T,d,1)
beta_fe_v1_S3 <- est.beta_fe(X_S3,Y_S3,N,G,T,d,1)



beta_fe_v2_S1 <- est.beta_fe(X_S1,Y_S1,N,G,T,d,2)
beta_fe_v2_S2 <- est.beta_fe(X_S2,Y_S2,N,G,T,d,2)
beta_fe_v2_S3 <- est.beta_fe(X_S3,Y_S3,N,G,T,d,2)


###### weights ######

weights_S1_v1 <- est.weights_v1(X_S1,Y_S1,N,G,T,d)
weights_S1_v2 <- est.weights_v2(X_S1,Y_S1,N,G,T,d)

weights_S2_v1 <- est.weights_v1(X_S2,Y_S2,N,G,T,d)
weights_S2_v2 <- est.weights_v2(X_S2,Y_S2,N,G,T,d)

weights_S3_v1 <- est.weights_v1(X_S3,Y_S3,N,G,T,d)
weights_S3_v2 <- est.weights_v2(X_S3,Y_S3,N,G,T,d)


###################################################
###### Heterogenous treatments effects setup ######
###################################################


###### DGP Modifications #######
treatmenteffect_HTE_S1 <- c(5,6,4,7.5,8,6,4,3,20)


Treatmentdummy_HTE_S1 <- matrix(NA, nrow = N, ncol = T)
for(i in 1:(1*(N/G))){ ##group 1 and 2 are control groups
  Treatmentdummy_HTE_S1[i,] = 0
}
for(i in (1*(N/G)+1):(2*(N/G))){
  Treatmentdummy_HTE_S1[i,1:9] = 0
  Treatmentdummy_HTE_S1[i,10] = treatmenteffect_HTE_S1[9]
}
for(i in (2*(N/G)+1):(3*(N/G))){
  Treatmentdummy_HTE_S1[i,1:8] = 0
  Treatmentdummy_HTE_S1[i,9:10] = treatmenteffect_HTE_S1[8:9]
}
for(i in (3*(N/G)+1):(4*(N/G))){
  Treatmentdummy_HTE_S1[i,1:8] = 0
  Treatmentdummy_HTE_S1[i,9:10] = treatmenteffect_HTE_S1[8:9]
}
for(i in (4*(N/G)+1):(5*(N/G))){
  Treatmentdummy_HTE_S1[i,1:7] = 0
  Treatmentdummy_HTE_S1[i,8:10] = treatmenteffect_HTE_S1[7:9]
}
for(i in (5*(N/G)+1):N){
  Treatmentdummy_HTE_S1[i,1] = 0
  Treatmentdummy_HTE_S1[i,2:10] = treatmenteffect_HTE_S1[]
}




Outcome_matrix_HTE_S1 <- matrix(Randomdummy, nrow = N, ncol = T)
Outcome_matrix_HTE_S1 = Outcome_matrix_HTE_S1 + Timedummy + Groupdummy + Treatmentdummy_HTE_S1
for (t in 1:T){
  Outcome_matrix_HTE_S1[,t] = Outcome_matrix_HTE_S1[,t] + Controlsdummy_1 + Controlsdummy_2
}

Y_HTE_S1 <- c(NA) ###Y = Output per worker, our dependent variable
for (i in 1:N){
  for (t in 1:T){
    Y_HTE_S1[(i-1)*T+t]=Outcome_matrix_HTE_S1[i,t]
  }
}

X_HTE_S1 = X_S1

X_HTE_S1[,d] = 0
for(i in 1:N){
  for(t in 1:T){
    if(Treatmentdummy_HTE_S1[i,t] > 0){
      X_HTE_S1[((i-1)*T+t),d] = 1
    }
  }
}


###### TWFE Estimations ######

beta_HTE_S1 <- est.beta(X_HTE_S1, Y_HTE_S1)


###calculate deviation from true values for beta_hat estimates for HTE sample 
TWFE_beta_HTE_S1_deviation <- c()
for(t in 1:(T-1)){
  TWFE_beta_HTE_S1_deviation[t] = (timeeffects[1] + beta_HTE_S1[t] - timeeffects[t+1])/timeeffects[t+1]
}
for(g in 1:(G-1)){
  TWFE_beta_HTE_S1_deviation[(T-1)+g] = (groupeffects[1] + beta_HTE_S1[(T-1)+g] - groupeffects[g+1])/groupeffects[g+1]
}

HTE_avg_TE_deviation_S1 = mean(TWFE_beta_HTE_S1_deviation[1:(T-1)]) 
HTE_avg_GE_deviation_S1 = mean(TWFE_beta_HTE_S1_deviation[T:length(TWFE_beta_HTE_S1_deviation)])

###### Decomposition estimates #######

beta_FE_HTE_V1_S1 = est.beta_fe(X_HTE_S1,Y_HTE_S1,N,G,T,d,1) 
beta_FE_HTE_V2_S1 = est.beta_fe(X_HTE_S1,Y_HTE_S1,N,G,T,d,2)


weights_S1_v1 <- est.weights_v1(X_HTE_S1,Y_HTE_S1,N,G,T,d)
weights_S1_v2 <- est.weights_v2(X_HTE_S1,Y_HTE_S1,N,G,T,d)

###### Inference ######


r_squared_S1 = est.r2(X_S1,Y_S1,beta_hat_S1)
r_squared_S2 = est.r2(X_S2,Y_S2,beta_hat_S2)
r_squared_S3 = est.r2(X_S3,Y_S3,beta_hat_S3)

r2_HTE_S1 = est.r2(X_HTE_S1,Y_HTE_S1,beta_HTE_S1)


### Test null hypothesis beta_hat_fe = 0, returns 1 for 5%, 2 for 1% and 3 for 0.1% confidence
test_t_stat(X_S1,Y_S1,beta_hat_S1, 0, d)
test_t_stat(X_S2,Y_S2,beta_hat_S2, 0, d)
test_t_stat(X_S3,Y_S3,beta_hat_S3, 0, d)

test_t_stat(X_HTE_S1,Y_HTE_S1,beta_HTE_S1,0,d)


### Tests Beta_fe estimates for HTE sample on robustness to treatment effect heterogeneity
inference(X_HTE_S1,Y_HTE_S1,N,G,T,d,2)


### Sigma(w) for the HTE sample
HTE_S1_Sigma_Weights_V1 <- est.sigma_weights(X_HTE_S1,Y_HTE_S1,N,G,T,d,1)
HTE_S1_Sigma_Weights_V2 <- est.sigma_weights(X_HTE_S1,Y_HTE_S1,N,G,T,d,2)


### Sigma(Delta) for the HTE sample
HTE_S1_Sigma_Delta <- est.sigma_delta(X_HTE_S1,Y_HTE_S1,N,G,T,d)


###############################################################
###### New Estimator by De Chaisemartin and Hautefoeille ######
###############################################################

DID_M_S1 = est.DID_M(X_S1,Y_S1,N,G,T,d)
DID_M_S2 = est.DID_M(X_S2,Y_S2,N,G,T,d)
DID_M_S3 = est.DID_M(X_S3,Y_S3,N,G,T,d)

DID_M_HTE_S1 = est.DID_M(X_HTE_S1, Y_HTE_S1, N, G, T, d)



DID_PL_S1 = est.DID_M_PL(X_S1,Y_S1,N,G,T,d)
DID_PL_S2 = est.DID_M_PL(X_S2,Y_S2,N,G,T,d)
DID_PL_S3 = est.DID_M_PL(X_S3,Y_S3,N,G,T,d)

DID_PL_HTE_S1 = est.DID_M_PL(X_HTE_S1,Y_HTE_S1,N,G,T,d)





