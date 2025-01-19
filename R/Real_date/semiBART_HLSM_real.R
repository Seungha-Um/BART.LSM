library(statnet)
library(MASS)
library(reshape2)
library(dplyr)
library(HLSM)
library(zeallot) 
library(SoftBart)
library(gdata)
library(rlist)
library(ggplot2)
library(ggpubr)
library(fastDummies)

options(dplyr.summarise.inform = FALSE)
'%<-%' <- zeallot::'%<-%' 
load("APMathBeliefs.Rdata")
source("~/quantile_normalize.R")
expit <- function(x) (1/(1 + exp(-x)))

X1.f1 = list()
X2.f1 = list()
X1.f2 = list()
X2.f2 = list()

for(k in 1:14){
  # 2013 scores for Q1 - Q5
  X1.f1[[k]] = ap.covs[[k]][,2]
  # 2013 scores for Q6 - Q10
  X1.f2[[k]] = ap.covs[[k]][,3]
  # 2012 scores for Q1 - Q5
  X2.f1[[k]] = ap.covs[[k]][,4]
  # 2012 scores for Q6 - Q10
  X2.f2[[k]] = ap.covs[[k]][,5]  
}


#sch_ix <- c(2, 5, 6, 7, 8, 10, 11, 12, 14)
#sch_ix <- c(2, 4, 5, 6, 7, 9, 11, 12, 13, 14)
sch_ix <- c(1:14)
sch_ix <- sch_ix[-c(7,13,14)]

# 2013 scores for Q1 - Q5 
YY1 <-  X1.f1[sch_ix] 
# 2012 scores for Q1 - Q5 
TT1 <- X2.f1[sch_ix]  #Treatment

# # 2013 scores for Q6 - Q10 
YY2 <-  X1.f2[sch_ix]
# # 2012 scores for Q6 - Q10
TT2 <- X2.f2[sch_ix]  #Treatment

A = Sociolist1[sch_ix] 

AA1 <- A #direct 

num_net = length(A)
A_undirect <- list()
for(k in 1:num_net) {
  temp <- A[[k]] + t(A[[k]])
  temp[temp == 2 ] <- 1 
  A_undirect[[k]] <- temp 
}

AA2 <- A_undirect

num_tree <- 200
num_iter <- num_burn <- 5000
sim_num <- num_iter + num_burn

dir_type <- c(1, 2)
or_type <- c(1, 2)
Q_type <- c(1, 2)
set.seed(2023)

d <- 2
q <- 2
o <-1 

for(q in 1:length(Q_type)){
  for(o in 1:length(or_type)){
    for(d in 1:length(dir_type)){
      
      if(d == 1){AA = AA1
      }else{AA = AA2}
      
      if(q == 1){YY = YY1;  TT = TT1
      }else{YY = YY2;  TT = TT2}
      
      # degree2 <- lapply(AA2, rowSums)
      # ix <- lapply(degree2, function(x) x!=0 )
      # 
      # AA = sapply(1:num_net, function(xx) AA[[xx]][ix[[xx]], ix[[xx]]] )
      # degree <- lapply(AA, rowSums)
      # YY = sapply(1:num_net, function(xx) YY[[xx]][ix[[xx]]] )
      # TT = sapply(1:num_net, function(xx) TT[[xx]][ix[[xx]]] )

      num_net = length(AA)
      degree <- lapply(AA, rowSums)
      num_nodes = sapply(YY, length)
      N <- sum(num_nodes)
      
      if(o == 1){
        nei_treat <- list()
        for(j in 1:num_net){
          T_temp <- TT[[j]]
          A_temp <- AA[[j]]
          nodes_temp <- num_nodes[j]
          degree_temp <- unlist(degree[j])
          nei_temp <- sapply(1:nodes_temp, function(k) sum(T_temp[A_temp[k,] == 1])/degree_temp[k])
          nei_temp[is.na(nei_temp)] <- 0
          nei_treat[[j]] <- nei_temp
        }
      }else{
        nei_treat_2nd <- list()
        for(j in 1:num_net){
          T_temp <- TT[[j]]
          A_temp <- AA[[j]]
          nodes_temp <- num_nodes[j]
          nei_vec <- rep(NA, nodes_temp)
          for(k in 1:nodes_temp){
            ind <- A_temp[k,] == 1   # 1st order friends
            ind2_mat <- A_temp[ind,]     # Extract 2nd order friends' adj matrix
            ind_all <- colSums(rbind(ind, ind2_mat)) != 0
            ind_all[k] <- FALSE
            nei_vec[k] <- mean(T_temp[ind_all])
            nei_vec[is.na(nei_vec)] <- 0
          }
          nei_treat_2nd[[j]] <- nei_vec
        }
        nei_treat <- nei_treat_2nd
      }
      
      ### Latent location estimation
      priors = list(MuBeta = c(0,0), VarBeta = c(100,100), MuZ = c(0,0),
                    VarZ = c(5,5), PriorA =  100, PriorB = 150)
      
      mod_complete = HLSMrandomEF(Y = AA, tuneIn = TRUE, dd = 2, niter = sim_num,
                                  priors = priors)
      
      ZZlist = getLS(mod_complete, burnin = 0, thin = 1)
      est_median <- sapply(1:num_net, function(k) apply(ZZlist[[k]], c(1,2), median))
      
      
      temp <- data.frame(nei = unlist(nei_treat)) %>% filter(nei != 0 ) %>%
        summarize(nei_min = quantile(nei, 0.025),
                  nei_max = quantile(nei, 0.975))
      nei_min <- temp$nei_min
      nei_max <- temp$nei_max
      
      result_list <- list()
      nei_vec <- seq(nei_min, nei_max, length.out = 20)
      
      test_length <- length(nei_vec) * N
      grp_vec <- rep(1:num_net, num_nodes)
      dummy_grp <- fastDummies::dummy_cols(grp_vec)[,-1]
      
      Y <- unlist(YY)
      Y_scaled <- scale(unlist(YY))
      c(Y_dims, center_Y, scale_Y) %<-% attributes(Y_scaled)
      temp_train <- as.data.frame(cbind(unlist(nei_treat), unlist(TT)))
      colnames(temp_train) <- c("nei_treat", "Treat")
      temp_train <- quantile_normalize_bart(temp_train)
      hypers <- Hypers(X = temp_train, Y = Y_scaled, normalize_Y = FALSE)
      hypers$num_tree <- num_tree
      opts <- Opts(update_s = FALSE)
    
      forest <- MakeForest(hypers, opts, warn = FALSE)
      step_size <- 1000
      mu_hat <- matrix(NA, sim_num, N)
      mu_test_hat <- matrix(NA, sim_num/step_size, length(nei_vec)*N)
      
      ave_mu_test <- list()
      ave_mu_test_all <- list()
  
      dim_beta <- 3
      beta_mat <- matrix(1, num_net, dim_beta)
      
      ################### LSM-BART
      df_train <- cbind(unlist(nei_treat), unlist(TT))
      colnames(df_train) <- c("nei_treat", "Treat")
      #df_train <- quantile_normalize_bart(df_train)
      
      df_test <- data.frame(rep(nei_vec, each = N), unlist(TT))
      colnames(df_test) <- c("nei_treat", "Treat")
      #df_test <- quantile_normalize_bart(df_test)
      
      XXtest <- quantile_normalize(df_train, df_test)
      df_train <- XXtest$X
      df_test <- XXtest$test_X
      
      df_test <- as.matrix(df_test)
      beta_mat_list <- list()
      Z_beta_list <- list()
      Sigma_vec <- rep(NA, sim_num)
    
      est_Z <-  sapply(1:num_net, function(k) ZZlist[[k]][,,1])
      est_Z <- lapply(est_Z, scale)
      
      Z_beta <- list()
      for(g in 1:num_net){
        Z_beta[[g]] <- cbind(1, est_Z[[g]]) %*% beta_mat[g,]
      }
  
      for(s in 1:sim_num){
        
        est_Z <- sapply(1:num_net, function(k) ZZlist[[k]][,,s])
        est_Z <- lapply(est_Z, scale)
        
        R = Y_scaled - unlist(Z_beta)
        mu_hat[s,] <- forest$do_gibbs(df_train, R, df_train, 1)
        mu_test_temp <- forest$do_predict(df_test) 
        Sigma <- forest$get_sigma()
        Sigma_vec[s] <- Sigma
        delta <- Y_scaled - mu_hat[s,] 
        
        beta_mat_new <- matrix(NA, num_net, 3)
        
        for(g in 1:num_net){
          ZZ_g <- cbind(1, est_Z[[g]])
          Cov <- solve((1/Sigma^2) * t(ZZ_g) %*% ZZ_g)
          mu_vec <- (1/Sigma^2) * Cov %*%  t(ZZ_g) %*% delta[grp_vec==g]
          beta_mat_new[g,] <- mvrnorm(1, mu_vec, Cov)  
        }
        
        beta_mat <- beta_mat_new
        beta_mat_list[[s]] <- data.frame(beta_mat = beta_mat,
                                         iter = s,
                                         grp = 1:num_net)
        
        Z_beta <- list()
        for(g in 1:num_net){
          Z_beta[[g]] <- cbind(1, est_Z[[g]]) %*% beta_mat[g,]
        }
        
        Z_beta_list[[s]] <- data.frame(Z_beta = list.rbind(Z_beta), grp = grp_vec, iter = s)
        
        EY_est <- data.frame(f_hat = mu_test_temp, linear = list.rbind(Z_beta) ) %>% 
          mutate(Y_hat = f_hat + linear) %>% pull(Y_hat)
        
        EY_est <- EY_est * scale_Y + center_Y  
        
        if(s %% step_size == 0) mu_test_hat[(s/step_size),] <- EY_est
        
        # To save estimates by group, creat a new df
        df_test2 <- data.frame(rep(nei_vec, each = N), cbind(unlist(TT), grp_vec))
        colnames(df_test2) <- c("nei_treat", "Treat", "grp")
        df_temp <- data.frame(df_test2, est_Y = EY_est)
        
        ave_mu_test[[s]] <- data.frame(df_temp %>% group_by(nei_treat, grp) %>% 
                                         summarize(mean = mean(est_Y)) %>%  mutate(iter = s))
        
        ave_mu_test_all[[s]] <- data.frame(df_temp %>% group_by(nei_treat) %>% 
                                             summarize(mean = mean(est_Y)) %>%  mutate(iter = s))
        
        if(s %% 100 == 0) cat("MCMC iter :", s, "or:", o,"dir:", d, "Q:", q, "\n") 
      }
      
      mu_hat_list <- mu_hat * scale_Y + center_Y
      ave_mu_test_list <- list.rbind(ave_mu_test)
      ave_mu_test_all_list <- list.rbind(ave_mu_test_all)
      mu_test_hat_list <- mu_test_hat
      
      out <- list()
      out$mu_hat_list <- mu_hat_list
      out$ave_mu_test_list <- ave_mu_test_list
      out$ave_mu_test_all_list <- ave_mu_test_all_list
      out$mu_test_hat_list <- mu_test_hat_list
      out$obs <- as.data.frame(cbind(treat = unlist(TT), nei = unlist(nei_treat), 
                                     Y = unlist(YY), grp = grp_vec))
      out$ZZ <-  list.rbind(Z_beta_list) %>% group_by(iter, grp) %>% 
                  summarize(meanZ = mean(Z_beta))
    
      saveRDS(out, paste0("semiBART", "Q", q, "Or", o, "dir", d, ".rds"))
    }  
  }
}
    
