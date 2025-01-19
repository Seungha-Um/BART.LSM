rm(list=ls())
library(statnet)
library(MASS)
library(reshape2)
library(dplyr)
library(HLSM)
library(zeallot) 
library(SoftBart)
library(stan4bart)
library(gdata)
library(rlist)
library(ggplot2)
library(ggpubr)
library(BART)
library(fastDummies)
source("~/quantile_normalize.R")

options(dplyr.summarise.inform = FALSE)
expit <- function(x) (1/(1 + exp(-x)))
'%<-%' <- zeallot::'%<-%' 

num_iter <- num_burn <- 5000
sim_num <- num_iter + num_burn
num_tree <- 200

net_vec <- seq(2, 14, 2)
rep_iter <- 10

set.seed(2023)
for(rr in 1:rep_iter){
  
  set.seed(rr)
  for(q in 1:length(net_vec)){
                              
    num_grp <- net_vec[q]
    num_ind <- rep(40, num_grp)
    N <- sum(num_ind)
    grp_vec <- rep(1:num_grp, num_ind)
    dummy_grp <- fastDummies::dummy_cols(grp_vec)[,-1]
    
    p <- 2
    
    X <- matrix(NA, N, p)
    AA <- list()
    XX <- list()
    
    homo_num <- 4
    mu_min <- -1.5      
    mu_max <- 1.5      
    mu_grp <- rbind(c(mu_max, mu_max), c(mu_min,mu_min), c(mu_max,mu_min), c(mu_min,mu_max))
    
    homo_temp <- rep(1:homo_num, each = num_ind[1]/homo_num)
    homo_grp <- sapply(1:num_grp, function(k) sample(homo_temp, num_ind[k], replace=F))
    for(a in 1:num_grp){
      degree_sum <- 1
      nn <- num_ind[a]
      homo_vec <- homo_grp[,a]
      while(degree_sum != 0) {  
        mu_mat <- mu_grp[homo_vec,]
        X <- t(apply(mu_mat, 1, function(x) mvrnorm(1, mu = x, Sigma = diag(0.3, p))))
        dist_X <- as.matrix(dist(X[,1:2], method = 'euclidean', diag = FALSE, upper = TRUE))
        temp <- 3 - 2 * dist_X
        A_p <- expit(temp)
        A <- sapply(1:nn, function(i) rbinom(nn, 1, A_p[i,]))
        lowerTriangle(A) <- upperTriangle(A, byrow=TRUE)
        diag(A) <- 0
        degree <- rowSums(A)
        degree_sum <- sum(degree == 0)
        if(sum(degree == 0)){print("isolated node exists")}
      }
      AA[[a]] <- A
      XX[[a]] <- scale(X)
    }
    
    X <- list.rbind(XX)
    degree_list <- lapply(AA, rowSums)
    
    Treat_fun <-  function(xx){5 + 0.5 * xx[,1] + 0.5 * xx[,2]  +
                  rnorm(length(xx[,1]), mean = 0, sd = 1)}
    
    set.seed(rr)
    Treat <- lapply(XX, Treat_fun)
    
    nei_treat <- list()
    for(j in 1:num_grp){
      T_temp <- Treat[[j]] 
      A_temp <- AA[[j]]
      nodes_temp <-  num_ind[j]
      degree_temp <- degree_list[[j]]
      temp_nei <- sapply(1:nodes_temp, function(k) sum(T_temp[A_temp[k,] == 1])/degree_temp[k]) 
      nei_treat[[j]] <- temp_nei
    }
    
    # c(nei_min, nei_max) %<-% quantile(unlist(nei_treat), c(0.025, 0.975))
    # nei_vec <- seq(nei_min, nei_max, length.out = 10)
    nei_vec <- seq(4, 6, length = 10)
    
    
    mu_fun <- function(Treat, G, X){
      Y_mu <- 10 + 3 * Treat + 10 * expit( 10 * (G - nei_vec[5]) ) + 10 * X[,1] + 10 * X[,2]
      return(Y_mu)
    }
    
    Y <-  rnorm(N, mu_fun(unlist(Treat), unlist(nei_treat), X), 1)
    
    true_potential <- 10 + 3 * 5 + 10 * expit( 10 * (nei_vec - nei_vec[5]) ) 
  
    Y_scaled <- scale(Y)
    c(Y_dims, center_Y, scale_Y) %<-% attributes(Y_scaled)
    temp_train <- as.data.frame(cbind(unlist(nei_treat), unlist(Treat)))
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
    ACPI_df <- list()
    
    dim_beta <- 1
    beta_mat <- matrix(1, num_grp, dim_beta)
    
    ################### LSM-BART
    df_train <- cbind(unlist(nei_treat), unlist(Treat))
    colnames(df_train) <- c("nei_treat", "Treat")
    
    df_test <- data.frame(rep(nei_vec, each = N), unlist(Treat))
    colnames(df_test) <- c("nei_treat", "Treat")
    
    XXtest <- quantile_normalize(df_train, df_test)
    df_train <- XXtest$X
    df_test <- XXtest$test_X
    df_test <- as.matrix(df_test)
    
    beta_mat_list <- list()
    Z_beta_list <- list()
    Sigma_vec <- rep(NA, sim_num)
    
    Z_beta <- as.matrix(dummy_grp) %*% beta_mat
    
    for(s in 1:sim_num){
      
      R = Y_scaled - unlist(Z_beta)
      mu_hat[s,] <- forest$do_gibbs(df_train, R, df_train, 1)
      mu_test_temp <- forest$do_predict(df_test) 
      Sigma <- forest$get_sigma()
      Sigma_vec[s] <- Sigma
      delta <- Y_scaled - mu_hat[s,] 
      
      beta_mat_new <- matrix(NA, num_grp, 1)
      
      for(g in 1:num_grp){
        Cov <- 1/(num_ind[g]/Sigma^2)
        mu_vec <-  sum(delta[grp_vec==g])/Sigma^2
        beta_mat_new[g,] <- rnorm(1, Cov*mu_vec, sd = Cov)  
      }
    
      beta_mat <- beta_mat_new
      beta_mat_list[[s]] <- data.frame(beta_mat = beta_mat,
                                       iter = s,
                                       grp = 1:num_grp)
      
      Z_beta <- as.matrix(dummy_grp) %*% beta_mat
      
      Z_beta_list[[s]] <- data.frame(Z_beta = list.rbind(Z_beta), grp = grp_vec, iter = s)
      
      EY_est <- data.frame(f_hat = mu_test_temp, linear = list.rbind(Z_beta) ) %>% 
        mutate(Y_hat = f_hat + linear) %>% pull(Y_hat)
      
      EY_est <- EY_est * scale_Y + center_Y  
      
      if(s %% step_size == 0) mu_test_hat[(s/step_size),] <- EY_est
      
      # To save estimates by group, creat a new df
      df_test2 <- data.frame(rep(nei_vec, each = N), cbind(unlist(Treat), grp_vec))
      colnames(df_test2) <- c("nei_treat", "Treat", "grp")
      df_temp <- data.frame(df_test2, est_Y = EY_est)
      
      ave_mu_test[[s]] <- data.frame(df_temp %>% group_by(nei_treat, grp) %>% 
                                       summarize(mean = mean(est_Y)) %>%  mutate(iter = s))
      
      ave_mu_test_all[[s]] <- data.frame(df_temp %>% group_by(nei_treat) %>% 
                                           summarize(mean = mean(est_Y)) %>%  mutate(iter = s))
      
      ACPI_df[[s]] <-  ave_mu_test_all[[s]]  %>%  
                        mutate(foo = lead(mean),
                               IE = foo - mean,
                               iter = s)
      if(s %% 100 == 0) cat("MCMC iter :", s, "net_num:", net_vec[q], "\n") 
    }
    
    mu_hat_list <- mu_hat
    ave_mu_test_list <- list.rbind(ave_mu_test)
    ave_mu_test_all_list <- list.rbind(ave_mu_test_all)
    mu_test_hat_list <- mu_test_hat
    
    # save all outputs
    # out <- list()
    # out$scale_Y <- scale_Y
    # out$center_Y <- center_Y 
    # out$Z_beta_list <- Z_beta_list
    # out$beta_mat_list <- beta_mat_list
    # out$Sigma_vec <- Sigma_vec
    # out$mu_hat_list <- mu_hat_list
    # out$homo_grp <- homo_grp
    # out$mu_test_hat_list <- mu_test_hat_list
    # out$ave_mu_test_list <- ave_mu_test_list 
    # out$true_potential <- true_potential
    # out$Y <- Y
    # out$Treat_list <-  Treat
    # out$nei_list <- nei_treat
    # out$num_grp <- num_grp 
    # out$N <- N
    # out$ave_mu_test_all_list <- ave_mu_test_all_list
    # out$X <- X
    # 
    out1 <- ave_mu_test_list %>% filter(iter > num_burn) %>% group_by(nei_treat, grp) %>%
            summarize(est_mean = mean(mean),
                      LI = quantile(mean, 0.025),
                      UI = quantile(mean, 0.975)) 
    
    df_True <- data.frame(nei_treat = nei_vec, true_po = true_potential) 
    
    out2 <- ave_mu_test_all_list %>% filter(iter > num_burn) %>% group_by(nei_treat) %>%
            summarize(est_mean = mean(mean),
                      LI = quantile(mean, 0.025),
                      UI = quantile(mean, 0.975)) %>% left_join(df_True, by = "nei_treat")
          
    true_PI <- df_True %>% 
                mutate(foo2 = lead(true_po), 
                       IE_true = foo2 - true_po) %>% 
                filter(!is.na(foo2))
    
    out3 <- list.rbind(ACPI_df) %>% 
            filter(iter > num_burn, !is.na(IE)) %>% group_by(nei_treat) %>% 
            summarize(PI_mean = mean(IE),
                      PI_LI = quantile(IE, 0.025),
                      PI_UI = quantile(IE, 0.975)) %>%
            left_join(true_PI, by = "nei_treat") %>%
            mutate(leng = PI_UI - PI_LI,
                   cov = PI_LI <= IE_true & IE_true <= PI_UI,
                   SE = (IE_true - PI_mean)^2)
    
    out <- list()
    out$out1 <- out1
    out$out2 <- out2
    out$out3 <- out3
    saveRDS(out, paste0("BART_incorrect_net", net_vec[q], "iter", rr, ".rds"))
  }
}

