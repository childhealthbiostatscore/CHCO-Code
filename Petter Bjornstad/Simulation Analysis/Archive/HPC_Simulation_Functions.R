sim.fxn <- function(i) {
  j <- sim.par$scenario[i]
  p_e_causal_value <- sim.par$P.e.causal[i]
  p_s_causal_value <- P.s.causal.all[[j]]
  current_n <- sim.par$N[i]
  OR_exposure <- sim.par$OR.exposure[i]
  
  # Simulate exposures with the correlation structure
  # Function to create a correlation matrix with a given correlation value
  create_corr_matrix <- function(P.e, Corr) {
    # Create a matrix filled with Corr
    corr_matrix <- matrix(Corr, nrow = P.e, ncol = P.e)
    
    # Set diagonal elements to 1
    diag(corr_matrix) <- 1
    
    return(corr_matrix)
  }
  
  # Example usage
  Corr <- sim.par$Corr[[1]]  # Desired correlation value
  corr_matrix <- create_corr_matrix(P.e, Corr)
  X <- mvrnorm(n = current_n, mu = rep(0, P.e), Sigma = corr_matrix)
  beta.X <- c(rep(log(OR_exposure), p_e_causal_value), rep(0, (P.e - p_e_causal_value)))
  profiles <- rbind(rep(-0.5, P.e), rep(0.5, P.e))  # for counterfactual contrasts for mixture effects
  
  Y <- do.call(cbind, lapply(1:P.s, FUN=function(v) {
    if(v %in% p_s_causal_value) { etaY <- log(Y.freq[v]/(1-Y.freq[v])) + scale(X, scale = TRUE, center = TRUE)%*%beta.X}
    if(!v %in% p_s_causal_value) { etaY <- log(Y.freq[v]/(1-Y.freq[v])) }
    ProbI_nonzero <- exp(etaY)/(1+exp(etaY))
    # ProbI_zero <- exp(etaY)/(1+exp(etaY))
    ProbI_zero=1-ProbI_nonzero
    I_zero <- rbinom(current_n, 1, ProbI_zero)
    # simulate the negative binomial part
    if(v %in% p_s_causal_value) { eta_mu <- log(Y.mean[v]) + scale(X, scale = TRUE, center = TRUE)%*%beta.X }
    if(!v %in% p_s_causal_value) { eta_mu <- log(Y.mean[v]) }
    #ProbY <- 1/(1+alpha.nb/exp(eta_mu))
    # Debugging print statements
    # print(paste("Species:", v, "etaY:", etaY, "eta_mu:", eta_mu))
    
    Y <- rnbinom(current_n, size=phi, mu=exp(eta_mu))
    # combine the two parts
    Y[I_zero==1] <- 0
    return(Y)
  }))  
  
  sim.output <- list(X = X, Y = Y, beta.X = beta.X, profiles = profiles,
                     p_e_causal_value = p_e_causal_value, p_s_causal_value = p_s_causal_value,
                     N = current_n, OR_exposure = OR_exposure)
  
  # # Calculate the sample correlation matrix
  # sample_corr_matrix <- cor(X)
  # 
  # # Print the correlation matrix
  # print(sample_corr_matrix)
  # # Check the correlation between the first and second column
  # cor(X[, 1], X[, 2])
  
  return(sim.output)
}

# 2. Zero-Inflated Negative Binomial Regression Model with G-Computation ----
ZING_Model <- function(r,dat,exposure) {
  ZING.results <- data.frame()
  poisson.results <- data.frame()
  if(min(dat[[r]])==0){
    tryCatch({
      m0 <- as.formula(paste0(r, "~", paste(exposure, collapse = " + ")))
      mod1 <- qgcomp.zi.noboot(f=m0, expnms = c(exposure),
                               data=dat, q=NULL, dist="negbin")
      #Means component
      #Mixture estimate
      mean.est <- mod1$coef[[1]][[2]]
      mean.pval <- mod1$pval[[1]][[2]]
      
      
      mean.se <- summary(mod1)$coeffients$count[2,2]
      
      
      #Individual estimates
      X_values <- data.frame(
        X.mean = summary(mod1$fit)[1]$coefficients$count[2:(length(exposure)+1)],
        X.se = summary(mod1$fit)[1]$coefficients$count[,2][2:(length(exposure)+1)],
        X.p = summary(mod1$fit)[1]$coefficients$count[,4][2:(length(exposure)+1)]
      )
      
      psi.mean <- data.frame(Taxa = r, Exposure = "Mixture", Mean = mean.est, SE = mean.se, P_Value = mean.pval,Component="Means")
      x.mean <- X_values
      x.mean <- x.mean %>% 
        mutate(Taxa = r,
               Exposure = rownames(x.mean),
               Component = "Means") %>% 
        rename(Mean=X.mean,
               SE=X.se,
               P_Value=X.p)
      x.mean <- x.mean %>% 
        dplyr::select(Taxa,Exposure,Mean,SE,P_Value,Component)
      rownames(x.mean) <- NULL
      mean.results <- rbind(psi.mean,x.mean)
      
      #Probability component
      #Mixture estimate
      mean.est <- mod1$coef[[2]][[2]]
      mean.pval <- mod1$pval[[2]][[2]]
      
      
      mean.se <- summary(mod1)$coeffients$zero[2,2]
      
      #Individual estimates
      X_values <- data.frame(
        X.mean = summary(mod1$fit)[1]$coefficients$zero[2:(length(exposure)+1)],
        X.se = summary(mod1$fit)[1]$coefficients$zero[,2][2:(length(exposure)+1)],
        X.p = summary(mod1$fit)[1]$coefficients$zero[,4][2:(length(exposure)+1)]
      )
      
      psi.prob <- data.frame(Taxa = r, Exposure = "Mixture", Mean = mean.est, SE = mean.se, P_Value = mean.pval,Component="Probability")
      x.prob <- X_values
      x.prob <- x.prob %>% 
        mutate(Taxa = r,
               Exposure = rownames(x.prob),
               Component = "Probability") %>% 
        rename(Mean=X.mean,
               SE=X.se,
               P_Value=X.p)
      x.prob <- x.prob %>% 
        dplyr::select(Taxa,Exposure,Mean,SE,P_Value,Component)
      rownames(x.prob) <- NULL
      prob.results <- rbind(psi.prob,x.prob)
      
      results <- rbind(mean.results,prob.results)
      ZING.results <- rbind(ZING.results,results)
      ZING.results$Model <- "ZINB"
    }, error = function(e) {
      # Check for the specific error condition
      if (grepl("glm.fit: algorithm did not converge", conditionMessage(e)) | grepl("glm.fit: fitted probabilities numerically 0 or 1 occurred", conditionMessage(e))) {
        psi.mean <- data.frame(Taxa = r, Exposure = "Mixture", Mean = NA, SE = NA, P_Value = NA,Component = "Means")
        psi.prob <- data.frame(Taxa = r, Exposure = "Mixture", Mean = NA, SE = NA, P_Value = NA,Component = "Probability")
        mean.x <- data.frame(Taxa= rep(r,times=length(exposure)),Exposure=exposure,Mean=rep(NA,times=length(exposure)),SE=rep(NA,times=length(exposure)),P_Value=rep(NA,times=length(exposure)),Component="Means")
        prob.x <- data.frame(Taxa= rep(r,times=length(exposure)),Exposure=exposure,Mean=rep(NA,times=length(exposure)),SE=rep(NA,times=length(exposure)),P_Value=rep(NA,times=length(exposure)),Component="Probability")
        results <- rbind(psi.mean,mean.x)
        results <- rbind(results,psi.prob)
        results <- rbind(results,prob.x)
        results$Model <- "ZINB"
        ZING.results <- rbind(ZING.results,results)
        ZING.results$Model <- "ZINB"
      } else {
        # For other errors, re-throw the error
        stop(e)
      }
    })
  } 
  
  if (min(dat[[r]])!=0) {
    #Perform poisson instead
    m0 <- as.formula(paste0(r, "~", paste(exposure, collapse = " + ")))
    mod1 <- qgcomp(f=m0, expnms = c(exposure),
                   data=dat, q=NULL, family=poisson())
   
   
    mean.est <- summary(mod1)$coefficients[2,1]
    mean.pval <- summary(mod1)$coefficients[2,6]
    mean.se <- summary(mod1)$coefficients[2,2]
    
    #Individual estimates
    X_values <- data.frame(
      X.mean = c(mod1$fit$coefficients[2:(length(exposure)+1)]))
    
    X.se <- c()
    for (i in 2:(length(exposure)+1)) {
      se <- summary(mod1$fit)$coefficients[i,2]
      X.se <- rbind(X.se,se)
    }
    rownames(X.se) <- NULL
    X.se <- as.vector(X.se)
    X.p <- c()
    for (i in 2:(length(exposure)+1)) {
      p <- summary(mod1$fit)$coefficients[i,4]
      X.p <- rbind(X.p,p)
    }
    rownames(X.p) <- NULL
    X.p <- as.vector(X.p)
    
    X_values$X.se <- X.se
    X_values$X.p <- X.p
    X_values$Taxa <- r
    X_values$Component <- "Means"
    X_values$Exposure <- exposure
    # colnames(X_values)
    X_values <- X_values %>% 
      rename(Mean=X.mean,
             SE=X.se,
             P_Value=X.p) %>% 
      dplyr::select(Taxa,Exposure,Mean,SE,P_Value,Component)
    rownames(X_values) <- NULL
    psi.mean <- data.frame(Taxa = r, Exposure = "Mixture", Mean = mean.est, SE = mean.se, P_Value = mean.pval,Component = "Means")
    poisson.results <- rbind(psi.mean,X_values)
    poisson.results$Model <- "Poisson"
  }
  model.results <- rbind(ZING.results,poisson.results)
}


# 3. Bayesian Hierarchical Zero-Inflated Negative Binomial Regression Model with G-Computation ----
BHRM.microbiome <-
  "model {
  for(r in 1:P.s) {
    for(i in 1:N) {
      Y[i,r] ~ dnegbin(mu[i,r], disp[r])
      mu[i,r] <- disp[r]/(disp[r]+(1-zero[i,r])*lambda[i,r]) - 0.000001*zero[i,r]
      log(lambda[i,r]) <- alpha[r] + inprod(species.beta[r,1:P.e], X[i,1:P.e]) 

      # zero-inflation
      zero[i,r] ~ dbern(pi[i,r])
      logit(pi[i,r]) <- alpha.zero[r] + inprod(species.beta.zero[r,1:P.e], X[i,1:P.e]) 
    }
    # prior on dispersion parameter
    disp[r] ~ dunif(0,50)

    # prior on intercept
    alpha[r] ~ dnorm(0, 1.0E-02)
    alpha.zero[r] ~ dnorm(0, 1.0E-02)

    # prior on proportion of non-zeros
    omega[r] ~ dunif(0,1)

    # prior on exposure effects
    for(p in 1:P.e) {
      # species.beta.zero[r,p] ~ dnorm(0, 1.0E-02)
      species.beta[r,p] ~ dnorm(mu.species[r,p], tau[r])
      mu.species[r,p] <- inprod(genus.beta[1:P.g,p], Z.s.g[r,1:P.g]) 
      #Zero inflation component
      species.beta.zero[r,p] ~ dnorm(mu.species.zero[r,p], tau[r])
      mu.species.zero[r,p] <- inprod(genus.beta.zero[1:P.g,p], Z.s.g[r,1:P.g])
    }

    # prior on precision
    tau[r] <- 1/(sigma[r]*sigma[r])
    sigma[r] ~ dunif(0,3)

    # g-estimation
    species.eta.low[r] <- inprod(species.beta[r,1:P.e], profiles[1,1:P.e])
    species.eta.high[r] <- inprod(species.beta[r,1:P.e], profiles[2,1:P.e])
    species.psi[r] <- species.eta.high[r]-species.eta.low[r]
    # zero-inflation
    species.eta.low.zero[r] <- inprod(species.beta.zero[r,1:P.e], profiles[1,1:P.e])
    species.eta.high.zero[r] <- inprod(species.beta.zero[r,1:P.e], profiles[2,1:P.e])
    species.psi.zero[r] <- species.eta.high.zero[r]-species.eta.low.zero[r]
  }

  # Genus level
  for(g.r in 1:P.g) {
    for(p in 1:P.e) {
      genus.beta[g.r,p] ~ dnorm(mu.family[g.r,p],genus.tau[g.r]) 
      mu.family[g.r,p] <- inprod(family.beta[1:P.f,p], Z.g.f[g.r,1:P.f])
      #Zero inflation component
      genus.beta.zero[g.r,p] ~ dnorm(mu.family.zero[g.r,p],genus.tau[g.r]) 
      mu.family.zero[g.r,p] <- inprod(family.beta.zero[1:P.f,p], Z.g.f[g.r,1:P.f])
    }
    # prior on precision
    genus.tau[g.r] <- 1/(genus.sigma[g.r]*genus.sigma[g.r])
    genus.sigma[g.r] ~ dunif(0,3)

    # g-estimation
    genus.eta.low[g.r] <- inprod(genus.beta[g.r,1:P.e], profiles[1,1:P.e])
    genus.eta.high[g.r] <- inprod(genus.beta[g.r,1:P.e], profiles[2,1:P.e])
    genus.psi[g.r] <- genus.eta.high[g.r]-genus.eta.low[g.r]
    #zero inflation
    genus.eta.low.zero[g.r] <- inprod(genus.beta.zero[g.r,1:P.e], profiles[1,1:P.e])
    genus.eta.high.zero[g.r] <- inprod(genus.beta.zero[g.r,1:P.e], profiles[2,1:P.e])
    genus.psi.zero[g.r] <- genus.eta.high.zero[g.r]-genus.eta.low.zero[g.r]
  }

  # Family level
  for(f.r in 1:P.f) {
    for(p in 1:P.e) {
      family.beta[f.r,p] ~ dnorm(mu.order[f.r,p], family.tau[f.r])
      mu.order[f.r,p] <- inprod(order.beta[1:P.o,p], Z.f.o[f.r,1:P.o])
      #Zero inflation component
      family.beta.zero[f.r,p] ~ dnorm(mu.order.zero[f.r,p], family.tau[f.r])
      mu.order.zero[f.r,p] <- inprod(order.beta.zero[1:P.o,p], Z.f.o[f.r,1:P.o])

    }
    # prior on precision
    family.tau[f.r] <- 1/(family.sigma[f.r]*family.sigma[f.r])
    family.sigma[f.r] ~ dunif(0,3)

    # g-estimation
    family.eta.low[f.r] <- inprod(family.beta[f.r,1:P.e], profiles[1,1:P.e])
    family.eta.high[f.r] <- inprod(family.beta[f.r,1:P.e], profiles[2,1:P.e])
    family.psi[f.r] <- family.eta.high[f.r]-family.eta.low[f.r]
    #zero inflation
    family.eta.low.zero[f.r] <- inprod(family.beta.zero[f.r,1:P.e], profiles[1,1:P.e])
    family.eta.high.zero[f.r] <- inprod(family.beta.zero[f.r,1:P.e], profiles[2,1:P.e])
    family.psi.zero[f.r] <- family.eta.high.zero[f.r]-family.eta.low.zero[f.r]
  }

  # Order level
  for(o.r in 1:P.o) {
    for(p in 1:P.e) {
      order.beta[o.r,p] ~ dnorm(mu.class[o.r,p], order.tau[o.r])
      mu.class[o.r,p] <- inprod(class.beta[1:P.c,p], Z.o.c[o.r,1:P.c])
      #Zero inflation component
      order.beta.zero[o.r,p] ~ dnorm(mu.class.zero[o.r,p], order.tau[o.r])
      mu.class.zero[o.r,p] <- inprod(class.beta.zero[1:P.c,p], Z.o.c[o.r,1:P.c])
    }
    # prior on precision
    order.tau[o.r] <- 1/(order.sigma[o.r]*order.sigma[o.r])
    order.sigma[o.r] ~ dunif(0,3)

    # g-estimation
    order.eta.low[o.r] <- inprod(order.beta[o.r,1:P.e], profiles[1,1:P.e])
    order.eta.high[o.r] <- inprod(order.beta[o.r,1:P.e], profiles[2,1:P.e])
    order.psi[o.r] <- order.eta.high[o.r]-order.eta.low[o.r]
    #zero infl
    order.eta.low.zero[o.r] <- inprod(order.beta.zero[o.r,1:P.e], profiles[1,1:P.e])
    order.eta.high.zero[o.r] <- inprod(order.beta.zero[o.r,1:P.e], profiles[2,1:P.e])
    order.psi.zero[o.r] <- order.eta.high.zero[o.r]-order.eta.low.zero[o.r]
  }

  # Class level
  for(c.r in 1:P.c) {
    for(p in 1:P.e) {
      class.beta[c.r,p] ~ dnorm(mu.phylum[c.r,p], class.tau[c.r])
      mu.phylum[c.r,p] <- inprod(phylum.beta[1:P.p,p], Z.c.p[c.r,1:P.p])
      #Zero inflation component
      class.beta.zero[c.r,p] ~ dnorm(mu.phylum.zero[c.r,p], class.tau[c.r])
      mu.phylum.zero[c.r,p] <- inprod(phylum.beta.zero[1:P.p,p], Z.c.p[c.r,1:P.p])
    }
    # prior on precision
    class.tau[c.r] <- 1/(class.sigma[c.r]*class.sigma[c.r])
    class.sigma[c.r] ~ dunif(0,3)

    # g-estimation
    class.eta.low[c.r] <- inprod(class.beta[c.r,1:P.e], profiles[1,1:P.e])
    class.eta.high[c.r] <- inprod(class.beta[c.r,1:P.e], profiles[2,1:P.e])
    class.psi[c.r] <- class.eta.high[c.r]-class.eta.low[c.r]

    #zero component
    class.eta.low.zero[c.r] <- inprod(class.beta.zero[c.r,1:P.e], profiles[1,1:P.e])
    class.eta.high.zero[c.r] <- inprod(class.beta.zero[c.r,1:P.e], profiles[2,1:P.e])
    class.psi.zero[c.r] <- class.eta.high.zero[c.r]-class.eta.low.zero[c.r]
  }

  # Phylum level
  for(p.r in 1:P.p) {
    for(p in 1:P.e) {
      phylum.beta[p.r,p] ~ dnorm(0, phylum.tau[p.r])
      #Zero inflation component
      phylum.beta.zero[p.r,p] ~ dnorm(0, phylum.tau[p.r])
    }
    # prior on precision
    phylum.tau[p.r] <- 1/(phylum.sigma[p.r]*phylum.sigma[p.r])
    phylum.sigma[p.r] ~ dunif(0,3)

    # g-estimation
    phylum.eta.low[p.r] <- inprod(phylum.beta[p.r,1:P.e], profiles[1,1:P.e])
    phylum.eta.high[p.r] <- inprod(phylum.beta[p.r,1:P.e], profiles[2,1:P.e])
    phylum.psi[p.r] <- phylum.eta.high[p.r]-phylum.eta.low[p.r]

    #Zero inflation
    phylum.eta.low.zero[p.r] <- inprod(phylum.beta.zero[p.r,1:P.e], profiles[1,1:P.e])
    phylum.eta.high.zero[p.r] <- inprod(phylum.beta.zero[p.r,1:P.e], profiles[2,1:P.e])
    phylum.psi.zero[p.r] <- phylum.eta.high.zero[p.r]-phylum.eta.low.zero[p.r]
  }

}"

#4. Bayesian Ridge Regression with G-Computation -----
RBHRM.microbiome <-
  "model {
  for(r in 1:P.s) {
    for(i in 1:N) {
      Y[i,r] ~ dnegbin(mu[i,r], disp[r])
      mu[i,r] <- disp[r]/(disp[r]+(1-zero[i,r])*lambda[i,r]) - 0.000001*zero[i,r]
      log(lambda[i,r]) <- alpha[r] + inprod(species.beta[r,1:P.e], X[i,1:P.e]) 

      # zero-inflation
      zero[i,r] ~ dbern(pi[i,r])
      logit(pi[i,r]) <- alpha.zero[r] + inprod(species.beta.zero[r,1:P.e], X[i,1:P.e]) 
    }
    # prior on dispersion parameter
    disp[r] ~ dunif(0,50)

    # prior on intercept
    alpha[r] ~ dnorm(0, 1.0E-02)
    alpha.zero[r] ~ dnorm(0, 1.0E-02)

    # prior on proportion of non-zeros
    omega[r] ~ dunif(0,1)

    # prior on exposure effects
    for(p in 1:P.e) {
      # species.beta.zero[r,p] ~ dnorm(0, 1.0E-02)
      species.beta[r,p] ~ dnorm(mu.species[r,p], tau[r])
      # species.beta[r,p] ~ dnorm(0, tau[r])
      # mu.species[r,p] <- inprod(genus.beta[1:P.g,p], Z.s.g[r,1:P.g]) #or here
      mu.species[r,p] <- 0
      #Zero inflation component
      species.beta.zero[r,p] ~ dnorm(mu.species.zero[r,p], tau[r])
      # species.beta.zero[r,p] ~ dnorm(0, tau[r])
      # mu.species.zero[r,p] <- inprod(genus.beta.zero[1:P.g,p], Z.s.g[r,1:P.g])
      mu.species.zero[r,p] <- 0
    }

    # prior on precision
    tau[r] <- 1/(sigma[r]*sigma[r])
    sigma[r] ~ dunif(0,3)

    # g-estimation
    species.eta.low[r] <- inprod(species.beta[r,1:P.e], profiles[1,1:P.e])
    species.eta.high[r] <- inprod(species.beta[r,1:P.e], profiles[2,1:P.e])
    species.psi[r] <- species.eta.high[r]-species.eta.low[r]
    # zero-inflation
    species.eta.low.zero[r] <- inprod(species.beta.zero[r,1:P.e], profiles[1,1:P.e])
    species.eta.high.zero[r] <- inprod(species.beta.zero[r,1:P.e], profiles[2,1:P.e])
    species.psi.zero[r] <- species.eta.high.zero[r]-species.eta.low.zero[r]
  }

  # Genus level
  for(g.r in 1:P.g) {
    for(p in 1:P.e) {
      genus.beta[g.r,p] ~ dnorm(mu.family[g.r,p],genus.tau[g.r])
      # genus.beta[g.r,p] ~ dnorm(0,genus.tau[g.r]) 
      # mu.family[g.r,p] <- inprod(family.beta[1:P.f,p], Z.g.f[g.r,1:P.f])
      mu.family[g.r,p] <- 0
      #Zero inflation component
      genus.beta.zero[g.r,p] ~ dnorm(mu.family.zero[g.r,p],genus.tau[g.r]) 
      # genus.beta.zero[g.r,p] ~ dnorm(0,genus.tau[g.r]) 
      # mu.family.zero[g.r,p] <- inprod(family.beta.zero[1:P.f,p], Z.g.f[g.r,1:P.f])
      mu.family.zero[g.r,p] <- 0
    }
    # prior on precision
    genus.tau[g.r] <- 1/(genus.sigma[g.r]*genus.sigma[g.r])
    genus.sigma[g.r] ~ dunif(0,3)

    # g-estimation
    genus.eta.low[g.r] <- inprod(genus.beta[g.r,1:P.e], profiles[1,1:P.e])
    genus.eta.high[g.r] <- inprod(genus.beta[g.r,1:P.e], profiles[2,1:P.e])
    genus.psi[g.r] <- genus.eta.high[g.r]-genus.eta.low[g.r]
    #zero inflation
    genus.eta.low.zero[g.r] <- inprod(genus.beta.zero[g.r,1:P.e], profiles[1,1:P.e])
    genus.eta.high.zero[g.r] <- inprod(genus.beta.zero[g.r,1:P.e], profiles[2,1:P.e])
    genus.psi.zero[g.r] <- genus.eta.high.zero[g.r]-genus.eta.low.zero[g.r]
  }

  # Family level
  for(f.r in 1:P.f) {
    for(p in 1:P.e) {
      family.beta[f.r,p] ~ dnorm(mu.order[f.r,p], family.tau[f.r])
      # family.beta[f.r,p] ~ dnorm(0, family.tau[f.r])
      # mu.order[f.r,p] <- inprod(order.beta[1:P.o,p], Z.f.o[f.r,1:P.o])
      mu.order[f.r,p] <- 0
      #Zero inflation component
      family.beta.zero[f.r,p] ~ dnorm(mu.order.zero[f.r,p], family.tau[f.r])
      # family.beta.zero[f.r,p] ~ dnorm(0, family.tau[f.r])
      # mu.order.zero[f.r,p] <- inprod(order.beta.zero[1:P.o,p], Z.f.o[f.r,1:P.o])
      mu.order.zero[f.r,p] <- 0

    }
    # prior on precision
    family.tau[f.r] <- 1/(family.sigma[f.r]*family.sigma[f.r])
    family.sigma[f.r] ~ dunif(0,3)

    # g-estimation
    family.eta.low[f.r] <- inprod(family.beta[f.r,1:P.e], profiles[1,1:P.e])
    family.eta.high[f.r] <- inprod(family.beta[f.r,1:P.e], profiles[2,1:P.e])
    family.psi[f.r] <- family.eta.high[f.r]-family.eta.low[f.r]
    #zero inflation
    family.eta.low.zero[f.r] <- inprod(family.beta.zero[f.r,1:P.e], profiles[1,1:P.e])
    family.eta.high.zero[f.r] <- inprod(family.beta.zero[f.r,1:P.e], profiles[2,1:P.e])
    family.psi.zero[f.r] <- family.eta.high.zero[f.r]-family.eta.low.zero[f.r]
  }

  # Order level
  for(o.r in 1:P.o) {
    for(p in 1:P.e) {
      order.beta[o.r,p] ~ dnorm(mu.class[o.r,p], order.tau[o.r])
      # order.beta[o.r,p] ~ dnorm(0, order.tau[o.r])
      # mu.class[o.r,p] <- inprod(class.beta[1:P.c,p], Z.o.c[o.r,1:P.c])
      mu.class[o.r,p] <- 0
      #Zero inflation component
      order.beta.zero[o.r,p] ~ dnorm(mu.class.zero[o.r,p], order.tau[o.r])
      # order.beta.zero[o.r,p] ~ dnorm(0, order.tau[o.r])
      # mu.class.zero[o.r,p] <- inprod(class.beta.zero[1:P.c,p], Z.o.c[o.r,1:P.c])
      mu.class.zero[o.r,p] <- 0
    }
    # prior on precision
    order.tau[o.r] <- 1/(order.sigma[o.r]*order.sigma[o.r])
    order.sigma[o.r] ~ dunif(0,3)

    # g-estimation
    order.eta.low[o.r] <- inprod(order.beta[o.r,1:P.e], profiles[1,1:P.e])
    order.eta.high[o.r] <- inprod(order.beta[o.r,1:P.e], profiles[2,1:P.e])
    order.psi[o.r] <- order.eta.high[o.r]-order.eta.low[o.r]
    #zero infl
    order.eta.low.zero[o.r] <- inprod(order.beta.zero[o.r,1:P.e], profiles[1,1:P.e])
    order.eta.high.zero[o.r] <- inprod(order.beta.zero[o.r,1:P.e], profiles[2,1:P.e])
    order.psi.zero[o.r] <- order.eta.high.zero[o.r]-order.eta.low.zero[o.r]
  }

  # Class level
  for(c.r in 1:P.c) {
    for(p in 1:P.e) {
      class.beta[c.r,p] ~ dnorm(mu.phylum[c.r,p], class.tau[c.r])
      # class.beta[c.r,p] ~ dnorm(0, class.tau[c.r])
      # mu.phylum[c.r,p] <- inprod(phylum.beta[1:P.p,p], Z.c.p[c.r,1:P.p])
      mu.phylum[c.r,p] <- 0
      #Zero inflation component
      class.beta.zero[c.r,p] ~ dnorm(mu.phylum.zero[c.r,p], class.tau[c.r])
      # class.beta.zero[c.r,p] ~ dnorm(0, class.tau[c.r])
      # mu.phylum.zero[c.r,p] <- inprod(phylum.beta.zero[1:P.p,p], Z.c.p[c.r,1:P.p])
      mu.phylum.zero[c.r,p] <- 0
    }
    # prior on precision
    class.tau[c.r] <- 1/(class.sigma[c.r]*class.sigma[c.r])
    class.sigma[c.r] ~ dunif(0,3)

    # g-estimation
    class.eta.low[c.r] <- inprod(class.beta[c.r,1:P.e], profiles[1,1:P.e])
    class.eta.high[c.r] <- inprod(class.beta[c.r,1:P.e], profiles[2,1:P.e])
    class.psi[c.r] <- class.eta.high[c.r]-class.eta.low[c.r]

    #zero component
    class.eta.low.zero[c.r] <- inprod(class.beta.zero[c.r,1:P.e], profiles[1,1:P.e])
    class.eta.high.zero[c.r] <- inprod(class.beta.zero[c.r,1:P.e], profiles[2,1:P.e])
    class.psi.zero[c.r] <- class.eta.high.zero[c.r]-class.eta.low.zero[c.r]
  }

  # Phylum level
  for(p.r in 1:P.p) {
    for(p in 1:P.e) {
      phylum.beta[p.r,p] ~ dnorm(0, phylum.tau[p.r])
      #Zero inflation component
      phylum.beta.zero[p.r,p] ~ dnorm(0, phylum.tau[p.r])
    }
    # prior on precision
    phylum.tau[p.r] <- 1/(phylum.sigma[p.r]*phylum.sigma[p.r])
    phylum.sigma[p.r] ~ dunif(0,3)

    # g-estimation
    phylum.eta.low[p.r] <- inprod(phylum.beta[p.r,1:P.e], profiles[1,1:P.e])
    phylum.eta.high[p.r] <- inprod(phylum.beta[p.r,1:P.e], profiles[2,1:P.e])
    phylum.psi[p.r] <- phylum.eta.high[p.r]-phylum.eta.low[p.r]

    #Zero inflation
    phylum.eta.low.zero[p.r] <- inprod(phylum.beta.zero[p.r,1:P.e], profiles[1,1:P.e])
    phylum.eta.high.zero[p.r] <- inprod(phylum.beta.zero[p.r,1:P.e], profiles[2,1:P.e])
    phylum.psi.zero[p.r] <- phylum.eta.high.zero[p.r]-phylum.eta.low.zero[p.r]
  }

}"
