kGDF <- function(k=1, y, fm, data, distr = "normal", mod = "GLM", nfit=1, noise = 0.25, perf=1, ndata=NULL,  pass, method="jackknife"){
  
  #' 
  #' @y observed output data
  #' @fm formula
  #' @data input data 
  #' @distr ditribution of the added noise
  #' @model model type
  #' @noise sd of the noise (whe it is a nomrmal distribution)
  #' @perf n.trees to predict BRT
  #' @pass number of loops to re run the model
  #' @method jacknife or bootstrap
  #' @nfit number of times that each model is fitted with the same data
  
  if (!is.null(ndata)) {data = data[1:ndata,]; y=y[1:ndata]}
  
  # sample size n 
  n <- length(y)
  
  #set noise
  if(distr == "normal" ) noise <- noise*sd(y)
  
  # dataframes for ye and yehat
  ye <- yehat <- data.frame(rep(NA,n))
  
  if(method == "jackknife"){ # makes sure that each data point is perturbed per pass
    
    # loop to re-run Model until each data point has been perturbed
    i <- 1 # subscription  
    for(t in 1:pass){
      # check for NAs
      check.NA <- rep(NA, n)
      
      while (sum(is.na(check.NA)) > 0){  
        y2 <- y # overwrite y2 each loop with original y
        
        to.fill <- which(is.na(check.NA)) # check which data point still needs to be perturbed
        
        # randomise subset of size k to perturb during this loop
        if ((length(to.fill) == 1) & (k == 1)){
          # if there is only one NA to be disturbed and sample size = k = 1
          y2[to.fill] <- y[to.fill]+rnorm(k, 0, noise) # perturb last remaining point
          check.NA[to.fill] <- 1 # mark perturbed subset
        }else{
          if ((length(to.fill) >= k ) & (length(to.fill) > 1)){
            to.use <- sample(to.fill, size= k, replace= F) # randomly select k data points
          }
          if ((length(to.fill) < k)){
            # if there are fewer remaining NAs than data points to be disturbed
            fill.up <- sample(1:n, size=k-length(to.fill), replace=F) # fill up sample of size k
            to.use <- c(to.fill, fill.up)
          }
          if (distr == "normal") y2[to.use] <- y2[to.use]+rnorm(k, 0, noise) #
          if(distr == "binary") y2[to.use] <- !y2[to.use]
          check.NA[to.use] <- 1 # mark perturbed subset 
        }
        
        assign("y2", y2, envir = .GlobalEnv) # provide "y2" in the global environment 
        # store intermediates
        ye[,i] <- y2
        
        yehat.temp <- matrix(NA, nrow=n, ncol=nfit)
        
        for (t in 1:nfit){ # fit the model many times with the same data and average allong resuls
          # accommodate model with perturbed data
          if(mod == "GLM" | mod == "GAM" | mod == "ANN") fm.updated <- update(fm, y2 ~., data=data)
          if(mod == "rF"){
            if(fm$type == "regression") fm.updated <- update(fm, y2 ~., data=data)
            if(fm$type == "classification") fm.updated <- update(fm, as.factor(y2) ~., data=data)
          }
          if(mod == "BRT") fm.updated <- update(fm, as.formula(y2~.), data=data)
          if(mod == "ranger"){ 
            ran <- data.frame(cbind(y2, data))
            fm.updated <- ranger( y2~  .,  data=ran)
          }
          
          if(mod == "rF"){
            if(fm$type == "regression") yehat.temp[,t] <- predict(fm.updated, newdata=data, type="response")
            #if(fm$type == "regression") yehat[,i] <- predict(fm.updated, type="response")
            if(fm$type == "classification") yehat.temp[,t] <- predict(fm.updated, newdata=data, type="prob")[, 2]           
          }
          if(mod == "ANN")  yehat.temp[,t] <- predict(fm.updated, type="raw")
          if(mod == "GLM" | mod == "GAM") yehat.temp[,t] <- predict(fm.updated, type="response")
          if(mod == "BRT") yehat.temp[,t] <- predict(fm.updated, type="response", n.trees=perf)
          if(mod == "ranger")  yehat.temp[,t] <- predict(fm.updated, data=data)$predictions
          
        }
        
        yehat[,i] <- rowMeans(yehat.temp)
        i <- i+1 # move one column further in the data frames
        # Remove y2 from global environment
        rm(y2, envir= .GlobalEnv)
      }
    }
    cat("i =",i-1)
  }
  
  if(method=="bootstrap"){ # meaning that the k perturbed data points are selected in every pass from the complete dataset, possibly leaving out a few
    for(i in 1:pass){
      # perturbation
      y2 <- y # perturb copy of y
      to.perturb <- sample(n,size=k,replace=FALSE)
      if (distr == "normal") y2[to.perturb] <- y2[to.perturb]+rnorm(k, 0, noise)
      if(distr == "binary") y2[to.perturb] <- !y2[to.perturb]
      assign("y2", y2, envir = .GlobalEnv) # for model update needed, apparently
      
      yehat.temp <- matrix(NA, nrow=n, ncol=nfit)
      for (t in 1:nfit){ # fit the model many times with the same data and average allong resuls
        
        # accommodate model with perturbed data
        if(mod == "GLM" | mod == "GAM" | mod == "ANN") fm.updated <- update(fm, y2 ~., data=data)
        #In this line update should maybe designed to fm.update
        if(mod == "rF" & distr == "normal") fm.updated <- update(fm, y2 ~., data=data)
        if(mod == "rF" & distr == "binary") fm.updated <- update(fm, as.factor(y2) ~., data=data)
        if(mod == "BRT") fm.updated <- update(fm, as.formula(y2~.), data = data)
        if(mod == "ranger"){ 
          ran <- data.frame(cbind(y2, data))
          fm.updated <- ranger( y2~  .,  data=data)
        }
        # store intermediates
        ye[,i] <- y2
        if(mod == "rF"){
          # in this line newdata should be detailed to avoid out of bag prediction
          if(fm$type == "regression") yehat.temp[,t] <- predict(fm.updated, newdata=data, type="response")
          if(fm$type == "classification") yehat.temp[,t] <- predict(fm.updated, newdata=data, type="prob")[, 2]           
        }
        if(mod == "ANN")  yehat.temp[,t] <- predict(fm.updated, type="raw")
        if(mod == "GLM" | mod == "GAM") yehat.temp[,t] <- predict(fm.updated, type="response")
        if(mod == "BRT") yehat.temp[,t] <- predict(fm.updated, type="response", n.trees=perf)
        if(mod == "ranger")  yehat.temp[,t] <- predict(fm.updated, data=data)$predictions
      }
      
      yehat[,i] <- rowMeans(yehat.temp)
      rm(y2, envir= .GlobalEnv) 
    }
  }
  
  # LR: Horizontal Method:
  ye <- as.matrix(ye); yehat <- as.matrix(yehat)
  #Fit a LR to yehat[i,] vs ye[,i]
  fm.LR <- lapply(seq_along(ye[,1]), function(x)lm(yehat[x,]~ye[x,]))
  #Calculate sum of slopes m[i] as GDF
  GDF.hor <- sum(sapply(fm.LR, function(i)coefficients(i)[2]), na.rm=T)
  
  return(GDF.hor)
}
