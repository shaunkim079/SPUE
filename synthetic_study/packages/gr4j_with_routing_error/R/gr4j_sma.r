

gr4j.sma<-function(param, initial_state, state_error, input,
                      etmult = 1, return_state = FALSE, 
                      transformed = FALSE, run_compiled = T){
  tryerror<-try({
    
    x1<-param[1]
    if(length(param)>1){
      param2<-param[2]
    } else {
      param2<-2
    }
    if (transformed) 
      x1 <- exp(x1)
    stopifnot(x1 >= 0)
    stopifnot(etmult >= 0)
    stopifnot(param2 > 0)
    
    if(is.loaded("sma_gr4j") & run_compiled){
      P <- as.double(input$P)
      E <- as.double(input$E * etmult)
      initial_state <- as.double(initial_state)
      x1 <- as.double(x1)
      state_error <- as.double(state_error)
      param2 <- as.double(param2)
      n <- as.integer(length(P))
      U <- as.double(rep(0,n))
      S <- as.double(rep(0,n))
      ET <- as.double(rep(0,n))
      
      out <- .C("sma_gr4j",
                P=P,
                E=E,
                n=n,
                x1=x1,
                param2=param2,
                initial_state=initial_state,
                state_error=state_error,
                U=U,
                S=S,
                ET=ET,
                success=as.integer(1))
      
      if(out$success!=1) stop("negative value in simulation when there shouldn't be")
      ans <- out$U
      if (return_state) {
        attributes(out$S) <- attributes(out$ET) <- attributes(out$U)
        ans <- cbind(U = out$U, S = out$S, ET = out$ET)
      }
      
    } else {
      P <- input$P
      E <- input$E * etmult
      bad <- is.na(P) | is.na(E)
      P[bad] <- 0
      E[bad] <- 0
      U <- S <- ET <- P
      S_prev <- initial_state
      for (t in seq(1, length(P))) {
        if(t>1){
          S_prev <- S_prev + state_error[t-1]
        }
        Pn <- max(P[t] - E[t], 0)
        En <- max(E[t] - P[t], 0)
        St_x1 <- S_prev/x1
        Ps <- 0
        ET[t] <- 0
        if (Pn > 0) {
          Ps <- (x1 * (1 - St_x1^param2) * tanh(Pn/x1)/(1 + 
                                                          St_x1 * tanh(Pn/x1)))
        }
        else {
          ET[t] <- (S_prev * (2 - St_x1) * tanh(En/x1)/(1 + 
                                                          (1 - St_x1) * tanh(En/x1)))
        }
        S[t] <- S_prev - ET[t] + Ps
        perc <- S[t] * (1 - (1 + ((4/9) * (S[t]/x1))^4)^(-0.25))
        S[t] <- S[t] - perc
        U[t] <- perc + (Pn - Ps)
        S_prev <- S[t]
        if(S_prev<0) stop("negative state value")
        if(U[t]<0) stop("negative discharge value")
        if(ET[t]<0) stop("negative ET value")
        if(Ps<0) stop("negative Ps value")
        if(perc<0) stop("negative perc value")
      }
      U[bad] <- NA
      ans <- U
      if (return_state) {
        attributes(S) <- attributes(ET) <- attributes(U)
        ans <- cbind(U = U, S = S, ET = ET)
      }
    }
  },silent=T)
  if(length(grep("Error",tryerror))>0) ans<-NA
  
  return(ans)
}



gr4j.sma.debug<-function(param, initial_state, state_error, input,
                   etmult = 1, return_state = FALSE, 
                   transformed = FALSE){
  
  layout(1:2)
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  
  x1<-param[1]
  if(length(param)>1){
    param2<-param[2]
  } else {
    param2<-2
  }
  
  if (transformed) 
    x1 <- exp(x1)
  stopifnot(x1 >= 0)
  stopifnot(etmult >= 0)
  stopifnot(param2 > 0)
  
  P <- input$P
  E <- input$E * etmult
  bad <- is.na(P) | is.na(E)
  P[bad] <- 0
  E[bad] <- 0
  All_Ps <- U <- S <- ET <- P
  S_prev <- initial_state
  for (t in seq(1, length(P))) {
    #     if(t>1){
    #      S_prev <- S_prev + state_error[t-1]
    #     }
    Pn <- max(P[t] - E[t], 0)
    En <- max(E[t] - P[t], 0)
    St_x1 <- S_prev/x1
    Ps <- 0
    ET[t] <- 0
    if (Pn > 0) {
      Ps <- (x1 * (1 - St_x1^param2) * tanh(Pn/x1)/(1 + 
                                                      St_x1 * tanh(Pn/x1)))
    }
    else {
      ET[t] <- (S_prev * (2 - St_x1) * tanh(En/x1)/(1 + 
                                                      (1 - St_x1) * tanh(En/x1)))
    }
    All_Ps[t] <- Ps
    S[t] <- S_prev - ET[t] + Ps
    perc <- S[t] * (1 - (1 + ((4/9) * (S[t]/x1))^4)^(-0.25))
    S[t] <- S[t] - perc
    U[t] <- perc + (Pn - Ps)
    S_prev <- S[t]
    if(S_prev<0) stop("negative state value")
    if(U[t]<0) stop("negative discharge value")
  }
  
  
  plot(S,type="l")
  var(S)
  sd(S)
  
  S1<-S
  perc1<-perc
  U1<-U
  ET1<-ET
  All_Ps1<-All_Ps
  

  x1<-param[1]
  if(length(param)>1){
    param2<-param[2]
  } else {
    param2<-2
  }
  
  if (transformed) 
    x1 <- exp(x1)
  stopifnot(x1 >= 0)
  stopifnot(etmult >= 0)
  stopifnot(param2 > 0)
  
  P <- input$P
  E <- input$E * etmult
  bad <- is.na(P) | is.na(E)
  P[bad] <- 0
  E[bad] <- 0
  All_Ps <- U <- S <- ET <- P
  S_prev <- initial_state
  for (t in seq(1, length(P))) {
    if(t>1){
     S_prev <- S_prev + state_error[t-1]
    }
    Pn <- max(P[t] - E[t], 0)
    En <- max(E[t] - P[t], 0)
    St_x1 <- S_prev/x1
    Ps <- 0
    ET[t] <- 0
    if (Pn > 0) {
      Ps <- (x1 * (1 - St_x1^param2) * tanh(Pn/x1)/(1 + 
                                                      St_x1 * tanh(Pn/x1)))
    }
    else {
      ET[t] <- (S_prev * (2 - St_x1) * tanh(En/x1)/(1 + 
                                                      (1 - St_x1) * tanh(En/x1)))
    }
    All_Ps[t] <- Ps
    S[t] <- S_prev - ET[t] + Ps
    perc <- S[t] * (1 - (1 + ((4/9) * (S[t]/x1))^4)^(-0.25))
    S[t] <- S[t] - perc
    U[t] <- perc + (Pn - Ps)
    S_prev <- S[t]
    if(S_prev<0) stop("negative state value")
    if(U[t]<0) stop("negative discharge value")
  }
  
  S2<-S
  perc2<-perc
  U2<-U
  ET2<-ET
  All_Ps2<-All_Ps
  
  
  plot(S,type="l")
  var(S)
  sd(S)
  
  
  plot(S1,type="l")
  lines(S2,col=2,lty=2)
  
  sd(S2-S1)
  sd(state_error)
  
  plot(U1,type="l")
  lines(U2,col=2,lty=2)
  
  plot(U1-U2,type="l")
  var(U1-U2)
  
  tt<-gr4j.sma(1000,700,actual_state_error,true_input,return_state = T)
  
  
  plot(tt[,2],type="l")
  lines(S2,col=2,lty=2)
  
  sd(tt[,2]-S2)
  
  var(tt[,1]-actual_discharge_error-U2)
  sd(tt[,1]-actual_discharge_error-U2)
  
  
  var_zero_mean(tt[,1]-actual_discharge_error-U2)
  
  plot(tt[,1]-actual_discharge_error,type="l")
  lines(U2,col=2,lty=2)
  
  plot(U2,type="l")
  plot(S2,type="l")
  
  
  
  
  
  
  
  
  
  
  
  browser()
  

  
  
  
  U[bad] <- NA
  ans <- U
  if (return_state) {
    attributes(S) <- attributes(ET) <- attributes(U)
    ans <- cbind(U = U, S = S, ET = ET)
  }
  return(ans)
}






gr4j.sma2<-function(param, initial_state, state_with_error, input,
                   etmult = 1, return_state = FALSE, 
                   transformed = FALSE){
  tryerror<-try({
    state_error<-rep(NA,length(state_with_error))
    x1<-param[1]
    if(length(param)>1){
      param2<-param[2]
    } else {
      param2<-2
    }
    
    if (transformed) 
      x1 <- exp(x1)
    stopifnot(x1 >= 0)
    stopifnot(etmult >= 0)
    stopifnot(param2 > 0)
    
    P <- input$P
    E <- input$E * etmult
    bad <- is.na(P) | is.na(E)
    P[bad] <- 0
    E[bad] <- 0
    U <- S <- ET <- P
    S_prev <- initial_state
    for (t in seq(1, length(P))) {
      if(t>1){
        state_error[t-1] <- state_with_error[t-1] - S_prev
        S_prev <- state_with_error[t-1]
      }
      Pn <- max(P[t] - E[t], 0)
      En <- max(E[t] - P[t], 0)
      St_x1 <- S_prev/x1
      Ps <- 0
      ET[t] <- 0
      if (Pn > 0) {
        Ps <- (x1 * (1 - St_x1^param2) * tanh(Pn/x1)/(1 + 
                                                        St_x1 * tanh(Pn/x1)))
      }
      else {
        ET[t] <- (S_prev * (2 - St_x1) * tanh(En/x1)/(1 + 
                                                        (1 - St_x1) * tanh(En/x1)))
      }
      S[t] <- S_prev - ET[t] + Ps
      perc <- S[t] * (1 - (1 + ((4/9) * (S[t]/x1))^4)^(-0.25))
      S[t] <- S[t] - perc
      U[t] <- perc + (Pn - Ps)
      S_prev <- S[t]
      if(S_prev<0) stop("negative state value")
      if(U[t]<0) stop("negative discharge value")
    }
    U[bad] <- NA
    ans <- U
    if (return_state) {
      attributes(S) <- attributes(ET) <- attributes(U)
      ans <- cbind(U = U, S = S, ET = ET, state_error = c(NA,state_error))
    }
  },silent=T)
  if(length(grep("Error",tryerror))>0) ans<-NA
  
  return(ans)
}





