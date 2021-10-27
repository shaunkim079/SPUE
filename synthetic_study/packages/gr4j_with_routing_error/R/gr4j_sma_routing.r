

gr4j.sma<-function(param, initial_state, state_error, input,
                      etmult = 1, return_state = FALSE, 
                      transformed = FALSE, run_compiled = T,
                      state_S_ts=NA){
  tryerror<-try({
    x1<-param[1]
    if(length(param)>1){
      param2<-param[2]
    } else {
      param2<-2
    }
    if (transformed){
      x1 <- exp(x1)
      initial_state<-exp(initial_state)
    }
    stopifnot(x1 >= 0)
    stopifnot(etmult >= 0)
    stopifnot(param2 > 0)
    stopifnot(length(state_error) == (length(input$P)-1))
    
    if(is.loaded("sma_gr4j_sk") & run_compiled){

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
      
      if(is.na(state_S_ts[1])){
        state_S_ts<-as.double(rep(-1e6,n-1))
        use_state_S_ts<-as.integer(0)
      } else {
        state_S_ts<-as.double(state_S_ts)
        use_state_S_ts<-as.integer(1)
      }
      
      out <- .C("sma_gr4j_sk",
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
                state_S_ts=state_S_ts,
                use_state_S_ts=use_state_S_ts,
                success=as.integer(1))
      
      if(out$success!=1) stop("negative value in SMA simulation when there shouldn't be")
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
          if(!is.na(state_S_ts[1])){
            if(state_S_ts[t-1]!=-1e6){
              S_prev<-state_S_ts[t-1]
            }
          }
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
          # if(Ps<0) browser()
        }
        else {
          ET[t] <- (S_prev * (2 - St_x1) * tanh(En/x1)/(1 + 
                                                          (1 - St_x1) * tanh(En/x1)))
        }
        S[t] <- S_prev - ET[t] + Ps
        perc <- S[t] * (1 - (1 + ((4/9) * (S[t]/x1))^4)^(-0.25))
        # perc<-0
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

# gr4jrouting.sim<-function (U, x2, x3, x4, R_0 = 0, split = 0.9, return_components = FALSE, 
#           epsilon = hydromad.getOption("sim.epsilon"), transformed = FALSE) 
# {
#   if (transformed) {
#     x2 <- sinh(x2)
#     x3 <- exp(x3)
#     x4 <- exp(x4) + 0.5
#   }
#   stopifnot(is.numeric(x2))
#   stopifnot(x3 >= 0)
#   stopifnot(x4 >= 0.5)
#   stopifnot(R_0 >= 0)
#   stopifnot(R_0 <= 1)
#   inAttr <- attributes(U)
#   U <- as.ts(U)
#   n <- ceiling(x4)
#   m <- ceiling(x4 * 2)
#   n2 <- floor(x4)
#   SH1 <- pmin((1:n/x4)^(5/2), 1)
#   SH2 <- pmin(c(0 + 0.5 * (1:n2/x4)^(5/2), 1 - 0.5 * (2 - n:m/x4)^(5/2)), 
#               1)
#   SH2[1:m/x4 > 2] <- 1
#   UH1 <- diff(c(0, SH1))
#   UH2 <- diff(c(0, SH2))
#   bad <- is.na(U)
#   U[bad] <- 0
#   filter.pad0 <- function(x, f) {
#     y <- x
#     y[] <- filter(c(rep(0, length(f)), x), filter = f, sides = 1)[-(1:length(f))]
#     y
#   }
#   Q9 <- filter.pad0(split * U, UH1)
#   Q1 <- filter.pad0((1 - split) * U, UH2)
#   COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
#   if (COMPILED) {
#     ans <- .C(routing_gr4j, as.double(Q9), as.double(Q1), 
#               as.integer(length(U)), as.double(x2), as.double(x3), 
#               as.double(R_0), Qr = double(length(U)), Qd = double(length(U)), 
#               R = double(length(U)), NAOK = FALSE, DUP = FALSE, 
#               PACKAGE = "hydromad")
#     Qr <- ans$Qr
#     Qd <- ans$Qd
#     R <- ans$R
#   }
#   else {
#     Qd <- Qr <- R <- U
#     R_prev <- R_0 * x3
#     for (t in seq(1, length(U))) {
#       Rt_x3 <- R_prev/x3
#       F <- x2 * Rt_x3^(7/2)
#       R[t] <- max(0, R_prev + Q9[t] + F)
#       Qr[t] <- R[t] * (1 - (1 + (R[t]/x3)^4)^(-0.25))
#       R[t] <- R[t] - Qr[t]
#       Qd[t] <- max(0, Q1[t] + F)
#       R_prev <- R[t]
#     }
#   }
#   Qr[Qr < epsilon] <- 0
#   Qd[Qd < epsilon] <- 0
#   attributes(Qr) <- attributes(Qd) <- attributes(R) <- inAttr
#   Qr[bad] <- NA
#   Qd[bad] <- NA
#   if (return_components) {
#     return(cbind(Xr = Qr, Xd = Qd, R = R))
#   }
#   else {
#     return(Qr + Qd)
#   }
# }


gr4jrouting.sim.sk<-function (U, x2, x3, x4, initial_state_R = 0, split = 0.9, return_components = FALSE, 
                              epsilon = 1e-05, transformed = FALSE, 
                              initial_condition_UH1=rep(0,ceiling(x4)-1),
                              initial_condition_UH2=rep(0,floor(x4)+length(ceiling(x4):ceiling(x4*2))-1),
                              run_compiled=T,
                              state_error_R,
                              state_R_ts=NA){
  # browser()
  tryerror<-try({
    if (transformed) {
      x2 <- sinh(x2)
      x3 <- exp(x3)
      initial_state_R<-exp(initial_state_R)
      x4 <- exp(x4) + 0.5
      
      initial_condition_UH1=rep(0,ceiling(x4)-1)
      initial_condition_UH2=rep(0,floor(x4)+length(ceiling(x4):ceiling(x4*2))-1)
    }
    stopifnot(is.numeric(x2))
    stopifnot(x3 >= 0)
    stopifnot(x4 >= 0.5)
    R_0<-initial_state_R/x3
    stopifnot(R_0 >= 0)
    stopifnot(R_0 <= 1)
    inAttr <- attributes(U)
    U <- as.ts(U)
    n <- ceiling(x4)
    m <- ceiling(x4 * 2)
    n2 <- floor(x4)
    SH1 <- pmin((1:n/x4)^(5/2), 1)
    SH2 <- pmin(c(0 + 0.5 * (1:n2/x4)^(5/2), 1 - 0.5 * (2 - n:m/x4)^(5/2)), 
                1)
    SH2[1:m/x4 > 2] <- 1
    UH1 <- diff(c(0, SH1))
    UH2 <- diff(c(0, SH2))
    
    stopifnot(length(initial_condition_UH1)+1 == length(UH1))
    stopifnot(length(initial_condition_UH2)+1 == length(UH2))
    
    bad <- is.na(U)
    U[bad] <- 0
    filter.pad0 <- function(x, f, initial_condition_UH=NA) {
      if(!return_components){
        y <- x
        x_with_zeros_at_beginning<-c(rep(0, length(f)), x)
        y[] <- filter(x_with_zeros_at_beginning, filter = f, sides = 1)[-(1:length(f))]
        if(!is.na(initial_condition_UH[1])){
          y[1:(length(initial_condition_UH))]<-y[1:(length(initial_condition_UH))]+initial_condition_UH
        }
        # f is the unit hydrograph, sums to 1
        return(list(flow=y,recorded_UH=NA))
      } else {
        # manual calculation
        # initial condition in unit hydrograph
        if(is.na(initial_condition_UH[1])){
          initial_condition_UH<-rep(0,length(f)-1)
          # initial_condition_UH<-rep(0,max(length(f)-1,1))
        }
        x_with_buffers<-c(x,rep(0,length(f)-1))
        y_manual<-rep(0,length(x_with_buffers))
        if(length(initial_condition_UH)>0){
          y_manual[1:length(initial_condition_UH)]<-initial_condition_UH
        }
        
        recorded_UH<-NULL
        for(i in 1:(length(y_manual)-length(f)+1)){
          current_UH<-y_manual[i:(i+length(f)-1)]+x_with_buffers[i]*f
          recorded_UH<-rbind(recorded_UH,current_UH,deparse.level = 0)
          y_manual[i:(i+length(f)-1)]<-current_UH
        }
        
        # get rid of buffers
        if(length(x_with_buffers)>length(x)){
          y_manual<-y_manual[-((length(y_manual)-length(f)+2):length(y_manual))]
        }
        
        
        #y
        return(list(flow=y_manual,recorded_UH=recorded_UH))
      }
    }
    Q9.run <- filter.pad0(split * U, UH1, initial_condition_UH1)
    Q1.run <- filter.pad0((1 - split) * U, UH2, initial_condition_UH2)
    
    Q9 <- Q9.run$flow
    Q1 <- Q1.run$flow
    #COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
    if (is.loaded("routing_gr4j_sk") & run_compiled) {
      
      if(is.na(state_R_ts[1])){
        state_R_ts<-as.double(rep(-1e6,n-1))
        use_state_R_ts<-as.integer(0)
      } else {
        state_R_ts<-as.double(state_R_ts)
        use_state_R_ts<-as.integer(1)
      }
      
      ans <- .C("routing_gr4j_sk", as.double(Q9), as.double(Q1), 
                as.integer(length(U)), as.double(x2), as.double(x3), 
                as.double(R_0), Qr = double(length(U)), Qd = double(length(U)), 
                R = double(length(U)), state_error_R = as.double(state_error_R),
                state_R_ts=state_R_ts,
                use_state_R_ts=use_state_R_ts,
                success=as.integer(1))
      if(ans$success!=1) stop("negative value in routing simulation when there shouldn't be")
      Qr <- ans$Qr
      Qd <- ans$Qd
      R <- ans$R
    }
    else {
      Qd <- Qr <- R <- U
      R_prev <- R_0 * x3
      for (t in seq(1, length(U))) {
        if(t>1){
          if(!is.na(state_R_ts[1])){
            if(state_R_ts[t-1]!=-1e6){
              R_prev<-state_R_ts[t-1]
            }
          }
          R_prev <- R_prev + state_error_R[t-1]
        }
        Rt_x3 <- R_prev/x3
        F <- x2 * Rt_x3^(7/2)
        R[t] <- max(0, R_prev + Q9[t] + F)
        Qr[t] <- R[t] * (1 - (1 + (R[t]/x3)^4)^(-0.25))
        R[t] <- R[t] - Qr[t]
        Qd[t] <- max(0, Q1[t] + F)
        R_prev <- R[t]
        if(R_prev<0) stop("negative routing state value")
        if(Rt_x3<0) stop("negative Rt_x3 value")
      }
    }
    Qr[Qr < epsilon] <- 0
    Qd[Qd < epsilon] <- 0
    attributes(Qr) <- attributes(Qd) <- attributes(R) <- inAttr
    Qr[bad] <- NA
    Qd[bad] <- NA
  },silent=T)
  
  if(length(grep("Error",tryerror))>0){
    return(NA)
  } else {
    if (return_components) {
      ans<-list(states_nonUH=cbind(Xr = Qr, Xd = Qd, R = R),states_UH1=Q9.run$recorded_UH,states_UH2=Q1.run$recorded_UH)
    } else {
      ans<-Qr + Qd
    }
    return(ans)
  }
}


gr4j.run<-function(param, initial_state_S, initial_state_R, state_error, input,
                   etmult = 1, return_state = FALSE, 
                   transformed = FALSE, run_compiled = T,
                   initial_condition_UH1=rep(0,ceiling(param[4])-1),
                   initial_condition_UH2=rep(0,floor(param[4])+length(ceiling(param[4]):ceiling(param[4]*2))-1),
                   state_error_R,
                   state_S_ts=NA,
                   state_R_ts=NA){
  output_sma<-gr4j.sma(param=param[1],
                       initial_state=initial_state_S, 
                       state_error=state_error, 
                       input=input, 
                       etmult = etmult, 
                       return_state = return_state, 
                       transformed = transformed, 
                       run_compiled = run_compiled,
                       state_S_ts=state_S_ts)
  if(is.na(output_sma[1])){
    return(NA)
  } else {
    if(return_state){
      routing_input<-output_sma[,1]
    } else {
      routing_input<-output_sma
    }
    output_rout<-gr4jrouting.sim.sk(routing_input,x2=param[2],x3=param[3],x4=param[4],
                                    initial_condition_UH1=initial_condition_UH1,
                                    initial_condition_UH2=initial_condition_UH2,
                                    return_components=return_state,initial_state_R=initial_state_R,
                                    run_compiled=run_compiled,
                                    state_error_R=state_error_R,
                                    transformed=transformed,
                                    state_R_ts=state_R_ts)
    if(any(is.na(output_rout))){
      return(NA)
    } else {
      if(return_state){
        output_rout_nonUH<-output_rout$states_nonUH
        output_flow<-output_rout_nonUH[,which(colnames(output_rout_nonUH)=="Xr")]+
          output_rout_nonUH[,which(colnames(output_rout_nonUH)=="Xd")]
        
        output_states_nonUH<-cbind(output_sma,output_rout_nonUH)
        
        all_output<-list(output_flow=output_flow,output_states_nonUH=output_states_nonUH,
                         output_states_UH1=output_rout$states_UH1,
                         output_states_UH2=output_rout$states_UH2)
        return(all_output)
      } else {
        return(output_rout)
      }
    }
  }
}







