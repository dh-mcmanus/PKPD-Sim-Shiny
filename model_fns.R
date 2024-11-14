
## Library
library(shiny)
library(dplyr)
library(ggplot2)
library(tidyr)
library(deSolve)
library(gridExtra)
library(EnvStats)

##
## Sample dropout time
##
fn_simweibull <- function(scale, shape) {
  A <- scale
  B <- 1/shape
  U <- runif(1,0,1)
  T <- (-1*log(U)/A)^B
  return(T)
}

##
## Beta Sampling
##
fn_samplebeta <- function(mu,cv) {
  sd <- mu*cv
  beta.params  <- getBetaParams(mu,sd^2)
  result <- rbeta(1,beta.params$alpha,beta.params$beta)
  return(result)
}
getBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

##
## Recruit the desired number of subjects based on trial design and arm settings
##
fn_recruitment <- function(nsubjects, 
                           trialdata, 
                           iiv) {
  
  ## Gout population TVs
  sua_mu  <- trialdata$mean[trialdata$para=="sua_mu"]
  wt_mu   <- trialdata$mean[trialdata$para=="wt_mu"]
  age_mu  <- trialdata$mean[trialdata$para=="age_mu"]
  vua_mu  <- trialdata$mean[trialdata$para=="vua_mu"]
  
  ## Gout population IIVs
  sua_eta <- trialdata$mean[trialdata$para=="sua_eta"]
  wt_eta  <- trialdata$mean[trialdata$para=="wt_eta"]
  age_eta <- trialdata$mean[trialdata$para=="age_eta"]
  vua_eta <- trialdata$mean[trialdata$para=="vua_eta"]
  
  ## IIV
  if (iiv==1) {
    
    ## SUA
    sua_sd <- sua_mu*sua_eta
    sua_logmean  <- log(sua_mu / (sqrt(1 + (sua_sd/sua_mu)^2)))
    sua_logsd    <- sqrt(log(1 + (sua_sd/sua_mu)^2))
    sua <- rlnormTrunc(nsubjects, meanlog = sua_logmean, sdlog = sua_logsd, min = 0, max = 150)  
    
    ## Weight
    wt_sd <- wt_mu*wt_eta
    wt_logmean <- log(wt_mu / (sqrt(1 + (wt_sd/wt_mu)^2)))
    wt_logsd <- sqrt(log(1 + (wt_sd/wt_mu)^2))
    wt <- rlnormTrunc(nsubjects, meanlog = wt_logmean, sdlog = wt_logsd, min = 50, max = 150)
    
    ## Age
    age_sd <- age_mu*age_eta
    age_logmean <- log(age_mu / (sqrt(1 + (age_sd/age_mu)^2)))
    age_logsd <- sqrt(log(1 + (age_sd/age_mu)^2))
    age <- rlnormTrunc(nsubjects, meanlog = age_logmean, sdlog = age_logsd, min = 40, max = 95)
    
    ## Volume of uric acid distribution
    vua  <- vua_mu*exp(rnorm(nsubjects,0,vua_eta))  
    
  }
  
  else {
    ## 
    sua <- rep(sua_mu,nsubjects)  
    wt  <- rep(wt_mu,nsubjects)  
    age <- rep(age_mu,nsubjects)  
    vua <- rep(vua_mu,nsubjects)  
    
  }
  
  ## Other inputs
  bua  <- sua*vua
  bxan = 8.94
  clua <- 4.11*(703/bua)
  
  ## Creatinine clearance calculation
  crcl = (((140-age)*wt*1.23)/110)-15
  crcl <- ifelse(crcl<45,45,ifelse(crcl>160,160,crcl))

  ## Combine in dataframe
  subset.patients <- data.frame(sua = sua,
                                bxan = bxan,
                                wt = wt,
                                age = age,
                                crcl = crcl,
                                vua = vua,
                                clua = clua)
                      
  return(subset.patients)
}

##
## Create dose history
##
mk.doseevents <- function(nsubj,
                          ndays,
                          dose,
                          imp_mu, 
                          dropout_mu,
                          do_scale,
                          do_shape,
                          timing_yn,
                          iiv) {
  
  ## Preliminaries
  popimp.eta <- imp_mu*0.1
  
  ## Result vector
  fullhist <- NULL
  
  ## Loop over each patient
  for (inti in 1:nsubj) {
    
    ## Inter-individual sampling
    if (iiv==1) {
      imp.mu_tmp <- ifelse(imp_mu>=0.98,0.98,imp_mu)
      i.imp  <- fn_samplebeta(imp.mu_tmp,popimp.eta)
      i.dropout <- ifelse(seq(0,(ndays-1),by=1) > fn_simweibull(do_scale,do_shape), 1, 0)
    }
    else {
      i.imp <- imp_mu
      i.dropout <- ifelse(seq(0,(ndays-1),by=1) > dropout_mu, 1, 0)
    }
    
    ## Dose events for implementation phase
    i.take  <- floor(runif(ndays,i.imp,i.imp+1))

    ## Create dosing history: doses of 80mg are for doses taken, 0mg are doses missed
    tmphist <- data.frame(id = inti,
                          var = "A1",
                          time = seq(8,(ndays)*24,by=24) + if(timing_yn==TRUE) {runif(ndays, min = -2, max = 6)} else {rep(0,ndays)},
                          dropout = i.dropout,
                          value = i.take*dose*(1-i.dropout),
                          method = "add") %>%
      mutate(day = floor(time/24)) %>%
      mutate(tod = time - day*24)
    
    ## Add subject data
    fullhist <- rbind(fullhist,tmphist)
  }
  return(fullhist)
}

##
## Get inital values for PKPD model compartments
##
fn_pkpdinits <- function(patient) {
  return(c(A1 = 0,
           A2 = 0,
           A3 = 0,
           AX = patient$bxan,
           AUA = patient$sua*patient$vua))
}

##
## Get Febuxostat PKPD model parameters
##
fbx_pkpdparams <- function(patient,
                           pkparams,
                           pdparams,
                           iiv) {
  
  ## Clearance
  cl_tv   <- pkparams$mean[pkparams$para=="cl_mu"]
  cl_crcl <- pkparams$mean[pkparams$para=="cl_crcl"]
  cl_wt   <- pkparams$mean[pkparams$para=="cl_wt"]
  
  ## Volumes, absorptionand inter-compartmental clearance
  v2_tv   <- pkparams$mean[pkparams$para=="v2_mu"]
  ka_tv   <- pkparams$mean[pkparams$para=="ka_mu"]
  q_tv    <- pkparams$mean[pkparams$para=="q_mu"]
  v3_tv   <- pkparams$mean[pkparams$para=="v3_mu"]
  
  ## Etas for random effects
  cl_eta  <- pkparams$mean[pkparams$para=="cl_eta"]
  ka_eta  <- pkparams$mean[pkparams$para=="ka_eta"]
  
  ## Drug PD parameters
  emax_mu    <- pdparams$mean[pdparams$para=="emax1"]
  ec50_mu    <- pdparams$mean[pdparams$para=="ec50"]
  ic50.1_mu  <- pdparams$mean[pdparams$para=="ic501"]
  ic50.2_mu  <- pdparams$mean[pdparams$para=="ic502"]
  ic50.1_eta <- pdparams$mean[pdparams$para=="ic501_eta"]
  ic50.2_eta <- pdparams$mean[pdparams$para=="ic502_eta"]
  
  ## Systems parameters
  bxan_mu  <- patient$bxan
  vxan_mu  <- 333
  clxan_mu <- 10.6
  vua_mu   <- patient$vua
  clua_mu  <- patient$clua
  
  ## Patient covariates
  wt <- patient$wt
  age <- patient$age
  crcl = patient$crcl # adjusted cockcroft-gault eqn
  
  ## TVs
  V2 = v2_tv   
  Q  = q_tv
  V3 = v3_tv 
  
  ## IIV
  if (iiv==1) {
    CL = (cl_tv + crcl*cl_crcl + wt*cl_wt)*exp(rnorm(1,0,cl_eta))
    KA = ka_tv*exp(rnorm(1,0,ka_eta))
    ic50.1 = ic50.1_mu*exp(rnorm(1,0,ic50.1_eta))
    ic50.2 = ic50.2_mu*exp(rnorm(1,0,ic50.2_eta))
  }
  else {
    CL = cl_tv + crcl*cl_crcl + wt*cl_wt
    KA = ka_tv
    ic50.1 = ic50.1_mu
    ic50.2 = ic50.2_mu
  }
  
  ## xanthine
  vx  <- vxan_mu # dl
  clx <- clxan_mu # mg/dl
  x.0 <- bxan_mu # mg
  
  ##
  sua <- patient$sua
  vua <- vua_mu
  bua <- sua * vua
  clua.0 <- clua_mu
  bua.0 <- patient$sua*vua 
  clua <- clua.0*(bua.0/bua)
  
  ## Ratio of uric acid to xanthine molar masses
  r <- 168.11/152.11
  
  # rate parameters
  k2=clx/vx # xanthine excretion
  k3=clua/vua # uric acid excretion
  k1=k3*bua/(x.0*r)
  k0=k1*x.0+k2*x.0
  
  return(c(ka = KA,
           k20 = CL/V2,
           k23 = Q/V2,
           k32 = Q/V3,
           V2 = V2,
           emax.1 = emax_mu,
           ec50.1 = ec50_mu,
           ic50.1 = ic50.1,
           ic50.2 = ic50.2,
           k0 = k0,
           k1 = k1,
           k2 = k2,
           k3 = k3,
           vua = vua))
}

##
## Function to implement the 2-compartment PKPD model
##
fbx_pkpdmodel <- function(t, inits, params) {
  
  with(as.list(c(inits, params)),{
    
    ## Ratio of uric acid to xanthine molar masses
    r <- 168.11/152.11
    
    ## rate of change
    dA1 <- -ka*A1
    dA2 <- ka*A1 - k20*A2 - k23*A2 + k32*A3
    dA3 <- k23*A2 - k32*A3
    
    ## Convert to drug concentration
    A2 <- ifelse(A2<0,0,A2)
    CP <- A2/V2
    
    ## Compute the PD models
    stim1 <- 1 + (emax.1*CP)/(ec50.1+CP)
    inh3  <- (ic50.1)/(ic50.1+CP)
    inh4  <- (ic50.2)/(ic50.2+CP)      
    
    ## Rate of change
    dAX  <- k0*inh4 - k1*AX*inh3 - k2*AX*stim1
    dAUA <- k1*AX*r*inh3 - k3*AUA
    
    ## Retrun the rate of change
    list(c(dA1, dA2, dA3, dAX, dAUA))
    
  }) 
}

##
## Simulate subjects in trial PKPD
##
TrialSim <- function(#name, # all or fbx
                     nsubj,
                     ndays,
                     pkinputs,
                     pdinputs,
                     cohort,
                     dhistories,
                     iiv) {

  ## Store outputs
  allout <- NULL
  
  ## Simulate each trial participant
  for (ipat in 1:nsubj) {
    
    ## Extract patient details
    tmp.patient <- cohort[ipat,]
    
    ## Extract subject dosing history    
    ihistory <- dhistories[dhistories$id==ipat,]
    times <- sort(ihistory$time,decreasing=F)
    
    ##
    iparams <- fbx_pkpdparams(tmp.patient, pkinputs, pdinputs, iiv)
    model <- fbx_pkpdmodel

    ## Seperate VUA from other params
    ivua    <- dplyr::last(iparams)
    iparams <- iparams[-length(iparams)]
    
    ## Setup PKPD model inital values
    iinits  <- fn_pkpdinits(tmp.patient)
    
    ## Compute sUA time course
    SimOut <- as.data.frame(ode(func   = model,
                                times  = times,
                                y      = iinits,
                                parms  = iparams,
                                method = "impAdams_d",
                                events = list(data = ihistory)))
    
    ## Convert from amount to concentration
    SimOut$SUA <- SimOut$AUA/ivua
    
    ## Extract final sUA measurement
    allout <- rbind(allout,
                    data.frame(id = ipat,
                               SimOut$time,
                               SimOut$A2,
                               SimOut$A2/iparams[5],
                               SimOut$SUA))
  }
  ##
  names(allout) <- c("id","time","A2","CP","sua")
  return(allout)
}

##
## Simulate subjects
##
fn_runsim1 <-  function(nsubj,
                        ndays,
                        dose,
                        imp_frac, 
                        n_do,
                        do_scale,
                        do_shape,
                        timing_yn,
                        df_subj,
                        df_pk,
                        df_pd,
                        iiv) {
  
  ##
  df_pk <- df_pk %>% rename("para" = "Parameter",
                            "mean" = "Mean",
                            "unit" = "Unit")
  df_pd <- df_pd %>% rename("para" = "Parameter",
                            "mean" = "Mean",
                            "unit" = "Unit")
  df_subj <- df_subj %>% rename("para" = "Parameter",
                                "mean" = "Mean",
                                "unit" = "Unit")
  
  ##
  cohort <- fn_recruitment(nsubj,
                           df_subj,
                           iiv)
  
  ##
  data <- mk.doseevents(nsubj,
                        ndays,
                        dose,
                        imp_frac, 
                        n_do, 
                        do_scale,
                        do_shape,
                        timing_yn,
                        iiv)

  ##
  return(Sim1 <- TrialSim(nsubj,
                          ndays,
                          df_pk,
                          df_pd,
                          cohort,
                          data,
                          iiv) %>%
           cbind(data %>% dplyr::select(-time,-id)))

}


##
##
##
fn_plot_adherence <- function(data, ndays) {
  
  plotdata <- data %>%
    mutate(day2 = floor(time/24)) %>%
    mutate(time2 = time - day*24+5) %>%
    mutate(take = ifelse(value==0,NA,1))
  
  
  plotdata %>% ggplot() +
    geom_point(aes(x=day2,y=time2*take,color="blue"),size=3,color="darkblue") +
    geom_vline(data = filter(plotdata,value==0), aes(xintercept = as.numeric(day2)),size=1,color="grey") +
    ylim(0,23.8) + ylab("Timnig of dose") + xlab("Day of dose") +
    theme(axis.line = element_line(colour = "black"),
          axis.text.y=element_text(size=0),
          axis.ticks.y=element_blank(),
          axis.text.x=element_text(size=14),
          text = element_text(size=14))
  
  
}


##
##
##
fn_plot_pkpd1 <- function(data, ndays) {
  
  ## Facet names
  fnames <- c(
    "CP" = "Drug Plasma Concentration",
    "sua" = "Serum Uric Acid Concentration"
  )
  
  ##
  plotdata <- data %>%
    mutate(day = floor(time/24)) %>%
    select(day,CP,sua) %>%
    gather(key,value,-day)
  
  ##
  plotdata %>%  ggplot() +
    geom_line(aes(x = day, y = value), color="#0000FF") +
    geom_hline(data = plotdata %>% filter(key == "sua"),
               aes(yintercept = 6), color="green") +
    scale_x_continuous(limits = c(0,ndays)) +
    expand_limits(y=0) +
    xlab("Time (days)") + 
    ylab("mg/dL") +  
    theme_bw()  +
    facet_wrap(~key, 
               nrow=2, 
               scales="free",
               labeller = as_labeller(fnames)) +
    theme(axis.line = element_line(colour = "black"),
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14),
          text = element_text(size=14))
  
}

##
##
##
fn_plot_pkpd2 <- function(data, ndays) {
  
  ## Facet names
  fnames <- c(
    "CP" = "Drug Plasma Concentration",
    "sua" = "Serum Uric Acid Concentration"
  )

  ##
  plotdata <- data %>%
    mutate(day = floor(time/24)) %>%
    select(id,day,CP,sua) %>%
    gather(key,value,-day, -id)
  
  ##
  plotdata %>%  ggplot() +
    geom_line(aes(x = day, y = value, color=as.factor(id))) +
    geom_hline(data = plotdata %>% filter(key == "sua"),
               aes(yintercept = 6), color="green") +
    scale_x_continuous(limits = c(0,ndays)) +
    scale_color_discrete(name = "Subject ID") +
    expand_limits(y=0) +
    xlab("Time (days)") + 
    ylab("mg/dL") +  
    theme_bw()  +
    facet_wrap(~key, 
               nrow=2, 
               scales="free",
               labeller = as_labeller(fnames)) +
    theme(legend.position = "none",
      axis.line = element_line(colour = "black"),
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14),
          text = element_text(size=14))
  
  
}

##
##
##
fn_plot_hist <- function(data, ndays) {
  
  data %>%
    group_by(id) %>%
    summarise(last = last(sua),
              first = first(sua)) %>%
    gather(key,value,-id) %>%
    ggplot(aes(x=value, fill=key)) +
    geom_histogram(aes(y=..density..), alpha=0.5, 
                   position="identity") +
    geom_density(alpha=.2) +
    xlab("Final sUA Conc. (mg/dL)") + 
    ylab("Density") +  
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          axis.text.y=element_text(size=14),
          axis.text.x=element_text(size=14),
          text = element_text(size=14))
}

##
##
##
fn_meansUA <- function(data, cutoff) {
  data %>%
    mutate(day = floor(time/24)) %>%
    filter(day > 14) %>%
    summarise(mean(sua))
}

##
##
##
fn_respondersUA <- function(data, cutoff, nsubj) {
  data %>%
    group_by(id) %>%
    summarise(final_sua = last(sua)) %>%
    ungroup() %>%
    summarise(respond = sum(final_sua < cutoff)/nsubj)
}
