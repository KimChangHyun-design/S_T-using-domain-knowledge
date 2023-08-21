library("fdaoutlier")

tvdmss2 <- function(dts,
                    emp_factor_mss = 1.5,
                    emp_factor_tvd = 1.5,
                    central_region_tvd = 0.5){

  depths_mss <- total_variation_depth(dts = dts)

  tvd <- tvd_old <- depths_mss$tvd
  mss <- depths_mss$mss


  dta_dim <- dim(dts)
  n_curves <- dta_dim[1]
  n_points <- dta_dim[2]
  index <- (1:n_curves)
  n_central_obs <- ceiling(n_curves/2)# set fbplot to .5

  # shape outliers
  shape_boxstats <- boxplot(mss, range = emp_factor_mss, plot=F);
  shape_outliers <- NULL
    
  if(length(shape_boxstats$out) != 0){
#segun:
#shape_outliers <- sapply(shape_boxstats$out, function(x) which(mss == x))
    #if (shape_boxstats$out < mean(mss)){
    shape_outliers <- which(mss %in% shape_boxstats$out[shape_boxstats$out < mean(mss)])
    if (length(shape_outliers) > 0){
      dts <- dts[-shape_outliers, ]
      tvd <- tvd[-shape_outliers]
      index <- index[-shape_outliers]   
    } 
  }
  magnitude_outliers  <- NULL
  
  # central region wrt to original number of curves
  outliers <- functional_boxplot(dts, depth_values = tvd,
                                 emp_factor = emp_factor_tvd,
                                 central_region = central_region_tvd)$outliers

  if(length(outliers) != 0) magnitude_outliers = index[outliers]
  return(list(outliers = sort(c(magnitude_outliers, shape_outliers)),
              shape_outliers = shape_outliers,
              magnitude_outliers = magnitude_outliers,
              tvd = tvd_old,
              mss = mss))
}

evaluation<-function(outlier_list,outlying,N){
  correct_anomaly <- intersect(outlying, outlier_list)
  TP <- length(correct_anomaly)
  FP <- length(outlier_list) - TP
  TN <- (N-length(outlying)) - FP
  FN <- length(outlying) - TP
  p.c <- round(TP/(TP+FN), 4)
  p.f <- round(FP/(FP+TN), 4)
  return(c(p.c,p.f))
}

plotting<-function(dd,outl=NULL,xl="",yl="",xaxis=NULL){
  if(is.null(xaxis)) xaxis<-c(1:ncol(dd))
  cc=adjustcolor("grey",alpha=1);lty=1
  if(1%in%outl){cc=2;lty=2}
  plot(xaxis,dd[1,],ylim=range(dd,na.rm=T),col=cc,lwd=2,type="l",ylab=yl,xlab=xl,cex.axis=1.8,cex.lab=2,lty=lty)
  for(i in setdiff(2:nrow(dd),outl)){
    lines(xaxis,dd[i,],col=cc,lwd=2,lty=1)
  }
  for(i in outl){
    lines(xaxis,dd[i,],col=adjustcolor("red",alpha=1),lwd=2,lty=2)
  }
}



## (transformation == "T0" || transformation == "D0")
functional_boxplot

## (transformation == "T1")

center_curves <- function(dt){
  return(dt - rowMeans(dt))
}

## (transformation == "T2"){

normalize_curves <- function(dt){
  return(dt/sqrt(rowSums(dt^2)))
}
                        
## (transformation == "T3"){
                        
nor_minmax = function(x){
  result = (x - min(x)) / (max(x) - min(x))
  return(result)
}
                        
center_curves2 <- function(dt){
  dt<-apply(X = dt, MARGIN = 2, FUN = "nor_minmax")  
  return(dt * rowMeans(dt))
}

ff_t<-function(dt){
  dd<-apply(X = dt, MARGIN = 1, FUN = "fft")  
  return(t(Mod(dd)))
}


## (transformation == "D1"|| transformation == "D2")
difference_curves <- function(dt){
  p <- dim(dt)[2]
  dt[,2:p] - dt[, 1:(p-1)]
}

## transformation == "O"
outlyingness_curves <- function(dt, n_projections = 500L, seed = NULL){
  dm <- dim(dt)
  repnp <- rep(dm[1], dm[2])
  if(length(dm) == 2){# data is univariate
    median_vec <- apply(dt, 2, median)
    mad_vec <- apply(dt, 2, mad)
    return(abs(dt- rep(median_vec, repnp))/rep(mad_vec, repnp))
  } else if (length(dm) == 3) { # data is multivariate # 100,50(t),dimension
    outlyingness <- apply(dt, 2, function(x){
      (1/projection_depth(x, n_projections = n_projections, seed = seed)) - 1
    })
    return(outlyingness)
  } else{
    stop("Argument \"dt\" must be a 2 or 3 dimensional array.")
  }
}


temporal_knowledge <- function(dt,sec){
  return(dt[,sec])
}




#######################

functional_boxplot2<-function (dts, depth_method = c("mbd", "tvd", "extremal", "dirout", 
                                                     "linfinity", "bd", "erld", "dq"),
                               depth_values = NULL, emp_factor = 1.5, direction = NULL, outlier_rate = NULL,
                               central_region = 0.5, erld_type = NULL, dq_quantiles = NULL){
  dm <- dim(dts)
  n <- dm[1]
  p <- dm[2]
  if (is.data.frame(dts)) {
    dt <- as.matrix(dts)
  }
  if (any(!is.finite(dts))) {
    stop("Missing or infinite values are not allowed in argument \"dts\"")
  }
  if (!is.array(dts) || !is.numeric(dts)) 
    stop("Argument \"dts\" must be a numeric matrix or dataframe.")
  if (length(dm) != 2) 
    stop("Dimension of 'dts' must be of length 2. Only univariate functional data is supported.")
  if (is.null(depth_values)) {
    depth_method <- match.arg(depth_method)
    if (depth_method == "mbd") {
      depth_values <- modified_band_depth(dts)
    }
    else if (depth_method == "tvd") {
      depth_values <- total_variation_depth(dts)$tvd
    }
    else if (depth_method == "extremal") {
      depth_values <- extremal_depth(dts)
    }
    else if (depth_method == "dirout") {
      depth_values <- -dir_out(dts, return_distance = T)$distance
    }
    else if (depth_method == "linfinity") {
      depth_values <- linfinity_depth(dts)
    }
    else if (depth_method == "bd") {
      depth_values <- band_depth(dts)
    }
    else if (depth_method == "erld") {
      if (is.null(erld_type)) {
          
        if (direction == "lower"){
            erld_type = "one_sided_left"
            depth_values <- extreme_rank_length(dts, type = erld_type)
        }
        else if (direction == "upper"){
            erld_type = "one_sided_right"
            depth_values <- extreme_rank_length(dts, type = erld_type)
        }
        else if (direction == "both"){
            erld_type = "two_sided"
            depth_values <- extreme_rank_length(dts, type = erld_type)
        }        
      }
      else {
        warning("The 'type' argument for extreme rank length depth not specified. Using the default type of 'two_sided'. ")
        depth_values <- extreme_rank_length(dts)
        }
    }
    else if (depth_method == "dq") {
      if (is.null(dq_quantiles)) {
        warning("Using the default quantile probabilites of 0.025 and 0.975 for directional quantile.")
        depth_values <- -directional_quantile(dts)
      }
      else {
        depth_values <- -directional_quantile(dts, quantiles = dq_quantiles)
      }
    }
  }
  else {
    if (length(depth_values) != n) {
      stop("Length of argument 'depth_values' must be equal to the number of rows in 'dts'.")
    }
  }
  if (central_region >= 1 || central_region <= 0) {
    stop("Argument 'central_region' must be greater than 0 and less than 1.")
  }
  
  sorted_depths <- sort(depth_values, decreasing = T, index.r = T)
  index_sorted_depth <- sorted_depths$ix
  sorted_depths <- sorted_depths$x
  median_curve <- index_sorted_depth[1]
  
  if (depth_method == "erld"){
    threshold <- quantile(depth_values, outlier_rate)
    outliers <- which(depth_values < threshold)
  }else{
    n_obs_central <- ceiling(n * central_region)
    center <- dts[index_sorted_depth[1:n_obs_central], ]
    inf <- apply(center, 2, min)
    sup <- apply(center, 2, max)
    distt <- emp_factor * (sup - inf)
    upper <- sup + distt
    lower <- inf - distt
    dts <- t(dts)
    outlier_test <- (dts <= lower) + (dts >= upper)
    outliers <- which(colSums(outlier_test) > 0)
  }
  
  
  return(list(outliers = unname(outliers), depth_values = depth_values, 
              median_curve = median_curve))
}


seq_transform2<-function (dts, dts1=NULL, sequence = c("T0", "T1", "T2"),
                          depth_method = c("mbd", "tvd", "extremal", "dirout", "linfinity", "bd", "erld", "dq"), 
                          save_data = FALSE, emp_factor = 1.5, central_region = 0.5, t_knlg,d_knlg=1, outlier_rate, direction,
                          erld_type = NULL, dq_quantiles = NULL, n_projections = 200L, 
                          seed = NULL) {
  outliers <- list()
  if (save_data) {
    transformed_data <- list()
  }
  for (i in seq_along(sequence)) {
    transformation = sequence[i]
    if (transformation == "T0" || transformation == "D0") {
      t0_outliers <- functional_boxplot2(dts, depth_method = depth_method, 
                                        central_region = central_region, emp_factor = emp_factor, direction = direction,
                                        outlier_rate = outlier_rate, erld_type = erld_type, dq_quantiles = dq_quantiles)$outliers
      outliers[[i]] <- t0_outliers
      if (save_data) 
        transformed_data[[i]] <- dts
    }
    else if (transformation == "TK") {
      dts <- temporal_knowledge(dts,t_knlg)
      dts1<-temporal_knowledge(dts1,t_knlg)
      tk_outliers <- functional_boxplot2(dts, depth_method = depth_method, 
                                        emp_factor = emp_factor, central_region = central_region, direction = direction,
                                        outlier_rate = outlier_rate, erld_type = erld_type, dq_quantiles = dq_quantiles)$outliers
      outliers[[i]] <- tk_outliers
      if (save_data) 
        transformed_data[[i]] <- dts
    }
    else if (transformation == "DK") {
      dts <- dts*d_knlg
      dk_outliers <- functional_boxplot2(dts, depth_method = depth_method, 
                                         emp_factor = emp_factor, central_region = central_region, direction = direction,
                                         outlier_rate = outlier_rate, erld_type = erld_type, dq_quantiles = dq_quantiles)$outliers
      outliers[[i]] <- dk_outliers
      if (save_data) 
        transformed_data[[i]] <- dts
    }
    else if (transformation == "SK") {
      dts <- dts1/dts
      sk_outliers <- functional_boxplot2(dts, depth_method = depth_method, 
                                         emp_factor = emp_factor, central_region = central_region, direction = direction,
                                         outlier_rate = outlier_rate, erld_type = erld_type, dq_quantiles = dq_quantiles)$outliers
      outliers[[i]] <- sk_outliers
      if (save_data) 
        transformed_data[[i]] <- dts
    }    
    else if (transformation == "T1") {
      dts <- center_curves(dts)
      t1_outliers <- functional_boxplot2(dts, depth_method = depth_method, 
                                        emp_factor = emp_factor, central_region = central_region, direction = direction, 
                                        outlier_rate = outlier_rate, erld_type = erld_type, dq_quantiles = dq_quantiles)$outliers
      outliers[[i]] <- t1_outliers
      if (save_data) 
        transformed_data[[i]] <- dts
    }
    else if (transformation == "T2") {
      dts <- normalize_curves(dts)
      t2_outliers <- functional_boxplot2(dts, depth_method = depth_method, 
                                        emp_factor = emp_factor, central_region = central_region, direction = direction, 
                                        outlier_rate = outlier_rate, erld_type = erld_type, dq_quantiles = dq_quantiles)$outliers
      outliers[[i]] <- t2_outliers
      if (save_data) 
        transformed_data[[i]] <- dts
    }
    else if(transformation == "T3"){
      dts <- center_curves2(dts)
      t3_outliers <- functional_boxplot2(dts, depth_method = depth_method,
                                        emp_factor = emp_factor, central_region = central_region, direction=direction,
                                        outlier_rate = outlier_rate, erld_type = erld_type, dq_quantiles = dq_quantiles)$outliers
      outliers[[i]] <- t3_outliers
      if (save_data) 
        transformed_data[[i]] <- dts
    }
    else if(transformation == "T5"){
      dts <- ff_t(dts)
      t5_outliers <- functional_boxplot2(dts, depth_method = depth_method,
                                         emp_factor = emp_factor, central_region = central_region, direction=direction,
                                         outlier_rate = outlier_rate, erld_type = erld_type, dq_quantiles = dq_quantiles)$outliers
      outliers[[i]] <- t5_outliers
      if (save_data) 
        transformed_data[[i]] <- dts
    }
    else if (transformation == "D1" || transformation ==  "D2") {
      dts <- difference_curves(dts)
      d1_outliers <- functional_boxplot2(dts, depth_method = depth_method, 
                                        emp_factor = emp_factor, central_region = central_region, direction = direction,  
                                        outlier_rate = outlier_rate, erld_type = erld_type, dq_quantiles = dq_quantiles)$outliers
      outliers[[i]] <- d1_outliers
      if (save_data) 
        transformed_data[[i]] <- dts
    }
    else if (transformation == "O") {
      dts <- outlyingness_curves(dts, n_projections = n_projections, 
                                 seed = seed)
      o_outliers <- functional_boxplot2(dts, depth_method = depth_method, 
                                       emp_factor = emp_factor, central_region = central_region, direction = direction, 
                                       outlier_rate = outlier_rate, erld_type = erld_type, dq_quantiles = dq_quantiles)$outliers
      outliers[[i]] <- o_outliers
      if (save_data) 
        transformed_data[[i]] <- dts
    }
    else {
      stop("Transformation ", transformation, " not supported. \n")
    }
  }
  if (any(duplicated(sequence))) {
    tt <- sequence[duplicated(sequence)][1]
    sequence[sequence == tt] <- paste0(tt, "_", 1:length(sequence[sequence == 
                                                                    tt]))
    warning("Duplicated transforms found in argument 'sequence',\n            changing output labels to: ", 
            paste(sequence, collapse = " "), ".")
  }
  names(outliers) <- sequence
  if (save_data) 
    names(transformed_data) <- sequence
  if (save_data) {
    return(list(outliers = outliers, transformed_data = transformed_data, transformation = transformation))
  }
  else {
    return(list(outliers = outliers, transformed_data = NULL, transformation = transformation))
  }
}

# Domain knowledge-informed transformation (R 막힘)
DT_R <- function(dt = NULL, t_knlg = R_sec, outlier_rate = outlier_rate){
    
    dt0 <- dt[,,3]
    dt1 <- dt[,,2]
    seq_R <- seq_transform2(dts = dt0, dts1 = dt1, sequence = c("SK", "TK", "DK"), t_knlg = t_knlg, depth_method = "erld", direction = "upper", outlier_rate = outlier_rate)
    outliers_R <- seq_R$outliers$DK
    return(outliers_R)
}

# Domain knowledge-informed transformation (F 막힘)

DT_F <- function(dt = NULL, t_knlg = F_sec, outlier_rate = outlier_rate){
    
    dt0 <- dt[,,3]
    dt1 <- dt[,,1]
    seq_F <- seq_transform2(dts = dt0, dts1 = dt1, sequence = c("SK", "T1", "T2", "TK", "DK"), t_knlg = t_knlg, depth_method = "erld", direction = "upper", outlier_rate = outlier_rate)
    outliers_F <- seq_F$outliers$DK
    return(outliers_F)    
    
}

# Framework

Domain_transformation_framework <- function(dt = NULL, DT = c("DT_R", "DT_F"), outlier_rate = c(0.01, 0.03, 0.05, 0.1)){
    
    total_outliers <- c()
    for (i in DT){
        if (i == "DT_R"){
            outliers_R <- DT_R(dt, outlier_rate = outlier_rate)
            total_outliers <- union(total_outliers, outliers_R)
        }
        else if (i == "DT_F"){
            outliers_F <- DT_F(dt, outlier_rate = outlier_rate)
            total_outliers <- union(total_outliers, outliers_F)
        }
    }
    return(total_outliers)
    
}

# Comparison with baseline
Experiments <- function(dt = NULL, outlier_rate = c(0.01, 0.03, 0.05, 0.1), target_num = target_num){
    
    evl_table <- c()
    N <- dim(dt)[1]
    power <- dt[,,3]
    outlying <- c(target_num)
    outlier_rate <- 0.01
    
    # Baseline
    msplot_object <- msplot(dts = power, return_mvdir=FALSE, plot=F)
    msplot_rst <- evaluation(msplot_object$outliers,outlying,N=N)

    TVD_object <- tvdmss2(power)
    TVD_rst <- evaluation(TVD_object$outliers,outlying,N=N)

    MUOD_object <- muod(power, cut_method = c("boxplot", "tangent"))
    MUOD_rst <- evaluation(MUOD_object$outliers$magnitude,outlying,N=N)

    seq0 <- seq_transform(dts = power, sequence = c("T1", "T2"),  depth_method = "erld", erld_type = "two_sided")
    ST_rst <- evaluation(seq0$outliers$T2,outlying,N=N)
    
    # Domain ST
    DT_outs <- Domain_transformation_framework(dt, DT = c("DT_R", "DT_F"), outlier_rate = outlier_rate)
    DST_rst <- evaluation(DT_outs,outlying,N=N)
    
    # Final results
    evl_table<-data.frame(rbind(evl_table,c(msplot_rst,TVD_rst,MUOD_rst,ST_rst,DST_rst)))
    names(evl_table)<-paste0(c("msplot_pc","msplot_pf","TVD_pc","TVD_pf","MUOD_pc","MUOD_pf","ST_pc","ST_pf","DST_pc","DST_pf"))
    evl_table
    
}

# Generate_simulation dataset
generate_simulation <- function(simulation = c("simulation1", "simulation2", "simulation3",
                                               "simulation4", "simulation5", "simulation6",
                                               "auxiliary1", "auxiliary2"), n = 1000, p = 100, outlier_rate = 0.05, seed = NULL){
    if (simulation == "simulation1"){
        dt0 <- simulation_model1(n = n, p = p, outlier_rate = outlier_rate,
                                 mu = 4, q = 2, kprob = 1.0,
                                 cov_alpha = 0.5, cov_beta = 1, cov_nu = 1,
                                 deterministic = TRUE, seed = seed)
        #dt0$data <- dt0$data + 8
    
    }else if (simulation == "simulation2"){
        dt0 <- simulation_model1(n = n, p = p, outlier_rate = outlier_rate,
                                 mu = 4, q = 2, kprob = 0,
                                 cov_alpha = 0.5, cov_beta = 1, cov_nu = 1,
                                 deterministic = TRUE, seed = seed)
        #dt0$data <- dt0$data + 8
        
    }else if (simulation == "simulation3"){
        dt0 <- simulation_model6(n = n, p = p, outlier_rate = outlier_rate,
                              mu = 4, q = 1, a = 0.2, b = 0.25, kprob = 1.0,
                              pi_coeff = 0.1, exp_pow = 2, exp_coeff = 20,
                              cov_alpha = 0.5, cov_beta = 1, cov_nu = 1,   
                              deterministic = TRUE, seed = seed)
        #dt0$data <- dt0$data + 8
        
    }else if (simulation == "simulation4"){
        dt0 <- simulation_model6(n = n, p = p, outlier_rate = outlier_rate,
                              mu = 4, q = 1, a = 0.2, b = 0.25, kprob = 0,
                              pi_coeff = 0.1, exp_pow = 2, exp_coeff = 20,
                              cov_alpha = 0.5, cov_beta = 1, cov_nu = 1,   
                              deterministic = TRUE, seed = seed)
        #dt0$data <- dt0$data + 8
        
    }else if (simulation == "simulation5"){
        dt0 <- simulation_model4(n = n, p = p, outlier_rate = outlier_rate,
                                 mu = 30, m = 3/2, # default mu=30, m=3/2
                                 cov_alpha = 0.3, cov_beta = (1/0.3), cov_nu = 1,
                                 deterministic = TRUE, seed = seed)
     
    }else if (simulation == "simulation6"){
        dt0 <- simulation_model7(n = n, p = p, outlier_rate = outlier_rate,
                                 a = 0.2, b = 0.25, sin_coeff = 1, pi_coeff = 5,
                                 deterministic = TRUE, seed = seed)
        dt0$data <- dt0$data + 8        
    #sin(2*pi_coeff/2/100 * pi * (tt + theta))  f=2.5, N =100
    }else if (simulation == "auxiliary1"){
        dt0 <- simulation_model1(n = n, p = p, outlier_rate = outlier_rate,
                                 mu = 4, q = 2, kprob = 0,
                                 cov_alpha = 1, cov_beta = 1, cov_nu = 1,
                                 deterministic = TRUE, seed = seed)
        dt0$data <- dt0$data + 6
    }else if (simulation == "auxiliary2"){
        dt0 <- simulation_model1(n = n, p = p, outlier_rate = outlier_rate,
                                 mu = 4, q = 2, kprob = 1,
                                 cov_alpha = 1, cov_beta = 1, cov_nu = 1,
                                 deterministic = TRUE, seed = seed)
        dt0$data <- dt0$data + 6
    }
    return (dt0)
}

# Plotting function
plotting2<-function(dd,outl=NULL,abnml=NULL,title=NULL){
  cc=1
  if(1%in%outl){cc=2}
  rescaleX<-seq(0,1.0,length.out=dim(dd)[2])
  plot(rescaleX, dd[1,],ylim=range(dd,na.rm=T),col=cc,lwd=cc,main=title,type="l",xaxt="n",yaxt="n",ann=FALSE)
  axis(1,cex.axis=1.5)
  axis(2,cex.axis=1.5)
  
  for(i in 1:nrow(dd)){
    col = "darkgrey"
    lty = "solid"
    if(i%in%outl & !(i%in%abnml))
    {col="Light Coral"
    cc=2
    lty="dotted"}
    lines(rescaleX,dd[i,],col=col,lwd=cc,lty=lty)
    cc=1
    lines(rescaleX,target_curve,col=2,lwd=2,lty="dashed")
  }
}
