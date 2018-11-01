## --- Utility functions --- ##

# -- for performance improvement -- #
rev_cumsum <- function(x){sum(x)-c(0L,cumsum(x[-length(x)]))}

# -- simulation wrapper -- #
sim_fcn <- function(n0,n1,mean, ...){
  arguments <- list(sapply(mean,function(x){rnorm(n0)}),sapply(mean,function(x){rnorm(n1,mean=x)}),...)
  do.call(pauc_fcn,arguments)
}

# -- Holm-Bonferroni implementation -- #
HolmBonferroni <- function(values, alpha, FUN, ...){
  arguments <- list(values, ...)
  p_vals <- 1-do.call(FUN,arguments)
  p_vals_sort <- sort(p_vals,index.return=TRUE)
  k <- suppressWarnings(min(which(p_vals_sort[[1]] > alpha/(length(values)+1-1:length(values)))))
  if(!is.finite(k)) return(p_vals_sort[[2]]) else {if(k > 1) return(p_vals_sort[[2]][1:(k-1)]) else return(NULL)}
}

# -- AUC computation -- #
pauc_fcn <- function(controls, cases, partial = TRUE, plot = NULL){
  
  Pop_bind <- rbind(cbind(controls,0),cbind(cases,1))
  nv <- dim(Pop_bind)[2]-1
  n0 <- match(1,Pop_bind[,nv+1])-1
  n1 <- dim(Pop_bind)[1]-n0
  
  if(nv==1){ Pop_bind_ind <- sort(sort(Pop_bind[,1],index.return=TRUE, method="radix")[[2]],index.return=TRUE, method="radix")[[2]]
  } else stop("As implemented, pauc_fcn is limited to univariate data")
  
  Pop_frame <- array(NA_integer_,dim=rep(dim(Pop_bind)[1],nv))
    
  Pop_frame[Pop_bind_ind] <- as.integer(Pop_bind[,dim(Pop_bind)[2]])
    
  sens <- c(rev_cumsum(Pop_frame),0)
  spec <- c(rev_cumsum(1L-Pop_frame),0)
  
  if(!is.null(plot)) plot(spec/n0,sens/n1, type="S", xlab="1-Specificity", ylab="Sensitivity", xlim=c(0,1),ylim=c(0,1))
  
  if(anyDuplicated(spec) > 0){
    sens <- sens[-which(duplicated(spec))]
  } 
  
  sens <- sens[-1]
  
  if(isTRUE(partial)){
    fl <- floor(length(sens)*.8)
    return((sum(sens[(1+fl):length(sens)])-sens[1+fl]*(length(sens)*.8-fl))/(n0*n1))
  } else {
    return(sum(sens)/(n0*n1))
  }
}

## /// Utility functions \\\ ##

## --- Simulation procedure --- ##

# -- Object generation -- #

sample_sizes <- (1:5)*50          # As implemented, same number of cases and controls
object_size_number <- 1000        # Ten million is a feasible number for many memory systems
partial <- TRUE                   # If TRUE, simulates partial AUC from 0 to 0.2
mixing_coefficients <- c(0.8,0.2) # Mix of non-qualitative biomarkers, e.g. 80% BC-phi, 20% NMP22

BC_sim_sort         <- new("list")
NMP_sim_sort        <- new("list")
FIT_sim_sort        <- new("list")
ColoGuard_sim_sort  <- new("list")
NMP_mixedCDF        <- new("list")
FIT_mixedCDF        <- new("list")
ColoGuard_mixedCDF  <- new("list")

for(i in 1:length(sample_sizes)){

  BC_sim_sort[[i]]        <- sort(sapply(rep(NA,object_size_number),function(m){sim_fcn(sample_sizes[i],sample_sizes[i],0.772,partial)}),decreasing=TRUE)
  NMP_sim_sort[[i]]       <- sort(sapply(rep(NA,object_size_number),function(m){sim_fcn(sample_sizes[i],sample_sizes[i],1.140,partial)}),decreasing=TRUE)
  FIT_sim_sort[[i]]       <- sort(sapply(rep(NA,object_size_number),function(m){sim_fcn(sample_sizes[i],sample_sizes[i],1.658,partial)}),decreasing=TRUE)
  ColoGuard_sim_sort[[i]] <- sort(sapply(rep(NA,object_size_number),function(m){sim_fcn(sample_sizes[i],sample_sizes[i],2.083,partial)}),decreasing=TRUE)

}

mixedCDF <- function(x,p,n){
  
  # x: (AUC-) value, p: mixing coefficients, n: sample_size (as implemented, same number of cases and controls)
  
  i <- match(n,sample_sizes)
  if(sum(p)!= 1){ p <- p/sum(p); warning("normalizing mixing coefficients") }
  
  p[1]*(1-sapply(x,function(z){ (match(TRUE,BC_sim_sort[[i]]<=z,nomatch=(length(BC_sim_sort[[i]])+1))-1)/length(BC_sim_sort[[i]])})) +
    p[2]*(1-sapply(x,function(z){ (match(TRUE,NMP_sim_sort[[i]]<=z,nomatch=(length(NMP_sim_sort[[i]])+1))-1)/length(NMP_sim_sort[[i]])}))
  
}

for(i in 1:length(sample_sizes)){
  NMP_mixedCDF[[i]]       <- mixedCDF(NMP_sim_sort[[i]],      p=mixing_coefficients,sample_sizes[i])
  FIT_mixedCDF[[i]]       <- mixedCDF(FIT_sim_sort[[i]],      p=mixing_coefficients,sample_sizes[i])
  ColoGuard_mixedCDF[[i]] <- mixedCDF(ColoGuard_sim_sort[[i]],p=mixing_coefficients,sample_sizes[i])
}

# -- Simulation utilities -- #

probWin <- function(size,p,n,k=1,m=167668500){
  
  # size: no. random numbers used, p: mixing coefficients, n: sample size, k: placement from top, m: no. competing entrants
  
  if(!identical(p,mixing_coefficients)) stop("as implemented, must use pre set mixing coefficients")
  
  i <- match(n,sample_sizes)
  apply(sapply(rbeta(size,m+1-k,k),function(x){c(mean(NMP_mixedCDF[[i]]>x),mean(FIT_mixedCDF[[i]]>x),mean(ColoGuard_mixedCDF[[i]]>x))}),MARGIN=1,FUN=mean)
  
}

trainingRoutine <- function(size,p,qualitative_biomarkers,n,no_advanced=1,m=167668500){
  
  # size: no. random numbers used, p: mixing coefficients, qualitative_biomarkers: no putative qualitative biomarkers (NMP,FIT,ColoGuard),
  # n: sample size, no_advanced: no. of biomarkers to advance, m: total no. of biomarkers entering training step
  
  if(max(qualitative_biomarkers > no_advanced)){ stop("more qualitative biomarkers than no. biomarkers advanced")}
  
  simple_probs <- array(0,dim=c(0,3))
  
  for(k in no_advanced+1-(1:max(qualitative_biomarkers))){
    
    simple_probs <- rbind(simple_probs,probWin(size,p,n,k,m))
    
  }
  
  composite_probs <- array(0,dim=c(0,3))
  
  for(j in 1:max(qualitative_biomarkers)){
    composite_probs <- rbind(composite_probs,dbinom(j,qualitative_biomarkers,simple_probs[j,]))
    
  }
  return(composite_probs)
}

validationRoutine_NMP <- function(composition, alpha, n, replicates=1){
  
  # composition: biomarkers entering the validation step (NMP22, FIT, ColoGuard), alpha: desired type-one error probability, n: sample size,
  # replicates: no. replicated validation steps
  
  i <- match(n,sample_sizes)
  
  internal_fcn <- function(w){
    
    NMP_pvals <- runif(composition[1])
    qual_pvals <- c(FIT_mixedCDF[[i]][ceiling(runif(composition[2])*length(FIT_mixedCDF[[i]]))],ColoGuard_mixedCDF[[i]][ceiling(runif(composition[3])*length(ColoGuard_mixedCDF[[i]]))])
    if(sum(composition[2:3])>0) NMP_pvals[which(NMP_pvals < min(qual_pvals))] <- 0
    p_vals <- c(NMP_pvals,qual_pvals)
    indices <- HolmBonferroni(p_vals,alpha,identity)
    
    sapply(1:3,function(i){sum((indices<=cumsum(composition)[i]) * (indices>c(0,cumsum(composition))[i]))})
  }
  Reduce(function(a,b){a+(unlist(b)>0)},c(list(c(0,0,0)),lapply(rep(NA,replicates),internal_fcn)))/replicates
}

simRoutine <- function(size,p=mixing_coefficients,qualitative_biomarkers=c(1,1,1),n_t=sample_sizes[1], n_v=sample_sizes[1] ,no_advanced=1,m=167668500,alpha=0.01){
  
  # size: no. replicates used, p: mixing coefficients, qualitative_biomarkers: no putative qualitative biomarkers (NMP,FIT,ColoGuard),
  # n_t: sample size training step, n_v: sample size validation step, no_advanced: no. of biomarkers to advance to validation,
  # m: total no. of biomarkers entering training step, alpha: desired type-one error probability
  
  t_probs <- new("list")
  sign_probs <- new("list")
  power_probs <- new("list")
  
  for(u in 1:length(n_t)){
    t_probs[[u]] <- trainingRoutine(size,p,qualitative_biomarkers,n_t[u],no_advanced,m)
  }
  
  for(v in 1:length(n_v)){
    sign_probs[[v]] <- array(0,dim=c(0,3))
    for(i in 1:max(qualitative_biomarkers)){
      sign_probs[[v]] <- rbind(sign_probs[[v]],c(
        
        if(qualitative_biomarkers[1] >= i){validationRoutine_NMP(c(no_advanced,0,0), alpha, n_v[v], size)[1]}else{0},
        
        if(qualitative_biomarkers[2] >= i){validationRoutine_NMP(c(no_advanced-i,i,0), alpha, n_v[v], size)[2]}else{0},
        
        if(qualitative_biomarkers[3] >= i){validationRoutine_NMP(c(no_advanced-i,0,i), alpha, n_v[v], size)[3]}else{0}))
    }
  }
  
  for(u in 1:length(n_t)){
    for(v in 1:length(n_v)){
      power_probs <- c(power_probs,list(colSums(t_probs[[u]]*sign_probs[[v]])))
      attributes(power_probs[[length(power_probs)]]) <- list(n_info = paste0("n_t=",n_t[u],", n_v=",n_v[v]))
    }
  }
  
  return(power_probs)
  
}

# -- test -- #
simRoutine(100)

## /// Simulation procedure \\\ ##
