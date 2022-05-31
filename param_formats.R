
# Book-keeping functions to expand from a vector of shared parameters to
# unit-specific copies, and to reverse this operation.

# Only estimated shared parameters are expanded. Fixed shared parameters,
# and unit-specific parameters that already have a copy for each unit, are
# left unchanged.

# The contracted form has a single entry for each shared parameter, and only
# unit-specific parameters have a terminal number

# The expanded form has U copies with appended 1:U for each shared or
# unit-specific parameter. Fixed parameters (i.e., parameters that are not
# estimated) have an appended 1. 

# FixedParNames,sharedParNames and unitParNames are the names of the
# corresponding parameters without any appended numbers.

# Another form of the parameters, known as the basic form, is a reduced version
# of the contracted form having just a single value for each parameter. It
# corresponds to the parameter vector for a spatPomp model with all shared
# parameters. expand_params() can also be used to generate an expanded form
# from the basic form, if both shared and unit-specific parameters are passed
# via the sharedParNames argument. From there, one can obtain the contracted
# form.

expand_params <- function(contracted_params,
  fixedParNames,sharedParNames,unitParNames,U){
  expand_shared_params <- foreach(par=sharedParNames,.combine=c)%do%{
    x <- rep(contracted_params[par],U)
    names(x) <- paste0(par,1:U)
    x
    }
  expand_fixed_params <- contracted_params[fixedParNames]
  names(expand_fixed_params) <- paste0(fixedParNames,"1")
  expand_unit_params <- foreach(par=unitParNames,.combine=c)%do%{
    contracted_params[paste0(par,1:U)]
  }
  c(expand_fixed_params,expand_shared_params,expand_unit_params)
}

contract_params <- function(params,
  fixedParNames,sharedParNames,unitParNames,
  average=FALSE,U){
  expandedParNames <- names(params)
  contract_fixed_params <- params[paste0(fixedParNames,"1")]
  names(contract_fixed_params) <- fixedParNames
  contract_unit_params <- foreach(par=unitParNames,.combine=c)%do%{
    params[paste0(par,1:U)]
  }
  contract_shared_params <- foreach(par=sharedParNames,.combine=c)%do%{
    x <- params[paste0(par,1:U)]
    if(sd(x)>0 & !average) stop("parameters must be the same for each unit in contract_params with averag=FALSE")
    x <- mean(x)
    names(x) <- par
    x
  }
  c(contract_fixed_params,contract_shared_params,contract_unit_params)
}


mean_by_unit <- function(params,sharedParNames,U){
  for(par in sharedParNames){
    params[paste0(par,1:U)] <- mean(params[paste0(par,1:U)])
  }
  params
}

check <- FALSE
if(check){
  contracted <- c(a=0.1,b=0.2,c=0.3,d=0.4,e1=0.3,e2=0.4,e3=0.5,e4=0.6,e5=0.7)
  fixedPars <- c("a","b")
  sharedPars <- c("c","d")
  unitPars <- "e"
  U <- 5

  expanded <- expand_params(contracted,U=U,
    fixedParNames=fixedPars,sharedParNames=sharedPars,unitParNames=unitPars)

  contracted2 <- contract_params(expanded,average=FALSE,
    fixedParNames=fixedPars,sharedParNames=sharedPars,unitParNames=unitPars)

  all(contracted==contracted2)

#  full_names <- c("a",paste0("b",1:12),paste0("R0",1:12))
#  gsub("[1-9]+[0-9]*$","",full_names)
  
}