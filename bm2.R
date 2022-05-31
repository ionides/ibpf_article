#' Brownian motion spatPomp generator with shared or unit-specific parameters
#'
#' Generate a class \sQuote{spatPomp} object representing a \code{U}-dimensional
#' Brownian motion with spatial correlation decaying geometrically with
#' distance around a circle. The model is defined in continuous time
#' though in this case an Euler approximation is exact at the evaluation
#' times.
#'
#' @param U A length-one numeric signifying dimension of the process.
#' @param N A length-one numeric signifying the number of observation time steps to evolve the process.
#' @param delta_t Process simulations are performed every \code{delta_t} time units
#' whereas observations occur every one time unit
#' @param shared_names identifies parameters that have common shared value for all units, which by default is all parameters.
#' @param unit_specific_names determines which parameters take a different value for each unit. Cannot be specified if shared_names is specified.
#' each unit. Other parameters are considered shared between all units.
#' @importFrom utils data
#' @return An object of class \sQuote{spatPomp} representing a simulation from a \code{U}-dimensional
#' Brownian motion
#' @examples
#' b <- bm2(U=4, N=20)
#' # See all the model specifications of the object
#' spy(b)
#' @export

bm2 <- function(U=5,N=100,delta_t=0.1,
  shared_names,unit_specific_names){
  dist <- function(u,v,n=U) min(abs(u-v),abs(u-v+U),abs(u-v-U))
  dmat <- matrix(0,U,U)
  for(u in 1:U) {
    for(v in 1:U) {
      dmat[u,v] <- dist(u,v)
    }
  }
  to_C_array <- function(v)paste0("{",paste0(v,collapse=","),"}")
  dist_C_rows <- apply(dmat,1,to_C_array)
  dist_C_array <- to_C_array(dist_C_rows)
  dist_C <- paste0("const int dist[",U,"][",U,"] = ",dist_C_array,"; ")

  obs_names <- paste0("U",1:U)
  bm_data <- data.frame(time=rep(1:N,U),unit=rep(obs_names,each=N),Y=rep(NA,U*N),stringsAsFactors=F)
  bm_unitnames <- unique(bm_data[["unit"]])
  bm_unitnames_level <- paste("U",sort(as.numeric(stringr::str_remove(bm_unitnames, "U"))),sep='')

# bm is written using spatPomp_Csnippet whereas measles is written via Csnippet
# spatPomp_Csnippet adds an assumption that initial states have a name of the
# form Xu_O whereas measles2 uses X_0u to be consistent with other unit-specific
# parameters. This should be reconciled. Meanwhile, we follow the measles2
# convention, and do not use the spatPomp_Csnippet capability for IVPs

  bm_unit_statenames <- c("X")
  bm_statenames <- paste0(bm_unit_statenames,1:U)

  bm_unit_IVPnames <- paste0(bm_unit_statenames,"_0")
  bm_unit_RPnames <- c("rho","sigma","tau")
  bm_unit_paramnames <- c(bm_unit_RPnames,bm_unit_IVPnames)


  if(!missing(shared_names)) {
    if(!missing(unit_specific_names)) {
      stop ("both shared_names and unit_specific names cannot be given to bm2")
    } else
      unit_specific_names <- bm_unit_paramnames[!bm_unit_paramnames%in%shared_names]
  }
  
  if(missing(shared_names)) {
    if(missing(unit_specific_names)) {
      shared_names <- bm_unit_paramnames
      unit_specific_names <- NULL
    } else shared_names <- bm_unit_paramnames[!bm_unit_paramnames%in%unit_specific_names]
  }

  set_unit_specific <- Csnippet(paste0("const int ", unit_specific_names,
    "_unit = 1;\n", collapse=" "))
  set_shared <- Csnippet(paste0("const int ", shared_names,
    "_unit = 0;\n", collapse=" "))

  bm_globals <- Csnippet(
    paste(dist_C,set_unit_specific, set_shared, sep = "\n")
  )

  # add a "1" for shared parameter names to make the pointers work
  bm_paramnames <- c(
    if(length(shared_names)>0){
      paste0(shared_names, "1")
    },
    if(length(unit_specific_names)>0){
      paste0(rep(unit_specific_names, each=U), 1:U)
    }
  )

  bm_rprocess <- spatPomp_Csnippet("
    const double *rho=&rho1;
    const double *sigma=&sigma1;
    double dW[U];
    double pow_rho[U];
    int u,v;

    pow_rho[0] = 1;
    for (u=1 ; u < U ; u++) {
      pow_rho[u] = pow_rho[u-1]*rho[u*rho_unit];
    }

    for (u = 0 ; u < U ; u++) {
      dW[u] = rnorm(0,sigma[u*sigma_unit]*sqrt(dt));
    }
    for (u = 0 ; u < U ; u++) {
      for (v=0; v < U ; v++) {
        X[u] += dW[v]*pow_rho[dist[u][v]];
      }
    }
  ", unit_statenames = c("X"))

  bm_skel <- spatPomp_Csnippet(
    unit_statenames = c("X"),
    unit_vfnames = c("X"),
    code = "
      for (int u = 0 ; u < U ; u++) {
        DX[u] = 0;
      }
    "
  )


  bm_rinit <- spatPomp_Csnippet(
    unit_statenames = c("X"),
    code = "
      const double *X_0 = &X_01;
      for (int u = 0; u < U; u++) {
        X[u]=X_0[u];
      }
    "
  )

  bm_dmeasure <- Csnippet("
    const double *tau = &tau1;
    const double *X = &X1;
    const double *Y = &Y1;
    double tol = pow(1.0e-18,U);
    int u;
    lik=0;
    for (u=0; u<U; u++) lik += dnorm(Y[u],X[u],tau[u*tau_unit],1);
    if(!give_log) lik = exp(lik) + tol;
  ")

  bm_eunit_measure <- Csnippet("
    ey = X;
  ")


## this munit_measure and vunit_measure make sense only if tau is shared
## should be removed or fixed

  bm_munit_measure <- Csnippet("
    M_tau1 = sqrt(vc);
  ")

  bm_vunit_measure <- Csnippet("
    vc = tau1*tau1;
  ")

  bm_rmeasure <- Csnippet("
    const double *tau = &tau1;
    const double *X = &X1;
    double *Y = &Y1;
    double tol = pow(1.0e-18,U);
    int u;
    for (u=0; u<U; u++) Y[u] = rnorm(X[u],tau[u*tau_unit]+tol);
  ")

  bm_dunit_measure <- Csnippet("
    const double *tau = &tau1;
    //double tol = 1.0e-18;
    lik = dnorm(Y,X,tau[u*tau_unit],1);
    if(!give_log) lik = exp(lik);
  ")

  bm_runit_measure <- Csnippet("
    const double *tau = &tau1;
    double tol = pow(1.0e-18,U);
    double Y;
    Y = rnorm(X,tau[u*tau_unit]+tol);
  ")

log_unit_names <- c("sigma", "tau")
logit_unit_names <- c("rho")
log_names <- unlist(lapply(log_unit_names,
  function(x,y,U){if(x%in%y)paste0(x,"1") else paste0(x,1:U)},
  y=shared_names,U=U))
logit_names <- unlist(lapply(logit_unit_names,
  function(x,y,U){if(x%in%y)paste0(x,"1") else paste0(x,1:U)},
  y=shared_names,U=U))

bm_partrans <- parameter_trans(log=log_names,logit=logit_names)


  bm_spatPomp <- spatPomp(bm_data %>% dplyr::arrange(time, factor(.data$unit, levels = bm_unitnames_level)),
                 times="time",
                 t0=0,
                 units="unit",
                 unit_statenames = bm_unit_statenames,
                 rprocess=euler(bm_rprocess,delta.t = delta_t),
                 skeleton=vectorfield(bm_skel),
                 paramnames=bm_paramnames,
                 globals=bm_globals,
                 rmeasure=bm_rmeasure,
                 dmeasure=bm_dmeasure,
                 eunit_measure=bm_eunit_measure,
                 munit_measure=bm_munit_measure,
                 vunit_measure=bm_vunit_measure,
                 dunit_measure=bm_dunit_measure,
                 runit_measure=bm_runit_measure,
                 partrans = bm_partrans,
                 rinit=bm_rinit
    )

  ## We need a parameter vector.
  bm_params <- rep(0,length=length(bm_paramnames))
  names(bm_params) <- bm_paramnames
  # parameter match the default for previous spatPomp bm() 
  bm_unit_params <- c(
    rho=0.4,
    sigma=1,  
    tau=1,
    X_0=0)
  for(p in unit_specific_names) {
    bm_params[paste0(p,1:U)] <- bm_unit_params[p]
  }
  for(p in shared_names) {
    bm_params[paste0(p,1)] <- bm_unit_params[p]
  }
  coef(bm_spatPomp) <- bm_params

  sim_bm <- simulate(bm_spatPomp,params=bm_params)
}
