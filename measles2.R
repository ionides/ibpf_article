#' Measles in UK: spatPomp generator with shared or unit-specific parameters
#'
#' Generate a spatPomp object for measles in the top-\code{U} most populous cities in England and Wales.
#' The model is adapted from He et al. (2010) with gravity transport following Park and Ionides (2019).
#' The data in the object is simulated using the process and measurement models of He et al. (2010).
#' The structure of this spatPomp is designed to be compatible with
#' spatiotemporal iterated filtering algorithms. Thus, shared parameters can
#' be represented with a value for each unit, which should be the same for each
#' unit in a valid model instance but may vary between units while optimizing.
#'
#' @param U A length-one numeric signifying the number of cities to be represented in the spatPomp object.
#' @importFrom utils data read.csv write.table
#' @param dt a numeric (in unit of years) that is used as the Euler time-increment for simulating measles data.
#' @param fixedParNames specifies parameters that have a shared value across
#' units and are not being optimized so can be represented by a single value.
#' @param sharedParNames identifies parameters that have common shared value
#' for all units which is represented separately for each unit for optimization
#' purposes.
#' @param unitParNames determines which parameters take a different value for each unit. Cannot be specified if sharedParNames is specified.
#' each unit. Other parameters are considered shared between all units.
#' @param simulated determines whether to return a simulation from the model or the
#' UK measles data
#' @return An object of class \sQuote{spatPomp} representing a \code{U}-dimensional spatially coupled measles POMP model.
#' @references
#'
#' \geosphere
#'
#' @note This function goes through a typical workflow of constructing
#' a typical spatPomp object (1-4 below). This allows the user to have a
#' file that replicates the exercise of model building as well as function
#' that creates a typical nonlinear model in epidemiology in case they want
#' to test a new inference methodology. We purposely do not modularize this
#' function because it is not an operational piece of the package and is
#' instead useful as an example.\cr
#' 1. Getting a measurements data.frame with columns for times,
#'    spatial units and measurements.\cr
#' 2. Getting a covariates data.frame with columns for times,
#'    spatial units and covariate data.\cr
#' 3. Constructing model components (latent state initializer,
#'    latent state transition simulator and measurement model). Depending
#'    on the methods used, the user may have to supply a vectorfield to
#'    be integrated that represents the deterministic skeleton of the latent
#'    process.\cr
#' 4. Bringing all the data and model components together to form a
#'    spatPomp object via a call to spatPomp().
#' @examples
#' m <- measles2(U = 5)
#' # See all the model specifications of the object
#' spy(m)
#' @export

# NOTE: As indicated in the Note section of the documentation, this
# this function goes through a typical workflow of constructing
# a spatPomp object. It is not meant to be operational, but
# instead an example of how one goes about going from getting data to creating
# a spatPomp object.

# NOTE: Code was written assuming that there are no fixed unit-specific parameters. What
# happens in that case?

measles2 <- function(U=6,dt=2/365,N=391,
  fixedParNames, sharedParNames, unitParNames, simulated=FALSE,
  basic_params =c(
    alpha = 1,
    iota = 0,  
    R0 = 30,
    cohort = 0,
    amplitude = 0.5,
    gamma = 52,
    sigma = 52,
    mu = 0.02,
    sigmaSE = 0.15, 
    rho = 0.5,
    psi = 0.15,
    g = 400,
    S_0 = 0.032, 
    E_0 = 0.00005, 
    I_0 = 0.00004
  )
){

# default basic parameter values based on Ionides et al (2021)
# "Bagged filters for partially observed interacting systems. 

  if(missing(sharedParNames)) sharedParNames <- NULL
  if(missing(unitParNames)) unitParNames <- NULL
  expandedParNames = c(sharedParNames,unitParNames)
  if(length(intersect(sharedParNames,unitParNames))>0)
    stop("sharedParNames and unitParNames should be disjoint")
  if(U>40) stop("U <= 40")
  if(N>391) stop("N <= 391")
  birth_lag <- 3*26  # delay until births hit susceptibles, in biweeks
  # 4yr was the estimate by He et al (2010)

  # pre-vaccine biweekly measles reports for the largest 40 UK cities, sorted by size
  measlesUK <- spatPomp::measlesUK
  measlesUK$city<-as.character(measlesUK$city)
  city_data_UK <- spatPomp::city_data_UK

  cities <- unique(measlesUK$city)[1:U]
  selected <- (measlesUK$city %in% cities) & (measlesUK$year>1949.99) &
    (measlesUK$year<1950.01+(N-1)/26)    
  measles_cases <- measlesUK[selected,c("year","city","cases")]
  covar_selected <- (measlesUK$city %in% cities) & (measlesUK$year<1950.01+(N-1)/26) 
  measles_covar <- measlesUK[covar_selected,c("year","city","pop","births")]
  u <- split(measles_covar$births,measles_covar$city)
  v <- sapply(u,function(x){c(rep(NA,birth_lag),x[1:(length(x)-birth_lag)])})
  measles_covar$lag_birthrate <- as.vector(v[,cities])*26
  measles_covar$births<- NULL
  measles_covarnames <- paste0(rep(c("pop","lag_birthrate"),each=U),1:U)
  measles_unit_covarnames <- c("pop","lag_birthrate")

  # Distance between two points on a sphere radius R
  # Adapted from geosphere package, which has been cited in the package
  distHaversine <- function (p1, p2, r = 6378137)
  {
      toRad <- pi/180
      p1 <- p1 * toRad
      p2 <- p2 * toRad
      p = cbind(p1[, 1], p1[, 2], p2[, 1], p2[, 2], as.vector(r))
      dLat <- p[, 4] - p[, 2]
      dLon <- p[, 3] - p[, 1]
      a <- sin(dLat/2) * sin(dLat/2) + cos(p[, 2]) * cos(p[, 4]) *
          sin(dLon/2) * sin(dLon/2)
      a <- pmin(a, 1)
      dist <- 2 * atan2(sqrt(a), sqrt(1 - a)) * p[, 5]
      return(as.vector(dist))
  }

  long_lat <- city_data_UK[1:U,c("lon","lat")]
  dmat <- matrix(0,U,U)
  for(u1 in 1:U) {
    for(u2 in 1:U) {
      dmat[u1,u2] <- round(distHaversine(long_lat[u1,],long_lat[u2,]) / 1609.344,1)
    }
  }

  p <- city_data_UK$meanPop[1:U]
  v_by_g <- matrix(0,U,U)
  dist_mean <- sum(dmat)/(U*(U-1))
  p_mean <- mean(p)
  for(u1 in 2:U){
    for(u2 in 1:(u1-1)){
      v_by_g[u1,u2] <- (dist_mean*p[u1]*p[u2]) / (dmat[u1,u2] * p_mean^2)
      v_by_g[u2,u1] <- v_by_g[u1,u2]
    }
  }
  to_C_array <- function(v)paste0("{",paste0(v,collapse=","),"}")
  v_by_g_C_rows <- apply(v_by_g,1,to_C_array)
  v_by_g_C_array <- to_C_array(v_by_g_C_rows)
  v_by_g_C <- Csnippet(paste0("const double v_by_g[",U,"][",U,"] = ",v_by_g_C_array,"; "))

  # here, basic parameters are as described in param_formats.R
  basic_statenames <- c('S','E','I','R','C')
  basic_RPnames <- c("alpha","iota","psi","R0","gamma","sigma","sigmaSE","cohort","amplitude","mu","rho","g")
  basic_IVPnames <- c("S_0", "E_0", "I_0")
  basicParNames <- c(basic_RPnames,basic_IVPnames)

  expandedParNames <- setdiff(basicParNames,fixedParNames)

  set_unit_specific <- Csnippet(paste0("const int ", expandedParNames,
    "_unit = 1;\n", collapse=" "))
  set_shared <- Csnippet(paste0("const int ", fixedParNames,
    "_unit = 0;\n", collapse=" "))

  measles_globals <- Csnippet(
    paste(v_by_g_C, set_unit_specific, set_shared, sep = "\n")
  )

  # add a "1" for shared parameter names to make the pointers work
  measles_paramnames <- c(
    if(length(fixedParNames)>0){
      paste0(fixedParNames, "1")
    },
    if(length(expandedParNames)>0){
      paste0(rep(expandedParNames, each=U), 1:U)
    }
  )

  measles_statenames <- paste0(rep(basic_statenames,each=U),1:U)

  measles_rprocess <- Csnippet('
    const double *amplitude=&amplitude1;
    const double *sigmaSE=&sigmaSE1;
    const double *mu=&mu1;
    const double *g=&g1;
    const double *cohort=&cohort1;
    const double *gamma=&gamma1;
    const double *sigma=&sigma1;
    const double *R0=&R01;
    const double *alpha=&alpha1;
    const double *iota=&iota1;
    double br, beta, seas, foi, dw, births;
    double rate[6], trans[6];
    double *S = &S1;
    double *E = &E1;
    double *I = &I1;
    double *R = &R1;
    double *C = &C1;
    //double *W = &W1;
    double powVec[U];
    //double *Acc = &Acc1;
    const double *pop = &pop1;
    const double *lag_birthrate = &lag_birthrate1;
    int u,v;
    // term-time seasonality
    t = (t-floor(t))*365.25;


    for (u = 0 ; u < U ; u++) {


      // needed for the Ensemble Kalman filter
      // or other methods making real-valued perturbations to the state
      // reulermultinom requires integer-valued double type for states
      S[u] = S[u]>0 ? floor(S[u]) : 0;
      E[u] = E[u]>0 ? floor(E[u]) : 0;
      I[u] = I[u]>0 ? floor(I[u]) : 0;
      R[u] = R[u]>0 ? floor(R[u]) : 0;

      // pre-computing this saves substantial time
      //powVec[u] = pow(I[u]/pop[u],alpha);
      powVec[u] = I[u]/pop[u];
      // IS THIS INTENDED TO BE FIXED TO ALPHA=1?
    }


    for (u = 0 ; u < U ; u++) {

       if ((t>=7&&t<=100) || (t>=115&&t<=199) || (t>=252&&t<=300) || (t>=308&&t<=356))
        seas = 1.0+amplitude[u*amplitude_unit]*0.2411/0.7589;
      else
        seas = 1.0-amplitude[u*amplitude_unit];

     // transmission rate
      beta = R0[u*R0_unit]*(gamma[u*gamma_unit]+mu[u*mu_unit])*seas;

      rate[1] = mu[u*mu_unit];		// natural S death
      rate[2] = sigma[u*sigma_unit];	// rate of ending of latent stage
      rate[3] = mu[u*mu_unit];		// natural E death
      rate[4] = gamma[u*mu_unit];	// recovery
      rate[5] = mu[u*mu_unit];		// natural I death

      // cohort effect
      if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt)
        br = cohort[u*cohort_unit]*lag_birthrate[u]/dt + (1-cohort[u*cohort_unit])*lag_birthrate[u];
      else
        br = (1.0-cohort[u*cohort_unit])*lag_birthrate[u];

      // expected force of infection
      if(alpha[u*alpha_unit]==1.0 && iota[u*iota_unit]==0.0)
        foi = I[u]/pop[u];
      else
        foi = pow( (I[u]+iota[u*iota_unit])/pop[u],alpha[u*alpha_unit]);
      // we follow Park and Ionides (2019) and raise pop to the alpha power
      // He et al (2010) did not do this.

      for (v=0; v < U ; v++) {
        if(v != u)
          foi += g[u*g_unit] * v_by_g[u][v] * (powVec[v] - powVec[u]) / pop[u];
      }

      // white noise (extrademographic stochasticity)
      dw = rgammawn(sigmaSE[u*sigmaSE_unit],dt);
      rate[0] = beta*foi*dw/dt;  // stochastic force of infection

      // Poisson births
      births = rpois(br*dt);

      // transitions between classes
      reulermultinom(2,S[u],&rate[0],dt,&trans[0]);
      reulermultinom(2,E[u],&rate[2],dt,&trans[2]);
      reulermultinom(2,I[u],&rate[4],dt,&trans[4]);

      S[u] += births   - trans[0] - trans[1];
      E[u] += trans[0] - trans[2] - trans[3];
      I[u] += trans[2] - trans[4] - trans[5];
      R[u] = pop[u] - S[u] - E[u] - I[u];
      C[u] += trans[4];           // true incidence
     }
  ')

  measles_dmeasure <- Csnippet("
    const double *C = &C1;
    const double *cases = &cases1;
    const double *rho=&rho1;
    const double *psi=&psi1;
    double m,v;
    double tol = 1e-300;
    double mytol = 1e-5;
    int u;

    lik = 0;
    for (u = 0; u < U; u++) {
      m = rho[u*rho_unit]*(C[u]+mytol);
      v = m*(1.0-rho[u*rho_unit]+psi[u*psi_unit]*psi[u*psi_unit]*m);
      // C < 0 can happen in bootstrap methods such as bootgirf
      if (C < 0) {lik += log(tol);} else {
        if (cases[u] > tol) {
          lik += log(pnorm(cases[u]+0.5,m,sqrt(v)+tol,1,0)-
            pnorm(cases[u]-0.5,m,sqrt(v)+tol,1,0)+tol);
        } else {
            lik += log(pnorm(cases[u]+0.5,m,sqrt(v)+tol,1,0)+tol);
        }
      }
    }
    if(!give_log) lik = (lik > log(tol)) ? exp(lik) : tol;
  ")

  measles_rmeasure <- Csnippet("
    const double *C = &C1;
    double *cases = &cases1;
    const double *rho = &rho1;
    const double *psi = &psi1;
    double m,v;
    double tol = 1.0e-300;
    int u;

    for (u = 0; u < U; u++) {
      m = rho[u*rho_unit]*(C[u]+tol);
      v = m*(1.0-rho[u*rho_unit]+psi[u*psi_unit]*psi[u*psi_unit]*m);
      cases[u] = rnorm(m,sqrt(v)+tol);
      if (cases[u] > 0.0) {
        cases[u] = nearbyint(cases[u]);
      } else {
        cases[u] = 0.0;
      }
    }
  ")

  measles_dunit_measure <- Csnippet('
    const double *rho = &rho1;
    const double *psi = &psi1;
    double mytol = 1e-5;
    double m = rho[(u-1)*rho_unit]*(C+mytol);
    double v = m*(1.0-rho[(u-1)*rho_unit]+psi[(u-1)*psi_unit]*psi[(u-1)*psi_unit]*m);
    double tol = 1e-300;
    // C < 0 can happen in bootstrap methods such as bootgirf
    if (C < 0) {lik = 0;} else {
      if (cases > tol) {
        lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-
          pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
      } else {
        lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
      }
    }
    if(give_log) lik = log(lik);
  ')

  measles_eunit_measure <- Csnippet("
    const double *rho = &rho1;
    ey = rho[u*rho_unit]*C;
  ")

  measles_vunit_measure <- Csnippet("
    //consider adding 1 to the variance for the case C = 0
    const double *rho = &rho1;
    const double *psi = &psi1;
    double mytol = 1e-5;
    double m;
    m = rho[u*rho_unit]*(C+mytol);
    vc = m*(1.0-rho[u*rho_unit]+psi[u*psi_unit]*psi[u*psi_unit]*m);
  ")

  measles_munit_measure <- Csnippet("
    const double *rho = &rho1;
    double binomial_var;
    double m;
    double mytol = 1e-5;
    m = rho[u*rho_unit]*(C+mytol);
    binomial_var = rho[u*rho_unit]*(1-rho[u*rho_unit])*C;
    if(vc > binomial_var) {
      M_psi1 = sqrt(vc - binomial_var)/m;
    }
  ")

  measles_rinit <- Csnippet("
    double *S = &S1;
    double *E = &E1;
    double *I = &I1;
    double *R = &R1;
    double *C = &C1;
    const double *S_0 = &S_01;
    const double *E_0 = &E_01;
    const double *I_0 = &I_01;
    double m;
    const double *pop = &pop1;
    int u;
    for(u = 0; u < U; u++) {
      m = (float)(pop[u]);
      S[u] = nearbyint(m*S_0[u*S_0_unit]);
      I[u] = nearbyint(m*I_0[u*I_0_unit]);
      E[u] = nearbyint(m*E_0[u*E_0_unit]);
      R[u] = pop[u]-S[u]-E[u]-I[u];
      C[u] = 0;
      // in any practical model fit, we expect R>0
      // though the model does not strictly enforce that 
    }
  ")

  measles_skel <- Csnippet('
    double beta, br, seas, foi;
    const double *amplitude=&amplitude1;
    const double *mu=&mu1;
    const double *g=&g1;
    const double *cohort=&cohort1;
    const double *gamma=&gamma1;
    const double *sigma=&sigma1;
    const double *R0=&R01;
    const double *alpha=&alpha1;
    const double *iota=&iota1;
    double *S = &S1;
    double *E = &E1;
    double *I = &I1;
    double *R = &R1;
    double *C = &C1;
    double *DS = &DS1;
    double *DE = &DE1;
    double *DI = &DI1;
    double *DR = &DR1;
    double *DC = &DC1;
    double powVec[U];
    const double *pop = &pop1;
    const double *lag_birthrate = &lag_birthrate1;
    int u,v;

    // pre-computing this saves substantial time
    for (u = 0 ; u < U ; u++) {
      powVec[u] = pow(I[u]/pop[u],alpha[u*alpha_unit]);
    }

    for (u = 0 ; u < U ; u++) {

      // term-time seasonality
      t = (t-floor(t))*365.25;
      if ((t>=7&&t<=100) || (t>=115&&t<=199) || (t>=252&&t<=300) || (t>=308&&t<=356))
        seas = 1.0+amplitude[u*amplitude_unit]*0.2411/0.7589;
      else
        seas = 1.0-amplitude[u*amplitude_unit];

    // transmission rate
    beta = R0[u*R0_unit]*(gamma[u*gamma_unit]+mu[u*amplitude_unit])*seas;

      // cannot readily put the cohort effect into a vectorfield for the skeleton
      // therefore, we ignore it here.
      // this is okay as long as the skeleton is being used for short-term forecasts
      //    br = lag_birthrate[u];

      // cohort effect, added back in with cohort arriving over a time interval 0.05yr
      if (fabs(t-floor(t)-251.0/365.0) < 0.5*0.05)
        br = cohort[u*cohort_unit]*lag_birthrate[u]/0.05 + (1-cohort[u*cohort_unit])*lag_birthrate[u];
      else
        br = (1.0-cohort[u*cohort_unit])*lag_birthrate[u];

      foi = I[u]/pop[u];
      for (v=0; v < U ; v++) {
        if(v != u)
          foi += g[u*g_unit] * v_by_g[u][v] * (I[v]/pop[v] - I[u]/pop[u]) / pop[u];
      }

      DS[u] = br - (beta*foi + mu[u*mu_unit])*S[u];
      DE[u] = beta*foi*S[u] - (sigma[u*sigma_unit]+mu[u*mu_unit])*E[u];
      DI[u] = sigma[u*sigma_unit]*E[u] - (gamma[u*gamma_unit]+mu[u*mu_unit])*I[u];
      DR[u] = gamma[u*gamma_unit]*I[u] - mu[u*mu_unit]*R[u];
      DC[u] = gamma[u*gamma_unit]*I[u];
    }
  ')

basic_log_names <- c("sigma", "gamma", "sigmaSE", "psi", "R0", "g","iota","mu")
basic_log_names <- setdiff(basic_log_names,fixedParNames)

basic_logit_names <- c("amplitude", "alpha","cohort","rho","S_0", "E_0", "I_0")
basic_logit_names <- setdiff(basic_logit_names,fixedParNames)

log_names <- unlist(lapply(basic_log_names, function(x,U) paste0(x,1:U),U))
logit_names <- unlist(lapply(basic_logit_names, function(x,U) paste0(x,1:U),U))
  
## it is possible for S+E+I to be greater than P,
## in which case R is negative, but that is not necessarily a critical problem.

measles_partrans <- parameter_trans(log=log_names,logit=logit_names)

m1 <-  spatPomp(measles_cases,
          units = "city",
          times = "year",
          t0 = min(measles_cases$year)-1/26,
          unit_statenames = basic_statenames,
          covar = measles_covar,
          rprocess=euler(measles_rprocess, delta.t=dt),
          skeleton=vectorfield(measles_skel),
          unit_accumvars = c("C"),
          paramnames=measles_paramnames,
          partrans=measles_partrans,
          globals=measles_globals,
          rinit=measles_rinit,
          dmeasure=measles_dmeasure,
          eunit_measure=measles_eunit_measure,
          munit_measure=measles_munit_measure,
          vunit_measure=measles_vunit_measure,
          rmeasure=measles_rmeasure,
          dunit_measure=measles_dunit_measure
  )

measles_params <- rep(0,length=length(measles_paramnames))
names(measles_params) <- measles_paramnames


for(p in fixedParNames) measles_params[paste0(p,1)] <- basic_params[p]
for(p in sharedParNames) measles_params[paste0(p,1:U)] <- basic_params[p]

# unit-specific parameters are perturbed  to make them different
for(p in unitParNames) {
  measles_params[paste0(p,1:U)] <- rnorm(
    U,
    mean=basic_params[p],
    sd=basic_params[p]*0.1
  )
}

coef(m1) <- measles_params
if(simulated) simulate(m1) else m1
}
