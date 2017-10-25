# never monkey-patch. never say never.
stan_model_builder <- rstan::stan_model
body(stan_model_builder) <- (function(){eval(substitute(substitute(stan_model,list(stanc=quote(function(..., model_name=NULL, model_code=NULL) stanc_builder(...)))),list(stan_model=body(rstan::stan_model))))})()
#' @export stan_model_builder
stan_model_builder
