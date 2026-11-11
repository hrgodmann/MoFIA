# Read-only constructor and exact conditional-predictor checks for commit b5f7ef5.
# No optimization, likelihood evaluations, fitted models, or data simulations.
path <- '/Users/henrikgodmann/Desktop/GitHub/MoFIA'
e <- new.env(); sys.source(file.path(path, 'R/model.R'), e)
old <- new.env(); src <- system2('git', c('-C', shQuote(path), 'show', '08b9acb:R/model.R'), stdout=TRUE)
source(textConnection(src), local=old)
cat(R.version.string, '\n')
cat('Current source commit:', system2('git',c('-C',shQuote(path),'rev-parse','HEAD'),stdout=TRUE),'\n')
its <- paste0('y',1:3); n <- 120L
base <- data.frame(y1=rep(0:2,40),y2=rep(c(1,2,0),40),y3=rep(c(2,0,1),40),
                   z=rep(0:1,60),a=rep(c(0,0,1,1),30),w=seq(-1,1,length.out=n))
run <- function(env, args) {
  warnings <- list(); error <- NULL
  m <- tryCatch(withCallingHandlers(do.call(env$mofia_model,args),warning=function(w) {
    warnings[[length(warnings)+1L]] <<- w; invokeRestart('muffleWarning')
  }),error=function(e) {error <<- e;NULL})
  list(model=m,warnings=warnings,error=error)
}
args <- list(data=base,items=its,moderators='z',threshold_mod=TRUE,loading_mod=TRUE,var_mod=FALSE,sparse=0)
checks <- list()
r <- run(e,args);checks$binary_exact <- inherits(r$error,'mofia_identification_error') && identical(r$error$issues$status,'not_identified')
args$moderators <- ~ z * a
r <- run(e,args); checks$binary_interactions_exact <- inherits(r$error,'mofia_identification_error')
args$moderators <- ~ w + I(w^2);args$mean_mod<-'w';args$loading_mod<-list(y1='w')
r <- run(e,args);checks$quadratic_absorbed_exact <- inherits(r$error,'mofia_identification_error')
args$moderators <- 'w';args$mean_mod<-TRUE;args$loading_mod<-TRUE
r <- run(e,args);checks$continuous_warning <- is.null(r$error) && length(r$warnings)==1L && inherits(r$warnings[[1]],'mofia_identification_warning') && identical(r$model$identification_issues$status,'identified_by_loading_moderation')
checks$warning_payload_matches_model <- identical(r$model$identification_issues,r$warnings[[1]]$issues)
args$threshold_mod <- c('y2','y3')
r <- run(e,args);checks$anchor_passes <- is.null(r$error) && !length(r$warnings) && nrow(r$model$identification_issues)==0
args$threshold_mod<-TRUE;args$mean_mod<-FALSE
r <- run(e,args);checks$no_mean_passes <- is.null(r$error) && !length(r$warnings)
cat('\nTargeted existing contracts:\n');print(unlist(checks));stopifnot(all(unlist(checks)))

# Regression: identical supported model family under an origin change.
a <- 10000;d <- base;d$x <- a+d$w
args<-list(data=d,items=its,moderators='x',threshold_mod=TRUE,loading_mod=TRUE,var_mod=FALSE,sparse=0)
before <- run(old,args);after <- run(e,args)
args$moderators <- 'w';centered<-run(e,args)
cat('\nOrigin-change regression:\n')
cat('Old raw source: error=',!is.null(before$error),'warnings=',length(before$warnings),'\n')
cat('New raw source: error=',!is.null(after$error),'class=',class(after$error)[1],'\n')
print(after$error$issues)
cat('New centered source: error=',!is.null(centered$error),'warnings=',length(centered$warnings),'\n')
print(centered$model$identification_issues)
cat('Product x^2 residual in span{1,x}, absolute max=',max(abs(qr.resid(qr(cbind(1,d$x)),d$x^2))),'\n')
cat('Centered polynomial design rank=',qr(cbind(1,d$w,d$w^2))$rank,'\n')
# Reparameterization check at fixed generic coefficients; eta residual is arbitrary,
# hence the conditional distributions match and so do their integrated versions.
gamma<-.3;lambda0c<-1.2;lambda1<-.1;betac<-.2;nu0c<-.5
lambda0<-lambda0c-lambda1*a;beta<-betac-lambda1*gamma*a
nu0<-nu0c-beta*a-lambda0c*gamma*a
residual<-rep(c(-1,0,1),40)
pred_raw<-nu0+beta*d$x+(lambda0+lambda1*d$x)*(gamma*d$x+residual)
pred_centered<-nu0c+betac*d$w+(lambda0c+lambda1*d$w)*(gamma*d$w+residual)
cat('Max conditional-predictor discrepancy after exact parameter transformation=',max(abs(pred_raw-pred_centered)),'\n')
stopifnot(is.null(before$error), is.null(after$error), is.null(centered$error), identical(after$model$identification_issues$status, centered$model$identification_issues$status), max(abs(pred_raw-pred_centered))<1e-8)  # follow-up: raw and centred agree
cat('\nFollow-up: seven contract checks passed; raw and centred models get the same classification.\n')
