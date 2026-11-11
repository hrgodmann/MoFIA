# Deterministic algebra review of b5f7ef5. No MoFIA fit, real likelihood,
# simulated data, package rebuild, or shared library loading is performed.
repo <- "/Users/henrikgodmann/Desktop/GitHub/MoFIA"
e <- new.env(parent = globalenv())
sys.source(file.path(repo, "R/fit.R"), e)
sys.source(file.path(repo, "R/tests.R"), e)
sys.source(file.path(repo, "R/scoretest.R"), e)
n <- 0L; nf <- 0L
check <- function(ok, label) {
  if (isTRUE(ok)) { n <<- n + 1L; cat(sprintf("PASS %02d %s\n", n, label)) }
  else { nf <<- nf + 1L; cat(sprintf("FAIL %02d %s\n", nf, label)) }
}
# 1e-6 relative tolerance accommodates cancellation in the deliberately extreme
# fixed-step quadratic-score fixture (largest observed covariance error ~2.5e-8).
# The initial 1e-8 run and diagnostic are retained separately, not treated as an
# engine fit failure or silently accepted snapshot.
close <- function(a, b, tol = 1e-6) isTRUE(all.equal(unname(a), unname(b), tolerance = tol, check.attributes = FALSE))

A <- matrix(c(4, 1, .5, 1, 3, .2, .5, .2, 2), 3)
L <- matrix(c(1, .1, .2, .3, 1.5, .1, .1, .3, 1.2), 3)
Sbase <- rbind(L, -L)
B <- crossprod(Sbase)
G <- diag(c(.2, -.1, .3))
g <- c(2e-6, -1e-6, 3e-6)
pnames <- c("p1", "p2", "p3")
rownames(A) <- colnames(A) <- pnames
rownames(B) <- colnames(B) <- pnames

# Exercise actual .mofia_se body with a deterministic quadratic score provider.
# The provider is an algebra fixture, not an engine likelihood.
e$.mofia_nodes <- function(...) NULL
e$.mofia_ll <- function(theta, model, ...) {
  S <- model$S0 - matrix(drop(model$A %*% theta) / nrow(model$S0), nrow(model$S0), length(theta), byrow = TRUE)
  colnames(S) <- names(theta)
  structure(0, scores = S, phi_scores = matrix(0, nrow(S), ncol(S)), resvar_scores = matrix(0, nrow(S), ncol(S)))
}
e$.constraint_curvature <- function(theta, model, W) model$G
run_se <- function(k, Gcase = G) {
  D <- diag(k)
  model <- list(A = D %*% A %*% D, G = D %*% Gcase %*% D,
    S0 = (Sbase + matrix(g / nrow(Sbase), nrow(Sbase), 3, byrow = TRUE)) %*% D,
    Q = 0, D = 1, par = data.frame(mod = rep(NA_integer_, 3)))
  f <- list(model = model, coefficients = setNames(rep(0, 3), pnames), nq = 1, link = "logit", adapt = NULL)
  e$.mofia_se(f)
}
base <- run_se(rep(1, 3))
Bactual <- B + tcrossprod(g) / nrow(Sbase)
for (k in list(c(1, 1, 1), c(1, 1e4, 1e-4), c(-1e2, 1e-2, -1e3))) {
  f <- run_se(k)
  check(identical(f$status, "ok"), paste("information status under", paste(k, collapse = ",")))
  for (kind in c("exact", "reduced")) {
    Ak <- if (kind == "exact") A else A + G
    Ai <- solve(Ak)
    check(close(f$vcov[[kind]]$ML * outer(k, k), Ai), paste(kind, "ML back transform"))
    check(close(f$vcov[[kind]]$MLR * outer(k, k), Ai %*% Bactual %*% Ai), paste(kind, "MLR back transform"))
    check(close(f$scaling[[kind]], sum(diag(Ai %*% Bactual)) / 3), paste(kind, "trace scaling"))
  }
  check(close(f$vcov$MLF * outer(k, k), solve(Bactual)), "MLF back transform")
  check(close(f$newton_decrement, .5 * drop(t(g) %*% solve(A, g))), "Newton decrement")
  check(f$converged, "Newton convergence")
}
f <- run_se(rep(1, 3), diag(c(0, 1, 1)) - A)
check(is.null(f$vcov$reduced), "singular reduced information stays unavailable")
check(is.na(f$scaling[["reduced"]]), "singular reduced scaling stays unavailable")
check(!is.null(f$vcov$exact), "exact covariance remains available when reduced unavailable")
for (mat in list(matrix(c(1,1,1,1),2), matrix(0,2,2), matrix(c(1,0,0,0),2)))
  check(identical(e$.info_inverse(mat)$status, "singular"), "true singular matrix rejected")
for (mat in list(diag(c(-1,1)), matrix(c(0,1,1,0),2)))
  check(identical(e$.info_inverse(mat)$status, "saddle"), "indefinite matrix rejected")
check(is.null(e$.info_inverse(matrix(NaN,2,2))$inverse), "nonfinite matrix rejected")

# Actual Wald code with a covariance stub and supplied fit coefficients.
e$vcov <- function(object, ...) object$V
theta <- setNames(c(.2, -.4, .1), pnames)
V <- solve(A)
Wref <- drop(t(theta) %*% solve(V, theta))
for (k in list(c(1,1,1), c(1,1e6,1e-5), c(-1e2,1e-2,-1e3))) {
  fit <- list(coefficients = theta / k, converged = TRUE, V = V / outer(k,k))
  check(close(e$mofia_wald(fit, pnames)$statistic, Wref), "actual Wald function invariant under units/sign")
}

# Actual score-test body on deterministic quadratic scores. These stubs do not
# call mofia_model(), mofia(), or the real .mofia_ll().
e$.add_params <- function(model, candidates) model$extended
e$.candidate_group <- function(candidates) rep("group", length(candidates))
AA <- matrix(c(5,.5,.4,.2, .5,4,.3,.5, .4,.3,3,.4, .2,.5,.4,2), 4)
SS <- rbind(diag(c(1,2,1.5,1)), -diag(c(1,2,1.5,1)))
ug <- c(0,0,.4,-.3)
SS <- SS + matrix(ug/nrow(SS), nrow(SS),4,byrow=TRUE)
pn <- c("n1","n2","c1","c2")
dimnames(AA) <- list(pn,pn)
score_run <- function(k, mat = AA) {
  D <- diag(k); S <- SS %*% D; colnames(S) <- pn
  ext <- list(A = D %*% mat %*% D, S0 = S, par = data.frame(name = pn, mod = rep(NA_integer_,4)), Q=0,D=1)
  A0 <- solve(mat[1:2,1:2]) / outer(k[1:2], k[1:2]); dimnames(A0)<-list(pn[1:2],pn[1:2])
  fit <- structure(list(coefficients = setNames(c(0,0),pn[1:2]), model=list(par=data.frame(name=pn[1:2]),extended=ext),
    vcov=list(exact=list(ML=A0)),converged=TRUE,nq=1,link="logit",adapt=NULL),class="mofia")
  list(single=e$mofia_score_test(fit,pn[3:4]),group=e$mofia_score_test(fit,pn[3:4],by="item"))
}
base_st <- score_run(rep(1,4))
for (k in list(c(1e2,1e-2,1e4,1e-3),c(-1e2,1e-2,-1e3,1e3))) {
  st <- score_run(k)
  for (kind in c("single","group")) {
    x<-st[[kind]]; y<-base_st[[kind]]; x<-x[match(y$test,x$test),]
    check(close(x$statistic,y$statistic),paste("actual score robust",kind,"unit invariance"))
    check(close(x$statistic_ML,y$statistic_ML),paste("actual score ML",kind,"unit invariance"))
    if (kind=="single") check(close(x$epc*k[match(x$test,pn)],y$epc),"score EPC back transform")
  }
}
# Exact conditional alias: candidate c1 duplicates nuisance n1's information.
sing <- AA; sing[3,] <- sing[1,]; sing[,3] <- sing[,1]
for (k in list(rep(1,4),c(1e2,1e-2,1e4,1e-3))) {
  st<-score_run(k,sing)$single
  check(is.na(st$statistic[st$test=="c1"]),"conditional alias score test unavailable under rescaling")
}

# Nonblocking limitations retained explicitly, rather than asserted away.
Z <- matrix(c(1, 1e-6, 1e-6, 0),2)
D <- diag(c(1,1e6))
cat("LIMIT zero-diagonal indefinite matrix status:",e$.info_inverse(Z)$status,
    "rescaled:",e$.info_inverse(D%*%Z%*%D)$status,"(both withhold inverse)\n")
fd <- function(k) e$.fd_steps(0,list(Q=1,par=data.frame(mod=1L),X=matrix(c(-1,1)*k,2)))
cat("LIMIT finite difference predictor displacement at X scales 1 and 1e-4:",fd(1),fd(1e-4)*1e-4,"\n")
cat("TOTAL",n,"passed;",nf,"failed; no engine fits or simulations run\n")
cat("R",as.character(getRversion()),R.version$platform,"\n")
