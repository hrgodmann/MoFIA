# Minimal deterministic regression of the numerical reduced-information boundary.
# No package is loaded and no actual MoFIA likelihood or model fit is executed.
root <- "/Users/henrikgodmann/Desktop/GitHub/jaspSem/tempAI/ReworkMNLFA-engine-fix-review-20261111/information"
lines <- readLines(file.path(root,"probe.R"))
end <- grep("^base <- run_se", lines) - 1L
eval(parse(text=lines[seq_len(end)]))
G0 <- diag(c(0,1,1)) - A
f <- run_se(rep(1,3),G0)
theta<-f$coefficients; model<-f$model
h<-e$.fd_steps(theta,model)
gf<-function(t) colSums(attr(e$.mofia_ll(t,model),"scores"))
H<-sapply(seq_along(theta),function(j) {
  d<-replace(numeric(length(theta)),j,h[j]); (gf(theta+d)-gf(theta-d))/(2*h[j])
})
H<-(H+t(H))/2; Ir <- -H+G0
cat("Analytical reduced information:\n");print(A+G0)
cat("Numerical reduced information entering .info_inverse:\n");print(Ir,digits=16)
cat("Equilibrated status/eigenvalues:\n");print(e$.info_inverse(Ir)[c("status","eigen")])
cat("Resulting reduced covariance and scaling:\n");print(f$vcov$reduced);print(f$scaling)
# Same mocked inputs and derivative step in original code. Only fitting functions
# are replaced; algebra score provider and curvature provider remain identical.
sys.source(file.path(root,"base-fit.R"),e)
old<-run_se(rep(1,3),G0)
cat("Original commit reduced covariance:\n");print(old$vcov$reduced)
cat("Original commit scaling:\n");print(old$scaling)
stopifnot(is.null(old$vcov$reduced), is.null(f$vcov$reduced), is.na(f$scaling[["reduced"]]))  # follow-up: corrected behaviour asserted
cat("Follow-up: reduced information withheld by both old and corrected code.\n")
# Candidate mitigation demonstration only; not an engine patch or adopted rule.
# Exact information supplies a scale that survives cancellation in the reduced
# matrix. Under I -> D I D and A -> D A D its spectrum remains invariant (up to
# sign congruence for negative units).
cat("Common scaling using exact-information diagonal (demonstration):\n")
for (k in list(rep(1,3),c(1,1e4,1e-4),c(-1e2,1e-2,-1e3))) {
  D<-diag(k); Ix<-D%*%A%*%D; Irt<-D%*%Ir%*%D
  sr<-sqrt(diag(Ix)); ev<-eigen(Irt/outer(sr,sr),symmetric=TRUE,only.values=TRUE)$values
  avail<-min(ev)>1e-8*max(abs(ev))
  cat("units",paste(k,collapse=","),"min",min(ev),"max",max(ev),"available",avail,"\n")
  stopifnot(!avail)
}
