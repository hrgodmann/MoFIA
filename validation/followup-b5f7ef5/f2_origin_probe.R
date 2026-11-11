.libPaths(c(commandArgs(TRUE)[1], .libPaths())); library(MoFIA)
its <- paste0('y',1:3); n <- 120
base <- data.frame(y1=rep(0:2,40),y2=rep(c(1,2,0),40),y3=rep(c(2,0,1),40), z=rep(0:1,60),a=rep(c(0,0,1,1),30),w=seq(-1,1,length.out=n))
cl <- function(expr) tryCatch(withCallingHandlers({expr; 'no issue'}, warning=function(w) if (inherits(w,'mofia_identification_warning')) {cat(' [warning:', w$issues$status, ']'); invokeRestart('muffleWarning')}), mofia_identification_error=function(e) paste('ERROR', e$issues$status, e$issues$moderators), error=function(e) paste('other error:', conditionMessage(e)))
for (o in c(0, 1e4, 1e9)) { d <- base; d$x <- o + d$w; cat('continuous x, origin', o, ':', cl(mofia_model(d, items=its, moderators='x', threshold_mod=TRUE, loading_mod=TRUE, var_mod=FALSE, sparse=0)), '\n') }
for (o in c(0, 1e4)) { d <- base; d$x <- o + d$w; cat('x + I(x^2), origin', o, ':', cl(mofia_model(d, items=its, moderators=~x+I(x^2), mean_mod='x', threshold_mod=TRUE, loading_mod=list(y1='x'), var_mod=FALSE, sparse=0)), '\n') }
for (cod in list(c(0,1), c(1,2), c(0,1000), c(-5e4, -5e4+1), c(1e9, 1e9+1))) { d <- base; d$zz <- cod[d$z+1]; cat('binary coded', cod, ':', cl(mofia_model(d, items=its, moderators='zz', threshold_mod=TRUE, loading_mod=TRUE, var_mod=FALSE, sparse=0)), '\n') }
d <- base; d$a2 <- 1e4 + d$a; d$z2 <- 1e4 + d$z; cat('binary interaction, offset 1e4:', cl(mofia_model(d, items=its, moderators=~z2*a2, threshold_mod=TRUE, loading_mod=TRUE, var_mod=FALSE, sparse=0)), '\n')
d <- base; d$x <- 1e4 + d$w; cat('binary mean, loading by far-origin continuous x:', cl(mofia_model(d, items=its, moderators=c('z','x'), mean_mod='z', threshold_mod=list(y1='z',y2='z',y3='z'), loading_mod=list(y1='x'), var_mod=FALSE, sparse=0)), '\n')
cat('binary mean+DIF, loadings by z and a (a not in thresholds):', cl(mofia_model(base, items=its, moderators=c('z','a'), mean_mod='z', threshold_mod=list(y1='z',y2='z',y3='z'), loading_mod=TRUE, var_mod=FALSE, sparse=0)), '\n')
cat('anchor:', cl(mofia_model(base, items=its, moderators='z', threshold_mod=c('y2','y3'), loading_mod=TRUE, var_mod=FALSE, sparse=0)), '\n')
cat('no mean:', cl(mofia_model(base, items=its, moderators='z', threshold_mod=TRUE, loading_mod=TRUE, mean_mod=FALSE, var_mod=FALSE, sparse=0)), '\n')
