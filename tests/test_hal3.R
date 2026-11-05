# Run after Rcpp::sourceCpp and source() of the R/hal_token_irls.R file.
stopifnot(requireNamespace('Matrix',quietly=TRUE), requireNamespace('glmnet',quietly=TRUE))
set.seed(883)
for (degree in 1:3) {
 n <- 300L
 dt <- data.table::data.table(x1=runif(n),x2=runif(n),x3=runif(n),
   cat=factor(sample(letters[1:3],n,TRUE)),node=factor(sample(paste0('n',1:3),n,TRUE),levels=paste0('n',1:3)))
 variables <- c('x1','x2','x3','cat','node')
 ref <- Learner_hal(covariates=variables[-5],cross_validation=FALSE,
   max_degree=degree,num_knots=rep(4L,degree))$hal_basis(
     vars=variables,DT=dt,max_interaction=degree,knots_per_order=rep(4L,degree))
 tok <- hal3_tokens(dt,variables,degree,rep(4L,degree))
 X <- ref$X
 stopifnot(ncol(X)==tok$ncol)
 # Compare columns individually rather than relying on column names.
 for (j in seq_len(ncol(X))) {
   support <- which(hal3_token_indicator(tok, dt, j))
   truth <- which(as.numeric(X[,j])!=0)
   stopifnot(identical(sort(as.integer(support)),as.integer(truth)))
 }
 y <- rpois(n,exp(-1.6+0.5*(dt$x1>0.5)))
 offset <- log(runif(n,0.5,2))
 pf <- tok$default_penalty_factor
 lambda <- c(0.15,0.08,0.03)
 ours <- hal3_fit_tokens(tok,y,offset,lambda,pf,tol=1e-8)
 refglm <- glmnet::glmnet(X,y,family='poisson',alpha=1,lambda=lambda,
                        offset=offset,penalty.factor=pf,intercept=TRUE,
                        standardize=FALSE,thresh=1e-10)
 for(j in seq_along(lambda)){
   theta <- as.matrix(stats::coef(refglm,s=lambda[j]))
   eta_ref <- as.numeric(offset+theta[1]+X%*%theta[-1])
   eta_ours <- hal3_predict_tokens(ours,dt,s=lambda[j],type='link',newoffset=offset)
   delta <- max(abs(eta_ref-eta_ours))
   cat('degree=',degree,'lambda=',lambda[j], 'predictor error=',delta,
       'KKT=',ours$max_kkt[j], 'converged=',ours$converged[j],'\n')
   stopifnot(delta<2e-3,ours$converged[j],ours$max_kkt[j]<1e-5)
 }
}
cat('PASS: token basis supports and Poisson fits match sparse references.\n')
