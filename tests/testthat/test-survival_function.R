# Import libraries ----
library(tmlensemble)
library(survival)
library(data.table)
library(reticulate)
library(keras)
library(xgboost)
library(glmnet)
library(rpart)
library(simevent)
# Define useful functions ----

sf <- function(t,lambda, beta, X1,beta1,X2){

  out <- exp(-lambda*exp(beta*X1+beta1*(X2))*t)
  return(out)
}

sf_simevent <- function(t,eta,nu){

  out <- exp(-eta*(t^{nu}))
  return(out)
}

inverse_sf<- function(lambda, beta, X1,beta1,X2, u){

  out <- -log(u)/(lambda*exp(beta*X1+beta1*(X2)))

  return(out)
}



# Define useful variables ----
seed=42
n=50

# Generate some data from a simple proportional model ----

{
  set.seed(seed)

  sample <- 1-runif(as.integer(n), min=0, max= 1)
  X <- rnorm(n)
  X2 <- sample(c(0,1), n, replace = TRUE)

  sample2 <- 1-runif(as.integer(n/3), min=0, max= 1)
  X.val <- rnorm(as.integer(n/3))
  X2.val <- sample(c(0,1), as.integer(n/3), replace = TRUE)

}

time<- inverse_sf(u=sample,
                  lambda=0.5,
                  beta=0.1,
                  beta1=0.2,
                  X2=X2,
                  X1=X)

time.val<- inverse_sf(u=sample2,
                  lambda=0.5,
                  beta=0.1,
                  beta1=0.2,
                  X2=X2.val,
                  X1=X.val)
{
  set.seed(1)
  cens <- runif(n=length(time),
                min=0,
                max=5.5)

  cens2 <- runif(n=length(time.val),
                min=0,
                max=5.5)

}

tmp <- time
tmp2 <- time.val
tmp[time > cens] <- cens[time > cens]
tmp2[time.val > cens2] <- cens2[time.val > cens2]


dt <- data.table(
  id=1:length(tmp),
  time = tmp,
  covariate = X,
  covariate2 = as.character(X2),
  status = 1-as.numeric(time > cens)


)

dt_val <- data.table(
  id=1:length(tmp2),
  time = tmp2,
  covariate = X.val,
  covariate2 = as.character(X2.val),
  status = 1-as.numeric(time.val > cens2)


)

# Fit ----
## Learners ----

l1 <- Learner_glm(covariates = c("covariate","covariate2"))

l2 <- Learner_glm(covariates = c("covariate2")) # there is no CV for a glm, so I return a warning

## Fit Superlearner ----

learners <- list(l1,l2)

# I make different examples below. One can either use a glm or a glmnet and
# either choose to use or not to use the nodes in the meta learner.

sl <- Superlearner(data=dt,
                   learners=learners,
                   id="id",
                   status="status",
                   # nodes=seq(0,6,.5),
                   nfold = 3,
                   meta_learner_algorithm = "glm",
                   add_nodes_metalearner = TRUE,
                   add_intercept_metalearner = TRUE,
                   event_time = "time")

sl1 <- Superlearner(data=dt,
                   learners=list(Learner_glm(covariates = c("covariate")),l2),
                   id="id",
                   status="status",
                   # nodes=seq(0,6,.5),
                   nfold = 3,
                   meta_learner_algorithm = "glm",
                   add_nodes_metalearner = TRUE,
                   add_intercept_metalearner = TRUE,
                   event_time = "time")

coef(sl1$superlearner$`1`$meta_learner_fit)

coef(sl$superlearner$`1`$meta_learner_fit)

log(min(dt$time))

# Integrated functionalities into RiskRegression: new predict method consistent with coxph and stuff ----
out <- predict(sl,dt_val,times=c(0),cause=1)

cox <- coxph(formula("Surv(time, status) ~  covariate  +   covariate2"), data=dt,x=T)

baseline_hazard <- basehaz(cox)

surv_fit <- survfit(cox,
                    data.frame(covariate = 5,
                               covariate2=factor("0",levels = c("0","1"))))


## When predictRisk takes in times and newdata, in the prediction phase, does the model disregard the time
## and simply computes feature based predictions based on the given times? aka dt.test <- data.frame(x1 = .., x2=..., time= ...) is
## time disregarded when predict.somemodel(dt.test,times)? it seems like it.


# we can now compare different superlearners in terms of some scoring metrics ----
library(prodlim)
library(survival)
Score(
  object = list("SL (repl. cox)" = sl,
                "SL(weak learners)" = sl1),
  formula = Surv(time,status)~1,
  times = quantile(dt$time),
  data=dt,
  metrics="Brier",
  conf.int=FALSE)

# Interestingly, one can also compare standard literature benchmark ----

Score(
  object = list("SL (repl. cox)" = sl,
                "COX" = cox),
  formula = Surv(time,status)~1,
  times = quantile(dt$time),
  data=dt,
  metrics="Brier",
  conf.int=FALSE)



#here it seems fine
1-predictRisk(cox,newdata=data.frame(covariate = 5,
                                   covariate2=factor("0",levels = c("0","1"))),
              times= c(3.389589344 ) )


1- predictRisk(sl,newdata=data.frame(covariate = 5,
                                     covariate2=factor("0",levels = c("0","1"))),
               cause=1,
               times= c( 3.389589344))


# cox acts weird?
1-predictRisk(cox,newdata=data.frame(covariate = 5,
                                     covariate2=factor("0",levels = c("0","1"))),times= c(2.274364488 ) )


1- predictRisk(sl,newdata=data.frame(covariate = 5,
                                     covariate2=factor("0",levels = c("0","1"))),
               cause=1,
               times= c(2.274364488))


surv_fit$time[37]
surv_fit$surv[37]

# When predicting, should we get the whole interval like cox?

out_cox <- 1-predictRisk(cox,data.frame(covariate = 5,
                             covariate2=factor("0",levels = c("0","1"))),
               times=sort(unique(dt$time)))

out <- predict(sl,data.frame(covariate = 5,
                             covariate2=factor("0",levels = c("0","1"))),
                                               times=sort(unique(dt$time)),
                                               cause=1)

out

# example
median(dt$time)

out[time>= 1.1 & time < 1.2, ]

1-predictRisk(cox,data.frame(covariate = 5,
                             covariate2=factor("0",levels = c("0","1"))),
              times=median(dt$time))

1-predictRisk(sl,data.frame(covariate = 5,
                             covariate2=factor("0",levels = c("0","1"))),
              times=median(dt$time))

full_out_cox <- predictRisk(cox, dt_val, times = sort(unique(dt$time)))

full_out_sl <- predictRisk(sl, dt_val, times = sort(unique(dt$time)))

predictRisk(sl,data.frame(covariate = 5,
                          covariate2=factor("0",levels = c("0","1"))),
            times=sort(unique(dt$time)),
            cause=1)

predict(sl,data.frame(covariate = 5,
                      # time=.2,
                      covariate2=factor("0",levels = c("0","1"))),
        times=c(min(dt$time)),
        cause=1)


1-predictRisk(cox,newdata=data.frame(covariate = 5,

                                     covariate2=factor("0",levels = c("0","1"))),cause=1,times=.25 )


predict(sl,data.frame(covariate = 5,
                      time=.2,
                      covariate2=factor("0",levels = c("0","1"))),
        times=0,
        cause=1)


predictRisk(cox,dt_val,times=c(2.486492854),type = "survival")

1-predictRisk(sl,dt_val,times=c(2.486492854),type = "survival")




# Combinations we want to predict for
## Right now the predict method is mostly something that I used for the plots I showed.
## Once should think how to make it userfriendly in the future.
tmp_05 <- data.frame(node = factor(sort(unique(seq(0,5,0.5)
))),
tij = .5,
covariate = .5,
deltaij=0,
covariate2=factor("0",levels = c("0","1")))


tmp_1 <- data.frame(node = factor(sort(unique(seq(0,5,1)
))),
tij = 1,
covariate = .5,
deltaij=0,
covariate2=factor("0",levels = c("0","1")))


# Benchmark cox model ----
cox <- coxph(formula("Surv(time, status) ~  covariate  +   covariate2"), data=dt)

baseline_hazard <- basehaz(cox)

surv_fit <- survfit(cox,
                    data.frame(covariate = .5,
                               covariate2=factor("0",levels = c("0","1"))))



preds <- predict(sl,tmp_05)
preds_1 <- predict(sl1,tmp_1)

xtrue <- seq(0, 6, 0.001)
ytrue <- sf(
  xtrue,
  lambda = 0.5,
  beta = 0.1,
  beta1 = 0.2,
  X1 = .5,
  X2 = 0
)

plot(xtrue,
     ytrue,
     type="l",
     xlab="Time", ylab="Survival Function")
lines(surv_fit$time,
      surv_fit$surv,
      col="yellow",
      type="S")
lines(seq(0,5.5,by=.5),c(1,preds$survival_function), col="red",type="S")
lines(seq(0,6,by=1),c(1,preds_1$survival_function), col="green",type="S")

# Test on Micaelas code ----

{
  set.seed(7)
  dat <- simSurvData(1e2)
}


l2 <- Learner_glm(covariates = c("L0"),
                  intercept= TRUE)

l3 <- Learner_glm(covariates = c("L0","A0"),
                  intercept= TRUE)

tmp_5 <- data.frame(node = factor(sort(unique(seq(0,27,.5)
))),
tij = 0.5,
L0 = dat[,median(L0)],
A0=factor("0",levels = c("0","1")),
deltaij=0
)


tmp_1 <- data.frame(node = factor(sort(unique(seq(0,27,1)
))),
tij = 1,
L0 = dat[,median(L0)],
A0=factor("0",levels = c("0","1")),
deltaij=0
)

observed_nodes <- sort(unique(dat$Time))

tmp_observed <- data.frame(node = factor(c(0,observed_nodes[observed_nodes!=max(observed_nodes)])),
                           tij = diff(c(0,observed_nodes)),
                           L0 = dat[,median(L0)],
                           A0=factor("0",levels = c("0","1")),
                           deltaij=0
)

dat[,A0:=factor(A0,levels=c("0","1"))]

learners <- list(l3,l2)

sl <- Superlearner(data=dat,
                   learners=learners,
                   stratified_k_fold = TRUE,
                   id="ID",
                   status="Delta",
                   # nodes=seq(0,35,.5),
                   event_time =  "Time",
                   nfold=10,
                   meta_learner_algorithm = "glm",
                   add_nodes_metalearner = TRUE,
                   add_intercept_metalearner = TRUE)

sl1 <- Superlearner(data=dat,
                   learners=learners,
                   stratified_k_fold = FALSE,
                   id="ID",
                   status="Delta",
                   nodes=seq(0,35,.5),
                   event_time =  "Time",
                   nfold=10,
                   meta_learner_algorithm = "glm",
                   add_nodes_metalearner = TRUE,
                   add_intercept_metalearner = TRUE)

preds <- predict(sl,tmp_observed)
preds1 <- predict(sl1,tmp_5)

cox <- coxph(formula("Surv(Time, Delta) ~  L0  +   A0"), data=dat)

baseline_hazard <- basehaz(cox)

surv_fit <- survfit(cox,
                    data.frame(L0 = dat[,median(L0)],
                               A0=factor("0",levels = c("0","1"))))

x_true <- seq(0,28,.01)
y_true <- sf_simevent(x_true,.1,1.1)

{
  plot(c(0,observed_nodes),
       c(1, preds$survival_function),
       col = "green",
       type = "S")
  lines(surv_fit$time, surv_fit$surv, type = "S")
  lines(x_true, y_true, col = "red")

}

## Thomas' code

dt <- synthesize_td1(100)





