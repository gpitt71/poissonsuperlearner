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
n=1000

# Generate some data from a simple proportional model ----

{
  set.seed(seed)

  sample <- 1-runif(as.integer(n), min=0, max= 1)
  X <- rnorm(n)
  X2 <- sample(c(0,1), n, replace = TRUE)

}

time<- inverse_sf(u=sample,
                  lambda=0.5,
                  beta=0.1,
                  beta1=0.2,
                  X2=X2,
                  X1=X)
max(time)

{
  set.seed(1)
  cens <- runif(n=length(time),
                min=0,
                max=5.5)

}

tmp <- time
tmp[time > cens] <- cens[time > cens]


dt <- data.table(
  id=1:length(tmp),
  time = tmp,
  covariate = X,
  covariate2 = as.character(X2),
  status = 1-as.numeric(time > cens)


)

# Fit ----
## Learners ----

l1 <- Learner_glm(covariates = c("covariate2"),
                     cross_validation=FALSE,
                     lambda=.1)

l2 <- Learner_glm(covariates = c("covariate","covariate2"),
                     cross_validation=TRUE) # there is no CV for a glm, so I return a warning

## Fit Superlearner ----

learners <- list(l1,l2)

# I make different examples below. One can either use a glm or a glmnet and
# either choose to use or not to use the nodes in the meta learner.

sl <- Superlearner(data=dt,
                   learners=learners,
                   id="id",
                   status="status",
                   nodes=seq(0,6,.5),
                   nfold = 5,
                   meta_learner_algorithm = "glm",
                   add_nodes_metalearner = FALSE,
                   event_time = "time")

sl1 <- Superlearner(data=dt,
                   learners=learners,
                   id="id",
                   status="status",
                   nodes=seq(0,6,1),
                   nfold = 5,
                   meta_learner_algorithm = "glmnet",
                   add_nodes_metalearner = T,
                   event_time = "time")

sl2 <- Superlearner(data=dt,
                    learners=list(Learner_glmnet(covariates = c("covariate2"),
                                                 cross_validation=FALSE,
                                                 lambda=.1)),
                    id="id",
                    status="status",
                    nodes=seq(0,6,1),
                    nfold = 5,
                    meta_learner_algorithm = "glmnet",
                    add_nodes_metalearner = TRUE,
                    event_time = "time")


# Combinations we want to predict for
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



preds <- predict(sl,tmp_05,step = .5)
preds_1 <- predict(sl1,tmp_1,step = 1)
preds_2 <- predict(sl2,tmp_1,step = 1)

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
lines(seq(0,6,by=1),c(1,preds_2$survival_function), col="blue",type="S")

# Test on Micaelas code ----

{
  set.seed(7)
  dat <- simSurvData(1e4)
}


l2 <- Learner_glm(covariates = c("L0"))

l3 <- Learner_glm(covariates = c("L0","A0"))

tmp_5 <- data.frame(node = factor(sort(unique(seq(0,27,.5)
))),
tij = .5,
L0 = dat[,median(L0)],
A0=factor("0",levels = c("0","1")),
deltaij=0
)


dat[,A0:=factor(A0,levels=c("0","1"))]

learners <- list(l3,l2)

sl <- Superlearner(data=dat,
                   learners=learners,
                   id="ID",
                   status="Delta",
                   nodes=seq(0,35,.5),
                   event_time =  "Time",
                   nfold=10,
                   meta_learner_algorithm = "glmnet")

preds <- predict(sl,tmp_5,step=.5)
preds1 <- predict(sl1,tmp_5,step=.5)

cox <- coxph(formula("Surv(Time, Delta) ~  L0  +   A0"), data=dat)

baseline_hazard <- basehaz(cox)

surv_fit <- survfit(cox,
                    data.frame(L0 = dat[,median(L0)],
                               A0=factor("0",levels = c("0","1"))))

x_true <- seq(0,28,.01)
y_true <- sf_simevent(x_true,.1,1.1)

{
  plot(c(seq(0, 27.5, .5)),
       c(1, preds$survival_function),
       col = "green",
       type = "S")
  lines(surv_fit$time, surv_fit$surv, type = "S")
  lines(x_true, y_true, col = "red")
}




