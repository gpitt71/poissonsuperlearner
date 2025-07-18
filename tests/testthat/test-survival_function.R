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

# Check xgboost
data_pre_processing <- function(data,
                                id,
                                status,
                                event_time,
                                nodes=NULL
){


  # browser()

  setDT(data)

  # Handle competing risks ----
  ## for each of the competing risks (CR) we need to create a table
  n_crisks <- pmax(length(unique(data[[status]])) - 1,1)
  ## the CR tables are stuck on top of each other to allow for possible interactions
  dt_fit <- do.call(rbind, replicate(n_crisks, data, simplify = FALSE))
  ## we create an artificial k index. Table specific.
  dt_fit <- dt_fit[, k := rep(1:n_crisks, each = dim(data)[1])]


  # Data Transformation ----
  tmp <- c(id, "k")

  dt_fit <- eval(parse(text = paste("dt_fit[, .(node = create_offset_variable(nodes, time_to_event = ",
                                    event_time,
                                    ")[, 1]",
                                    ", tij = create_offset_variable(nodes, time_to_event = ",
                                    event_time,
                                    ")[,2]",
                                    ", deltaij = create_response_variable_c_risks(nodes,time_to_event = ",
                                    event_time,
                                    ", delta=",
                                    status,
                                    ", event_type = k)",
                                    ")",
                                    ", by = .(",
                                    id,
                                    ", k)",
                                    "]")))

  ## Retrieve covariates

  dt_fit <- merge(dt_fit, data, by = id, all.x = TRUE)

  setnames(dt_fit, c(id),c("id"))

  dt_fit[,c("node",
            "k"):=list(factor(node,levels=as.character(nodes)),
                       as.factor(k))]


  # browser()
  dt_fit[,node:=relevel(node,ref=as.character(last(nodes)))]

  return(dt_fit)

}


create_response_variable_c_risks <- function(nodes, time_to_event, delta, event_type){
  # browser()
  p_holder <- ifelse(delta == event_type, 1, 0)

  l <- sum(nodes < time_to_event)

  out <- c(rep(0, max(0, l - 1)),
           p_holder)

  return(out)
}

create_offset_variable <- function(nodes, delta, time_to_event){

  # browser()
  tmp <- c(nodes[nodes < time_to_event],
           first(nodes[nodes >= time_to_event]))



  if (all(nodes < time_to_event)) {
    tmp <- c(tmp, time_to_event)
  } else{
    tmp[length(tmp)] <- time_to_event
  }


  tij <- diff(c(tmp))

  grid_nodes <- c(nodes[nodes < time_to_event])

  return(cbind(grid_nodes,tij))
}



data_pp_1 <- data_pre_processing(data=dt,
                                 id="id",
                                 status="status",
                                 # setting = "survival",
                                 # covariates="covariate",
                                 # treatment=NULL,
                                 nodes = seq(0,5,1),
                                 event_time="time")


xgtrain.mf  <- model.frame(as.formula("deltaij ~ covariate+ covariate2-1+node+offset(log(tij))"),data_pp_1)
xgtrain.m  <- model.matrix(attr(xgtrain.mf,"terms"),data = data_pp_1)
xgtrain  <- xgb.DMatrix(xgtrain.m,label = data_pp_1$deltaij)
setinfo(xgtrain, "base_margin", log(data_pp_1$tij))

xgb = xgb.train(
  nrounds = 800
  , params = list(objective="count:poisson",
                  eta=.01,
                  max_depth =2,
                  alpha=.5,
                  lambda=.5)
  , data = xgtrain
)

tmp_1 <- data.frame(node = factor(sort(unique(
  data_pp_1$node
))),
tij = 1,
covariate = .5,
covariate2=factor(0,levels=c("0","1")))

tmp_1_xgb <-tmp_1%>% dplyr::mutate(deltaij=1)
xgtrain.mf  <- model.frame(as.formula("deltaij ~ covariate+ covariate2-1+node+offset(log(tij))"),tmp_1_xgb)
xgtrain.m  <- model.matrix(attr(xgtrain.mf,"terms"),data = tmp_1_xgb)
xgtest  <- xgb.DMatrix(xgtrain.m,label = tmp_1_xgb$deltaij)
setinfo(xgtest, "base_margin", log(tmp_1_xgb$tij))

pw_constant_hazard_1_xgb <- predict(xgb, newdata = xgtest)

tmp_1['sf_xgb'] <- exp(-cumsum(pw_constant_hazard_1_xgb*1))
tmp_1

sf <- function(t,lambda, beta, X1,beta1,X2){

  out <- exp(-lambda*exp(beta*X1+beta1*(X2))*t)
  return(out)
}


inverse_sf<- function(lambda, beta, X1,beta1,X2, u){

  out <- -log(u)/(lambda*exp(beta*X1+beta1*(X2)))

  return(out)
}

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

lines(seq(0,6,by=1),c(1,tmp_1$sf_xgb), col="blue",type="S")

# Fit ----
## Learners ----

l1 <- Learner_glmnet(covariates = c("covariate","covariate2"),
                     cross_validation=FALSE,
                     lambda=0)

l2 <- Learner_glmnet(covariates = c("covariate"),
                     cross_validation=T,
                     penalise_nodes=T)

## Fit Superlearner ----

learners <- list(l1,l2)

# I make different examples below. One can either use a glm or a glmnet and
# either choose to use or not to use the nodes in the meta learner.

sl <- Superlearner(data=dt,
                   learners=list(l1,l2),
                   stratified_k_fold = T,
                   id="id",
                   status="status",
                   nodes=seq(0,6,.5),
                   nfold = 3,
                   meta_learner_algorithm = "glm",
                   add_nodes_metalearner = T,
                   add_intercept_metalearner = F,
                   event_time = "time")


coef(sl$superlearner$`1`$meta_learner_fit)

log(min(dt$time))



sl <- Superlearner(data=dt,
                   learners=list(l1,l2),
                   stratified_k_fold = F,
                   id="id",
                   status="status",
                   nodes=seq(0,6,1),
                   nfold = 4,
                   meta_learner_algorithm = "glm",
                   matrix_transformation = T,
                   variable_transformation="covariate ~ sl_cut(covariate,breaks=seq(-4,4,.5))",
                   add_nodes_metalearner = T,
                   add_intercept_metalearner = F,
                   event_time = "time")



# Integrated functionalities into RiskRegression: new predict method consistent with coxph and stuff ----
out <- predict(sl,dt_val,times=c(0.01),cause=1)

library(survival)

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
                "COX" = cox),
  formula = Surv(time,status)~1,
  times = quantile(dt$time),
  data=dt_val,
  metrics="Brier",
  conf.int=FALSE)



#here it seems fine
1-predictRisk(cox,newdata=data.frame(covariate = 5,
                                   covariate2=factor("0",levels = c("0","1"))),
              times= c(.2
                       ) )


1- predictRisk(sl,newdata=data.frame(covariate = 5,
                                     covariate2=factor("0",levels = c("0","1"))),
               cause=1,
               times= c(.2))


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
  dat <- simSurvData(1000000)
}


l2 <- Learner_glmnet(covariates = c("L0"),
                  intercept= TRUE,
                  lambda=0)

l3 <- Learner_glmnet(covariates = c("L0","A0"),
                  intercept= TRUE,
                  lambda=0)

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
                   stratified_k_fold = F,
                   id="ID",
                   status="Delta",
                   nodes=seq(0,35,.5),
                   event_time =  "Time",
                   nfold=3,
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





