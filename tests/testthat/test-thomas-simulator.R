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

synthesize_td1 <- function(n){
  ## studypop <- select_studypop(studypop,how = "first")
  ## cc = synthesize(Hist(time_cvd,status_cvd)~sex+age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_Albuminuria+value_Albuminuria+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion,
  ## data = studypop,return_code = TRUE,recursive = TRUE)
  sim_model <- lava::lvm()
  sim_model <- lava::categorical(sim_model,~value_Albuminuria,labels=c('Normal','Micro','Macro'),K=3,p=c(0.839037544077992,0.121344119477287))
  lava::distribution(sim_model,~value_Motion) <- lava::binomial.lvm(p=0.68823895457374)
  lava::transform(sim_model, value_AlbuminuriaMicro~value_Albuminuria) <- function(x){1*(c(x)=='Micro')}
  lava::transform(sim_model, value_AlbuminuriaMacro~value_Albuminuria) <- function(x){1*(c(x)=='Macro')}
  lava::distribution(sim_model,~time.event.0) <- lava::coxWeibull.lvm(scale=0.0105086718102262,shape=2.2578792372401)
  lava::regression(sim_model,time.event.0~sex+age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.0296337593270627,-0.0130146695831196,-0.0125850223581346,-0.00757123095797519,0.0509424048213361,-0.00617576911315076,-0.0330093143101983,0.0763711547050251,-0.118521604228188,-0.137101392329059,-0.106325304898881,-0.038926777962038)
  lava::distribution(sim_model,~time.event.1) <- lava::coxWeibull.lvm(scale=0.00344127508109558,shape=1.1383691898171)
  lava::regression(sim_model,time.event.1~sex+age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.224377188536925,0.0412496430494417,0.00980515944727106,0.00607734573133048,0.0701435790868005,0.0143694959616522,0.458766977417934,0.760940706393913,0.480004039250349,0.442184803104845,0.311354666798778,-0.176043482185139)
  lava::distribution(sim_model,~time.event.2) <- lava::coxWeibull.lvm(scale=8.83388856676543e-05,shape=1.49710154034773)
  lava::regression(sim_model,time.event.2~sex+age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.276186036048047,0.0803874049498844,0.00678536173681053,0.00156423178426313,0.0287051547482761,0.0102472549995989,0.225523885897842,0.604004891420424,0.313118375461698,0.3289175620287,0.833815403131808,-0.410553409305998)
  sim_model <- lava::eventTime(sim_model,time_cvd ~ min(time.event.0=0, time.event.1=1, time.event.2=2),'status_cvd')
  lava::distribution(sim_model,~sex) <- lava::binomial.lvm(p=0.999496889372283)
  lava::regression(sim_model,sex~age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.00419945004262922,0.0167362987885574,-0.0234306227369814,-0.0359511647069019,0.0081536492208128,-0.303067265578663,-0.823374688355303,0.772344765719494,0.788866896108341,-0.128104780308398,0.0703784675251428)
  lava::distribution(sim_model,~age) <- lava::normal.lvm(mean=87.1360552525798,sd=7.85900822570684)
  lava::regression(sim_model,age~diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(0.13464807792832,0.108830396484175,-0.228701089912775,-0.0195655719093812,-0.509766842238701,-6.46432603931621,10.5444260864881,7.44269735464302,0.297454851105592,-0.329813521948485)
  lava::distribution(sim_model,~diabetes_duration) <- lava::normal.lvm(mean=49.1838155099675,sd=11.6488358803787)
  lava::regression(sim_model,diabetes_duration~value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(0.108793588422259,-0.573777909470319,0.0439160668983721,5.15075043729618,3.43606097921423,7.44982926983961,6.09938047892439,-0.895805379395887,-1.12258284322591)
  lava::distribution(sim_model,~value_SBP) <- lava::normal.lvm(mean=170.680834820675,sd=15.7859824102622)
  lava::regression(sim_model,value_SBP~value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(0.830042721503775,-0.0474646569265241,4.18139990647481,6.16920566367668,6.72926859500419,5.14364597888868,-2.31453562054907,1.25968644021521)
  lava::distribution(sim_model,~value_LDL) <- lava::normal.lvm(mean=1.22192152932358,sd=0.779215865602543)
  lava::regression(sim_model,value_LDL~value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(0.00573229770575866,-0.174067045226142,0.165936095360042,-0.131104050180492,-0.137783184696127,0.0731956223747194,-0.0804342224264081)
  lava::distribution(sim_model,~value_HBA1C) <- lava::normal.lvm(mean=36.2672701451412,sd=14.7522519423839)
  lava::regression(sim_model,value_HBA1C~value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(6.28998068372534,10.4761828959306,-4.25721881658791,-4.04060258145281,6.07448371233417,-2.52459417319143)
  lava::distribution(sim_model,~log2_eGFR_young) <- lava::normal.lvm(mean=-6.80484451537998,sd=0.359220470986263)
  lava::regression(sim_model,log2_eGFR_young~log2_eGFR_old+value_Smoking+value_Motion) <- c(-1.04557250788844,-0.0581689722337218,-0.027529777561561)
  lava::distribution(sim_model,~log2_eGFR_old) <- lava::normal.lvm(mean=-3.28270154306741,sd=3.26496141102303)
  lava::regression(sim_model,log2_eGFR_old~value_Smoking+value_Motion) <- c(-0.0505405444942584,-0.234580959408549)
  lava::distribution(sim_model,~value_Smoking) <- lava::binomial.lvm(p=0.334664005322723)
  lava::regression(sim_model,value_Smoking~value_Motion) <- c(-0.409037540062743)
  lava::transform(sim_model,eGFR~log2_eGFR_old+log2_eGFR_young+age) <- function(x){as.numeric(x[["age"]]<40)*100*2^{x[["log2_eGFR_young"]]}+as.numeric(x[["age"]]>40)*100*2^{x[["log2_eGFR_old"]]}}
  d <- lava::sim(sim_model,n)
  data.table::setDT(d)
  d[,log2_eGFR_old := NULL]
  d[,log2_eGFR_young := NULL]
  d[,value_AlbuminuriaMicro := NULL]
  d[,value_AlbuminuriaMacro := NULL]
  d[,time.event.0 := NULL]
  d[,time.event.1 := NULL]
  d[,time.event.2 := NULL]
  d[]
}

# load the data
d <- synthesize_td1(n = 1000)
d[,id:=1:dim(d)[1]]
d1 <- synthesize_td1(n = 300)
d1[,id:=1:dim(d1)[1]]

# define the learners

d[, (c("value_Motion", "sex", "value_Smoking")) := lapply(.SD, as.factor), .SDcols = c("value_Motion", "sex", "value_Smoking")]
d1[, (c("value_Motion", "sex", "value_Smoking")) := lapply(.SD, as.factor), .SDcols = c("value_Motion", "sex", "value_Smoking")]


l0 <- Learner_glmnet(covariates = c("value_Albuminuria",
                                 "value_Motion",
                                 "age",
                                 "sex",
                                 "value_SBP",
                                 "diabetes_duration",
                                 "value_LDL",
                                 "value_HBA1C",
                                 "value_Smoking"),
                  cross_validation=TRUE,
                  intercept=FALSE,
                  penalise_nodes=F)

l1 <- Learner_glmnet(covariates = c("value_Albuminuria",
                                 "value_Motion"),
                  cross_validation=TRUE,
                  intercept=FALSE,
                  penalise_nodes=F)

l2 <- Learner_glmnet(covariates = c("age",
                                 "sex"),
                     cross_validation=TRUE,
                     intercept=FALSE,
                     penalise_nodes=F)

l3 <- Learner_glmnet(covariates = c("diabetes_duration",
                                    "value_SBP"),
                     cross_validation=TRUE,
                     intercept=FALSE,
                     penalise_nodes=F)

l4 <- Learner_glmnet(covariates = c("value_LDL",
                                    "value_HBA1C",
                                    "value_Smoking"),
                     cross_validation=TRUE,
                     intercept=FALSE,
                     penalise_nodes=F)




learners <- list(l0,l1,l2,l3,l4)

# I make different examples below. One can either use a glm or a glmnet and
# either choose to use or not to use the nodes in the meta learner.

sl <- Superlearner(data=d,
                   learners=learners,
                   id="id",
                   status="status_cvd",
                   nfold = 3,
                   meta_learner_algorithm = "glmnet",
                   nodes=seq(1,35,1),
                   add_nodes_metalearner = TRUE,
                   add_intercept_metalearner = FALSE,
                   penalise_nodes_metalearner=TRUE,
                   event_time = "time_cvd")


sl1 <- Superlearner(data=d,
                   learners=learners,
                   id="id",
                   status="status_cvd",
                   nfold = 3,
                   meta_learner_algorithm = "glmnet",
                   nodes=seq(1,35,.5),
                   add_nodes_metalearner = TRUE,
                   add_intercept_metalearner = FALSE,
                   penalise_nodes_metalearner=TRUE,
                   event_time = "time_cvd")

sl2 <- Superlearner(data=d,
                   learners=learners,
                   id="id",
                   status="status_cvd",
                   nfold = 5,
                   meta_learner_algorithm = "glmnet",
                   nodes=seq(1,35,1),
                   add_nodes_metalearner = TRUE,
                   add_intercept_metalearner = FALSE,
                   penalise_nodes_metalearner=TRUE,
                   event_time = "time_cvd")


l0 <- Learner_glmnet(covariates = c("value_Albuminuria",
                                 "value_Motion",
                                 "age",
                                 "sex",
                                 "value_SBP",
                                 "diabetes_duration",
                                 "value_LDL",
                                 "value_HBA1C",
                                 "value_Smoking"),
                  lambda=0,
                  intercept=FALSE,
                  penalise_nodes=T)

l1 <- Learner_glmnet(covariates = c("value_Albuminuria",
                                 "value_Motion",
                                 "age",
                                 "sex",
                                 "value_SBP",
                                 "diabetes_duration",
                                 "value_LDL",
                                 "value_HBA1C",
                                 "value_Smoking"),
                  cross_validation=TRUE,
                  intercept=FALSE,
                  penalise_nodes=T)

learners <- list(l0,l1)

sl3 <- Superlearner(data=d,
                    learners=learners,
                    id="id",
                    status="status_cvd",
                    nfold = 3,
                    meta_learner_algorithm = "glm",
                    nodes=seq(1,35,1),
                    add_nodes_metalearner = TRUE,
                    add_intercept_metalearner = FALSE,
                    penalise_nodes_metalearner=TRUE,
                    event_time = "time_cvd")




l0 <- Learner_glmnet(covariates = c("value_Albuminuria",
                                    "value_Motion",
                                    "age",
                                    "sex",
                                    "value_SBP",
                                    "diabetes_duration",
                                    "value_LDL",
                                    "value_HBA1C",
                                    "value_Smoking"),
                     lambda=0,

                     intercept=FALSE)

l1 <- Learner_glmnet(covariates = c("value_Albuminuria",
                                    "value_Motion"),
                     lambda=0,
                     intercept=FALSE)

l2 <- Learner_glmnet(covariates = c("age",
                                    "sex"),
                     lambda=0,
                     intercept=FALSE)

l3 <- Learner_glmnet(covariates = c("diabetes_duration",
                                    "value_SBP"),
                     lambda=0,
                     intercept=FALSE)

l4 <- Learner_glmnet(covariates = c("value_LDL",
                                    "value_HBA1C",
                                    "value_Smoking"),
                     lambda=0,
                     intercept=FALSE)




learners <- list(l0,l1,l2,l3,l4)

sl4 <- Superlearner(data=d,
                    learners=learners,
                    id="id",
                    status="status_cvd",
                    nfold = 5,
                    meta_learner_algorithm = "glm",
                    nodes=seq(1,35,1),
                    add_nodes_metalearner = TRUE,
                    add_intercept_metalearner = FALSE,
                    penalise_nodes_metalearner=FALSE,
                    event_time = "time_cvd")

library(survival)
library(prodlim)

f=CSC(Hist(time_cvd,status_cvd)~.,data=d[, !"id"])
# f1=CSC(Hist(time_cvd,status_cvd)~diabetes_duration,data=d)
# x=Score(list("3 variables"=f1,"duration"=f),data=d,split.method = "cv10",formula=Hist(time_cvd,status_cvd)~1,times=5,se.fit = FALSE,contrasts = FALSE)
# x
# predict(sl,
#         newdata = d[1,],
#         times=1,
#         cause = 1)



numbers <- quantile(d$time_cvd)
vector <- d$time_cvd  # Example vector

closest_numbers <- sapply(numbers, function(x) vector[which.min(abs(vector - x))])
# print(closest_numbers)
x=Score(
  object = list("SL (3V, step 1)" = sl,
                "SL (3V, step 0.5)" = sl1,
                "SL (5V, step 1)" = sl2,
                "SL (3V, glm+glmnet, step 1)"=sl3,
                "SL (3V, lambas 0, step 1)"=sl4,
                "f"=f),
  formula = Hist(time_cvd,status_cvd)~1,
  times = closest_numbers,
  data=d1,
  metrics="Brier",
  cause=1,
  conf.int=FALSE)


