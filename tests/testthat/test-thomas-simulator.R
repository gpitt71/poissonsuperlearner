# Import libraries ----
library(tmlensemble)
library(survival)
library(data.table)
library(reticulate)
library(keras)
library(xgboost)
library(glmnet)
library(rpart)
library(ParBayesianOptimization)


find_true_age_curve2 <- function(n,keep, value_age=50){

  sim_model <- lava::lvm()
  sim_model <- lava::categorical(sim_model,~value_Albuminuria,labels=c('Normal','Micro','Macro'),K=3,p=c(0.839037544077992,0.121344119477287))
  lava::distribution(sim_model,~value_Motion) <- lava::binomial.lvm(p=0.68823895457374)
  lava::transform(sim_model, value_AlbuminuriaMicro~value_Albuminuria) <- function(x){1*(c(x)=='Micro')}
  lava::transform(sim_model, value_AlbuminuriaMacro~value_Albuminuria) <- function(x){1*(c(x)=='Macro')}
  lava::distribution(sim_model,~time.event.0) <- lava::coxWeibull.lvm(scale=0.0105086718102262,shape=2.2578792372401)
  lava::regression(sim_model,time.event.0~sex+age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.0296337593270627,-0.0130146695831196,-0.0125850223581346,-0.00757123095797519,0.0509424048213361,-0.00617576911315076,-0.0330093143101983,0.0763711547050251,-0.118521604228188,-0.137101392329059,-0.106325304898881,-0.038926777962038)
  lava::distribution(sim_model,~time.event.1) <- lava::coxWeibull.lvm(scale=0.00344127508109558,shape=1.1383691898171)


  lava::regression(sim_model,time.event.1~sex+age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.224377188536925,0.0412496430494417,0.00980515944727106,0.00607734573133048,0.0701435790868005,0.0143694959616522,0.458766977417934,0.760940706393913,0.480004039250349,0.442184803104845,0.311354666798778,-0.176043482185139)
  # }
  lava::distribution(sim_model,~time.event.2) <- lava::coxWeibull.lvm(scale=8.83388856676543e-05,shape=1.49710154034773)



  lava::regression(sim_model,time.event.2~sex+age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.276186036048047,0.0803874049498844,0.00678536173681053,0.00156423178426313,0.0287051547482761,0.0102472549995989,0.225523885897842,0.604004891420424,0.313118375461698,0.3289175620287,0.833815403131808,-0.410553409305998)
  # }
  sim_model <- lava::eventTime(sim_model,time_cvd ~ min(time.event.0=0, time.event.1=1, time.event.2=2),'status_cvd')
  sim_model <- lava::eventTime(sim_model,uncensored_time_cvd ~ min(time.event.1=1, time.event.2=2),'uncensored_status_cvd')
  lava::distribution(sim_model,~sex) <- lava::binomial.lvm(p=0.999496889372283)
  lava::regression(sim_model,sex~age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.00419945004262922,0.0167362987885574,-0.0234306227369814,-0.0359511647069019,0.0081536492208128,-0.303067265578663,-0.823374688355303,0.772344765719494,0.788866896108341,-0.128104780308398,0.0703784675251428)

  lava::distribution(sim_model, ~age) <- value_age

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
  # if (only_age) {
  #   d[, age:=value_age]
  # }
  d[,log2_eGFR_old := NULL]
  d[,log2_eGFR_young := NULL]
  d[,value_AlbuminuriaMicro := NULL]
  d[,value_AlbuminuriaMacro := NULL]
  d[,time.event.0 := NULL]
  d[,time.event.1 := NULL]
  d[,time.event.2 := NULL]
  if (length(keep)>0)
    d[,keep,with=FALSE]
  else d[]
}

# true_curve <- sapply(seq(20, 90, length=100), function(age) {
#   find_true_age_curve2(n=5e5, value_age=age, keep=c("uncensored_time_cvd","uncensored_status_cvd","age"))[, mean(1*(uncensored_time_cvd<=5&uncensored_status_cvd==1))]
# })

synthesize_td1 <- function(n,keep){
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
  # if (only_age) {
  #   lava::regression(sim_model,time.event.1~1) <- c(-0.224377188536925*value_age,0.0412496430494417,0.00980515944727106,0.00607734573133048,0.0701435790868005,0.0143694959616522,0.458766977417934,0.760940706393913,0.480004039250349,0.442184803104845,0.311354666798778,-0.176043482185139)[1]
  # } else {
    lava::regression(sim_model,time.event.1~sex+age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.224377188536925,0.0412496430494417,0.00980515944727106,0.00607734573133048,0.0701435790868005,0.0143694959616522,0.458766977417934,0.760940706393913,0.480004039250349,0.442184803104845,0.311354666798778,-0.176043482185139)
  # }
  lava::distribution(sim_model,~time.event.2) <- lava::coxWeibull.lvm(scale=8.83388856676543e-05,shape=1.49710154034773)

  # if (only_age) {
  #   lava::regression(sim_model,time.event.2~1) <- c(-0.276186036048047*value_age,0.0803874049498844,0.00678536173681053,0.00156423178426313,0.0287051547482761,0.0102472549995989,0.225523885897842,0.604004891420424,0.313118375461698,0.3289175620287,0.833815403131808,-0.410553409305998)[1]
  # } else {
    lava::regression(sim_model,time.event.2~sex+age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.276186036048047,0.0803874049498844,0.00678536173681053,0.00156423178426313,0.0287051547482761,0.0102472549995989,0.225523885897842,0.604004891420424,0.313118375461698,0.3289175620287,0.833815403131808,-0.410553409305998)
  # }
  sim_model <- lava::eventTime(sim_model,time_cvd ~ min(time.event.0=0, time.event.1=1, time.event.2=2),'status_cvd')
  sim_model <- lava::eventTime(sim_model,uncensored_time_cvd ~ min(time.event.1=1, time.event.2=2),'uncensored_status_cvd')
  lava::distribution(sim_model,~sex) <- lava::binomial.lvm(p=0.999496889372283)
  lava::regression(sim_model,sex~age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.00419945004262922,0.0167362987885574,-0.0234306227369814,-0.0359511647069019,0.0081536492208128,-0.303067265578663,-0.823374688355303,0.772344765719494,0.788866896108341,-0.128104780308398,0.0703784675251428)
  # if (only_age) {
  #   lava::distribution(sim_model,~age) <- value_age
  # } else {
  lava::distribution(sim_model,~age) <- lava::normal.lvm(mean=87.1360552525798,sd=7.85900822570684)
  # }
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
  # if (only_age) {
  #   d[, age:=value_age]
  # }
  d[,log2_eGFR_old := NULL]
  d[,log2_eGFR_young := NULL]
  d[,value_AlbuminuriaMicro := NULL]
  d[,value_AlbuminuriaMacro := NULL]
  d[,time.event.0 := NULL]
  d[,time.event.1 := NULL]
  d[,time.event.2 := NULL]
  if (length(keep)>0)
  d[,keep,with=FALSE]
  else d[]
}


true_curve <-fread("C:\\Users\\pwt887\\Downloads\\true_curve_predicted_risk1e6.txt")

d <- synthesize_td1(n = 100000,keep=c("uncensored_time_cvd","uncensored_status_cvd","age"))
d[,event_cvd_5:=1*(uncensored_time_cvd<=5&uncensored_status_cvd==1)]
library(rms)
library(Publish)

d[,id:=1:nrow(d)]
f=lrm(event_cvd_5~age,data=d)
splinePlot.lrm(f,xvar="age",xvalues=seq(20,90,1))
library(riskRegression)
predictRisk(f,newdata=data.frame(age=c(30,100)))


l1 <- Learner_gam(covariates = c("s(age)"),
                  intercept=FALSE)



out <- Superlearner(d,
                    id="id",
                    status="uncensored_status_cvd",
                    stratified_k_fold=FALSE,
                    event_time = "uncensored_time_cvd",
                    learners=list(l1),
                    # nodes= seq(1,35,1.749019),
                    nodes = seq(1,1500,10), #1.779275
                    meta_learner_algorithm = "glmnet",
                    add_nodes_metalearner=TRUE,
                    add_intercept_metalearner=FALSE,
                    nested_cross_validation_meta_learner = TRUE,
                    penalise_nodes_metalearner=TRUE,
                    # lambda=0.6772677,
                    # alpha=0.9274484,
                    nfold = 3
)


l2 <- Learner_glmnet(covariates = c("age"),
                     cross_validation=FALSE,
                     lambda=0,
                     # alpha=0.2302469,
                     intercept=FALSE,
                     penalise_nodes=T)

out1 <- Superlearner(d,
                    id="id",
                    status="uncensored_status_cvd",
                    stratified_k_fold=FALSE,
                    event_time = "uncensored_time_cvd",
                    learners=list(l2),
                    # nodes= seq(1,35,1.749019),
                    nodes = seq(1,1500,1), #1.779275
                    meta_learner_algorithm = "glmnet",
                    add_nodes_metalearner=TRUE,
                    add_intercept_metalearner=FALSE,
                    nested_cross_validation_meta_learner = TRUE,
                    penalise_nodes_metalearner=TRUE,
                    # lambda=0.6772677,
                    # alpha=0.9274484,
                    nfold = 10
)

# l3 <- Learner_glmnet(covariates = c("age"),
#                      cross_validation=TRUE,
#                      # lambda=c(0,.1,.5,1),
#                      # alpha=0.2302469,
#                      intercept=FALSE,
#                      penalise_nodes=T)
#
# out2 <- Superlearner(d,
#                      id="id",
#                      status="uncensored_status_cvd",
#                      stratified_k_fold=FALSE,
#                      event_time = "uncensored_time_cvd",
#                      learners=list(l3),
#                      # nodes= seq(1,35,1.749019),
#                      nodes = seq(1,1000,1), #1.779275
#                      meta_learner_algorithm = "glmnet",
#                      add_nodes_metalearner=TRUE,
#                      add_intercept_metalearner=FALSE,
#                      nested_cross_validation_meta_learner = TRUE,
#                      penalise_nodes_metalearner=TRUE,
#                      # lambda=0.6772677,
#                      # alpha=0.9274484,
#                      nfold = 3
# )
#
#
#
# l4 <- Learner_glmnet(covariates = c("age"),
#                      cross_validation=F,
#                      lambda=c(0),
#                      # alpha=0.2302469,
#                      intercept=FALSE,
#                      penalise_nodes=T)
#
# out3 <- Superlearner(d,
#                      id="id",
#                      status="uncensored_status_cvd",
#                      stratified_k_fold=FALSE,
#                      event_time = "uncensored_time_cvd",
#                      learners=list(l4),
#                      # nodes= seq(1,35,1.749019),
#                      nodes = seq(1,1000,.5), #1.779275
#                      meta_learner_algorithm = "glmnet",
#                      add_nodes_metalearner=TRUE,
#                      add_intercept_metalearner=FALSE,
#                      nested_cross_validation_meta_learner = TRUE,
#                      penalise_nodes_metalearner=TRUE,
#                      # lambda=0.6772677,
#                      # alpha=0.9274484,
#                      nfold = 3
# )


sl_out <- predict(out,,newdata=data.frame(age=c(30,90)),
                  times=5)

sl_out1 <- predict(out1,,newdata=data.frame(age=c(30,60)),
                  times=5)

# sl_out2 <- predict(out2,,newdata=data.frame(age=c(30,60)),
#                    times=5)




pr<-predictRisk(out,,newdata=data.frame(age=seq(20, 75, 1)),
        times=5,cause = 1)

pr2<-predictRisk(out1,,newdata=data.frame(age=seq(20, 75, 1)),
                times=5,cause = 1)

# pr3<-predictRisk(out2,newdata=data.frame(age=seq(20,90,1)),
#                  times=5,cause = 1)
# pr4<-predictRisk(out3,newdata=data.frame(age=seq(20,90,1)),
#                  times=5,cause = 1)



# coef(out$superlearner$learners_fit[[1]])

Publish::splinePlot.lrm(f,xvar="age",xvalues=seq(20, 75, length=100))
lines(seq(20, 90, length=100),
      true_curve[['V1']],
      col="black")
lines(seq(20, 75,1),
      pr,
      col="blue")
lines(seq(20, 75, 1),
      pr2,
      col="red")
# lines(seq(20,90,1),
#       pr3,
#       col="green")
# lines(seq(20,90,1),
#       pr4,
#       col="yellow")

legend("topleft",                       # position of the legend
       legend = c("Gam, k=1", "Glmnet - lambda = 0","true risk"),   # labels
       col = c("blue","red", "black","yellow"),  # matching colors
       lty = 1,                          # line type (solid)
       bty = "n")                        # box type: "n" means no border


# compute true va -----

calculate_true_va <- function(n,keep,plugin_probs=c(1e-12,1e-12)){

  sim_model <- lava::lvm()
  sim_model <- lava::categorical(sim_model,~value_Albuminuria,labels=c('Normal','Micro','Macro'),K=3,p=plugin_probs)
  lava::distribution(sim_model,~value_Motion) <- lava::binomial.lvm(p=0.68823895457374)
  lava::transform(sim_model, value_AlbuminuriaMicro~value_Albuminuria) <- function(x){1*(c(x)=='Micro')}
  lava::transform(sim_model, value_AlbuminuriaMacro~value_Albuminuria) <- function(x){1*(c(x)=='Macro')}
  lava::distribution(sim_model,~time.event.0) <- lava::coxWeibull.lvm(scale=0.0105086718102262,shape=2.2578792372401)
  lava::regression(sim_model,time.event.0~sex+age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.0296337593270627,-0.0130146695831196,-0.0125850223581346,-0.00757123095797519,0.0509424048213361,-0.00617576911315076,-0.0330093143101983,0.0763711547050251,-0.118521604228188,-0.137101392329059,-0.106325304898881,-0.038926777962038)
  lava::distribution(sim_model,~time.event.1) <- lava::coxWeibull.lvm(scale=0.00344127508109558,shape=1.1383691898171)
  # if (only_age) {
  #   lava::regression(sim_model,time.event.1~1) <- c(-0.224377188536925*value_age,0.0412496430494417,0.00980515944727106,0.00607734573133048,0.0701435790868005,0.0143694959616522,0.458766977417934,0.760940706393913,0.480004039250349,0.442184803104845,0.311354666798778,-0.176043482185139)[1]
  # } else {
  lava::regression(sim_model,time.event.1~sex+age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.224377188536925,0.0412496430494417,0.00980515944727106,0.00607734573133048,0.0701435790868005,0.0143694959616522,0.458766977417934,0.760940706393913,0.480004039250349,0.442184803104845,0.311354666798778,-0.176043482185139)
  # }
  lava::distribution(sim_model,~time.event.2) <- lava::coxWeibull.lvm(scale=8.83388856676543e-05,shape=1.49710154034773)

  # if (only_age) {
  #   lava::regression(sim_model,time.event.2~1) <- c(-0.276186036048047*value_age,0.0803874049498844,0.00678536173681053,0.00156423178426313,0.0287051547482761,0.0102472549995989,0.225523885897842,0.604004891420424,0.313118375461698,0.3289175620287,0.833815403131808,-0.410553409305998)[1]
  # } else {
  lava::regression(sim_model,time.event.2~sex+age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.276186036048047,0.0803874049498844,0.00678536173681053,0.00156423178426313,0.0287051547482761,0.0102472549995989,0.225523885897842,0.604004891420424,0.313118375461698,0.3289175620287,0.833815403131808,-0.410553409305998)
  # }
  sim_model <- lava::eventTime(sim_model,time_cvd ~ min(time.event.0=0, time.event.1=1, time.event.2=2),'status_cvd')
  sim_model <- lava::eventTime(sim_model,uncensored_time_cvd ~ min(time.event.1=1, time.event.2=2),'uncensored_status_cvd')
  lava::distribution(sim_model,~sex) <- lava::binomial.lvm(p=0.999496889372283)
  lava::regression(sim_model,sex~age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.00419945004262922,0.0167362987885574,-0.0234306227369814,-0.0359511647069019,0.0081536492208128,-0.303067265578663,-0.823374688355303,0.772344765719494,0.788866896108341,-0.128104780308398,0.0703784675251428)
  # if (only_age) {
  #   lava::distribution(sim_model,~age) <- value_age
  # } else {
  lava::distribution(sim_model,~age) <- lava::normal.lvm(mean=87.1360552525798,sd=7.85900822570684)
  # }
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
  # if (only_age) {
  #   d[, age:=value_age]
  # }
  d[,log2_eGFR_old := NULL]
  d[,log2_eGFR_young := NULL]
  d[,value_AlbuminuriaMicro := NULL]
  d[,value_AlbuminuriaMacro := NULL]
  d[,time.event.0 := NULL]
  d[,time.event.1 := NULL]
  d[,time.event.2 := NULL]
  if (length(keep)>0)
    d[,keep,with=FALSE]
  else d[]
}

true_d_macro <- calculate_true_va(n=1e6, keep=c("uncensored_time_cvd","uncensored_status_cvd","value_Albuminuria"))[, mean(1*(uncensored_time_cvd<=5&uncensored_status_cvd==1))]

#'Normal','Micro','Macro'

true_d_micro <- calculate_true_va(n=1e6, plugin_probs = c(1e-12,0.999999), keep=c("uncensored_time_cvd","uncensored_status_cvd","value_Albuminuria"))[, mean(1*(uncensored_time_cvd<=5&uncensored_status_cvd==1))]
true_d_normal <- calculate_true_va(n=1e6, plugin_probs = c(0.999999,1e-12), keep=c("uncensored_time_cvd","uncensored_status_cvd","value_Albuminuria"))[, mean(1*(uncensored_time_cvd<=5&uncensored_status_cvd==1))]



tmp <- calculate_true_va(n=1e2, keep=c("uncensored_time_cvd","uncensored_status_cvd","value_Albuminuria"))


synthesize_td1 <- function(n,keep){
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
  # if (only_age) {
  #   lava::regression(sim_model,time.event.1~1) <- c(-0.224377188536925*value_age,0.0412496430494417,0.00980515944727106,0.00607734573133048,0.0701435790868005,0.0143694959616522,0.458766977417934,0.760940706393913,0.480004039250349,0.442184803104845,0.311354666798778,-0.176043482185139)[1]
  # } else {
  lava::regression(sim_model,time.event.1~sex+age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.224377188536925,0.0412496430494417,0.00980515944727106,0.00607734573133048,0.0701435790868005,0.0143694959616522,0.458766977417934,0.760940706393913,0.480004039250349,0.442184803104845,0.311354666798778,-0.176043482185139)
  # }
  lava::distribution(sim_model,~time.event.2) <- lava::coxWeibull.lvm(scale=8.83388856676543e-05,shape=1.49710154034773)

  # if (only_age) {
  #   lava::regression(sim_model,time.event.2~1) <- c(-0.276186036048047*value_age,0.0803874049498844,0.00678536173681053,0.00156423178426313,0.0287051547482761,0.0102472549995989,0.225523885897842,0.604004891420424,0.313118375461698,0.3289175620287,0.833815403131808,-0.410553409305998)[1]
  # } else {
  lava::regression(sim_model,time.event.2~sex+age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.276186036048047,0.0803874049498844,0.00678536173681053,0.00156423178426313,0.0287051547482761,0.0102472549995989,0.225523885897842,0.604004891420424,0.313118375461698,0.3289175620287,0.833815403131808,-0.410553409305998)
  # }
  sim_model <- lava::eventTime(sim_model,time_cvd ~ min(time.event.0=0, time.event.1=1, time.event.2=2),'status_cvd')
  sim_model <- lava::eventTime(sim_model,uncensored_time_cvd ~ min(time.event.1=1, time.event.2=2),'uncensored_status_cvd')
  lava::distribution(sim_model,~sex) <- lava::binomial.lvm(p=0.999496889372283)
  lava::regression(sim_model,sex~age+diabetes_duration+value_SBP+value_LDL+value_HBA1C+value_AlbuminuriaMicro+value_AlbuminuriaMacro+log2_eGFR_young+log2_eGFR_old+value_Smoking+value_Motion) <- c(-0.00419945004262922,0.0167362987885574,-0.0234306227369814,-0.0359511647069019,0.0081536492208128,-0.303067265578663,-0.823374688355303,0.772344765719494,0.788866896108341,-0.128104780308398,0.0703784675251428)
  # if (only_age) {
  #   lava::distribution(sim_model,~age) <- value_age
  # } else {
  lava::distribution(sim_model,~age) <- lava::normal.lvm(mean=87.1360552525798,sd=7.85900822570684)
  # }
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
  # if (only_age) {
  #   d[, age:=value_age]
  # }
  d[,log2_eGFR_old := NULL]
  d[,log2_eGFR_young := NULL]
  d[,value_AlbuminuriaMicro := NULL]
  d[,value_AlbuminuriaMacro := NULL]
  d[,time.event.0 := NULL]
  d[,time.event.1 := NULL]
  d[,time.event.2 := NULL]
  if (length(keep)>0)
    d[,keep,with=FALSE]
  else d[]
}

d <- synthesize_td1(n = 10000,keep=c("uncensored_time_cvd","uncensored_status_cvd","value_Albuminuria"))
d[,event_cvd_5:=1*(uncensored_time_cvd<=5&uncensored_status_cvd==1)]

library(rms)
library(Publish)

d[,id:=1:nrow(d)]
f=lrm(event_cvd_5~value_Albuminuria,data=d)

prop_lrm <- predict(f,data.frame(value_Albuminuria=c("Macro","Micro","Normal")), type = "fitted")



l1 <- Learner_glmnet(covariates = c("value_Albuminuria"),
                     cross_validation=TRUE,
                     # lambda=0,
                     # alpha=0.2302469,
                     intercept=FALSE,
                     penalise_nodes=T)



out <- Superlearner(d,
                    id="id",
                    status="uncensored_status_cvd",
                    stratified_k_fold=FALSE,
                    event_time = "uncensored_time_cvd",
                    learners=list(l1),
                    # nodes= seq(1,35,1.749019),
                    nodes = seq(1,1000,1), #1.779275
                    meta_learner_algorithm = "glmnet",
                    add_nodes_metalearner=TRUE,
                    add_intercept_metalearner=FALSE,
                    nested_cross_validation_meta_learner = TRUE,
                    penalise_nodes_metalearner=TRUE,
                    # lambda=0.6772677,
                    # alpha=0.9274484,
                    nfold = 3
)


l2 <- Learner_glmnet(covariates = c("value_Albuminuria"),
                     cross_validation=FALSE,
                     lambda=0,
                     # alpha=0.2302469,
                     intercept=FALSE,
                     penalise_nodes=T)

out1 <- Superlearner(d,
                     id="id",
                     status="uncensored_status_cvd",
                     stratified_k_fold=FALSE,
                     event_time = "uncensored_time_cvd",
                     learners=list(l2),
                     # nodes= seq(1,35,1.749019),
                     nodes = seq(1,1000,1), #1.779275
                     meta_learner_algorithm = "glmnet",
                     add_nodes_metalearner=TRUE,
                     add_intercept_metalearner=FALSE,
                     nested_cross_validation_meta_learner = TRUE,
                     penalise_nodes_metalearner=TRUE,
                     # lambda=0.6772677,
                     # alpha=0.9274484,
                     nfold = 10
)




sl_out <- predict(out,,newdata=data.frame(value_Albuminuria=c("Macro","Micro","Normal")),
                  times=5)

sl_out1 <- predict(out1,,newdata=data.frame(value_Albuminuria=c("Macro","Micro","Normal")),
                   times=5)

# sl_out2 <- predict(out2,,newdata=data.frame(age=c(30,60)),
#                    times=5)




pr<-predictRisk(out,,newdata=data.frame(value_Albuminuria=c("Macro","Micro","Normal")),
                times=5,cause = 1)

pr2<-predictRisk(out1,,newdata=data.frame(value_Albuminuria=c("Macro","Micro","Normal")),
                 times=5,cause = 1)



########### Comparing the glmnet models ----


# Lambda values to iterate over
lambdas <- c(0, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5, 1, 1.5, 2)
colors <- rainbow(length(lambdas))
predictions <- list()

# Create a plot based on the spline model
splinePlot.lrm(f, xvar="age", xvalues=seq(20, 90, 1))

# Loop over lambda values
for (i in seq_along(lambdas)) {
  l <- Learner_glmnet(covariates = c("age"),
                      cross_validation = FALSE,
                      lambda = lambdas[i],
                      intercept = FALSE,
                      penalise_nodes = TRUE)

  out <- Superlearner(d,
                      id = "id",
                      status = "uncensored_status_cvd",
                      stratified_k_fold = FALSE,
                      event_time = "uncensored_time_cvd",
                      learners = list(l),
                      nodes = seq(1, 1000, 10),
                      meta_learner_algorithm = "glmnet",
                      add_nodes_metalearner = TRUE,
                      add_intercept_metalearner = FALSE,
                      nested_cross_validation_meta_learner = TRUE,
                      penalise_nodes_metalearner = TRUE,
                      nfold = 3)

  pr_lambda <- predictRisk(out,
                           newdata = data.frame(age = seq(20, 90, 1)),
                           times = 5)

  predictions[[i]] <- pr_lambda

  lines(seq(20, 90, 1), pr_lambda, col = colors[i], type = "S")
}

# Add legend
legend("topleft",
       legend = paste("Glmnet - lambda =", lambdas),
       col = colors,
       lty = 1,
       bty = "n")


###############  Comparing the gam models ----

library(data.table)
library(rms)
library(Publish)
library(mgcv)
library(survival)

# Step 1: TRUE MODEL from large uncensored dataset
{set.seed(123)
d_true <- synthesize_td1(n = 100000, keep = c("uncensored_time_cvd", "uncensored_status_cvd", "age"))
d_true[, event_cvd_5 := 1 * (uncensored_time_cvd <= 5 & uncensored_status_cvd == 1)]
d_true <- d_true[uncensored_time_cvd < 1000]
f_true <- lrm(event_cvd_5 ~ age, data = d_true)}

# Step 2: Define common prediction grid and sample sizes
age_seq <- seq(20, 90, 1)
sample_sizes <- c(50, 500, 1000, 5000, 10000)
colors <- rainbow(length(sample_sizes))
sl_predictions <- list()

# Step 3: Base plot: True curve
splinePlot.lrm(f_true, xvar = "age", xvalues = age_seq)

# Step 4: Loop over sample sizes and fit GAM via Superlearner
for (i in seq_along(sample_sizes)) {
  n <- sample_sizes[i]
  set.seed(100 + n)

  # Simulate data
  d <- synthesize_td1(n = n, keep = c("uncensored_time_cvd", "uncensored_status_cvd", "age"))
  d <- d[uncensored_time_cvd < 1000]
  d[, event_cvd_5 := 1 * (uncensored_time_cvd <= 5 & uncensored_status_cvd == 1)]
  d[, id := 1:.N]

  # Define GAM learner
  l <- Learner_gam(covariates = c("s(age)"), intercept = FALSE)

  # Fit Superlearner
  sl <- Superlearner(d,
                     id = "id",
                     status = "uncensored_status_cvd",
                     event_time = "uncensored_time_cvd",
                     stratified_k_fold = FALSE,
                     learners = list(l),
                     nodes = seq(1, 1000, 10),
                     meta_learner_algorithm = "glmnet",
                     add_nodes_metalearner = TRUE,
                     add_intercept_metalearner = FALSE,
                     nested_cross_validation_meta_learner = TRUE,
                     penalise_nodes_metalearner = TRUE,
                     nfold = 3)

  # Predict 5-year risk
  pr <- predictRisk(sl, newdata = data.frame(age = age_seq), times = 5)
  sl_predictions[[i]] <- pr

  # Add to plot
  lines(age_seq, pr, col = colors[i], type = "S")
}

# Step 5: Legend
legend("topleft",
       legend = paste("GAM Superlearner, n =", sample_sizes),
       col = colors,
       lty = 1,
       bty = "n")





############### Systemic behaviour of gams ----

library(data.table)
library(rms)
library(Publish)
library(mgcv)  # for s()
library(survival)

# Step 1: TRUE MODEL from large uncensored dataset
set.seed(123)
d_true <- synthesize_td1(n = 100000, keep = c("uncensored_time_cvd", "uncensored_status_cvd", "age"))
d_true[, event_cvd_5 := 1 * (uncensored_time_cvd <= 5 & uncensored_status_cvd == 1)]
d_true <- d_true[uncensored_time_cvd < 1000]
f_true <- lrm(event_cvd_5 ~ age, data = d_true)

# Step 2: Define common prediction grid and sample sizes
age_seq <- seq(20, 90, 1)
sample_sizes <- c(1000, 5000, 10000)
colors <- rainbow(length(sample_sizes))
sl_predictions <- list()

# Step 3: Base plot: True curve
splinePlot.lrm(f_true, xvar = "age", xvalues = age_seq)

# Step 4: Loop over sample sizes and fit GAM via Superlearner
for (i in seq_along(sample_sizes)) {
  n <- sample_sizes[i]
  set.seed(100 + n)

  # Simulate data
  d <- synthesize_td1(n = n, keep = c("uncensored_time_cvd", "uncensored_status_cvd", "age"))
  d <- d[uncensored_time_cvd < 1000]
  d[, event_cvd_5 := 1 * (uncensored_time_cvd <= 5 & uncensored_status_cvd == 1)]
  d[, id := 1:.N]

  # Define GAM learner
  l <- Learner_gam(covariates = c("s(age,k=2,fx=TRUE)"), intercept = FALSE)

  # Fit Superlearner
  sl <- Superlearner(d,
                     id = "id",
                     status = "uncensored_status_cvd",
                     event_time = "uncensored_time_cvd",
                     stratified_k_fold = FALSE,
                     learners = list(l),
                     nodes = seq(1, 1000, 10),
                     meta_learner_algorithm = "glmnet",
                     add_nodes_metalearner = TRUE,
                     add_intercept_metalearner = FALSE,
                     nested_cross_validation_meta_learner = TRUE,
                     penalise_nodes_metalearner = TRUE,
                     nfold = 3)

  # Predict 5-year risk
  pr <- predictRisk(sl, newdata = data.frame(age = age_seq), times = 5)
  sl_predictions[[i]] <- pr

  # Add to plot
  lines(age_seq, pr, col = colors[i], type = "S")
}

# Step 5: Legend
legend("topleft",
       legend = paste("GAM Superlearner, n =", sample_sizes),
       col = colors,
       lty = 1,
       bty = "n")




############### Systemic behaviour of glmnets ----

library(data.table)
library(rms)
library(Publish)
library(mgcv)  # for s()
library(survival)

# Step 1: TRUE MODEL from large uncensored dataset
{set.seed(123)
d_true <- synthesize_td1(n = 500000, keep = c("uncensored_time_cvd", "uncensored_status_cvd", "age"))
d_true[, event_cvd_5 := 1 * (uncensored_time_cvd <= 5 & uncensored_status_cvd == 1)]
d_true <- d_true[uncensored_time_cvd < 1000]
f_true <- lrm(event_cvd_5 ~ age, data = d_true)}

# Step 2: Define common prediction grid and sample sizes
age_seq <- seq(20, 90, 1)
sample_sizes <- c(50,100,1000, 5000, 10000)
colors <- rainbow(length(sample_sizes))
sl_predictions <- list()

# Step 3: Base plot: True curve
splinePlot.lrm(f_true, xvar = "age", xvalues = age_seq)

# Step 4: Loop over sample sizes and fit GAM via Superlearner
for (i in seq_along(sample_sizes)) {
  n <- sample_sizes[i]
  set.seed(100 + n)

  # Simulate data
  d <- synthesize_td1(n = n, keep = c("uncensored_time_cvd", "uncensored_status_cvd", "age"))
  d <- d[uncensored_time_cvd < 1000]
  d[, event_cvd_5 := 1 * (uncensored_time_cvd <= 5 & uncensored_status_cvd == 1)]
  d[, id := 1:.N]

  l <- Learner_glmnet(covariates = c("age"),
                       cross_validation=T,
                       lambda=c(.1,.5,1),
                       alpha=0,
                       intercept=F,
                       penalise_nodes=T)

  # Fit Superlearner
  sl <- Superlearner(d,
                     id = "id",
                     status = "uncensored_status_cvd",
                     event_time = "uncensored_time_cvd",
                     stratified_k_fold = FALSE,
                     learners = list(l),
                     nodes = seq(1, 1000, 1),
                     meta_learner_algorithm = "glmnet",
                     add_nodes_metalearner = TRUE,
                     add_intercept_metalearner = FALSE,
                     nested_cross_validation_meta_learner = TRUE,
                     penalise_nodes_metalearner = TRUE,
                     nfold = 5)

  # Predict 5-year risk
  pr <- predictRisk(sl, newdata = data.frame(age = age_seq), times = 5)
  sl_predictions[[i]] <- pr

  # Add to plot
  lines(age_seq, pr, col = colors[i], type = "S")
}

# Step 5: Legend
legend("topleft",
       legend = paste("glmnet , n =", sample_sizes),
       col = colors,
       lty = 1,
       bty = "n")

###############

# load the data
{
  set.seed(1)
d <- synthesize_td1(n = 1000,keep=NULL)
d[,id:=1:dim(d)[1]]
d1 <- synthesize_td1(n = 300,keep=NULL)
d1[,id:=1:dim(d1)[1]]}

# define the learners

d[, (c("value_Motion", "sex", "value_Smoking")) := lapply(.SD, as.factor), .SDcols = c("value_Motion", "sex", "value_Smoking")]
d1[, (c("value_Motion", "sex", "value_Smoking")) := lapply(.SD, as.factor), .SDcols = c("value_Motion", "sex", "value_Smoking")]


# Select SL hyperparameters

bounds <- list(
  nfold = c(2L,10L),
  sequence_step=c(.1,2),
  possible_algorithms=c(1L,2L)
)


scoringFunction <- function(nfold,
                            sequence_step,
                            possible_algorithms) {


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
                       penalise_nodes=T#,
                       # lambda=lambda_l1,
                       # alpha=alpha_l1

                       )

  learners <- list(l0,l1)

  possible_algorithms <- c("glm","glmnet")



  out <- Superlearner.cv(d,
                         id="id",
                         status="status_cvd",
                         stratified_k_fold=FALSE,
                         event_time = "time_cvd",
                         learners=learners,
                         n_cross_validation_folds=5,
                         nodes = seq(1,35,sequence_step),
                         meta_learner_algorithm = "glmnet",
                         add_nodes_metalearner=TRUE,
                         add_intercept_metalearner=FALSE,
                         penalise_nodes_metalearner=TRUE,
                         nested_cross_validation_meta_learner=TRUE,
                         nfold = nfold
                         # lambda=lambda_sl,
                         # alpha=alpha_sl

                         )

  return(
    list(
      train_likelihood = out$train_likelihood,
      Score = out$validation_likelihood
    )
  )
}


optObj <- ParBayesianOptimization::bayesOpt(
  FUN = scoringFunction
  , bounds = bounds
  , initPoints = length(bounds)+1
  , iters.n = 8
  , iters.k = 6
)

optObj$scoreSummary[Score==min(Score)]


# Epoch Iteration nfold sequence_step possible_algorithms gpUtility acqOptimum
# <num>     <int> <num>         <num>               <num>     <num>     <lgcl>
#   1:     0         1     6      1.882441                   1        NA      FALSE
# inBounds Elapsed    Score train_likelihood errorMessage
# <lgcl>   <num>    <num>            <num>       <lgcl>
#   1:     TRUE     132 600.7633         2462.772           NA

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
                     # lambda=0.08483658,
                     # alpha=0.2302469,
                     intercept=FALSE,
                     penalise_nodes=T)



      # 1.749019 0.08483658 0.2302469 0.6772677 0.9274484

plot(optObj)
optObjSimp <- addIterations(optObj,4,verbose=FALSE)

{
  set.seed(1)
  out <- Superlearner(d,
                       id="id",
                       status="status_cvd",
                       stratified_k_fold=FALSE,
                       event_time = "time_cvd",
                       learners=list(l0,l1),
                    # nodes= seq(1,35,1.749019),
                       nodes = seq(1,35,1.882441), #1.779275
                       meta_learner_algorithm = "glmnet",
                       add_nodes_metalearner=TRUE,
                       add_intercept_metalearner=FALSE,
                      nested_cross_validation_meta_learner = TRUE,
                       penalise_nodes_metalearner=TRUE,
                    # lambda=0.6772677,
                    # alpha=0.9274484,
                       nfold = 6
                    )}




f=CSC(Hist(time_cvd,status_cvd)~1+age + sex + diabetes_duration + value_SBP + value_LDL +
        value_HBA1C + I(value_Albuminuria == "Micro") + I(value_Albuminuria == "Macro")+
        I(log2(eGFR) * (age < 40)) + I(log2(eGFR) * (age >= 40)) +
        value_Smoking + I(value_Motion==1),data=d[, !"id"])


numbers <- quantile(d$time_cvd)
vector <- d$time_cvd  # Example vector

closest_numbers <- sapply(numbers, function(x) vector[which.min(abs(vector - x))])
# print(closest_numbers)
x=Score(
  object = list("SL (Bayes CV)" = out,
                "COX"=f),
  formula = Hist(time_cvd,status_cvd)~1,
  times = closest_numbers,
  data=d1,
  metrics="Brier",
  cause=1,
  conf.int=FALSE)
x
optObjSimp <- addIterations(optObj,2,verbose=FALSE)

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
                "COX"=f),
  formula = Hist(time_cvd,status_cvd)~1,
  times = closest_numbers,
  data=d1,
  metrics="Brier",
  cause=1,
  conf.int=FALSE)

# fixed covariates experiment ----

find_true_age_curve <- function(n, keep, value_age = 50) {
  sim_model <- lava::lvm()

  sim_model <- lava::categorical(sim_model, ~value_Albuminuria,
                                 labels = c("Normal", "Micro", "Macro"),
                                 K = 3, p = c(1e-12, 1e-12))
  lava::transform(sim_model, value_AlbuminuriaMicro ~ value_Albuminuria) <- function(x) (x == "Micro") * 1
  lava::transform(sim_model, value_AlbuminuriaMacro ~ value_Albuminuria) <- function(x) (x == "Macro") * 1

  lava::distribution(sim_model, ~value_Motion) <- lava::binomial.lvm(p = 1)

  # Event models
  lava::distribution(sim_model, ~time.event.0) <- lava::coxWeibull.lvm(scale = 0.0105086718102262, shape = 2.2578792372401)
  lava::regression(sim_model, time.event.0 ~ sex + age + diabetes_duration + value_SBP + value_LDL + value_HBA1C +
                     value_AlbuminuriaMicro + value_AlbuminuriaMacro + log2_eGFR_young + log2_eGFR_old +
                     value_Smoking + value_Motion) <-
    c(-0.02963376, -0.01301467, -0.01258502, -0.00757123, 0.05094240, -0.00617577,
      -0.03300931, 0.07637115, -0.11852160, -0.13710139, -0.10632530, -0.03892678)

  lava::distribution(sim_model, ~time.event.1) <- lava::coxWeibull.lvm(scale = 0.00344127508109558, shape = 1.1383691898171)
  lava::regression(sim_model, time.event.1 ~ sex + age + diabetes_duration + value_SBP + value_LDL + value_HBA1C +
                     value_AlbuminuriaMicro + value_AlbuminuriaMacro + log2_eGFR_young + log2_eGFR_old +
                     value_Smoking + value_Motion) <-
    c(-0.22437719, 0.04124964, 0.00980516, 0.00607735, 0.07014358, 0.01436950,
      0.45876698, 0.76094071, 0.48000404, 0.44218480, 0.31135467, -0.17604348)

  lava::distribution(sim_model, ~time.event.2) <- lava::coxWeibull.lvm(scale = 8.83388856676543e-05, shape = 1.49710154034773)
  lava::regression(sim_model, time.event.2 ~ sex + age + diabetes_duration + value_SBP + value_LDL + value_HBA1C +
                     value_AlbuminuriaMicro + value_AlbuminuriaMacro + log2_eGFR_young + log2_eGFR_old +
                     value_Smoking + value_Motion) <-
    c(-0.27618604, 0.08038740, 0.00678536, 0.00156423, 0.02870515, 0.01024725,
      0.22552389, 0.60400489, 0.31311838, 0.32891756, 0.83381540, -0.41055341)

  sim_model <- lava::eventTime(sim_model, time_cvd ~ min(time.event.0 = 0, time.event.1 = 1, time.event.2 = 2), "status_cvd")
  sim_model <- lava::eventTime(sim_model, uncensored_time_cvd ~ min(time.event.1 = 1, time.event.2 = 2), "uncensored_status_cvd")

  # Fix age deterministically
  # lava::distribution(sim_model, ~age) <- function(n) rep(value_age, n)
  lava::distribution(sim_model, ~age) <- value_age#lava::numeric.lvm(value = value_age)

  # All other variables: set deterministic distributions (sd = 0) + regression + zero residual variance
  lava::distribution(sim_model, ~diabetes_duration) <- lava::normal.lvm(mean = 49.18, sd = 0)
  lava::regression(sim_model, diabetes_duration ~ value_SBP + value_LDL + value_HBA1C + value_AlbuminuriaMicro +
                     value_AlbuminuriaMacro + log2_eGFR_young + log2_eGFR_old + value_Smoking + value_Motion) <-
    c(0.10879359, -0.57377791, 0.04391607, 5.15075044, 3.43606098, 7.44982927, 6.09938048, -0.89580538, -1.12258284)
  lava::variance(sim_model, ~diabetes_duration) <- 0

  lava::distribution(sim_model, ~value_SBP) <- lava::normal.lvm(mean = 170.68, sd = 0)
  lava::regression(sim_model, value_SBP ~ value_LDL + value_HBA1C + value_AlbuminuriaMicro + value_AlbuminuriaMacro +
                     log2_eGFR_young + log2_eGFR_old + value_Smoking + value_Motion) <-
    c(0.83004272, -0.04746466, 4.18139991, 6.16920566, 6.72926860, 5.14364598, -2.31453562, 1.25968644)
  lava::variance(sim_model, ~value_SBP) <- 0

  lava::distribution(sim_model, ~value_LDL) <- lava::normal.lvm(mean = 1.22, sd = 0)
  lava::regression(sim_model, value_LDL ~ value_HBA1C + value_AlbuminuriaMicro + value_AlbuminuriaMacro +
                     log2_eGFR_young + log2_eGFR_old + value_Smoking + value_Motion) <-
    c(0.00573230, -0.17406705, 0.16593610, -0.13110405, -0.13778318, 0.07319562, -0.08043422)
  lava::variance(sim_model, ~value_LDL) <- 0

  lava::distribution(sim_model, ~value_HBA1C) <- lava::normal.lvm(mean = 36.27, sd = 0)
  lava::regression(sim_model, value_HBA1C ~ value_AlbuminuriaMicro + value_AlbuminuriaMacro +
                     log2_eGFR_young + log2_eGFR_old + value_Smoking + value_Motion) <-
    c(6.28998068, 10.47618290, -4.25721882, -4.04060258, 6.07448371, -2.52459417)
  lava::variance(sim_model, ~value_HBA1C) <- 0

  lava::distribution(sim_model, ~log2_eGFR_young) <- lava::normal.lvm(mean = -6.80, sd = 0)
  lava::regression(sim_model, log2_eGFR_young ~ log2_eGFR_old + value_Smoking + value_Motion) <-
    c(-1.04557251, -0.05816897, -0.02752978)
  lava::variance(sim_model, ~log2_eGFR_young) <- 0

  lava::distribution(sim_model, ~log2_eGFR_old) <- lava::normal.lvm(mean = -3.28, sd = 0)
  lava::regression(sim_model, log2_eGFR_old ~ value_Smoking + value_Motion) <-
    c(-0.05054054, -0.23458096)
  lava::variance(sim_model, ~log2_eGFR_old) <- 0

  lava::distribution(sim_model, ~value_Smoking) <- lava::binomial.lvm(p = 0)
  lava::regression(sim_model, value_Smoking ~ value_Motion) <- -0.40903754
  lava::variance(sim_model, ~value_Smoking) <- 0

  lava::distribution(sim_model, ~sex) <- lava::binomial.lvm(p = 1)
  lava::regression(sim_model, sex ~ age + diabetes_duration + value_SBP + value_LDL + value_HBA1C +
                     value_AlbuminuriaMicro + value_AlbuminuriaMacro + log2_eGFR_young + log2_eGFR_old +
                     value_Smoking + value_Motion) <-
    c(-0.00419945, 0.01673630, -0.02343062, -0.03595116, 0.00815365, -0.30306727,
      -0.82337469, 0.77234477, 0.78886690, -0.12810478, 0.07037847)
  lava::variance(sim_model, ~sex) <- 0

  lava::transform(sim_model, eGFR ~ log2_eGFR_old + log2_eGFR_young + age) <- function(x) {
    ifelse(x[["age"]] < 40,
           100 * 2^(x[["log2_eGFR_young"]]),
           100 * 2^(x[["log2_eGFR_old"]]))
  }

  d <- lava::sim(sim_model, n)

  data.table::setDT(d)
  d[, `:=`(log2_eGFR_old = NULL,
           log2_eGFR_young = NULL,
           value_AlbuminuriaMicro = NULL,
           value_AlbuminuriaMacro = NULL,
           time.event.0 = NULL,
           time.event.1 = NULL,
           time.event.2 = NULL)]

  if (length(keep) > 0)
    return(d[, ..keep])
  else
    return(d[])
}

