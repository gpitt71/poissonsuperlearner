library(data.table)
library(tmlensemble)
library(prodlim)
library(parallel)

synthesize_TD1_riskengine_data <- function(n,keep = NULL){
  #
  # empty model
  #
  sim_model <- lava::lvm()
  #
  # covariates
  #
  lava::distribution(sim_model,~sex) <- lava::binomial.lvm(p=0.462559634930512)
  lava::distribution(sim_model,~age) <- lava::normal.lvm(mean=42.2853710165578,sd=16.1358967048911)
  # this parameter seems to affect the age range:
  lava::covariance(sim_model,~age) <- 40
  lava::distribution(sim_model,~diabetes_duration) <- lava::normal.lvm(mean=-0.721917658912036,sd=11.8217218392295)
  lava::distribution(sim_model,~value_SBP) <- lava::normal.lvm(mean=113.015053986869,sd=15.1994730394504)
  lava::distribution(sim_model,~value_LDL) <- lava::normal.lvm(mean=2.34380835266364,sd=0.78907040283751)
  lava::distribution(sim_model,~value_HBA1C) <- lava::normal.lvm(mean=66.1642378865059,sd=15.1450866080122)
  lava::distribution(sim_model,~log2_eGFR_young) <- lava::normal.lvm(mean=-10.7725954841897,sd=1.85777815372659)
  lava::distribution(sim_model,~log2_eGFR_old) <- lava::normal.lvm(mean=-7.16808301118136,sd=0.295067430820889)
  lava::distribution(sim_model,~value_Smoking) <- lava::binomial.lvm(p=0.00406502756222992)
  lava::distribution(sim_model,~value_Motion) <- lava::binomial.lvm(p=0.647591503001192)
  # albuminuria
  sim_model <- lava::categorical(sim_model,~value_Albuminuria,labels=c('Normal','Micro','Macro'),K=3,p=c(0.839037544077992,0.121344119477287))
  lava::transform(sim_model, value_AlbuminuriaMicro~value_Albuminuria) <- function(x){1*(c(x)=='Micro')}
  lava::transform(sim_model, value_AlbuminuriaMacro~value_Albuminuria) <- function(x){1*(c(x)=='Macro')}
  #
  # outcome
  #
  # latent event times: .1 is cvd .2 death without cvd .0 is censoring
  lava::distribution(sim_model,~time.event.1) <- lava::coxWeibull.lvm(scale=0.00344127508109558,shape=1.1383691898171)
  lava::distribution(sim_model,~time.event.2) <- lava::coxWeibull.lvm(scale=8.8338885667654e-05,shape=1.49710154034773)
  lava::distribution(sim_model,~time.event.0) <- lava::coxWeibull.lvm(scale=0.0105086718102263,shape=2.2578792372401)
  # observed event time
  sim_model <- lava::eventTime(sim_model,time_cvd ~ min(time.event.0=0, time.event.1=1, time.event.2=2),'status_cvd')
  # uncensored event time
  sim_model <- lava::eventTime(sim_model,uncensored_time_cvd ~ min(time.event.1=1, time.event.2=2),'uncensored_status_cvd')
  #
  # structural equations
  #
  ## lava::regression(sim_model,age~sex) <- c(-0.605030442652732)
  lava::regression(sim_model,age~sex) <- c(-0.105030442652732)
  lava::regression(sim_model,diabetes_duration~age+sex) <- c(0.461453096059993,1.8834280372868)
  lava::regression(sim_model,value_SBP~diabetes_duration+age+sex) <- c(0.173711840258398,0.408090766933998,-5.26413576607228)
  lava::regression(sim_model,value_LDL~value_SBP+diabetes_duration+age+sex) <- c(0.0017997232707461,-0.00215485947423231,-0.00194439878409113,-0.0247941295129788)
  lava::regression(sim_model,value_HBA1C~value_LDL+value_SBP+diabetes_duration+age+sex) <- c(2.41425803932559,-0.0366520181671025,0.105793793763282,-0.113517351234731,1.14432850532257)
  lava::regression(sim_model,log2_eGFR_young~value_AlbuminuriaMicro+value_AlbuminuriaMacro+value_HBA1C+value_LDL+value_SBP+diabetes_duration+age+sex) <- c(-0.0378605757282762,0.369860071363398,-0.000381556996934381,0.0759866855173213,-0.000921667597434942,0.0109443806788427,0.173985375155216,-0.05148192662188)
  # Fixme: it is odd to let the log2_eGFR_old depend on log2_eGFR_young
  lava::regression(sim_model,log2_eGFR_old~log2_eGFR_young+value_AlbuminuriaMicro+value_AlbuminuriaMacro+value_HBA1C+value_LDL+value_SBP+diabetes_duration+age+sex) <- c(-0.998329448050345,0.13212871605835,0.627054487201391,-0.0017321066682523,-0.0168673394900242,0.000656653716889479,0.00231859330125351,0.0105241259989393,0.0666215469676934)
  lava::regression(sim_model,value_Smoking~log2_eGFR_old+log2_eGFR_young+value_AlbuminuriaMicro+value_AlbuminuriaMacro+value_HBA1C+value_LDL+value_SBP+diabetes_duration+age+sex) <- c(-0.573224437209381,-0.533358383827049,0.417778953216348,0.811962094950751,0.0258095967578417,0.128112260674147,-0.0102291702067243,-0.00575108475194562,0.00529875681026287,-0.130756787272673)
  lava::regression(sim_model,value_Motion~value_Smoking+log2_eGFR_old+log2_eGFR_young+value_AlbuminuriaMicro+value_AlbuminuriaMacro+value_HBA1C+value_LDL+value_SBP+diabetes_duration+age+sex) <- c(-0.31142087785739,-0.134144429182559,-0.0862952685994639,-0.121332306983678,-0.334529806176446,-0.0100110418345244,-0.138145721076664,0.00690508520692886,-0.0080292376964733,-0.00517733988644732,0.070748099701443)
  lava::regression(sim_model,time.event.0~value_Motion+value_Smoking+log2_eGFR_old+log2_eGFR_young+value_AlbuminuriaMicro+value_AlbuminuriaMacro+value_HBA1C+value_LDL+value_SBP+diabetes_duration+age+sex) <- c(-0.0389267779620385,-0.106325304898881,-0.137101392329058,-0.118521604228187,-0.0330093143101979,0.0763711547050245,-0.00617576911315076,0.0509424048213361,-0.00757123095797519,-0.0125850223581346,-0.0130146695831196,-0.029633759327063)
  lava::regression(sim_model,time.event.1~value_Motion+value_Smoking+log2_eGFR_old+log2_eGFR_young+value_AlbuminuriaMicro+value_AlbuminuriaMacro+value_HBA1C+value_LDL+value_SBP+diabetes_duration+age+sex) <- c(-0.176043482185139,0.311354666798778,0.442184803104845,0.480004039250349,0.458766977417934,0.760940706393913,0.0143694959616522,0.0701435790868005,0.00607734573133048,0.00980515944727106,0.0412496430494417,-0.224377188536924)
  lava::regression(sim_model,time.event.2~value_Motion+value_Smoking+log2_eGFR_old+log2_eGFR_young+value_AlbuminuriaMicro+value_AlbuminuriaMacro+value_HBA1C+value_LDL+value_SBP+diabetes_duration+age+sex) <- c(-0.410553409305998,0.833815403131808,0.3289175620287,0.313118375461698,0.225523885897843,0.604004891420425,0.0102472549995989,0.0287051547482761,0.00156423178426314,0.00678536173681053,0.0803874049498844,-0.276186036048047)
  # construction of eGFR
  lava::transform(sim_model,eGFR~log2_eGFR_old+log2_eGFR_young+age) <- function(x){as.numeric(x[["age"]]<40)*100*2^{x[["log2_eGFR_young"]]}+as.numeric(x[["age"]]>40)*100*2^{x[["log2_eGFR_old"]]}}
  # simulate from model
  d <- lava::sim(sim_model,n)
  data.table::setDT(d)
  d[,log2_eGFR_old := NULL]
  d[,log2_eGFR_young := NULL]
  d[,value_AlbuminuriaMicro := NULL]
  d[,value_AlbuminuriaMacro := NULL]
  d[,time.event.0 := NULL]
  d[,time.event.1 := NULL]
  d[,time.event.2 := NULL]
  if (length(keep)>0) d[,keep,with=FALSE]
  else d[]
}

gam_model_curve <- function(k_value, number_of_nodes,nfolds) {


  pr <- 0

  for(seed in 1:10){

    set.seed(seed)
    # I would repeat the following lines 100 times and then average
    d <- synthesize_TD1_riskengine_data(100000,
                                        keep = c("age", "uncensored_time_cvd", "uncensored_status_cvd"))

    d[['id']] <- 1:nrow(d)

    l1 <- Learner_gam(covariates = c(paste0("s(age,k=", k_value, ")")), intercept =
                        FALSE)

    out <- Superlearner(
      d,
      id = "id",
      status = "uncensored_status_cvd",
      stratified_k_fold = FALSE,
      event_time = "uncensored_time_cvd",
      learners = list(l1),
      nodes = seq(1, 1500, length.out = number_of_nodes),
      meta_learner_algorithm = "glmnet",
      add_nodes_metalearner = TRUE,
      add_intercept_metalearner = FALSE,
      nested_cross_validation_meta_learner = TRUE,
      penalise_nodes_metalearner = TRUE,
      nfold = nfolds
    )


    xseq <- seq(19, 66, 1)

    pr <- pr + predictRisk(out,
                           newdata = data.frame(age = xseq),
                           times = 5,
                           cause = 1)

  }

  model_out <- data.table(
    age=xseq,
    predicted_absr=pr/length(1:10),
    number_of_nodes=number_of_nodes,
    nfolds=nfolds,
    k_value=k_value
  )


  return(model_out)
}

# Define the parameter grid
params_list <- list(
  list(k_value = -1, number_of_nodes = 150, nfolds = 5),
  list(k_value = -1, number_of_nodes = 500, nfolds = 5),
  list(k_value = -1, number_of_nodes = 800, nfolds = 5),
  list(k_value = 1, number_of_nodes = 150, nfolds = 5),
  list(k_value = 1, number_of_nodes = 500, nfolds = 5),
  list(k_value = 1, number_of_nodes = 800, nfolds = 5),
  list(k_value = 2, number_of_nodes = 150, nfolds = 5),
  list(k_value = 2, number_of_nodes = 500, nfolds = 5),
  list(k_value = 2, number_of_nodes = 800, nfolds = 5)
)


run_gam_model <- function(params) {
  gam_model_curve(k_value = params$k_value,
                  number_of_nodes = params$number_of_nodes,
                  nfolds = params$nfolds)
}


# Set up cluster
n_cores <- 2
cl <- makeCluster(n_cores)

# Export necessary items to workers
clusterExport(cl, varlist = c("synthesize_TD1_riskengine_data",
                              "gam_model_curve",
                              "run_gam_model"))
clusterEvalQ(cl, {
  library(lava)
  library(data.table)
  library(prodlim)
  library(tmlensemble)
})


# Run in parallel
results_list <- parLapply(cl, params_list, run_gam_model)

# Stop cluster
stopCluster(cl)

# Combine all data.tables into one
results_dt <- rbindlist(results_list, use.names = TRUE, fill = TRUE)

# Save output to file
output_path <- "~/modi_mount/poisson_sl/gam_models_1807.txt"
write.table(result_dt, file = output_path, row.names = FALSE, col.names = FALSE)

# Stop the cluster
stopCluster(cl)

