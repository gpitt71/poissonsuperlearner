### simulateStenoT1.R ---
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: feb 12 2026 (06:22)
## Version:
## Last-Updated: feb 12 2026 (08:12)
##           By: Thomas Alexander Gerds
##     Update #: 18
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
##' Simulate data that are alike those behind the Steno type-1 risk engine
##'
##' The data are computer generated and do not belong to real people.
##' @title Simulate time to CVD of hypothetical type-1 diabetes patients
##' @param n Sample size
##' @param beta_age_rate_cvd Log-hazard ratio for age on the rate of CVD
##' @param beta_LDL_rate_cvd Log-hazard ratio for LDL on the rate of CVD
##' @param beta_diabetes_duration_rate_cvd Log-hazard ratio for diabetes duration on the rate of CVD
##' @param keep Optional vector of variables to keep in the data set
##' @return data.table with baseline covariates and
##' time to cvd (time,event) where event has values 0 for right censored
##'         1 for cvd and 2 for death without cvd.
##' @examples
##' simulateStenoT1(n=3)
##' @export
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
simulateStenoT1 <- function(
    n,
    coefficient_age = 0.05,
    coefficient_LDL = 0.10,
    value_diabetis = 0.02,
    keep = NULL,
    scenario = c("alpha", "beta"),
    competing_risks = FALSE
) {
  scenario <- match.arg(scenario)

  # ------------------------------------------------------------------
  # Shared model parts (covariates + structural equations)
  # ------------------------------------------------------------------
  sim_model <- lava::lvm()

  # covariates
  lava::distribution(sim_model, ~sex) <- lava::binomial.lvm(p = 0.462559634930512)
  lava::distribution(sim_model, ~age) <- lava::normal.lvm(mean = 42.2853710165578, sd = 16.1358967048911)
  # affects the age range
  lava::covariance(sim_model, ~age) <- 40
  lava::distribution(sim_model, ~diabetes_duration) <- lava::normal.lvm(mean = -0.721917658912036, sd = 11.8217218392295)
  lava::distribution(sim_model, ~value_SBP) <- lava::normal.lvm(mean = 113.015053986869, sd = 15.1994730394504)
  lava::distribution(sim_model, ~value_LDL) <- lava::normal.lvm(mean = 2.34380835266364, sd = 0.78907040283751)
  lava::distribution(sim_model, ~value_HBA1C) <- lava::normal.lvm(mean = 66.1642378865059, sd = 15.1450866080122)
  lava::distribution(sim_model, ~log2_eGFR_young) <- lava::normal.lvm(mean = -10.7725954841897, sd = 1.85777815372659)
  lava::distribution(sim_model, ~log2_eGFR_old) <- lava::normal.lvm(mean = -7.16808301118136, sd = 0.295067430820889)
  lava::distribution(sim_model, ~value_Smoking) <- lava::binomial.lvm(p = 0.00406502756222992)
  lava::distribution(sim_model, ~value_Motion) <- lava::binomial.lvm(p = 0.647591503001192)

  # albuminuria
  sim_model <- lava::categorical(
    sim_model, ~value_Albuminuria,
    labels = c("Normal", "Micro", "Macro"),
    K = 3, p = c(0.839037544077992, 0.121344119477287)
  )
  lava::transform(sim_model, value_AlbuminuriaMicro ~ value_Albuminuria) <- function(x) {
    1 * (c(x) == "Micro")
  }
  lava::transform(sim_model, value_AlbuminuriaMacro ~ value_Albuminuria) <- function(x) {
    1 * (c(x) == "Macro")
  }

  # structural equations (shared)
  lava::regression(sim_model, age ~ sex) <- c(-0.105030442652732)
  lava::regression(sim_model, diabetes_duration ~ age + sex) <- c(0.461453096059993, 1.8834280372868)
  lava::regression(sim_model, value_SBP ~ diabetes_duration + age + sex) <- c(0.173711840258398, 0.408090766933998, -5.26413576607228)
  lava::regression(sim_model, value_LDL ~ value_SBP + diabetes_duration + age + sex) <- c(
    0.0017997232707461, -0.00215485947423231, -0.00194439878409113, -0.0247941295129788
  )
  lava::regression(sim_model, value_HBA1C ~ value_LDL + value_SBP + diabetes_duration + age + sex) <- c(
    2.41425803932559, -0.0366520181671025, 0.105793793763282, -0.113517351234731, 1.14432850532257
  )

  lava::regression(
    sim_model,
    log2_eGFR_young ~ value_AlbuminuriaMicro + value_AlbuminuriaMacro + value_HBA1C +
      value_LDL + value_SBP + diabetes_duration + age + sex
  ) <- c(
    -0.0378605757282762, 0.369860071363398, -0.000381556996934381, 0.0759866855173213,
    -0.000921667597434942, 0.0109443806788427, 0.173985375155216, -0.05148192662188
  )

  lava::regression(
    sim_model,
    log2_eGFR_old ~ log2_eGFR_young + value_AlbuminuriaMicro + value_AlbuminuriaMacro +
      value_HBA1C + value_LDL + value_SBP + diabetes_duration + age + sex
  ) <- c(
    -0.998329448050345, 0.13212871605835, 0.627054487201391, -0.0017321066682523,
    -0.0168673394900242, 0.000656653716889479, 0.00231859330125351, 0.0105241259989393,
    0.0666215469676934
  )

  lava::regression(
    sim_model,
    value_Smoking ~ log2_eGFR_old + log2_eGFR_young + value_AlbuminuriaMicro +
      value_AlbuminuriaMacro + value_HBA1C + value_LDL + value_SBP +
      diabetes_duration + age + sex
  ) <- c(
    -0.573224437209381, -0.533358383827049, 0.417778953216348, 0.811962094950751,
    0.0258095967578417, 0.128112260674147, -0.0102291702067243, -0.00575108475194562,
    0.00529875681026287, -0.130756787272673
  )

  lava::regression(
    sim_model,
    value_Motion ~ value_Smoking + log2_eGFR_old + log2_eGFR_young +
      value_AlbuminuriaMicro + value_AlbuminuriaMacro + value_HBA1C +
      value_LDL + value_SBP + diabetes_duration + age + sex
  ) <- c(
    -0.31142087785739, -0.134144429182559, -0.0862952685994639, -0.121332306983678,
    -0.334529806176446, -0.0100110418345244, -0.138145721076664, 0.00690508520692886,
    -0.0080292376964733, -0.00517733988644732, 0.070748099701443
  )

  # construction of eGFR
  lava::transform(sim_model, eGFR ~ log2_eGFR_old + log2_eGFR_young + age) <- function(x) {
    as.numeric(x[["age"]] < 40) * 100 * 2^(x[["log2_eGFR_young"]]) +
      as.numeric(x[["age"]] > 40) * 100 * 2^(x[["log2_eGFR_old"]])
  }

  # ------------------------------------------------------------------
  # Outcome model (scenario-specific)
  # ------------------------------------------------------------------
  # latent event times: .1 is CVD, .2 is death without CVD (alpha only), .0 is censoring
  lava::distribution(sim_model, ~time.event.1) <- lava::coxWeibull.lvm(
    scale = 0.00344127508109558,
    shape = 1.1383691898171
  )
  lava::distribution(sim_model, ~time.event.0) <- lava::coxWeibull.lvm(
    scale = 0.0105086718102263,
    shape = 2.2578792372401
  )

  if (scenario == "alpha") {
    lava::distribution(sim_model, ~time.event.2) <- lava::coxWeibull.lvm(
      scale = 8.8338885667654e-05,
      shape = 1.49710154034773
    )
  }

  # censoring model (shared)
  lava::regression(
    sim_model,
    time.event.0 ~ value_Motion + value_Smoking + log2_eGFR_old + log2_eGFR_young +
      value_AlbuminuriaMicro + value_AlbuminuriaMacro + value_HBA1C + value_LDL +
      value_SBP + diabetes_duration + age + sex
  ) <- c(
    -0.0389267779620385, -0.106325304898881, -0.137101392329058, -0.118521604228187,
    -0.0330093143101979, 0.0763711547050245, -0.00617576911315076, 0.0509424048213361,
    -0.00757123095797519, -0.0125850223581346, -0.0130146695831196, -0.029633759327063
  )

  if (scenario == "alpha") {
    # event 1 hazard (linear)
    lava::regression(
      sim_model,
      time.event.1 ~ value_Motion + value_Smoking + log2_eGFR_old + log2_eGFR_young +
        value_AlbuminuriaMicro + value_AlbuminuriaMacro + value_HBA1C + value_LDL +
        value_SBP + diabetes_duration + age + sex
    ) <- c(
      -0.176043482185139, 0.311354666798778, 0.0442184803104845, 0.480004039250349,
      0.0458766977417934, 0.0760940706393913, 0.0143694959616522, coefficient_LDL,
      0.00607734573133048, value_diabetis, coefficient_age, -0.224377188536924
    )

    # event 2 hazard (death)
    lava::regression(
      sim_model,
      time.event.2 ~ value_Motion + value_Smoking + log2_eGFR_old + log2_eGFR_young +
        value_AlbuminuriaMicro + value_AlbuminuriaMacro + value_HBA1C + value_LDL +
        value_SBP + diabetes_duration + age + sex
    ) <- c(
      -0.410553409305998, 0.833815403131808, 0.3289175620287, 0.313118375461698,
      0.225523885897843, 0.604004891420425, 0.0102472549995989, 0.09,
      0.00156423178426314, 0.06, 0.09, -0.276186036048047
    )
  }

  if (scenario == "beta") {
    # Nonlinear terms for event 1 hazard
    lava::transform(sim_model, age_hinge55_sq ~ age) <- function(x) {
      a <- x[["age"]]
      (pmax(a - 55, 0)^2) / 100
    }
    lava::transform(sim_model, LDL_hinge3_sq ~ value_LDL) <- function(x) {
      ldl <- x[["value_LDL"]]
      (pmax(ldl - 3.0, 0)^2)
    }
    nl_age_strength <- 1
    nl_ldl_strength <- 1.5

    lava::regression(
      sim_model,
      time.event.1 ~ value_Motion + value_Smoking + log2_eGFR_old + log2_eGFR_young +
        value_AlbuminuriaMicro + value_AlbuminuriaMacro + value_HBA1C + value_LDL +
        value_SBP + diabetes_duration + age + sex +
        age_hinge55_sq + LDL_hinge3_sq
    ) <- c(
      -0.176043482185139, 0.311354666798778, 0.0442184803104845, 0.480004039250349,
      0.0458766977417934, 0.0760940706393913, 0.0143694959616522, coefficient_LDL,
      0.00607734573133048, value_diabetis, coefficient_age, -0.224377188536924,
      nl_age_strength, nl_ldl_strength
    )
  }

  # ------------------------------------------------------------------
  # Simulate and post-process
  # ------------------------------------------------------------------
  d <- lava::sim(sim_model, n)
  data.table::setDT(d)

  # drop internal variables
  d[, log2_eGFR_old := NULL]
  d[, log2_eGFR_young := NULL]
  d[, value_AlbuminuriaMicro := NULL]
  d[, value_AlbuminuriaMacro := NULL]

  if (scenario == "beta") {
    d[, age_hinge55_sq := NULL]
    d[, LDL_hinge3_sq := NULL]
  }

  # observed time/event
  if (scenario == "alpha" && isTRUE(competing_risks)) {
    d[, time_cvd := pmin(time.event.0, time.event.1, time.event.2)]
    d[, status_cvd := data.table::fifelse(time.event.1 < time.event.0 & time.event.1 <= time.event.2, 1,
                                          data.table::fifelse(time.event.2 < time.event.0 & time.event.2 < time.event.1, 2, 0)
    )]
  } else {
    d[, time_cvd := pmin(time.event.0, time.event.1)]
    d[, status_cvd := as.numeric(time.event.1 < time.event.0)]
  }

  # uncensored event time/status
  if (scenario == "alpha") {
    d[, uncensored_time_cvd := pmin(time.event.1, time.event.2)]
    d[, uncensored_status_cvd := as.numeric(time.event.2 < time.event.1) + 1L]
  } else {
    d[, uncensored_time_cvd := time.event.1]
    d[, uncensored_status_cvd := 1L]
  }

  # keep current simulator output names (time/event)
  d[, time := time_cvd]
  d[, event := status_cvd]
  d[, uncensored_time := uncensored_time_cvd]
  d[, uncensored_event := uncensored_status_cvd]

  # beta script drops time.event.*
  if (scenario == "beta") {
    d[, time.event.0 := NULL]
    d[, time.event.1 := NULL]
  }

  # id + ordering
  data.table::setorder(d, time, -event)
  d[, id := 1:.N]
  data.table::setcolorder(d, "id")

  # factors
  d[, c("sex", "value_Smoking", "value_Motion") := lapply(.SD, as.factor),
    .SDcols = c("sex", "value_Smoking", "value_Motion")
  ]

  if (length(keep) > 0) {
    return(d[, keep, with = FALSE][])
  }
  d[]
}
