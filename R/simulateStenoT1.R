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
                            # Coefficients for age, sex, diabetes duration, and SBP
                            beta_sex_age = 0.105,
                            beta_age_diabetes_duration = 0.461,
                            beta_sex_diabetes_duration = 1.883,
                            beta_diabetes_duration_value_SBP = 0.174,
                            beta_age_value_SBP = 0.408,
                            beta_sex_value_SBP = -5.264,
                            # Coefficients for LDL
                            beta_value_SBP_value_LDL = 0.0018,
                            beta_diabetes_duration_value_LDL = -0.00215,
                            beta_age_value_LDL = -0.00194,
                            beta_sex_value_LDL = -0.02479,
                            # Coefficients for HBA1C
                            beta_value_LDL_value_HBA1C = 2.414,
                            beta_value_SBP_value_HBA1C = -0.03665,
                            beta_diabetes_duration_value_HBA1C = 0.1058,
                            beta_age_value_HBA1C = -0.11352,
                            beta_sex_value_HBA1C = 1.144,
                            # Coefficients for log2_eGFR_young
                            beta_value_AlbuminuriaMicro_log2_eGFR_young = -0.038,
                            beta_value_AlbuminuriaMacro_log2_eGFR_young = 0.37,
                            beta_value_HBA1C_log2_eGFR_young = -0.00038,
                            beta_value_LDL_log2_eGFR_young = 0.076,
                            beta_value_SBP_log2_eGFR_young = -0.00092,
                            beta_diabetes_duration_log2_eGFR_young = 0.011,
                            beta_age_log2_eGFR_young = 0.174,
                            beta_sex_log2_eGFR_young = -0.051,
                            # Coefficients for log2_eGFR_old
                            beta_log2_eGFR_young_log2_eGFR_old = -0.998,
                            beta_value_AlbuminuriaMicro_log2_eGFR_old = 0.132,
                            beta_value_AlbuminuriaMacro_log2_eGFR_old = 0.627,
                            beta_value_HBA1C_log2_eGFR_old = -0.00173,
                            beta_value_LDL_log2_eGFR_old = -0.01687,
                            beta_value_SBP_log2_eGFR_old = 0.00066,
                            beta_diabetes_duration_log2_eGFR_old = 0.00232,
                            beta_age_log2_eGFR_old = 0.01052,
                            beta_sex_log2_eGFR_old = 0.06662,
                            # Coefficients for Smoking
                            beta_log2_eGFR_old_value_Smoking = -0.573,
                            beta_log2_eGFR_young_value_Smoking = -0.533,
                            beta_value_AlbuminuriaMicro_value_Smoking = 0.418,
                            beta_value_AlbuminuriaMacro_value_Smoking = 0.812,
                            beta_value_HBA1C_value_Smoking = 0.026,
                            beta_value_LDL_value_Smoking = 0.128,
                            beta_value_SBP_value_Smoking = -0.0102,
                            beta_diabetes_duration_value_Smoking = -0.00575,
                            beta_age_value_Smoking = 0.0053,
                            beta_sex_value_Smoking = -0.131,
                            # Coefficients for Motion
                            beta_value_Smoking_value_Motion = -0.311,
                            beta_log2_eGFR_old_value_Motion = -0.134,
                            beta_log2_eGFR_young_value_Motion = -0.086,
                            beta_value_AlbuminuriaMicro_value_Motion = -0.121,
                            beta_value_AlbuminuriaMacro_value_Motion = -0.335,
                            beta_value_HBA1C_value_Motion = -0.01001,
                            beta_value_LDL_value_Motion = -0.138,
                            beta_value_SBP_value_Motion = 0.0069,
                            beta_diabetes_duration_value_Motion = -0.00803,
                            beta_age_value_Motion = -0.00518,
                            beta_sex_value_Motion = 0.071,
                            # Coefficients for rate_censored
                            beta_value_Motion_rate_censored = -0.039,
                            beta_value_Smoking_rate_censored = -0.106,
                            beta_log2_eGFR_old_rate_censored = -0.137,
                            beta_log2_eGFR_young_rate_censored = -0.119,
                            beta_value_AlbuminuriaMicro_rate_censored = -0.033,
                            beta_value_AlbuminuriaMacro_rate_censored = 0.076,
                            beta_value_HBA1C_rate_censored = -0.00618,
                            beta_value_LDL_rate_censored = 0.051,
                            beta_value_SBP_rate_censored = -0.00757,
                            beta_diabetes_duration_rate_censored = -0.01259,
                            beta_age_rate_censored = -0.030,
                            beta_sex_rate_censored = -0.02963,
                            # Coefficients for rate_cvd
                            beta_value_Motion_rate_cvd = -0.176,
                            beta_value_Smoking_rate_cvd = 0.311,
                            beta_log2_eGFR_old_rate_cvd = 0.044,
                            beta_log2_eGFR_young_rate_cvd = 0.48,
                            beta_value_AlbuminuriaMicro_rate_cvd = 0.046,
                            beta_value_AlbuminuriaMacro_rate_cvd = 0.076,
                            beta_value_HBA1C_rate_cvd = 0.01437,
                            beta_value_LDL_rate_cvd = 0.1,
                            beta_value_SBP_rate_cvd = 0.00608,
                            beta_diabetes_duration_rate_cvd = 0.02,
                            beta_age_rate_cvd = 0.05,
                            beta_sex_rate_cvd = -0.224,
                            # Coefficients for rate_death
                            beta_value_Motion_rate_death = -0.411,
                            beta_value_Smoking_rate_death = 0.834,
                            beta_log2_eGFR_old_rate_death = 0.329,
                            beta_log2_eGFR_young_rate_death = 0.313,
                            beta_value_AlbuminuriaMicro_rate_death = 0.226,
                            beta_value_AlbuminuriaMacro_rate_death = 0.604,
                            beta_value_HBA1C_rate_death = 0.01025,
                            beta_value_LDL_rate_death = 0.09,
                            beta_value_SBP_rate_death = 0.00156,
                            beta_diabetes_duration_rate_death = 0.06,
                            beta_age_rate_death = 0.09,
                            beta_sex_rate_death = -0.276,
                            # Variables to keep
                            keep = NULL
                            ) {
    # Lava model
    sim_model <- lava::lvm()
    # Covariate distributions
    lava::distribution(sim_model, ~sex) <- lava::binomial.lvm(p = 0.4626)
    lava::distribution(sim_model, ~age) <- lava::normal.lvm(mean = 42.285, sd = 16.136)
    lava::covariance(sim_model, ~age) <- 40
    lava::distribution(sim_model, ~diabetes_duration) <- lava::normal.lvm(mean = -0.722, sd = 11.822)
    lava::distribution(sim_model, ~value_SBP) <- lava::normal.lvm(mean = 113.015, sd = 15.199)
    lava::distribution(sim_model, ~value_LDL) <- lava::normal.lvm(mean = 2.344, sd = 0.789)
    lava::distribution(sim_model, ~value_HBA1C) <- lava::normal.lvm(mean = 66.164, sd = 15.145)
    lava::distribution(sim_model, ~log2_eGFR_young) <- lava::normal.lvm(mean = -10.773, sd = 1.858)
    lava::distribution(sim_model, ~log2_eGFR_old) <- lava::normal.lvm(mean = -7.168, sd = 0.295)
    lava::distribution(sim_model, ~value_Smoking) <- lava::binomial.lvm(p = 0.0041)
    lava::distribution(sim_model, ~value_Motion) <- lava::binomial.lvm(p = 0.648)
    # latent event times: .1 is cvd .2 death without cvd .0 is censoring
    lava::distribution(sim_model,~time.event.1) <- lava::coxWeibull.lvm(scale=0.00344127508109558,shape=1.1383691898171)
    lava::distribution(sim_model,~time.event.2) <- lava::coxWeibull.lvm(scale=8.8338885667654e-05,shape=1.49710154034773)
    lava::distribution(sim_model,~time.event.0) <- lava::coxWeibull.lvm(scale=0.0105086718102263,shape=2.2578792372401)

    # Structural equations
    lava::regression(sim_model, age ~ sex) <- c(beta_sex_age)
    lava::regression(sim_model, diabetes_duration ~ age + sex) <- c(beta_age_diabetes_duration, beta_sex_diabetes_duration)
    lava::regression(sim_model, value_SBP ~ diabetes_duration + age + sex) <- c(beta_diabetes_duration_value_SBP, beta_age_value_SBP, beta_sex_value_SBP)
    lava::regression(sim_model, value_LDL ~ value_SBP + diabetes_duration + age + sex) <- c(beta_value_SBP_value_LDL, beta_diabetes_duration_value_LDL, beta_age_value_LDL, beta_sex_value_LDL)
    lava::regression(sim_model, value_HBA1C ~ value_LDL + value_SBP + diabetes_duration + age + sex) <- c(beta_value_LDL_value_HBA1C, beta_value_SBP_value_HBA1C, beta_diabetes_duration_value_HBA1C, beta_age_value_HBA1C, beta_sex_value_HBA1C)
    lava::regression(sim_model, log2_eGFR_young ~ value_AlbuminuriaMicro + value_AlbuminuriaMacro + value_HBA1C + value_LDL + value_SBP + diabetes_duration + age + sex) <- c(
        beta_value_AlbuminuriaMicro_log2_eGFR_young, beta_value_AlbuminuriaMacro_log2_eGFR_young, beta_value_HBA1C_log2_eGFR_young, beta_value_LDL_log2_eGFR_young, beta_value_SBP_log2_eGFR_young, beta_diabetes_duration_log2_eGFR_young, beta_age_log2_eGFR_young, beta_sex_log2_eGFR_young
    )
    lava::regression(sim_model, log2_eGFR_old ~ log2_eGFR_young + value_AlbuminuriaMicro + value_AlbuminuriaMacro + value_HBA1C + value_LDL + value_SBP + diabetes_duration + age + sex) <- c(
        beta_log2_eGFR_young_log2_eGFR_old, beta_value_AlbuminuriaMicro_log2_eGFR_old, beta_value_AlbuminuriaMacro_log2_eGFR_old, beta_value_HBA1C_log2_eGFR_old, beta_value_LDL_log2_eGFR_old, beta_value_SBP_log2_eGFR_old, beta_diabetes_duration_log2_eGFR_old, beta_age_log2_eGFR_old, beta_sex_log2_eGFR_old
    )
    lava::regression(sim_model, value_Smoking ~ log2_eGFR_old + log2_eGFR_young + value_AlbuminuriaMicro + value_AlbuminuriaMacro + value_HBA1C + value_LDL + value_SBP + diabetes_duration + age + sex) <- c(
        beta_log2_eGFR_old_value_Smoking, beta_log2_eGFR_young_value_Smoking, beta_value_AlbuminuriaMicro_value_Smoking, beta_value_AlbuminuriaMacro_value_Smoking, beta_value_HBA1C_value_Smoking, beta_value_LDL_value_Smoking, beta_value_SBP_value_Smoking, beta_diabetes_duration_value_Smoking, beta_age_value_Smoking, beta_sex_value_Smoking
    )
    lava::regression(sim_model, value_Motion ~ value_Smoking + log2_eGFR_old + log2_eGFR_young + value_AlbuminuriaMicro + value_AlbuminuriaMacro + value_HBA1C + value_LDL + value_SBP + diabetes_duration + age + sex) <- c(
        beta_value_Smoking_value_Motion, beta_log2_eGFR_old_value_Motion, beta_log2_eGFR_young_value_Motion, beta_value_AlbuminuriaMicro_value_Motion, beta_value_AlbuminuriaMacro_value_Motion, beta_value_HBA1C_value_Motion, beta_value_LDL_value_Motion, beta_value_SBP_value_Motion, beta_diabetes_duration_value_Motion, beta_age_value_Motion, beta_sex_value_Motion
    )
    lava::regression(sim_model, time.event.0 ~ value_Motion + value_Smoking + log2_eGFR_old + log2_eGFR_young + value_AlbuminuriaMicro + value_AlbuminuriaMacro + value_HBA1C + value_LDL + value_SBP + diabetes_duration + age + sex) <- c(
        beta_value_Motion_rate_censored, beta_value_Smoking_rate_censored, beta_log2_eGFR_old_rate_censored, beta_log2_eGFR_young_rate_censored, beta_value_AlbuminuriaMicro_rate_censored, beta_value_AlbuminuriaMacro_rate_censored, beta_value_HBA1C_rate_censored, beta_value_LDL_rate_censored, beta_value_SBP_rate_censored, beta_diabetes_duration_rate_censored, beta_age_rate_censored, beta_sex_rate_censored
    )
    lava::regression(sim_model, time.event.1 ~ value_Motion + value_Smoking + log2_eGFR_old + log2_eGFR_young + value_AlbuminuriaMicro + value_AlbuminuriaMacro + value_HBA1C + value_LDL + value_SBP + diabetes_duration + age + sex) <- c(
        beta_value_Motion_rate_cvd, beta_value_Smoking_rate_cvd, beta_log2_eGFR_old_rate_cvd, beta_log2_eGFR_young_rate_cvd, beta_value_AlbuminuriaMicro_rate_cvd, beta_value_AlbuminuriaMacro_rate_cvd, beta_value_HBA1C_rate_cvd, beta_value_LDL_rate_cvd, beta_value_SBP_rate_cvd, beta_diabetes_duration_rate_cvd, beta_age_rate_cvd, beta_sex_rate_cvd
    )
    lava::regression(sim_model, time.event.2 ~ value_Motion + value_Smoking + log2_eGFR_old + log2_eGFR_young + value_AlbuminuriaMicro + value_AlbuminuriaMacro + value_HBA1C + value_LDL + value_SBP + diabetes_duration + age + sex) <- c(
        beta_value_Motion_rate_death, beta_value_Smoking_rate_death, beta_log2_eGFR_old_rate_death, beta_log2_eGFR_young_rate_death, beta_value_AlbuminuriaMicro_rate_death, beta_value_AlbuminuriaMacro_rate_death, beta_value_HBA1C_rate_death, beta_value_LDL_rate_death, beta_value_SBP_rate_death, beta_diabetes_duration_rate_death, beta_age_rate_death, beta_sex_rate_death
    )

    
    # Transformation
    # Albuminuria 
    sim_model <- lava::categorical(sim_model, ~value_Albuminuria, labels = c("Normal", "Micro", "Macro"), K = 3, p = c(0.839, 0.121))
    lava::transform(sim_model, value_AlbuminuriaMicro ~ value_Albuminuria) <- function(x) 1 * (c(x) == "Micro")
    lava::transform(sim_model, value_AlbuminuriaMacro ~ value_Albuminuria) <- function(x) 1 * (c(x) == "Macro")

    # eGFR 
    lava::transform(sim_model, eGFR ~ log2_eGFR_old + log2_eGFR_young + age) <- function(x) {
        as.numeric(x[["age"]] < 40) * 100 * 2^x[["log2_eGFR_young"]] +
            as.numeric(x[["age"]] > 40) * 100 * 2^x[["log2_eGFR_old"]]
    }

    # observed event time
    sim_model <- lava::eventTime(sim_model,time ~ min(time.event.0=0, time.event.1=1, time.event.2=2),'event')

    # uncensored event time
    sim_model <- lava::eventTime(sim_model,uncensored_time ~ min(time.event.1=1, time.event.2=2),'uncensored_event')

    # Simulate data
    d <- lava::sim(sim_model, n)
    data.table::setDT(d)

    # Clean simulated data
    d[, log2_eGFR_old := NULL]
    d[, log2_eGFR_young := NULL]
    d[, value_AlbuminuriaMicro := NULL]
    d[, value_AlbuminuriaMacro := NULL]
    setorder(d, time, -event)
    d[, id := 1:.N]
    setcolorder(d, "id")

    # Convert specific variables to factor
    d[, c("sex", "value_Smoking", "value_Motion") := lapply(.SD, as.factor), .SDcols = c("sex", "value_Smoking", "value_Motion")]
    
    # Return filtered data or entire dataset
    if (length(keep) > 0) {
        return(d[, keep,with = FALSE][])
    }
    return(d[])
}



######################################################################
### simulateStenoT1.R ends here

