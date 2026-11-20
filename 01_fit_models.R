# ==============================================================================
# 01_fit_models.R
# Methadone and buprenorphine consumption vs. GDP, health spending, life expectancy
# ------------------------------------------------------------------------------
# What this script does, in order:
#   1. Reads the INCB / World Bank panel and puts it in long format
#      (one row per country x period x drug).
#   2. Codes the consumption symbols (dash, "<<", blanks), screens for unusual
#      values, and builds the primary analysis sample (the eight non-overlapping
#      3-year periods, 2000-2002 ... 2021-2023).
#   3. Fits three models separately for methadone and buprenorphine, in TWO PARTS
#      because many countries report nil (a dash) rather than a number:
#        - AMOUNT: how much is used, among country-periods with positive use
#          (log consumption; zeros cannot be logged)
#        - ANY USE: whether any use is reported (1) or nil (0) (logistic models with
#          country-clustered standard errors, plus a linear fixed-effects model)
#      Each part uses three model types:
#        - a mixed model (random intercept for country, fixed effects for period)
#        - a between-within ("hybrid") mixed model
#        - a two-way fixed-effects regression with country-clustered SEs
#   4. Re-fits those models under each sensitivity analysis (Table 3), and fits the
#      unadjusted one-predictor models (Supplementary Table 7).
#   5. Runs the other checks: single-equation models (Supplementary Table 4), a joint
#      two-part model (extra, not in the paper), the between/within test (Supplementary
#      Table 6), nil-entry counts (Supplementary Table 3) and the country-average checks
#      (Supplementary Table 1).
#   6. Cross-checks the amount models against an independent fit (glmmTMB).
#   7. Saves everything to results/ so 02_tables_figures.R can pick it up.
#
# How to run: put this file next to incb_worldbank_panel_excel.csv, set the working
# directory to that folder, and run the whole script (source() or Rscript).
# ==============================================================================


# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------
# Install anything that is missing (only does something the first time you run it).
pkgs <- c("tidyverse", "lme4", "lmerTest", "broom.mixed", "sandwich", "lmtest", "glmmTMB", "MASS", "numDeriv")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

library(tidyverse)    # data wrangling (dplyr, tidyr, purrr, readr, stringr, ...)
library(lme4)         # mixed-effects models
library(lmerTest)     # adds p-values and denominator degrees of freedom to lmer()
library(broom.mixed)  # tidy() turns mixed-model output into a data frame
library(sandwich)     # cluster-robust variance matrices (vcovCL)
library(lmtest)       # coeftest() applies a robust variance matrix to an lm


# ------------------------------------------------------------------------------
# Settings you might want to change
# ------------------------------------------------------------------------------
infile      <- "incb_worldbank_panel_excel.csv"
us          <- "US"                    # how the US is spelled in the `country` column
high_income <- "HIC"                   # income_type is HIC, UIC (upper middle), LMIC_LIC
drugs       <- c("methadone", "buprenorphine")

# The three predictors in the main model (after log-transforming GDP and health share).
main_terms  <- c("log_gdp", "log_health", "life_exp")

# Every association is reported per 10% higher GDP per person or health spending share (the predictors are
# logged, so a 10% rise is a step of log(1.10)) and per 1 year of life expectancy.
step <- 1.10

# Zeros (nil entries). An INCB dash means "nil", so it is always coded as 0. Because
# the log of zero is undefined, we do not force the zeros into one model. Instead:
#   AMOUNT part   uses only positive values (log consumption)
#   ANY USE part  uses every reported value: 1 = some use, 0 = nil
# ("<<" means use below one dose, so it counts as some use in the ANY USE part.)
# Single-equation checks that keep the zeros (log(1 + x) and Poisson) are run as
# sensitivity analyses (Supplementary Table 4).

# Peer high-income countries (for the US / Canada comparison): high-income countries
# other than the US and Canada that reported a numeric (positive) value for that drug
# in at least this many of the eight periods. Set to 1 to include every reporter.
peer_min_periods <- 4

dir.create("results", showWarnings = FALSE)   # folder for the CSVs we write


# ------------------------------------------------------------------------------
# 1. Read the data and reshape it
# ------------------------------------------------------------------------------
raw <- read_csv(infile, show_col_types = FALSE)

# The file is "wide": q_b is buprenorphine and q_m is methadone, each with its own
# *_symbol column that holds a dash or "<<" when there is no number.
# First, keep just the columns we need and give them short names.
panel <- raw %>%
  transmute(
    country,
    income_group    = income_type,
    region          = country_type,
    period          = time,                                   # e.g. "2009-2011"
    start_year      = as.integer(substr(time, 1, 4)),         # 2009
    end_year        = as.integer(substr(time, 6, 9)),         # 2011
    non_overlapping = as.logical(non_overlapping_period),     # TRUE for the 8 primary periods
    gdp             = gdp_per_capita_current_usd,
    life_exp        = life_expectancy,
    health_total    = health_expenditure_pct_gdp,             # total health spending, % of GDP
    health_govt     = govt_health_expenditure_pct_gdp,        # government health spending, % of GDP
    b_value = q_b, b_symbol = q_b_symbol,
    m_value = q_m, m_symbol = q_m_symbol
  ) %>%
  # Only countries with a World Bank income group are analysed (199 countries).
  filter(!is.na(income_group))

# Second, stack the two drugs so each row is one country x period x drug.
# (Same panel twice, once labelled buprenorphine and once methadone.)
long <- bind_rows(
  panel %>% mutate(drug = "buprenorphine", value = b_value, symbol = b_symbol),
  panel %>% mutate(drug = "methadone",     value = m_value, symbol = m_symbol)
) %>%
  select(-b_value, -b_symbol, -m_value, -m_symbol)


# ------------------------------------------------------------------------------
# 2. Code consumption and build the log variables
# ------------------------------------------------------------------------------
dash <- "\u2014"   # the em dash the INCB uses to mean "nil"

# `less_than_one` decides what "<<" (less than one dose) becomes in the AMOUNT part:
#   NA   -> treated as missing (main analysis)
#   0.5  -> used in one sensitivity analysis
code_consumption <- function(df, less_than_one = NA_real_) {
  df %>%
    mutate(
      # case_when() takes the first rule that fits, top to bottom
      consumption = case_when(
        !is.na(value)  ~ value,           # a reported number: keep it
        symbol == dash ~ 0,               # dash = nil = 0
        symbol == "<<" ~ less_than_one,   # less than one dose
        TRUE           ~ NA_real_         # anything else = not reported
      ),

      # AMOUNT outcome: log of consumption. Zeros can't be logged, so they become NA
      # here (and are handled by the ANY USE part instead).
      log_cons = if_else(consumption > 0, log(consumption), NA_real_),

      # ANY USE outcome: 1 = some use, 0 = nil, NA = not reported.
      # "<<" is use below one dose, so it counts as some use.
      any_use = case_when(
        !is.na(consumption) & consumption > 0  ~ 1,
        !is.na(consumption) & consumption == 0 ~ 0,
        symbol == "<<"                         ~ 1,
        TRUE                                   ~ NA_real_
      ),

      # Used only in the single-equation sensitivity check that keeps zeros in
      log1p_cons = log1p(consumption),

      # Predictors that get logged: GDP per person and health spending as % of GDP
      log_gdp    = log(gdp),
      log_health = log(health_total),
      log_govt   = log(health_govt)       # used only in the government-spending sensitivity run
    )
}

# The primary sample: the eight non-overlapping 3-year periods from 2000-2002 on.
# (The other windows in the file overlap, so they are used only in robustness checks.)
primary_periods <- function(df) filter(df, non_overlapping, start_year >= 2000)


# ------------------------------------------------------------------------------
# 3. Screen for unusual values
# ------------------------------------------------------------------------------
# A positive value is flagged if it is...
#   (a) above 100,000, or
#   (b) "isolated tenfold": at least 10x the country's median value AND at
#       least 10x the value just before it AND just after it (so it can't be the
#       first or last value in a country's series), or
#   (c) followed by a "persistent decline": at least 10x every later value
#       (we need at least 2 later values so a single dip doesn't count).
# Flagged values stay in the primary models and are dropped in a sensitivity run.
screen_values <- function(df) {
  df %>%
    filter(consumption > 0) %>%                 # screen only positive numeric values
    arrange(country, drug, start_year) %>%      # order each country's series by time
    group_by(country, drug) %>%                 # one series per country and drug
    mutate(
      i       = row_number(),                   # position within the series
      typical = median(consumption),            # the country's median value for this drug
      before  = lag(consumption),               # previous value in the series
      after   = lead(consumption),              # next value in the series

      isolated = !is.na(before) & !is.na(after) &
                 consumption >= 10 * typical &
                 consumption >= 10 * before &
                 consumption >= 10 * after,

      persistent = map_lgl(i, function(k) {
        later <- consumption[-seq_len(k)]       # everything after position k
        length(later) >= 2 && all(consumption[k] >= 10 * later)
      }),

      huge = consumption > 100000
    ) %>%
    ungroup() %>%
    filter(isolated | persistent | huge) %>%    # keep only the flagged rows
    # Build a readable label such as "isolated tenfold; persistent tenfold decline"
    mutate(flag = pmap_chr(list(isolated, persistent, huge), function(i, p, h) {
      str_c(c("isolated tenfold", "persistent tenfold decline", "above 100\u2009000")[c(i, p, h)],
            collapse = "; ")
    })) %>%
    select(country, period, drug, value = consumption, flag)
}

all_coded <- code_consumption(long)                      # every window, every row
flags     <- screen_values(primary_periods(all_coded))   # screen the primary periods
write_csv(flags, "results/flagged_values.csv")           # Supplementary Table 2

# Mark each row as a flagged value and/or in a flagged country.
# A "flagged country" is any country with a flagged value for EITHER drug.
add_flags <- function(df) {
  df %>%
    left_join(flags %>% transmute(country, period, drug, flag_value = TRUE),
              by = c("country", "period", "drug")) %>%
    mutate(flag_value   = coalesce(flag_value, FALSE),   # no match = not flagged
           flag_country = country %in% flags$country)
}

# `dat` is THE analysis dataset: primary periods only, with flags attached.
dat <- all_coded %>% primary_periods() %>% add_flags()
write_csv(dat, "results/analysis_data.csv")              # used for Table 1 and the figures


# ------------------------------------------------------------------------------
# 3b. US and Canada vs. peer high-income countries (Figure 2, Results text)
# ------------------------------------------------------------------------------
positive <- dat %>% filter(consumption > 0)       # numeric, positive country-period values

# Peers = high-income countries other than the US and Canada, kept only if they
# reported this drug in at least `peer_min_periods` periods. Requiring several
# periods keeps the peer group from changing a lot just because a country happened
# to report (or not report) in one period.
peers <- positive %>%
  filter(income_group == high_income, !country %in% c(us, "Canada")) %>%
  group_by(drug, country) %>%
  filter(n_distinct(period) >= peer_min_periods) %>%
  ungroup()

# Peer median and interquartile range in each period (the grey band and line in Figure 2)
peer_series <- peers %>%
  group_by(drug, period) %>%
  summarise(n_peers = n(), median = median(consumption),
            q1 = quantile(consumption, 0.25), q3 = quantile(consumption, 0.75),
            .groups = "drop")

# The US and Canada, period by period
us_canada <- positive %>%
  filter(country %in% c(us, "Canada")) %>%
  select(drug, period, country, consumption)

# The 2021-2023 comparison quoted in the Results: how many times the peer median?
peer_comparison <- us_canada %>%
  filter(period == "2021-2023") %>%
  left_join(peer_series %>% filter(period == "2021-2023") %>% select(drug, peer_median = median, n_peers),
            by = "drug") %>%
  mutate(times_peer_median = consumption / peer_median)

write_csv(peer_series,     "results/peer_series.csv")
write_csv(us_canada,       "results/us_canada_series.csv")
write_csv(peer_comparison, "results/peer_comparison.csv")


# ------------------------------------------------------------------------------
# 4. Model-fitting functions
# ------------------------------------------------------------------------------

# Complete-case analysis: keep one drug, and keep only rows that have every
# variable the model needs. Optionally keep only countries with >= `min_periods`
# complete periods.
complete_cases <- function(df, drug_name, predictors, extras = NULL, min_periods = NULL,
                           y = "log_cons") {
  df <- df %>%
    filter(drug == drug_name) %>%
    drop_na(all_of(c(y, predictors, extras)))

  if (!is.null(min_periods)) {
    df <- df %>%
      group_by(country) %>%
      filter(n_distinct(period) >= min_periods) %>%
      ungroup()
  }
  df
}

# Pull the fixed effects out of a mixed model, with Wald 95% CIs and normal (z) inference.
#
# The standard errors use the FULL observed information: the Hessian of the exact
# Gaussian log-likelihood with respect to the fixed effects AND the two variance
# parameters, evaluated at the maximum-likelihood fit. (lme4's own standard errors
# treat the variance parameters as known, so its intervals are slightly narrower.)
tidy_fixed <- function(m) {
  X    <- getME(m, "X")                      # fixed-effects design matrix
  y    <- getME(m, "y")                      # outcome
  g    <- getME(m, "flist")[[1]]             # country (the grouping factor)
  nc   <- as.numeric(rowsum(rep(1, length(y)), g))   # observations per country
  p    <- ncol(X)
  est  <- fixef(m)
  s_e  <- sigma(m)                                    # residual SD
  s_u  <- sqrt(as.numeric(VarCorr(m)[[1]]))           # country random-intercept SD

  # Negative log-likelihood of a random-intercept model, written out in closed form.
  # par = c(fixed effects, log residual SD, log country SD)
  neg_loglik <- function(par) {
    r   <- as.numeric(y - X %*% par[1:p])             # residuals given the fixed effects
    ve  <- exp(2 * par[p + 1])                        # residual variance
    vu  <- exp(2 * par[p + 2])                        # country variance
    rs  <- as.numeric(rowsum(r, g))                   # residual sums by country
    rr  <- as.numeric(rowsum(r^2, g))                 # residual sums of squares by country
    den <- ve + nc * vu
    0.5 * sum((nc - 1) * log(ve) + log(den) + (rr - vu * rs^2 / den) / ve) +
      0.5 * length(y) * log(2 * pi)
  }

  H  <- numDeriv::hessian(neg_loglik, c(est, log(s_e), log(s_u)))   # observed information
  se <- sqrt(diag(solve(H)))[1:p]                                   # SEs for the fixed effects

  tibble(term = names(est), estimate = as.numeric(est), std.error = se) %>%
    mutate(conf.low  = estimate - 1.96 * std.error,
           conf.high = estimate + 1.96 * std.error,
           p.value   = 2 * pnorm(-abs(estimate / std.error)))
}

# Pooled logistic regression with country-clustered standard errors (ANY USE part).
# There is no country random intercept here: standard errors are clustered by country
# instead, so repeated observations of the same country are not treated as independent.
# The 95% CI uses a t distribution with (number of countries - 1) degrees of freedom.
pooled_logit <- function(f, df, keep) {
  # (When a region such as North America ALWAYS reports use, its own coefficient cannot be
  # estimated and glm warns about it. That does not affect the other coefficients, so the
  # warning is silenced.)
  m  <- suppressWarnings(glm(f, family = binomial, data = df))
  ct <- coeftest(m, vcov = vcovCL(m, cluster = df$country))
  g  <- n_distinct(df$country)
  tibble(term = rownames(ct), estimate = ct[, 1], std.error = ct[, 2]) %>%
    filter(term %in% keep) %>%
    mutate(crit      = qt(0.975, g - 1),
           conf.low  = estimate - crit * std.error,
           conf.high = estimate + crit * std.error,
           p.value   = 2 * pt(-abs(estimate / std.error), g - 1)) %>%
    select(-crit)
}

# `y` is the outcome column and `binary` says whether it is a yes/no outcome:
#   AMOUNT part:   y = "log_cons",   binary = FALSE  (linear mixed models)
#   ANY USE part:  y = "any_use",    binary = TRUE   (pooled logistic, clustered SEs)
#   single equation check: y = "log1p_cons", binary = FALSE

# --- Primary model --------------------------------------------------------------
# AMOUNT: outcome ~ predictors + period (fixed effects) + random intercept for country.
#   The random intercept accounts for repeated measures of the same country; the
#   period effects absorb changes that are common to all countries. Fitted by maximum
#   likelihood (REML = FALSE), which reproduces the manuscript.
# ANY USE: the same without the random intercept (a pooled logistic model); standard
#   errors are clustered by country instead.
fit_mixed <- function(df, predictors, extras = NULL, y = "log_cons", binary = FALSE) {
  if (binary) {
    f <- reformulate(c(predictors, extras, "period"), y)
    return(pooled_logit(f, df, predictors) %>% mutate(model = "Primary logistic"))
  }
  f <- reformulate(c(predictors, extras, "period", "(1 | country)"), y)
  tidy_fixed(lmerTest::lmer(f, data = df, REML = FALSE)) %>%
    filter(term %in% predictors) %>%        # drop intercept and period dummies
    mutate(model = "Primary mixed")
}

# --- Between-within (hybrid) model --------------------------------------------
# Splits every predictor into
#   _between = the country's own average across periods (differences BETWEEN countries)
#   _within  = this period's value minus that average   (changes WITHIN a country)
# and puts both in the same model, so each gets its own coefficient.
fit_hybrid <- function(df, predictors, extras = NULL, y = "log_cons", binary = FALSE) {
  for (p in predictors) {
    df[[paste0(p, "_between")]] <- ave(df[[p]], df$country)                 # country mean
    df[[paste0(p, "_within")]]  <- df[[p]] - df[[paste0(p, "_between")]]    # deviation from it
  }
  split_terms <- c(paste0(predictors, "_between"), paste0(predictors, "_within"))

  tab <- if (binary) {
    pooled_logit(reformulate(c(split_terms, extras, "period"), y), df, split_terms)
  } else {
    f <- reformulate(c(split_terms, extras, "period", "(1 | country)"), y)
    tidy_fixed(lmerTest::lmer(f, data = df, REML = FALSE)) %>% filter(term %in% split_terms)
  }

  # Relabel so "log_health_between" becomes model = "Between", term = "log_health"
  tab %>%
    mutate(model = if_else(str_detect(term, "_between$"), "Between", "Within"),
           term  = str_remove(term, "_(between|within)$"))
}

# --- Fixed-effects regression -------------------------------------------------
# A country dummy and a period dummy for every country and period, so each country is
# compared only with itself over time. Standard errors are clustered by country. For the
# ANY USE part this is a linear probability model, so its coefficient is a change in the
# probability of any use.
fit_twfe <- function(df, predictors, y = "log_cons") {
  m  <- lm(reformulate(c(predictors, "factor(country)", "factor(period)"), y), data = df)
  # coeftest() swaps in the cluster-robust variance matrix
  ct <- coeftest(m, vcov = vcovCL(m, cluster = df$country))

  tibble(term = rownames(ct), estimate = ct[, 1], std.error = ct[, 2], p.value = ct[, 4]) %>%
    filter(term %in% predictors) %>%        # drop the dozens of dummy coefficients
    # 95% CI uses a t distribution with (number of clusters - 1) degrees of freedom
    mutate(crit      = qt(0.975, n_distinct(df$country) - 1),
           conf.low  = estimate - crit * std.error,
           conf.high = estimate + crit * std.error,
           model     = "Country FE") %>%
    select(-crit)
}

# Fit all three models for one drug under one specification, and return one tidy table.
# part = "amount" (how much, among positive values), "any" (any use vs nil) or
#        "single" (one equation, log(1 + consumption), zeros kept).
fit_all <- function(df, drug_name, predictors = main_terms, extras = NULL,
                    min_periods = NULL, spec = "Primary", part = "amount") {
  y      <- c(amount = "log_cons", any = "any_use", single = "log1p_cons")[[part]]
  binary <- part == "any"
  d <- complete_cases(df, drug_name, predictors, extras, min_periods, y)

  # If a subset leaves too few countries, skip it rather than fit a silly model
  if (n_distinct(d$country) < 5) return(NULL)

  out <- bind_rows(fit_mixed(d, predictors, extras, y, binary),
                   fit_hybrid(d, predictors, extras, y, binary),
                   fit_twfe(d, predictors, y)) %>%
    select(model, term, estimate, std.error, p.value, conf.low, conf.high) %>%
    mutate(drug = drug_name, spec = spec, part = part,
           n_obs = nrow(d), n_countries = n_distinct(d$country))   # sample size used

  # Put every association on one scale: per 10% higher GDP per person or health share, and per year of life
  # expectancy. `effect` = exp(coefficient x step): a multiplier on consumption (amount part), an odds ratio
  # (any-use part), or a ratio on 1 + consumption (single-equation check). The linear fixed-effects model
  # for ANY USE is shown instead as a change in percentage POINTS in the probability of any use.
  out <- out %>%
    mutate(per_unit    = if_else(str_detect(term, "^log_"), log(step), 1),
           is_linear   = part == "any" & model == "Country FE",
           effect      = if_else(is_linear, estimate  * per_unit * 100, exp(estimate  * per_unit)),
           effect_low  = if_else(is_linear, conf.low  * per_unit * 100, exp(conf.low  * per_unit)),
           effect_high = if_else(is_linear, conf.high * per_unit * 100, exp(conf.high * per_unit)),
           effect_type = case_when(is_linear ~ "percentage points",
                                   part == "any" ~ "odds ratio",
                                   TRUE ~ "ratio")) %>%
    select(-per_unit, -is_linear)
  out
}

# Poisson regression (PPML) on the level of consumption, zeros kept: a standard way to
# handle zeros without adding a constant. Period effects; SEs clustered by country.
fit_ppml <- function(df, drug_name, spec = "Poisson (PPML)") {
  d  <- complete_cases(df, drug_name, main_terms, y = "consumption")
  m  <- glm(consumption ~ log_gdp + log_health + life_exp + factor(period),
            family = quasipoisson, data = d)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = d$country))
  crit <- qt(0.975, n_distinct(d$country) - 1)
  tibble(term = rownames(ct), estimate = ct[, 1], std.error = ct[, 2], p.value = ct[, 4]) %>%
    filter(term %in% main_terms) %>%
    mutate(conf.low = estimate - crit * std.error, conf.high = estimate + crit * std.error,
           model = "Poisson pooled", drug = drug_name, spec = spec, part = "single",
           n_obs = nrow(d), n_countries = n_distinct(d$country),
           per_unit    = if_else(str_detect(term, "^log_"), log(step), 1),       # same scale as the other models
           effect      = exp(estimate * per_unit),
           effect_low  = exp(conf.low  * per_unit),
           effect_high = exp(conf.high * per_unit),
           effect_type = "ratio") %>%
    select(-per_unit)
}


# ------------------------------------------------------------------------------
# 5. Primary models (Table 2 and Table 2 continued, Figure 3)
# ------------------------------------------------------------------------------
# Both parts, once for methadone and once for buprenorphine
primary <- bind_rows(map_dfr(drugs, ~ fit_all(dat, .x, part = "amount")),
                     map_dfr(drugs, ~ fit_all(dat, .x, part = "any")))
write_csv(primary, "results/primary_models.csv")                    # Table 2, Figure 3 and Supplementary Table 5

# --- Unadjusted (single-predictor) models (Supplementary Table 7) ---------------
# The same model as the primary one (country random intercept and period effects for the
# amount part; pooled logistic with clustered SEs for the any-use part), but with ONE predictor
# at a time. They use exactly the observations of the primary model (complete for all three
# predictors), so any difference from the primary estimates reflects adjustment, not a
# different sample.
unadjusted <- map_dfr(c("amount", "any"), function(part) {
  y      <- c(amount = "log_cons", any = "any_use")[[part]]
  binary <- part == "any"
  map_dfr(drugs, function(drug_name) {
    d <- complete_cases(dat, drug_name, main_terms, y = y)
    map_dfr(main_terms, function(term) {
      fit_mixed(d, term, y = y, binary = binary) %>%
        mutate(drug = drug_name, spec = "Unadjusted", part = part,
               n_obs = nrow(d), n_countries = n_distinct(d$country))
    })
  })
}) %>%
  mutate(per_unit    = if_else(str_detect(term, "^log_"), log(step), 1),
         effect      = exp(estimate  * per_unit),
         effect_low  = exp(conf.low  * per_unit),
         effect_high = exp(conf.high * per_unit),
         effect_type = if_else(part == "any", "odds ratio", "ratio")) %>%
  select(-per_unit)
write_csv(unadjusted, "results/unadjusted_models.csv")      # Supplementary Table 7


# ------------------------------------------------------------------------------
# 6. Sensitivity analyses (Table 3)
# ------------------------------------------------------------------------------
# For "below-one symbols set to 0.5" we need to re-code consumption from scratch.
dat_lt1 <- code_consumption(long, less_than_one = 0.5) %>% primary_periods()

# For the ANY USE part, the matching check drops the "<<" entries instead of counting them as use.
dat_no_lt1 <- dat %>% mutate(any_use = if_else(symbol %in% "<<", NA_real_, any_use))

# One list entry per sensitivity analysis. Each says which data to use and, where
# relevant, a different set of predictors, extra covariates, or a minimum number of periods.
sensitivity <- list(
  list(spec = "Primary",                        data = dat),

  # Use government instead of total health spending
  list(spec = "Government health spending",     data = dat,
       predictors = c("log_gdp", "log_govt", "life_exp")),

  # Split by income
  list(spec = "High income only",               data = filter(dat, income_group == high_income)),
  list(spec = "Non-high income only",           data = filter(dat, income_group != high_income)),

  # Drop life expectancy (it is closely related to income)
  list(spec = "Without life expectancy",        data = dat,
       predictors = c("log_gdp", "log_health")),

  # Add income group or region as extra covariates
  list(spec = "Adding income group",            data = dat, extras = "income_group"),
  list(spec = "Adding region",                  data = dat, extras = "region"),

  # Drop the US
  list(spec = "Excluding US",                   data = filter(dat, country != us)),

  # Drop flagged values / countries from the screen in section 3
  list(spec = "Screened values omitted",        data = filter(dat, !flag_value)),
  list(spec = "Flagged countries omitted",      data = filter(dat, !flag_country)),

  # Keep countries with at least four complete periods
  list(spec = "At least four complete periods", data = dat, min_periods = 4),

  # AMOUNT part: set "<<" (less than one) to 0.5 instead of missing
  list(spec = "Below-one symbols set to 0.5",   data = dat_lt1,    parts = "amount"),
  # ANY USE part: drop "<<" entries instead of counting them as some use
  list(spec = "Below-one entries excluded",     data = dat_no_lt1, parts = "any")
)

# Fit every specification for both drugs. `%||%` means "use the right side if the
# left side is NULL", so entries without their own predictors fall back to main_terms.
sens <- map_dfr(sensitivity, function(s) {
  map_dfr(s$parts %||% c("amount", "any"), function(part) {      # both parts unless the entry says otherwise
    map_dfr(drugs, ~ fit_all(s$data, .x,
                             predictors  = s$predictors  %||% main_terms,
                             extras      = s$extras,
                             min_periods = s$min_periods,
                             spec        = s$spec,
                             part        = part))
  })
})
write_csv(sens, "results/sensitivity_models.csv")                  # Table 3


# ------------------------------------------------------------------------------
# 6b. Keeping the nil values in a single equation (Supplementary Table 4)
# ------------------------------------------------------------------------------
# Two common one-equation alternatives to the two-part approach:
#   log(1 + consumption) in the same three model types, and Poisson regression (PPML).
single <- bind_rows(
  map_dfr(drugs, ~ fit_all(dat, .x, part = "single", spec = "Single equation: log(1 + consumption)")),
  map_dfr(drugs, ~ fit_ppml(dat, .x))
)
write_csv(single, "results/single_equation_checks.csv")          # Supplementary Table 4

# ------------------------------------------------------------------------------
# 6c. Joint two-part model (extra check, not in the paper)
# ------------------------------------------------------------------------------
# The primary analysis fits the two parts separately. As a sensitivity analysis we fit them
# TOGETHER in one model with a country random intercept in BOTH parts:
#   - a logistic part for "no use" (a nil entry), and
#   - a Gamma part (log link) for the amount among country-periods with positive use.
# glmmTMB::ziGamma does exactly this. Because the logistic part models the chance of NO use,
# its sign is flipped below so that it reads as the odds of ANY use.
#
# Three results per predictor (per 10% higher GDP or health share, per year of life expectancy):
#   "Any use (odds ratio)"           odds of reporting any use
#   "Amount among users (ratio)"     multiplies mean consumption among users
#   "Expected consumption (ratio)"   multiplies expected consumption overall: (chance of any
#                                    use) x (amount among users); interval by simulation
fit_joint <- function(df, drug_name, n_sims = 1000) {
  d <- complete_cases(df, drug_name, main_terms, y = "consumption")   # a number or a dash; "<<" entries are excluded
  d$country <- factor(d$country)

  m <- glmmTMB::glmmTMB(
    consumption ~ log_gdp + log_health + life_exp + period + (1 | country),
    ziformula = ~ log_gdp + log_health + life_exp + period + (1 | country),
    family    = glmmTMB::ziGamma(link = "log"), data = d)
  converged <- isTRUE(m$fit$convergence == 0) && isTRUE(m$sdr$pdHess)

  b_amt  <- glmmTMB::fixef(m)$cond                       # amount among users (log scale)
  b_zero <- glmmTMB::fixef(m)$zi                         # chance of NO use (log-odds)
  co     <- summary(m)$coefficients
  per_unit <- c(log_gdp = log(step), log_health = log(step), life_exp = 1)

  # Part 1 and 2: Wald intervals from the model's own standard errors
  parts <- map_dfr(main_terms, function(tm) {
    k <- per_unit[[tm]]
    a <- co$cond[tm, ]; z <- co$zi[tm, ]
    bind_rows(
      tibble(component = "Any use (odds ratio)", term = tm,
             effect = exp(-z[1] * k), effect_low = exp(-(z[1] + 1.96 * z[2]) * k),
             effect_high = exp(-(z[1] - 1.96 * z[2]) * k)),
      tibble(component = "Amount among users (ratio)", term = tm,
             effect = exp(a[1] * k), effect_low = exp((a[1] - 1.96 * a[2]) * k),
             effect_high = exp((a[1] + 1.96 * a[2]) * k)))
  })

  # Combined effect on expected consumption. Expected consumption for each row is
  #   (1 - chance of no use) x (mean among users),
  # using the fitted country effects. Raise one predictor (+10%, or +1 year) and compare the
  # average. The interval comes from drawing the fixed effects from their joint distribution.
  X   <- model.matrix(~ log_gdp + log_health + life_exp + period, d)
  u_a <- glmmTMB::ranef(m)$cond$country[as.character(d$country), 1]
  u_z <- glmmTMB::ranef(m)$zi$country[as.character(d$country), 1]
  expected <- function(X, ba, bz) mean((1 - plogis(as.numeric(X %*% bz) + u_z)) * exp(as.numeric(X %*% ba) + u_a))
  ratio <- function(tm, ba, bz) {
    X2 <- X; X2[, tm] <- X2[, tm] + per_unit[[tm]]
    expected(X2, ba, bz) / expected(X, ba, bz)
  }
  n_a <- length(b_amt); n_z <- length(b_zero)
  V   <- vcov(m, full = TRUE)[1:(n_a + n_z), 1:(n_a + n_z)]       # amount block, then no-use block
  set.seed(2026)
  draws <- MASS::mvrnorm(n_sims, c(b_amt, b_zero), V)
  combined <- map_dfr(main_terms, function(tm) {
    r <- apply(draws, 1, function(b) ratio(tm, b[1:n_a], b[n_a + 1:n_z]))
    tibble(component = "Expected consumption (ratio)", term = tm,
           effect = ratio(tm, b_amt, b_zero),
           effect_low = unname(quantile(r, 0.025)), effect_high = unname(quantile(r, 0.975)))
  })

  bind_rows(parts, combined) %>%
    mutate(drug = drug_name, n_obs = nrow(d), n_countries = n_distinct(d$country), converged = converged)
}

joint <- map_dfr(drugs, ~ fit_joint(dat, .x))
write_csv(joint, "results/joint_model.csv")                          # extra, not in the paper


# ------------------------------------------------------------------------------
# 6d. Do the between-country and within-country associations differ? (Supplementary Table 6)
# ------------------------------------------------------------------------------
# A Mundlak-type test. Add each predictor's country average to the primary model. The
# coefficient on a country average is the BETWEEN coefficient minus the WITHIN coefficient,
# so testing whether it is zero tests whether the two are equal. If they differ, the primary
# model (which blends them) assumes something that does not hold: that a country's baseline
# is unrelated to its GDP, health spending and life expectancy.
#   AMOUNT part:  linear mixed model, standard errors from the full observed information
#   ANY USE part: pooled logistic model, country-clustered standard errors (log-odds scale)
fit_mundlak <- function(df, drug_name, part = "amount") {
  y      <- c(amount = "log_cons", any = "any_use")[[part]]
  d      <- complete_cases(df, drug_name, main_terms, y = y)
  bar    <- paste0(main_terms, "_bar")
  for (p in main_terms) d[[paste0(p, "_bar")]] <- ave(d[[p]], d$country)     # country average

  tab <- if (part == "any") {
    pooled_logit(reformulate(c(main_terms, bar, "period"), y), d, bar)
  } else {
    f <- reformulate(c(main_terms, bar, "period", "(1 | country)"), y)
    tidy_fixed(lmerTest::lmer(f, data = d, REML = FALSE)) %>% filter(term %in% bar)
  }
  tab %>%
    transmute(term = str_remove(term, "_bar$"),
              difference = estimate, conf.low, conf.high, p.value) %>%     # between minus within (raw scale)
    # The same contrast as a ratio on the reporting scale: (between multiplier) / (within multiplier)
    mutate(per_unit    = if_else(str_detect(term, "^log_"), log(step), 1),
           effect      = exp(difference * per_unit),
           effect_low  = exp(conf.low    * per_unit),
           effect_high = exp(conf.high   * per_unit)) %>%
    select(-per_unit) %>%
    mutate(drug = drug_name, part = part, n_obs = nrow(d), n_countries = n_distinct(d$country))
}
mundlak <- bind_rows(
  map_dfr(drugs, ~ fit_mundlak(dat, .x, "amount")),
  map_dfr(drugs, ~ fit_mundlak(dat, .x, "any"))
)
write_csv(mundlak, "results/between_within_test.csv")      # Supplementary Table 6


# ------------------------------------------------------------------------------
# 6e. Nil entries by income group (Supplementary Table 3)
# ------------------------------------------------------------------------------
# How common are nil entries?
nil_by_group <- dat %>%
  filter(!is.na(consumption)) %>%                       # reported entries only (number or dash)
  group_by(drug, income_group) %>%
  summarise(reported = n(), nil = sum(consumption == 0), .groups = "drop")
nil_all <- nil_by_group %>% group_by(drug) %>%
  summarise(income_group = "All", reported = sum(reported), nil = sum(nil), .groups = "drop")
nil_share <- bind_rows(nil_by_group, nil_all) %>% mutate(pct_nil = 100 * nil / reported)
write_csv(nil_share, "results/nil_share.csv")             # Supplementary Table 3


# ------------------------------------------------------------------------------
# 7. Country-average checks (Supplementary Table 1)
# ------------------------------------------------------------------------------
# These look only at differences BETWEEN countries: average each country's
# variables across the periods in a window, then run ordinary least squares (OLS) with
# one row per country. We keep only the health-share coefficient.
country_ols <- function(df, drug_name, exclude = NULL) {
  avg <- complete_cases(df, drug_name, main_terms) %>%
    filter(country != exclude %||% "") %>%         # optionally leave one country out
    group_by(country) %>%
    summarise(across(c(log_cons, all_of(main_terms)), mean), .groups = "drop")

  m <- lm(reformulate(main_terms, "log_cons"), data = avg)

  # Heteroskedasticity-robust (HC1) standard errors, with normal-based 95% CIs.
  # (With one row per country there are no repeated measures to cluster on, so
  # robust errors are the natural protection against unequal variances.)
  ct <- coeftest(m, vcov = vcovHC(m, type = "HC1"))

  tibble(estimate  = ct["log_health", 1],
         std.error = ct["log_health", 2]) %>%
    mutate(conf.low    = estimate - 1.96 * std.error,
           conf.high   = estimate + 1.96 * std.error,
           effect      = exp(estimate  * log(step)),          # multiplier per 10% higher health share
           effect_low  = exp(conf.low  * log(step)),
           effect_high = exp(conf.high * log(step)),
           n_countries = nrow(avg))
}

# Different ways of grouping the years (the file contains many overlapping windows):
all_periods <- all_coded %>% add_flags() %>% filter(start_year >= 1996)   # every window
three_year  <- filter(all_periods, end_year - start_year == 2, start_year >= 2000)

windows <- list(
  "Between-country OLS (primary periods)" = dat,    # the main country-average model
  # Three-year windows that tile the calendar starting one or two years later
  "Three-year periods starting 2001"      = filter(three_year, (start_year - 2001) %% 3 == 0),
  "Three-year periods starting 2002"      = filter(three_year, (start_year - 2002) %% 3 == 0),
  "All three-year windows"                = three_year,
  "All windows"                           = all_periods   # also the early 5-year windows
)

country_avg <- map_dfr(names(windows), function(w) {
  map_dfr(drugs, ~ country_ols(windows[[w]], .x) %>% mutate(drug = .x, spec = w))
})

# Leave one country out: refit the main country-average model once per country.
loo <- map_dfr(drugs, function(dr) {
  countries <- unique(complete_cases(dat, dr, main_terms)$country)
  map_dfr(countries, ~ country_ols(dat, dr, exclude = .x) %>% mutate(drug = dr, left_out = .x))
})

# Summarise: the range of estimates, and whether every interval excluded zero
loo_summary <- loo %>%
  group_by(drug) %>%
  summarise(min = min(estimate), max = max(estimate),
            effect_min = min(effect), effect_max = max(effect),
            all_exclude_zero = all(conf.low > 0 | conf.high < 0),
            .groups = "drop")

write_csv(country_avg, "results/country_average.csv")                     # Supplementary Table 1
write_csv(loo_summary, "results/country_average_leave_one_out.csv")       # Supplementary Table 1 (last row)


# ------------------------------------------------------------------------------
# 8. Cross-check: the amount models against an independent fit (glmmTMB)
# ------------------------------------------------------------------------------
# glmmTMB fits the same Gaussian random-intercept model by maximum likelihood, and its standard
# errors come from the Hessian over ALL parameters (fixed effects and variances), the same idea as
# tidy_fixed() above. If both fits are right, estimates and standard errors agree closely. The check
# stops the script if they do not.
check_glmmTMB <- map_dfr(drugs, function(dr) {
  d  <- complete_cases(dat, dr, main_terms, y = "log_cons")
  m  <- glmmTMB::glmmTMB(log_cons ~ log_gdp + log_health + life_exp + period + (1 | country),
                         data = d, family = gaussian(), REML = FALSE)
  co <- summary(m)$coefficients$cond
  tibble(drug = dr, term = main_terms,
         est_glmmTMB = co[main_terms, 1], se_glmmTMB = co[main_terms, 2])
}) %>%
  left_join(primary %>%
              filter(part == "amount", model == "Primary mixed") %>%
              select(drug, term, est_mine = estimate, se_mine = std.error),
            by = c("drug", "term")) %>%
  mutate(abs_diff_estimate = abs(est_glmmTMB - est_mine),
         rel_diff_se       = abs(se_glmmTMB - se_mine) / se_mine)

check_glmmTMB

stopifnot(all(check_glmmTMB$abs_diff_estimate < 0.001),   # estimates agree to 3 decimals
          all(check_glmmTMB$rel_diff_se       < 0.01))    # standard errors agree within 1%

message("Done. Results saved in results/")
