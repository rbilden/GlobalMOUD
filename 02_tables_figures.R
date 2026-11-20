# ==============================================================================
# 02_tables_figures.R
# Builds the manuscript tables (one Word file + a CSV each) and the figures
# (PNG and PDF) from the results written by 01_fit_models.R.
# ------------------------------------------------------------------------------
#   Tables:  Table 1               consumption by income group
#            Table 2               primary and between / within coefficients
#                                  for the AMOUNT part (log consumption among users)
#            Table 2 continued     the same models for the ANY USE part (any use vs nil)
#            Table 3               health-share sensitivity analyses (amount part)
#            Table 3 continued     health-share sensitivity analyses (any-use part)
#            Supplementary Table 1 between-country robustness checks
#            Supplementary Table 2 values flagged by the screen
#            Supplementary Table 3 how common nil entries are, by income group
#            Supplementary Table 4 single-equation checks that keep nil values
#            Supplementary Table 5 country fixed-effects models (both parts)
#            Supplementary Table 6 test of whether between- and within-country associations differ
#            Supplementary Table 7 unadjusted (one predictor at a time) versus adjusted estimates
#            Extra (CSV only)      joint two-part model, not in the paper
#   Figures: Figure 1              consumption by income group (boxplots + dots)
#            Figure 2              US, Canada and other high-income countries over time
#            Figure 3              coefficient plot (income and health share)
#
# Run 01_fit_models.R first, from the same folder, so results/ exists.
# ==============================================================================


# ------------------------------------------------------------------------------
# Packages
# ------------------------------------------------------------------------------
pkgs <- c("tidyverse", "flextable", "officer", "ggbeeswarm", "ragg")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

library(tidyverse)
library(flextable)   # builds nicely formatted tables
library(officer)     # writes those tables into a Word file
library(ggbeeswarm)  # spreads dots out so they don't overlap (Figure 1)

dir.create("tables",  showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)


# ------------------------------------------------------------------------------
# Read the results from script 01
# ------------------------------------------------------------------------------
dat        <- read_csv("results/analysis_data.csv",                 show_col_types = FALSE)  # the analysis dataset
primary_all <- read_csv("results/primary_models.csv",               show_col_types = FALSE)  # Table 2 / Figure 3
sens_all    <- read_csv("results/sensitivity_models.csv",           show_col_types = FALSE)  # Table 3
single_all  <- read_csv("results/single_equation_checks.csv",       show_col_types = FALSE)  # Supp. Table 4
nil_share   <- read_csv("results/nil_share.csv",                    show_col_types = FALSE)  # Supp. Table 3
joint_all   <- read_csv("results/joint_model.csv",                  show_col_types = FALSE)  # Supp. Table 5
bw_test     <- read_csv("results/between_within_test.csv",          show_col_types = FALSE)  # Supp. Table 6
unadjusted  <- read_csv("results/unadjusted_models.csv",            show_col_types = FALSE)  # Supp. Table 7

# Two parts: AMOUNT (how much, among positive values) and ANY USE (any use vs nil)
primary     <- filter(primary_all, part == "amount")
primary_any <- filter(primary_all, part == "any")
sens        <- filter(sens_all,    part == "amount")
sens_any    <- filter(sens_all,    part == "any")
country_av <- read_csv("results/country_average.csv",               show_col_types = FALSE)  # Supp. Table 1
loo        <- read_csv("results/country_average_leave_one_out.csv", show_col_types = FALSE)  # Supp. Table 1 (last row)
flags      <- read_csv("results/flagged_values.csv",                show_col_types = FALSE)  # Supp. Table 2


# ------------------------------------------------------------------------------
# Small helpers and labels
# ------------------------------------------------------------------------------
drugs <- c("methadone", "buprenorphine")

# Nicer names for the drug and income-group codes in the data
drug_label  <- c(methadone = "Methadone", buprenorphine = "Buprenorphine")
group_label <- c(HIC = "High income", UIC = "Upper middle income",
                 LMIC_LIC = "Low or lower middle income")

# Format a number with a fixed number of decimals and a true minus sign (U+2212),
# e.g. -0.39 -> "−0.39". The manuscript uses this character.
num <- function(x, digits = 2) {
  formatC(x, format = "f", digits = digits) %>%
    str_replace("^-", "\u2212") %>%
    str_replace("\\.", "\u00b7")          # Lancet style: mid-dot as the decimal mark
}

# Thousands: no separator for four-digit numbers, a thin space (U+2009) from five digits (SI / Lancet style)
thou <- function(x) {
  x   <- as.character(x)
  big <- !is.na(x) & str_detect(x, "^\\d{5,}$")
  x[big] <- formatC(as.numeric(x[big]), format = "d", big.mark = "\u2009")
  x
}

# Format an estimate and its confidence interval: "0·72 (0·46 to 0·98)"
est_ci <- function(est, lo, hi, digits = 2) {
  str_c(num(est, digits), " (", num(lo, digits), " to ", num(hi, digits), ")")
}

# Format big numbers for Table 1: one decimal, thin-space thousands from five digits, and no ".0"
# (so 51 stays "51", 657.5 becomes "657·5" and 2345.25 becomes "2345·2").
big <- function(x) {
  s   <- str_remove(formatC(x, format = "f", digits = 1), "\\.0$")      # one decimal, no ".0"
  int <- str_extract(s, "^\\d+")
  str_replace(s, "^\\d+", thou(int)) %>% str_replace("\\.", "\u00b7")  # thin-space thousands, mid-dot decimal
}

# Turn a data frame into a clean, journal-style Word table
make_flextable <- function(df) {
  flextable(df) %>%
    theme_booktabs() %>%            # horizontal rules only, no vertical lines
    bold(part = "header") %>%
    fontsize(size = 9, part = "all") %>%
    autofit()
}


# ==============================================================================
# TABLE 1: consumption by income group
# ==============================================================================
# Median (interquartile range) of the positive numeric country-period values,
# plus how many countries and observations contribute.
positive <- dat %>% filter(consumption > 0)

# One function so the by-group rows and the "All countries" row are computed identically
summarise_group <- function(df) {
  df %>% summarise(
    median      = median(consumption),
    q1          = quantile(consumption, 0.25),
    q3          = quantile(consumption, 0.75),
    n_countries = n_distinct(country),
    n_obs       = n(),
    .groups = "drop")
}

table1 <- bind_rows(
  positive %>% group_by(drug, income_group) %>% summarise_group(),                    # by income group
  positive %>% group_by(drug)               %>% summarise_group() %>%                 # all countries
    mutate(income_group = "All")
) %>%
  mutate(
    # Paste everything into one cell, e.g. "657·5 (53·5–2345·2) 63 countries / 428 observations"
    cell = str_c(big(median), " (", big(q1), "\u2013", big(q3), ") ",
                 n_countries, " countries / ", n_obs, " observations"),
    income_group = recode(income_group, !!!group_label, All = "All countries")
  ) %>%
  select(income_group, drug, cell) %>%
  pivot_wider(names_from = drug, values_from = cell) %>%           # one column per drug
  arrange(match(income_group, c(group_label, "All countries"))) %>% # keep the usual row order
  select(`Income group` = income_group,
         Methadone = methadone, Buprenorphine = buprenorphine)


# ==============================================================================
# TABLE 2: ratio of consumption (amount among users), one panel per drug
# ==============================================================================
model_order      <- c("Primary mixed", "Between", "Within", "Country FE")   # all four (used in Figure 3)
model_order_main <- c("Primary mixed", "Between", "Within")                  # main-text table columns
# (The country fixed-effects results are in Supplementary Table 5: they answer almost the same
# question as the "within" column, comparing each country with itself over time.)
term_label  <- c(log_gdp = "GDP per person",
                 log_health = "Health share",
                 life_exp = "Life expectancy")

make_table2 <- function(drug_name) {
  primary %>%
    filter(drug == drug_name, model %in% model_order_main) %>%
    mutate(
      # `effect` = ratio of consumption per 10% higher GDP / health share, or per year of life expectancy
      cell      = est_ci(effect, effect_low, effect_high),
      Predictor = factor(term_label[term], levels = term_label),
      model     = factor(model, levels = model_order_main)
    ) %>%
    select(Predictor, model, cell) %>%
    pivot_wider(names_from = model, values_from = cell) %>%   # one column per model type
    arrange(Predictor) %>%
    mutate(Predictor = as.character(Predictor))
}

# Panel heading, e.g. "Methadone — 661 observations / 110 countries"
n_label <- function(drug_name) {
  p <- primary %>% filter(drug == drug_name) %>% slice(1)
  str_c(drug_label[[drug_name]], " \u2014 ", p$n_obs, " observations / ", p$n_countries, " countries")
}

table2_meth <- make_table2("methadone")
table2_bup  <- make_table2("buprenorphine")


# ==============================================================================
# TABLE 3: health-share ratio (per 10% higher) under each sensitivity analysis
# ==============================================================================
# Order the rows the way they appear in the manuscript
spec_order <- c("Primary", "Government health spending", "High income only",
                "Non-high income only", "Without life expectancy", "Adding income group",
                "Adding region", "Excluding US", "Screened values omitted",
                "Flagged countries omitted", "At least four complete periods",
                "Below-one symbols set to 0.5")

table3 <- sens %>%
  # Take the primary mixed model, and the "health share" term. In the government-spending
  # run that term is called log_govt, so we keep both names.
  filter(model == "Primary mixed", term %in% c("log_health", "log_govt")) %>%
  mutate(cell = est_ci(effect, effect_low, effect_high),
         n    = str_c(n_obs, " / ", n_countries)) %>%
  select(spec, drug, cell, n) %>%
  pivot_wider(names_from = drug, values_from = c(cell, n)) %>%   # cell_methadone, n_methadone, ...
  arrange(match(spec, spec_order)) %>%
  select(Specification = spec,
         `Methadone ratio per 10% (95% CI)`      = cell_methadone,
         `Methadone N obs / countries`           = n_methadone,
         `Buprenorphine ratio per 10% (95% CI)`  = cell_buprenorphine,
         `Buprenorphine N obs / countries`       = n_buprenorphine)


# ==============================================================================
# TABLE 2 CONTINUED: the ANY USE part (did the country report any use?)
# ==============================================================================
# Columns are odds ratios: for GDP and health share, per 10% higher; for life expectancy, per YEAR.
model_order_any <- c("Primary logistic", "Between", "Within")

make_table2_any <- function(drug_name) {
  primary_any %>%
    filter(drug == drug_name, model %in% model_order_any) %>%
    mutate(
      digits = 2,
      cell   = map_chr(seq_along(effect), ~ est_ci(effect[.x], effect_low[.x], effect_high[.x], digits[.x])),
      Predictor = factor(term_label[term], levels = term_label),
      model     = factor(model, levels = model_order_any)
    ) %>%
    select(Predictor, model, cell) %>%
    pivot_wider(names_from = model, values_from = cell) %>%
    arrange(Predictor) %>%
    mutate(Predictor = as.character(Predictor))
}
n_label_any <- function(drug_name) {
  p <- primary_any %>% filter(drug == drug_name) %>% slice(1)
  str_c(drug_label[[drug_name]], " \u2014 any use: ", p$n_obs, " observations / ", p$n_countries, " countries")
}
table2_any_meth <- make_table2_any("methadone")
table2_any_bup  <- make_table2_any("buprenorphine")

# TABLE 3 CONTINUED: health-share sensitivity analyses for the ANY USE part
spec_order_any <- c(setdiff(spec_order, "Below-one symbols set to 0.5"), "Below-one entries excluded")
table3_any <- sens_any %>%
  filter(model == "Primary logistic", term %in% c("log_health", "log_govt")) %>%
  mutate(cell = est_ci(effect, effect_low, effect_high),
         n    = str_c(n_obs, " / ", n_countries)) %>%
  select(spec, drug, cell, n) %>%
  pivot_wider(names_from = drug, values_from = c(cell, n)) %>%
  arrange(match(spec, spec_order_any)) %>%
  select(Specification = spec,
         `Methadone odds ratio per 10% (95% CI)`     = cell_methadone,
         `Methadone N obs / countries`                    = n_methadone,
         `Buprenorphine odds ratio per 10% (95% CI)` = cell_buprenorphine,
         `Buprenorphine N obs / countries`                = n_buprenorphine)


# ==============================================================================
# SUPPLEMENTARY TABLE 1: between-country robustness
# ==============================================================================
supp1 <- country_av %>%
  mutate(cell = est_ci(effect, effect_low, effect_high)) %>%
  select(Specification = spec, drug, cell) %>%
  pivot_wider(names_from = drug, values_from = cell) %>%
  select(Specification, Methadone = methadone, Buprenorphine = buprenorphine) %>%
  # Add the leave-one-country-out row: the range of estimates across all the refits
  bind_rows(tibble(
    Specification = "Leave one country out",
    Methadone = loo %>% filter(drug == "methadone") %>%
      with(str_c(num(effect_min), "\u2013", num(effect_max),
                 if_else(all_exclude_zero, "; all intervals exclude 1", ""))),
    Buprenorphine = loo %>% filter(drug == "buprenorphine") %>%
      with(str_c(num(effect_min), "\u2013", num(effect_max),
                 if_else(all_exclude_zero, "; all intervals exclude 1", "")))
  ))


# ==============================================================================
# SUPPLEMENTARY TABLE 2: values flagged by the screen
# ==============================================================================
supp2 <- flags %>%
  arrange(match(drug, drugs), country) %>%        # methadone first, then alphabetical
  transmute(Country    = country,
            Period     = str_replace(period, "-", "\u2013"),    # en dash between years
            Medication = drug_label[drug],
            Value      = thou(format(value, scientific = FALSE, trim = TRUE)),   # thin-space thousands from five digits
            Flag       = flag)


# ==============================================================================
# SUPPLEMENTARY TABLE 3: how common are nil entries?
# ==============================================================================
# Nil entries (dashes) as a share of all reported entries (a number or a dash).
supp3 <- nil_share %>%
  mutate(income = recode(income_group, !!!group_label, All = "All countries"),
         cell   = str_c(nil, " of ", reported, " (", floor(pct_nil + 0.5), "%)")) %>%   # round half up (R rounds 80.5 to 80)
  select(income, drug, cell) %>%
  pivot_wider(names_from = drug, values_from = cell) %>%
  arrange(match(income, c(group_label, "All countries"))) %>%
  select(`Income group` = income, Methadone = methadone, Buprenorphine = buprenorphine)

# ==============================================================================
# SUPPLEMENTARY TABLE 4: single-equation checks that keep the nil values
# ==============================================================================
# Health-share coefficient (95% CI). The first row is the primary AMOUNT model for comparison.
supp4_long <- bind_rows(
  primary    %>% filter(model == "Primary mixed",  term == "log_health") %>% mutate(label = "Amount among users (primary model)"),
  single_all %>% filter(model == "Primary mixed",  term == "log_health") %>% mutate(label = "log(1 + consumption), mixed model"),
  single_all %>% filter(model == "Between",        term == "log_health") %>% mutate(label = "log(1 + consumption), between countries"),
  single_all %>% filter(model == "Within",         term == "log_health") %>% mutate(label = "log(1 + consumption), within countries"),
  single_all %>% filter(model == "Country FE",     term == "log_health") %>% mutate(label = "log(1 + consumption), country fixed effects"),
  single_all %>% filter(model == "Poisson pooled", term == "log_health") %>% mutate(label = "Poisson (PPML), clustered standard errors")
)
supp4 <- supp4_long %>%
  mutate(cell = est_ci(effect, effect_low, effect_high),
         n    = str_c(n_obs, " / ", n_countries)) %>%
  select(label, drug, cell, n) %>%
  pivot_wider(names_from = drug, values_from = c(cell, n)) %>%
  arrange(match(label, unique(supp4_long$label))) %>%
  select(Specification = label,
         `Methadone ratio per 10% (95% CI)`    = cell_methadone,
         `Methadone N obs / countries`        = n_methadone,
         `Buprenorphine ratio per 10% (95% CI)` = cell_buprenorphine,
         `Buprenorphine N obs / countries`    = n_buprenorphine)


# ==============================================================================
# EXTRA (not in the paper): joint two-part model
# ==============================================================================
# Kept for the analysis record; written to tables/extra_joint_model.csv but not to the Word file.
# One model with a country random intercept in both parts. Effects are per 10% higher GDP
# and health share, and per YEAR of life expectancy.
component_order <- c("Any use (odds ratio)", "Amount among users (ratio)", "Expected consumption (ratio)")
extra_joint <- joint_all %>%
  mutate(Specification = str_c(term_label[term], ": ", str_to_lower(str_sub(component, 1, 1)), str_sub(component, 2)),
         cell = est_ci(effect, effect_low, effect_high)) %>%
  select(Specification, term, component, drug, cell) %>%
  pivot_wider(names_from = drug, values_from = cell) %>%
  arrange(match(term, names(term_label)), match(component, component_order)) %>%
  select(Specification, Methadone = methadone, Buprenorphine = buprenorphine)
joint_n <- joint_all %>% distinct(drug, n_obs, n_countries)     # sample sizes for the caption


# ==============================================================================
# SUPPLEMENTARY TABLE 5: country fixed-effects models (both parts)
# ==============================================================================
# Every country is compared only with itself over time (country and period effects, standard
# errors clustered by country). AMOUNT part: ratio of consumption. ANY USE part: a
# linear probability model, shown as the change in percentage POINTS in the probability of any use
# per 10% higher GDP and health share, or per year of life expectancy.
fe_amount <- primary %>%
  filter(model == "Country FE") %>%
  mutate(cell   = est_ci(effect, effect_low, effect_high),
         column = str_c(drug_label[drug], ": amount")) %>%
  select(term, column, cell)
fe_any <- primary_any %>%
  filter(model == "Country FE") %>%
  mutate(cell   = est_ci(effect, effect_low, effect_high, 2),
         column = str_c(drug_label[drug], ": any use (percentage points)")) %>%
  select(term, column, cell)
supp5 <- bind_rows(fe_amount, fe_any) %>%
  pivot_wider(names_from = column, values_from = cell) %>%
  arrange(match(term, names(term_label))) %>%
  mutate(Predictor = term_label[term]) %>%
  select(Predictor, `Methadone: amount`, `Methadone: any use (percentage points)`,
         `Buprenorphine: amount`, `Buprenorphine: any use (percentage points)`)


# ==============================================================================
# SUPPLEMENTARY TABLE 6: do the between- and within-country associations differ?
# ==============================================================================
# RATIO of the between-country to the within-country association (95% CI), with the p-value of the
# test that they are equal. A ratio far from 1 means the primary model, which blends the two, is not a
# good summary. AMOUNT: ratio of the consumption multipliers. ANY USE: ratio of the odds ratios.
p_fmt <- function(p) if_else(p < 0.001, "p < 0\u00b7001", str_c("p = ", num(p, 3)))
supp6 <- bw_test %>%
  mutate(cell   = str_c(est_ci(effect, effect_low, effect_high), "; ", p_fmt(p.value)),
         column = str_c(drug_label[drug], ": ", if_else(part == "amount", "amount", "any use"))) %>%
  select(term, column, cell) %>%
  pivot_wider(names_from = column, values_from = cell) %>%
  arrange(match(term, names(term_label))) %>%
  mutate(Predictor = term_label[term]) %>%
  select(Predictor, `Methadone: amount`, `Buprenorphine: amount`,
         `Methadone: any use`, `Buprenorphine: any use`)


# ==============================================================================
# SUPPLEMENTARY TABLE 7: unadjusted (one predictor at a time) versus adjusted estimates
# ==============================================================================
# Same model and same observations as the primary model. "Adjusted" is the primary model
# (all three predictors together); "Unadjusted" fits each predictor on its own.
supp7 <- bind_rows(
  unadjusted %>% mutate(analysis = "Unadjusted"),
  bind_rows(primary %>% filter(model == "Primary mixed"),
            primary_any %>% filter(model == "Primary logistic")) %>% mutate(analysis = "Adjusted")
) %>%
  filter(term %in% names(term_label)) %>%
  mutate(cell   = est_ci(effect, effect_low, effect_high),
         column = str_c(drug_label[drug], ": ", if_else(part == "amount", "amount", "any use"))) %>%
  select(term, analysis, column, cell) %>%
  pivot_wider(names_from = column, values_from = cell) %>%
  arrange(match(term, names(term_label)), match(analysis, c("Unadjusted", "Adjusted"))) %>%
  mutate(Predictor = term_label[term]) %>%
  select(Predictor, Analysis = analysis, `Methadone: amount`, `Buprenorphine: amount`,
         `Methadone: any use`, `Buprenorphine: any use`)


# ==============================================================================
# WRITE THE TABLES: one Word document, plus a CSV of each table
# ==============================================================================
# Adds a heading, a caption and the table to the Word document
add_table <- function(doc, title, caption, ft) {
  doc %>%
    body_add_par(title, style = "heading 2") %>%
    body_add_par(caption, style = "Normal") %>%
    body_add_flextable(ft) %>%
    body_add_par("", style = "Normal")            # blank line between tables
}

doc <- read_docx() %>%
  add_table("Table 1",
            "Median (interquartile range) of consumption in S-DDD per million inhabitants per day, with the number of contributing countries and observations. Positive numeric country-period observations in the eight non-overlapping 3-year periods, 2000\u20132023.",
            make_flextable(table1)) %>%
  # The panel heading (drug + sample size) sits in a merged row above the column names
  add_table("Table 2 (methadone)",
            "Ratio of consumption (95% CI) among country-periods with positive use, per 10% higher GDP per person or health share and per year of life expectancy. Country fixed-effects results are in Supplementary Table 5.",
            make_flextable(table2_meth) %>% add_header_row(values = n_label("methadone"), colwidths = 4)) %>%
  add_table("Table 2 (buprenorphine)", "",
            make_flextable(table2_bup) %>% add_header_row(values = n_label("buprenorphine"), colwidths = 4)) %>%
  add_table("Table 2 continued (any use, methadone)",
            "Whether any use was reported (1) or nil (0). Odds ratios (95% CI): for GDP and health share, per 10% higher; for life expectancy, per year. Country fixed-effects results are in Supplementary Table 5.",
            make_flextable(table2_any_meth) %>% add_header_row(values = n_label_any("methadone"), colwidths = 4)) %>%
  add_table("Table 2 continued (any use, buprenorphine)", "",
            make_flextable(table2_any_bup) %>% add_header_row(values = n_label_any("buprenorphine"), colwidths = 4)) %>%
  add_table("Table 3",
            "Ratio of consumption per 10% higher health-spending share (95% CI) from the primary mixed model (amount among users) under each sensitivity analysis, with number of observations / countries.",
            make_flextable(table3)) %>%
  add_table("Table 3 continued (any use)",
            "Odds ratio (95% CI) for any reported use per 10% higher health-spending share, from the primary logistic model under each sensitivity analysis.",
            make_flextable(table3_any)) %>%
  add_table("Supplementary Table 1",
            "Between-country ratio of consumption per 10% higher health-spending share (95% CI) from OLS regressions of country averages (heteroskedasticity-robust HC1 standard errors), adjusted for GDP and life expectancy.",
            make_flextable(supp1)) %>%
  add_table("Supplementary Table 2",
            "Values flagged by the diagnostic screen and retained in the primary models. Candidates for source verification, not confirmed errors.",
            make_flextable(supp2)) %>%
  add_table("Supplementary Table 3",
            "Nil (dash) entries as a share of reported country-period entries, by income group.",
            make_flextable(supp3)) %>%
  add_table("Supplementary Table 4",
            "Ratio per 10% higher health-spending share (95% CI) when nil values are kept in a single equation, compared with the amount-among-users model. For log(1 + consumption) the ratio applies to 1 + consumption.",
            make_flextable(supp4)) %>%
  add_table("Supplementary Table 5",
            "Country fixed-effects models (country and period effects; standard errors clustered by country). Amount: ratio of consumption (95% CI) among country-periods with positive use, per 10% higher GDP and health share and per year of life expectancy. Any use: a linear probability model, in percentage points (95% CI) of the probability of any use per 10% higher GDP and health share, or per year of life expectancy.",
            make_flextable(supp5)) %>%
  add_table("Supplementary Table 6",
            "Test of whether the between-country and within-country associations differ: the ratio of the between-country to the within-country association (95% CI) and p-value, from the primary model with each predictor's country average added. Amount: ratio of consumption multipliers among country-periods with positive use. Any use: ratio of odds ratios.",
            make_flextable(supp6)) %>%
  add_table("Supplementary Table 7",
            "Unadjusted and adjusted associations. Unadjusted: each predictor fitted on its own in the primary model; adjusted: the primary model with GDP per person, health share and life expectancy together. Same observations in both. Amount: ratio of consumption (95% CI) among country-periods with positive use, per 10% higher GDP or health share or per year of life expectancy. Any use: odds ratio (95% CI) for any reported use, same units.",
            make_flextable(supp7))

print(doc, target = "tables/tables.docx")

# CSV copies, handy for checking numbers or pasting into other software
walk2(list(table1, table2_meth, table2_bup, table2_any_meth, table2_any_bup, table3, table3_any,
           supp1, supp2, supp3, supp4, extra_joint, supp5, supp6, supp7),
      c("table1", "table2_methadone", "table2_buprenorphine", "table2_any_methadone",
        "table2_any_buprenorphine", "table3", "table3_any", "supp_table1", "supp_table2",
        "supp_table3", "supp_table4", "extra_joint_model", "supp_table5", "supp_table6", "supp_table7"),
      ~ write_csv(.x, str_glue("tables/{.y}.csv")))


# ==============================================================================
# FIGURES
# ==============================================================================
# A clean, simple look shared by all figures
# Axis labels in Lancet style: mid-dot decimals, no separator for four-digit numbers, thin space from five digits
lancet_int <- function(x) thou(format(x, scientific = FALSE, trim = TRUE))
lancet_dec <- function(x) str_replace(format(x, nsmall = 1, trim = TRUE), "\\.", "\u00b7")

theme_set(
  theme_classic(base_size = 11) +
    theme(strip.background = element_blank(),               # no grey box around facet titles
          strip.text       = element_text(face = "bold"),
          legend.position  = "bottom")
)

# Save each figure twice: PNG for Word/slides, PDF (vector) for journals.
# cairo_pdf is used so the en dashes in period labels print correctly.
save_figure <- function(plot, name, width, height) {
  ggsave(str_glue("figures/{name}.png"), plot, width = width, height = height,
         dpi = 300, device = ragg::agg_png)
  ggsave(str_glue("figures/{name}.pdf"), plot, width = width, height = height,
         device = cairo_pdf)
}


# ---- Figure 1: consumption by income group -----------------------------------
# Each dot is one positive country-period observation; boxes show the median and
# interquartile range. The y axis is logarithmic. The number in parentheses is the
# number of contributing countries.
fig1_data <- positive %>%
  mutate(
    group = factor(income_group, levels = c("HIC", "UIC", "LMIC_LIC"),
                   labels = c("HIC", "UMIC", "LMIC/LIC")),
    drug  = factor(drug_label[drug], levels = drug_label)
  ) %>%
  group_by(drug, group) %>%
  mutate(label = str_c(group, "\n(", n_distinct(country), ")")) %>%   # e.g. "HIC\n(63)"
  ungroup() %>%
  arrange(drug, group) %>%
  mutate(label = factor(label, levels = unique(label)))               # keep HIC, UMIC, LMIC order

fig1 <- ggplot(fig1_data, aes(label, consumption)) +
  geom_quasirandom(colour = "grey55", alpha = 0.5, size = 0.8) +      # the individual dots
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5) +          # box on top; dots already show outliers
  scale_y_log10(labels = lancet_int) +
  facet_wrap(~ drug, scales = "free_x") +                              # free_x so each panel shows its own counts
  labs(x = NULL, y = "S-DDD per million inhabitants per day (log scale)")

save_figure(fig1, "figure1_income_groups", 8, 4.5)


# ---- Figure 2: US, Canada and other high-income countries ---------------------
# The peer median and interquartile range come from script 01 (peers = high-income
# countries other than the US and Canada with at least `peer_min_periods` numeric periods).
peers   <- read_csv("results/peer_series.csv",     show_col_types = FALSE) %>%
  mutate(drug = factor(drug_label[drug], levels = drug_label), period = str_replace(period, "-", "\u2013"))
us_can  <- read_csv("results/us_canada_series.csv", show_col_types = FALSE) %>%
  mutate(drug = factor(drug_label[drug], levels = drug_label), period = str_replace(period, "-", "\u2013"))

fig2 <- ggplot() +
  # grey band = interquartile range of the peers
  geom_ribbon(data = peers, aes(period, ymin = q1, ymax = q3, group = drug), fill = "grey80") +
  # dark line = peer median
  geom_line(data = peers,
            aes(period, median, group = drug, colour = "Other high-income median (IQR)"),
            linewidth = 0.9) +
  # coloured lines = the US and Canada
  geom_line(data = us_can, aes(period, consumption, group = country, colour = country), linewidth = 0.9) +
  scale_colour_manual(values = c(US = "#b2182b", Canada = "#2166ac",
                                 "Other high-income median (IQR)" = "grey25"),
                      name = NULL) +
  scale_y_continuous(labels = lancet_int) +
  facet_wrap(~ drug) +
  labs(x = NULL, y = "S-DDD per million inhabitants per day") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))            # tilt the long period labels

save_figure(fig2, "figure2_us_canada_peers", 9, 4.5)


# ---- Figure 3: coefficient plot ------------------------------------------------
# Points are estimates; bars are 95% CIs. Open points are not statistically
# significant (p >= 0.05). Shows the primary model and the between / within split; the country
# fixed-effects results are in Supplementary Table 5.
fig3_data <- primary %>%
  filter(term %in% c("log_gdp", "log_health"),            # income and health spending only
         model %in% model_order_main) %>%
  mutate(
    model       = factor(model, levels = rev(model_order_main)),   # rev() so "Primary mixed" is on top
    predictor   = factor(term_label[term], levels = term_label[c("log_gdp", "log_health")]),
    drug        = factor(drug_label[drug], levels = drug_label),
    significant = p.value < 0.05
  )

fig3 <- ggplot(fig3_data, aes(effect, model)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +        # "no association" line (ratio = 1)
  geom_errorbarh(aes(xmin = effect_low, xmax = effect_high), height = 0.2) +
  geom_point(aes(fill = significant), shape = 21, size = 2.8, colour = "black") +
  scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none") +
  scale_x_log10(breaks = c(0.9, 1, 1.1, 1.2, 1.3, 1.4), labels = lancet_dec) +                    # a ratio scale: equal distances are equal multiples
  facet_grid(predictor ~ drug) +                                              # rows = predictor, columns = drug
  labs(x = "Ratio of consumption per 10% higher GDP or health share (95% CI)", y = NULL)

save_figure(fig3, "figure3_coefficients", 8, 4.2)

message("Done. See tables/ and figures/")
