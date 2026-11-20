# ==============================================================================
# 03_supplementary_figure1.R
# Supplementary Figure 1: mean methadone and buprenorphine consumption in
# low- or lower-middle-income (LMIC/LIC) and upper-middle-income (UMIC) countries,
# using overlapping three-year windows (2000-2002 to 2021-2023).
# ------------------------------------------------------------------------------
# Reads the compiled panel directly, so it does not need 01_fit_models.R.
# The mean is taken over positive numeric country-period observations (the same
# observations that Table 1 and Figure 1 summarise): nil (dash) and "<<" entries
# are not in the mean. The five-year windows 1996-2000 and 1997-2001 are left out
# because they start before the 2000-2023 study period.
# Writes figures/supp_figure1_lower_income_means.pdf and .png
# ==============================================================================
pkgs <- c("tidyverse", "ragg")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)
library(tidyverse)

dir.create("figures", showWarnings = FALSE)

infile    <- "incb_worldbank_panel_excel.csv"
panel_raw <- read_csv(infile, show_col_types = FALSE)

supp_fig1_windows <- panel_raw %>%
  filter(!is.na(income_type), income_type %in% c("LMIC_LIC", "UIC")) %>%
  mutate(start_year = as.integer(substr(time, 1, 4)),
         end_year   = as.integer(substr(time, 6, 9))) %>%
  filter(start_year >= 2000, end_year - start_year == 2)           # overlapping 3-year windows

supp_fig1_data <- bind_rows(
  supp_fig1_windows %>% transmute(country, income_type, time, drug = "Methadone",     value = q_m),
  supp_fig1_windows %>% transmute(country, income_type, time, drug = "Buprenorphine", value = q_b)
) %>%
  filter(!is.na(value), value > 0) %>%                             # positive numeric entries only
  group_by(income_type, drug, time) %>%
  summarise(mean_sddd = mean(value), n_countries = n_distinct(country), .groups = "drop") %>%
  mutate(
    time        = str_replace(time, "-", "\u2013"),
    income_type = factor(income_type, levels = c("LMIC_LIC", "UIC"),
                         labels = c("LMIC/LIC", "UMIC")),
    drug        = factor(drug, levels = c("Methadone", "Buprenorphine"))
  )

write_csv(supp_fig1_data, "figures/supp_figure1_data.csv")                   # the plotted values

# The same formatting helpers and figure style as in 02_tables_figures.R
# Thousands: no separator for four-digit numbers, a thin space (U+2009) from five digits (SI / Lancet style)
thou <- function(x) {
  x   <- as.character(x)
  big <- !is.na(x) & str_detect(x, "^\\d{5,}$")
  x[big] <- formatC(as.numeric(x[big]), format = "d", big.mark = "\u2009")
  x
}

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


supp_fig1 <- ggplot(supp_fig1_data, aes(time, mean_sddd, colour = drug, group = drug)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c(Methadone = "#d95f02", Buprenorphine = "#2166ac"), name = NULL) +
  scale_y_continuous(labels = lancet_int) +
  facet_wrap(~ income_type, scales = "free_y") +
  labs(x = NULL, y = "Mean S-DDD per million inhabitants per day") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))

save_figure(supp_fig1, "supp_figure1_lower_income_means", 9, 4.8)
message("Done. See figures/")
