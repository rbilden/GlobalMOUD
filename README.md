# GlobalMOUD
Contains code and data needed to generate results for the paper "Global Trends in Methadone and Buprenorphine Consumption."

To reproduce
Run Quarto notebook

OR

  1. Put the files in one folder and set it as the R working directory.
  2. Run 01_fit_models.R (or render the .qmd), then 02_tables_figures.R and 03_supplementary_figure1.R.
     Missing packages install on the first run.

Requires R (4.1 or later) and the packages listed at the top of each script.

All code was written with the AI assistance (Claude.ai and ChatGPT) and reviewed by authors.

Please contact Rebecca Bilden (rbilden@stanford.edu) to request the analytic dataset and with questions.

DATA SOURCES

International Narcotics Control Board. Narcotic Drugs: Estimated World Requirements and Statistics [technical report series]. United Nations. Country-level methadone and buprenorphine consumption data for reporting periods spanning 2000–2023.
https://www.incb.org/incb/en/narcotic-drugs/Technical_Reports/narcotic_drugs_reports.html

World Bank. Life expectancy at birth, total (years) [dataset]. World Development Indicators. Indicator SP.DYN.LE00.IN. Underlying sources include United Nations World Population Prospects, national statistical offices, and Eurostat.
https://data.worldbank.org/indicator/SP.DYN.LE00.IN

World Bank. Current health expenditure (% of GDP) [dataset]. World Development Indicators. Indicator SH.XPD.CHEX.GD.ZS. Underlying source is the World Health Organization’s Global Health Expenditure Database.
https://data.worldbank.org/indicator/SH.XPD.CHEX.GD.ZS

World Bank. World Bank Country and Lending Groups [dataset]. FY2027 geographic region and income-group classifications. World Bank Data Help Desk.
https://datahelpdesk.worldbank.org/knowledgebase/articles/906519-world-bank-country-and-lending-groups

World Bank dataset terms
https://www.worldbank.org/ext/en/legal/terms-conditions/datasets
