The health department is interested in whether dengue can invade the city if it is introduced and our modelling group is interested in the response from the model with /app/data/model_specification.md and local data, rather than literature data.

Constant-temperature Aedes aegypti life-history and vector-competence experiments that are maintained across three laboratories have been pooled into /app/data/trait_experiments.csv; /app/data/lab_registry.csv explains how each laboratory has recorded measurements. The data for daily mean temperatures are contained in /app/data/station_daily_tmean.csv which is located in the city weather station. The density-dependence coefficient δ is not known; the only information on mosquito density is the indoor survey on 2021 (aspirator_survey_2021.csv /app/data/aspirator_survey_2021.csv).

Please write the following to /app/output/:

- `trait_fits.csv` with columns `trait`, `q`, `T0`, `Tm`: one row for each of the eight traits in the specification, parameters as defined there.
- `r0_constant_temperature.csv` with columns `temperature_c`, `females_per_person`, `r0` for constant temperatures from 16.0 to 34.0 °C in 0.5 °C steps: the equilibrium number of adult female mosquitoes per resident and the basic reproduction number if the climate were held at that temperature all year.
- `seasonal_summary.json` with `mean_temperature_c` (annual mean of the city's mean annual temperature cycle), `mean_females_per_person` and `peak_density_day_of_year` (annual mean, and day of year 1–365 of the maximum, of adult females per resident in the disease-free state under that cycle), `r0_at_mean_temperature` (the constant-temperature value at that annual mean), and `r0_seasonal`, the basic reproduction number of the model under the seasonal cycle.
