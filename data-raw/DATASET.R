pkgs <- c("tidyverse", "lubridate", "xts", "fbi")
void <- lapply(pkgs, library, character.only = TRUE)

path = "data-raw/"
pap = pap_factory(path)

# Monthly data ####
data_m = read_csv(pap("m2021-12_20220115.csv"))

# Extract transformations (first row)
di = data_m %>% slice(1) %>% select(-1)
tmp = di %>% rbind(c(5,5,5,5,5,5,5,5,5,5,    # This is an alternative
                     5,5,5,5,5,5,5,5,1,2,    # stationarizing scheme,
                     4,5,5,1,1,1,1,1,1,1,    # which corresponds to
                     5,5,5,5,5,5,5,5,5,5,    # the light transformations
                     5,5,5,5,1,1,1,4,4,4,    # discussed in the paper.
                     4,4,4,4,4,4,4,5,5,5,    # It is more close
                     5,5,2,5,5,5,5,5,7,5,    # to the differencing
                     5,5,2,5,5,1,5,1,1,1,    # strategy employed by Forni
                     1,1,1,1,1,1,1,1,1,1,    # and Gambetti (JME,2010)
                     1,1,1,1,1,4,4,4,4,5,    # such that many more
                     5,5,5,5,5,5,5,5,5,5,    # variables are left in
                     5,5,5,5,5,5,5,5,5,5,    # levels compared to
                     5,5,2,5,5,5,1)) %>% t() # McCracken and Ng (2016))
di = tibble(variable = rownames(tmp),
            transform = tmp[,1],
            transform1 = tmp[,2])

# Parse dates (Month/Day/Year => lubridate::mdy() !!! )
data_m = data_m %>% slice(-1)
data_m = data_m %>%
  mutate(date = lubridate::mdy(sasdate)) %>%
  select(-sasdate)

# Define transformations and adjust info tibble ####
define_trafo = function(number){
  if(number == 1){
    function(x){x}
  } else if(number == 2) {
    function(x){x - lag(x)}
  } else if(number == 3) {
    function(x){x - 2*lag(x) + lag(x, 2)}
  } else if(number == 4) {
    function(x){log(x)}
  } else if(number == 5) {
    function(x){log(x) - lag(log(x))}
  } else if(number == 6) {
    function(x){log(x) - 2*lag(log(x)) + lag(log(x), 2)}
  } else if(number == 7) {
    function(x){x/lag(x) - 1}
  }
}

di %>%
  mutate(fct = map(transform, ~ define_trafo(.x)),
         fct1 = map(transform1, ~ define_trafo(.x))) -> di

# Create nested df ####
dd = data_m %>%
  gather(-date, key = "variable", value = "value") %>%
  group_by(variable) %>%
  nest() %>%
  left_join(di, by = "variable") %>%
  mutate(fct_data = map2(fct, data, ~tibble(date = .y$date,
                                            value = .y$value,
                                            value_trafo = .x(.y$value))),
         fct_data1 = map2(fct1, data, ~tibble(date = .y$date,
                                              value = .y$value,
                                              value_trafo = .x(.y$value)))) %>%
  dplyr::select(-data)

# Prepare and save data with heavy transformations
FRED_heavy <- dd %>%
  dplyr::select(variable, fct_data) %>%
  unnest(fct_data) %>%
  pivot_wider(id_cols = "date", names_from = "variable", values_from = "value_trafo")

FRED_heavy <- prune_data(df = FRED_heavy,
                         start_date = ymd(19730401),
                         end_date = ymd(20071101),
                         impute = TRUE,
                         sdize = FALSE,
                         trans_ix = dd$transform)

usethis::use_data(FRED_heavy, overwrite = TRUE)

# Prepare and save data with light transformations
FRED_light <- dd %>%
  dplyr::select(variable, fct_data1) %>%
  unnest(fct_data1) %>%
  pivot_wider(id_cols = "date", names_from = "variable", values_from = "value_trafo")

FRED_light <- prune_data(df = FRED_light,
                         start_date = ymd(19730401),
                         end_date = ymd(20071101),
                         impute = TRUE,
                         sdize = FALSE,
                         trans_ix = dd$transform1)
usethis::use_data(FRED_light, overwrite = TRUE)

# Prepare and save Forni and Gambetti (2010) data
FG_data <- list(df = readRDS(pap("dataGF.rds")))
FG_data$date <- seq(ymd(19730401), ymd(20071101), by = "month")
FG_data$trans_ix <- as.vector(readRDS(pap("transGF.rds")))
usethis::use_data(FG_data, overwrite = TRUE)

