rm(list = ls())

library(cmdstanr)
library(posterior)
library(bayesplot)

library(ggplot2)
library(tidyverse)

set_cmdstan_path("C:/cmdstan_dir/cmdstan-2.25.0")

# set directory for where model is (whichever device you are on)
#setwd()
setwd("L:/Lab_JamesR/lachlanW/sleep model/PhillipsSkeldon2017SleepStan")

cmdstan_path()
cmdstan_version()

#check_cmdstan_toolchain(fix=TRUE)

# compile model
mod <- cmdstan_model("Phillips_Skeldon_sleep_2017.stan")

# get data 
source("Data_Phillips_Skeldon_sleep_2017_stan.R")
#data_sleep <- list(T_n = T_n, ts = ts, y0 = y0)
#data_sleep <- list(T_n = T_n, ts = ts, y0 = y0, state_obs = state_test)
data_sleep <- list(T_n = T_n, ts = ts, y0 = y0)#, prop_obs = prop_test)

fit <- mod$sample(
  data = data_sleep,
  seed = 123,
  chains = 4,
  parallel_chains = 2,
  iter_warmup = 150,
  iter_sampling = 300,
  refresh = 150,
)

ests_sum <- fit$summary()

#ests_manual_summary <- summarise_draws(fit$draws(), "mean", "median", "sd", "mad", "quantile2", .args = list(na.rm = TRUE))

mcmc_hist(fit$draws("chi"))
mcmc_hist(fit$draws("prop_sleep"))
mcmc_hist(fit$draws("midsleep[1]")+24)

draws_array <- fit$draws()
str(draws_array)


yest <- data.frame(Vm = ests_sum$mean[c(1:T_n)+2], Vv = ests_sum$mean[c(1:T_n)+2+T_n], H = ests_sum$mean[c(1:T_n)+2+2*T_n], time = ts, state = ests_sum$mean[c(1:T_n)+2+6*T_n])

yest_med <- data.frame(Vm = ests_sum$median[c(1:T_n)+2], Vv = ests_sum$median[c(1:T_n)+2+T_n], H = ests_sum$median[c(1:T_n)+2+2*T_n], time = ts, state = ests_sum$median[c(1:T_n)+2+6*T_n])

yest_med_man <- data.frame(Vm = ests_manual_summary$median[c(1:T_n)+2], Vv = ests_manual_summary$median[c(1:T_n)+2+T_n], H = ests_manual_summary$median[c(1:T_n)+2+2*T_n], time = ts, state = ests_manual_summary$median[c(1:T_n)+2+6*T_n])


yest %>% select(-state) %>% 
  pivot_longer(Vm:H, names_to = "Variable", values_to = "values") %>%
  ggplot(aes(x = time, y = values, group = Variable, colour = Variable)) + 
  #geom_point() + 
  geom_line() + 
  #geom_line(aes(x = time, y = state*2)) + 
  facet_wrap(~Variable, scales = "free", ncol = 1)

yest_med %>% select(-state) %>% 
  pivot_longer(Vm:H, names_to = "Variable", values_to = "values") %>%
  ggplot(aes(x = time, y = values, group = Variable, colour = Variable)) + 
  #geom_point() + 
  geom_line() + 
  #geom_line(aes(x = time, y = state*2)) + 
  facet_wrap(~Variable, scales = "free", ncol = 1)


yest %>% ggplot(aes(x = time, y = state)) +
  geom_line()
yest_med %>% ggplot(aes(x = time, y = state)) +
  geom_line()

saveRDS(yest$state, "state_tester.rds")
