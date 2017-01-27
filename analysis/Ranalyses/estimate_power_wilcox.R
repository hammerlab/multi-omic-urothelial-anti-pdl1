
library(tidyverse)

## generate simulated data 
sim_data <- function(index = 1, 
                     n = 30,
                     loc_benefit = 12.4,
                     loc_non = 6.4,
                     prop_benefit = 0.2
) {
  n_benefit = n*prop_benefit
  n_non = n*(1-prop_benefit)
  sim_f = function(n, loc) rexp(n = n, rate = 1/loc)
  x_benefit <- sim_f(n = n_benefit, loc = loc_benefit)
  x_non <- sim_f(n = n_non, loc = loc_non)
  df <- tbl_df(list(x = x_benefit, group = 'benefit')) %>%
    bind_rows(tbl_df(list(x = x_non, group = 'non')))
  df
}

## generate 1000 replicates of simulated data
sim_df <- seq_len(1000) %>%
  map(sim_data, loc_benefit = 12.4, loc_non = 6.4, n = 30, prop_benefit = 0.19) 

## summarize median value for each replicate
median_summary <- sim_df %>%
  map_df(.f = ~ .x %>% group_by(group) %>% summarize(median = median(x)) %>% ungroup() %>% spread(key = group, value = median))

## compute mann whitney u (wilcoxon ranksum) 
## statistical test for each replicate
res <- sim_df %>%
  map(~ wilcox.test(x ~ group, data = .)) %>%
  map_df(broom::tidy, .id = 'simulation') %>%
  bind_cols(median_summary)

## summarize power to detect effect at alpha 0.05
mean(res$p.value < 0.05)

## data-prep for plot of power (cum proportion of samples) by alpha (p-value threshold)
sum_res <- res %>%
  dplyr::group_by(`p.value`) %>%
  dplyr::summarize(n_pvalue = n()) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(`p.value`) %>%
  dplyr::mutate(total_n = sum(n_pvalue),
                cum_n = cumsum(n_pvalue),
                power = cum_n/total_n
  ) %>%
  dplyr::rename(alpha = `p.value`)

## plot of power by alpha threshold
ggplot(sum_res, aes(x = alpha, y = power)) + 
  geom_line() + 
  geom_vline(aes(xintercept=0.05), color = 'green', linetype = 'dashed') +
  theme_minimal() +
  ggtitle('Power to detect assoc of mutation load with DCB\nAssuming n = 30, prop = 0.3, and medians = 12.4 vs 6.4')

## plot of simulated median values, for each simulation
ggplot(res %>% 
         gather(benefit, non, value = 'simulated_median', key = 'group'),
       aes(y = simulated_median, x = group, colour=group, group = group)) + 
  geom_boxplot() +
  ggtitle('Distribution of simulated median values by DCB\nShowing results for n=1000 simulations')


## plot of simulated data from a single simulation, selected at random
ggplot(sim_df[runif(min=0, max=1000, n=6)] %>% 
         dplyr::bind_rows(.id='simulation') %>%
         dplyr::left_join(res, by='simulation') %>%
         dplyr::mutate(sim_label = stringr::str_c('Sim ', simulation, ' (p = ', format.pval(`p.value`, digits = 2), ')')),
       aes(y = x, x = group, colour = group, group = group)) + 
  geom_boxplot() +
  facet_wrap(~ sim_label) +
  ggtitle('Distribution of simulated data by DCB\n(showing results for 6 simulations selected at random)')

