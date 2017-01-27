Estimate power for mutational load &lt;&gt; PFS
================
Jacqueline Buros
September 20, 2016

Introduction
------------

In general, we will estimate the power to detect a "significant" effect using simulation.

The process to do this is fairly straightforward:

1.  Write a function to simulate data according to the alternative hypothesis (H1)
2.  Repeat the data simulation X times (we choose X=1000)
3.  Analyze each replicated dataset as you would the real data (ie do the same statistical test)
4.  Compute the power to detect an effect at p &lt; alpha (ie, the proportion of trials with p &lt; alpha)

The Plan
--------

In our case, we have a fixed sample size of ~29 patients, so we treat this as given.

We want to simulate `mutation_load` values for `benefit` and `non-benefit` patients, which will be similar to those observed in the Rosenberg et al 2016 cohort. The only difference in our simulations should be the reduced sample size -- they had a sample size of ~150 patients whereas we will use a sample size of ~29.

The question at hand is, if we saw the same association in our cohort as they saw in the larger cohort, would this be statistically significant?

They observed

1.  median mutation load among benefit patients of 12.4
2.  median mutation load among non-benefit patients of 6.4
3.  approx 20% of patients had a durable reponse

We will use these values in our simulations.

Data simulation
---------------

One challenge in simulating these data is that the Rosenberg et al. paper used a non-parametric statistical test, and so we do not have an easy parametric form from which to simulate these data.

As a starting point, we will try using an exponential distribution to simulate these data, since counts are often distributed exponentially. Note that this is a parameter that we can modify later if it doesn't work well; we will start here and see how reasonable the simulated data are.

``` r
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
```

Notice that we have put default values here which match those we observed in the Rosenberg et al publication.

We can use this function to simulate a dataset like so:

``` r
df <- sim_data()
df %>% head()
```

    ## # A tibble: 6 × 2
    ##           x   group
    ##       <dbl>   <chr>
    ## 1  2.180910 benefit
    ## 2 29.556676 benefit
    ## 3  5.523065 benefit
    ## 4  9.126047 benefit
    ## 5  3.209550 benefit
    ## 6  6.017778 benefit

Here is a plot of these simulated data:

``` r
ggplot(df, aes(y = x, x = group, group = group, colour = group)) +
  geom_boxplot()
```

![](estimate_power_wilcox_rexp_files/figure-markdown_github/sim-data-test-plot-1.png)

Simulating replicate samples
----------------------------

Next we will generate 1000 simulations, each using the same default inputs.

``` r
sim_df <- seq_len(1000) %>%
  map(sim_data, loc_benefit = 12.4, loc_non = 6.4, n = 30, prop_benefit = 0.19) 
```

This results in a list of 1000 datasets.

Computing results for each simulation
-------------------------------------

Next, we summarize the median by DCB-group for each simulation so that we can see how well our simulated data matches the observed values in the Rosenberg publication.

``` r
median_summary <- sim_df %>%
  map_df(.f = ~ .x %>% 
           group_by(group) %>% 
           summarize(median = median(x)) %>% 
           ungroup() %>% 
           spread(key = group, value = median)
         )
```

Finally, we apply the wilcoxon ranksum (aka Mann Whitney U) test to each sample, and join these results to our summary of medians by group.

``` r
## compute mann whitney u (wilcoxon ranksum) 
## statistical test for each replicate
res <- sim_df %>%
  map(~ wilcox.test(x ~ group, data = .)) %>%
  map_df(broom::tidy, .id = 'simulation') %>%
  bind_cols(median_summary)
```

The resulting data.frame contains the following pieces of information:

``` r
str(res)
```

    ## 'data.frame':    1000 obs. of  7 variables:
    ##  $ simulation : chr  "1" "2" "3" "4" ...
    ##  $ statistic  : num  68 99 88 71 106 103 87 62 64 102 ...
    ##  $ p.value    : num  0.67446 0.02267 0.11445 0.5557 0.00546 ...
    ##  $ method     : Factor w/ 1 level "Wilcoxon rank sum test": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ alternative: Factor w/ 1 level "two.sided": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ benefit    : num  10.37 18.09 8.04 6.23 16.59 ...
    ##  $ non        : num  5.98 6.23 4.18 5.39 5.88 ...

Summarizing results
-------------------

The power to detect an effect at *α* &lt;= 0.05 is approximated by the proportion of simulations yielding a p-value &lt; *α*.

It is fairly trivial to calculate this as:

``` r
mean(res$p.value < 0.05)
```

    ## [1] 0.191

We can also plot the power at different levels of alpha, by computing the cumulative proportion of samples with p&lt;=*α*.

![](estimate_power_wilcox_rexp_files/figure-markdown_github/plot-power-by-alpha-1.png)

Evaluating quality of simulated datasets
----------------------------------------

However, it's useful to keep in mind that our power calculation here is only as good as the data simulations we used to generate it.

It's thus equally as important to inspect the simulated data, both graphically and numerically.

Following are some plots that may aid in this process:

![](estimate_power_wilcox_rexp_files/figure-markdown_github/plot-simulated-medians-1.png)

![](estimate_power_wilcox_rexp_files/figure-markdown_github/plot-example-simulations-1.png)
