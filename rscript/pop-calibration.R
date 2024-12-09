if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5")
install.packages("tidyr")
library(rhdf5)
library(dplyr)
library(ggplot2)
library(colorspace)
library(RColorBrewer)
library(cowplot)
library(here)
library(tidyr)




load_age_data <- function(no_runs=10, duration=100) {
  df <- NULL
  max_seed <- no_runs - 1
  for (gr in growth_rates) {
    for (i in 0:max_seed) {
      path = sprintf("main/output/params_growth_rate=%s/seed_%s/population.hd5",gr, i)
      if ( file.exists( path ) ) {
        cur_df <- h5read(path, "population/age_dists", compoundAsDataFrame=FALSE)
        print(i)
        cur_df <- cur_df[["dist"]][,duration]
        pop <- sum(cur_df)
        
        cur_df <- data.frame(pop=cur_df/pop) |>
          mutate(age = row_number()-1, seed=i) |>
          mutate(age_group= cut(age, breaks=c(seq(0,80,5),105),
                                right=FALSE)) |>
          group_by(age_group) |>
          summarise(pop=sum(pop), growth_rate=gr, .groups="drop")
        
        if (is.null(df)) {
          df <- cur_df
        } else {
          df <- rbind(df, cur_df)
        }
      }
    }
  }
  df |>
    group_by(age_group, growth_rate) |>
    summarise(pop=mean(pop), .groups="drop")
}

load_hh_data <- function(no_runs=10, duration=100) {
  df <- NULL
  max_seed <- no_runs - 1
  for (c_prob in couple_probs) {
    for (l_prob in leaving_probs) {
      for (d_prob in divorce_probs) {
        for (i in 0:max_seed) {
          path = sprintf("main/output/params_couple_prob=%s_leaving_prob=%s_divorce_prob=%s/seed_%s/population.hd5",
                         c_prob, l_prob, d_prob, i)
          if ( file.exists( path ) ) {
            cur_df <- h5read(path, "population/hh_size_dists_count", compoundAsDataFrame=FALSE)
            cur_df <- cur_df[["dist"]][,duration]
            
            cur_df <- data.frame(hh_per=cur_df/sum(cur_df)) |>
              mutate(size = row_number()) |>
              mutate(size_ed= cut(size, breaks=c(0, 1.5, 3.5, 5.5, 21),
                                  right=FALSE, labels=c("1", "2-3", "4-5", "6+"))) |>
              group_by(size_ed) |>
              summarise(hh_per=sum(hh_per), .groups="drop") |>
              mutate(seed=i, cprob=as.factor(c_prob), lprob=as.factor(l_prob),
                     dprob=as.factor(d_prob))
            if (is.null(df)) {
              df <- cur_df
            } else {
              df <- rbind(df, cur_df)
            }
          }
        }
      }
    }
  }
  df
}

cal_mse <- function(df_hh, real_hh_df, lp, dp) {
  mean((filter(df_hh, lprob==lp, dprob==dp)$hh_per - real_hh_df$hh_per)^2)
}

h5ls("main/output/population.hd5")
h5ls("main/output/disease_100_120.hd5")
path = "main/output/params_growth_rate=0/seed_0/population.hd5"

cur_df <- h5read(path, "population/age_dists", compoundAsDataFrame=FALSE)
cur_df <- cur_df[["dist"]][,200]

growth_rates = c("0")
duration = 200
age_df <- load_age_data(no_runs=10, duration=duration) |>
  mutate(Source = "synthetic") |>
  subset(select=c("age_group", "pop", "Source", "growth_rate"))

real_age_df <- data.frame(pop=c(4.3, 5.1, 5.7, 6.1, 6.5, 7.2, 7.1, 7.1, 7.5, 7.4, 7.6, 7.3, 6.5, 5.2, 3.9, 2.5, 3.2)/100) |>
  mutate(age_group = row_number(), Source="data", growth_rate=growth_rates[1])
real_age_df$age_group <- age_df$age_group
real_age_df_dup <- real_age_df[rep(seq_len(nrow(real_age_df)), each = 1), ] |>
  mutate(growth_rate=growth_rates[2])

age_df_2 <- load_age_data(no_runs=10, duration=duration) |>
  mutate(Source = "synthetic") |>
  filter(growth_rate==growth_rates[2]) |>
  subset(select=c("age_group", "pop", "Source", "growth_rate"))
plot_age_df <- age_df |>
  rbind(real_age_df)
  #rbind(age_df_2, real_age_df_dup)

# print(mean((filter(age_df, growth_rate=="0.028")$pop - real_age_df$pop)^2))


ggplot(plot_age_df, aes(x=age_group, y=pop*100, group=Source, fill=Source)) +
  geom_col(position="identity", width=0.8, alpha=0.5) +
  # coord_flip() +
  scale_fill_manual(values=c("red", "grey35")) +
  theme_minimal() +
  xlab("Age group") +
  ylab("%") +
  #facet_grid(,vars(growth_rate)) +
  coord_flip()

couple_probs = c("0.03")
leaving_probs = c("0.01")
divorce_probs = c("0.01")

real_hh_df <- data.frame(hh_per=c(0.28,0.49,0.19,0.04)) |>
  mutate(size_ed = c("1", "2-3", "4-5", "6+"), Source="data")

df_hh <- load_hh_data(no_runs = 10, duration=190) |>
  group_by(cprob, lprob, dprob, size_ed) |>
  summarise(hh_per=mean(hh_per), Source="synthetic", .groups='drop')

cal_mse(df_hh, real_hh_df, "0.004", "0.0015")

plot_df <- df_hh |>
  subset(select=c("hh_per", "size_ed", "Source")) |>
  rbind(real_hh_df)

ggplot(plot_df, aes(x=size_ed, y=hh_per*100, group=Source, fill=Source)) +
  geom_col(position="identity", width=0.8, alpha=0.5) +
  # facet_grid(vars(lprob), vars(dprob)) +
  scale_fill_manual(values=c("red", "grey35")) +
  ylab("%") + xlab("Household size") +
  theme_minimal()
rfggplot(filter(df_hh, dprob=="0.0005"), aes(x=size, y=hh_per)) +
  geom_bar(stat="identity", position="dodge", width=0.8) +
  scale_x_continuous(limits=c(0.5,14), breaks=seq(1,14,1)) +
  facet_grid(vars(lprob), vars(cprob)) +
  geom_line(data=real_hh_df, aes(x=size, y=hh_per), color="red", size=0.8)
ggplot(filter(df_hh, cprob=="0.08"), aes(x=size, y=hh_per)) +
  geom_bar(stat="identity", position="dodge", width=0.8) +
  scale_x_continuous(limits=c(0.5,14), breaks=seq(1,14,1)) +
  facet_grid(vars(lprob), vars(dprob)) +
  geom_line(data=real_hh_df, aes(x=size, y=hh_per), color="red", size=0.8)

path = "main/output/params_mda_rounds=0_q_red=0/seed_0/disease_1_101.hd5"

group_df <- h5read(path, "prevalence/prev_group", compoundAsDataFrame=FALSE)
age_df <- group_df[["prev_age"]][,100]
