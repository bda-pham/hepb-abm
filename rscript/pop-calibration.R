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




load_age_data <- function(generation_threshold="1", max_migrant_age=101, separate_growth="False", no_runs=10, duration=200, origin=NULL) {
  df <- NULL
  max_seed <- no_runs - 1
  for (i in 0:max_seed) {
    path = sprintf("main/output/params_generation_threshold=%s_max_migrant_age=%s_separate_growth=%s/seed_%s/population.hd5",
                   generation_threshold, max_migrant_age, separate_growth, i)
    if ( file.exists( path ) ) {
      label = "population/age_dists"
      if (!is.null(origin)) {
        label = paste("population.",origin,"/age_dists",sep="")
      }
      cur_df <- h5read(path, label, compoundAsDataFrame=FALSE)
      print(i)
      cur_df <- cur_df[["dist"]][,duration]
      pop <- sum(cur_df)
      origin_text = origin
      if (is.null(origin)) {
        origin_text = "All"
      }
      cur_df <- data.frame(pop=cur_df/pop) |>
        mutate(age = row_number()-1, seed=i) |>
        mutate(age_group= cut(age, breaks=c(seq(0,80,5),105),
                              right=FALSE)) |>
        group_by(age_group, seed) |>
        summarise(pop=sum(pop), sg=separate_growth, gt=generation_threshold, mma=max_migrant_age, origin=origin_text, .groups="drop")
      
      if (is.null(df)) {
        df <- cur_df
      } else {
        df <- rbind(df, cur_df)
      }
    }
  }
  df
}

load_hh_data <- function(couple_probs, leaving_probs, divorce_probs, no_runs=10, duration=100) {
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
              mutate(size_ed= cut(size, breaks=c(0, 2.5, 4.5, 7.5, 21),
                                  right=FALSE, labels=c("1-2", "3-4", "5-7", "8+"))) |>
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

cal_mse <- function(df_hh, real_hh_df) {
  mean((df_hh$hh_per - real_hh_df$hh_per)^2)
}

h5ls("main/output/population.hd5")
h5ls("main/output/disease_100_120.hd5")
path = "main/output/params_growth_rate=0/seed_0/population.hd5"

cur_df <- h5read(path, "population/age_dists", compoundAsDataFrame=FALSE)
cur_df <- cur_df[["dist"]][,200]

growth_rates = c("0")
duration = 200


real_age_df <- data.frame(pop=c(4.3, 5.1, 5.7, 6.1, 6.5, 7.2, 7.1, 7.1, 7.5, 7.4, 7.6, 7.3, 6.5, 5.2, 3.9, 2.5, 3.2)/100) |>
  mutate(age_group = row_number(), Source="data", Scenario="")
real_age_df$age_group <- age_df$age_group

age_df5 <- load_age_data(generation_threshold="1", max_migrant_age="101", separate_growth="False",
                        origin="THAI", no_runs=5, duration=180) |>
  mutate(Source = "Base", Scenario="Base (Thai population)") |>
  subset(select=c("age_group", "pop", "Source", "Scenario", "origin"))

age_df <- load_age_data(generation_threshold="1", max_migrant_age="101", separate_growth="False",
                        origin="MIGRANT", no_runs=5, duration=180) |>
  mutate(Source = "1", Scenario="Base") |>
  subset(select=c("age_group", "pop", "Source", "Scenario", "origin"))

# age_df2 <- load_age_data(generation_threshold="2", max_migrant_age="101", separate_growth="False",
#                         origin="MIGRANT", no_runs=5, duration=180) |>
#   mutate(Source = "Thailand's growth,\n generation threshold 2", Scenario="Thailand's growth,\n generation threshold 2") |>
#   subset(select=c("age_group", "pop", "Source", "Scenario", "origin"))


age_df3 <- load_age_data(generation_threshold="1", max_migrant_age="101", separate_growth="True",
                         origin="MIGRANT", no_runs=5, duration=180) |>
  mutate(Source = "Separate growth,\n generation threshold 1", Scenario="Separate growth,\n generation threshold 1") |>
  subset(select=c("age_group", "pop", "Source", "Scenario", "origin"))

age_df4 <- load_age_data(generation_threshold="2", max_migrant_age="101", separate_growth="True",
                         origin="MIGRANT", no_runs=5, duration=180) |>
  mutate(Source = "Separate growth,\n generation threshold 2", Scenario="Separate growth,\n generation threshold 2") |>
  subset(select=c("age_group", "pop", "Source", "Scenario", "origin"))


# plot_age_df <- age_df |>
#   rbind(real_age_df)
  #rbind(age_df_2, real_age_df_dup)

# print(mean((filter(age_df, growth_rate=="0.028")$pop - real_age_df$pop)^2))


# ggplot(plot_age_df, aes(x=age_group, y=pop*100, group=Source, fill=Source)) +
#   geom_col(position="identity", width=0.8, alpha=0.5) +
#   # coord_flip() +
#   scale_fill_manual(values=c("red", "grey35")) +
#   theme_minimal() +
#   xlab("Age group") +
#   ylab("%") +
#   #facet_grid(,vars(growth_rate)) +
#   coord_flip()

# age_df2base <- age_df
# age_df2base$Source = age_df2$Source
age_df3base <- age_df
age_df3base$Source = age_df3$Source
age_df4base <- age_df
age_df4base$Source = age_df4$Source
age_df5base <- age_df
age_df5base$Source = age_df5$Source

plot_origin_df <- rbind(
  # age_df2, age_df2base,
                        age_df3, age_df3base,
                        age_df4, age_df4base,
                        age_df5, age_df5base)
plot_origin_df <- plot_origin_df |>
  group_by(age_group, Scenario, origin, Source) |>
  summarise(pop=mean(pop), q05=quantile(pop, 0.025), q95=quantile(pop, 0.975), .groups="drop")

ggplot(filter(plot_origin_df,Scenario!="Base"), aes(x=age_group, y=pop*100, group=Scenario, fill=Scenario)) +
  #geom_point(data=filter(plot_origin_df,Scenario=="Base"), shape="|", size=3) +
  geom_col(data=filter(plot_origin_df, Scenario=="Base"), position="identity", width=0.7, alpha=0.5, colour="black") +
  geom_col(position="identity", width=0.8, alpha=0.5) +
  # coord_flip() +
  scale_fill_manual(values=c("white", "brown", "red", "blue", "purple")) +
  theme_minimal() +
  xlab("Age group") +
  ylab("%") +
  theme(legend.position="bottom", axis.text.y = element_text(size=7)) +
  guides(fill=guide_legend(ncol=2,bycol=TRUE)) +
  coord_flip() +
  facet_grid(,vars(Source))
  

couple_probs = 
leaving_probs = 
divorce_probs = 

real_hh_df <- data.frame(hh_per=c(0.609, 0.288, 0.098, 0.005)) |>
  mutate(size_ed = c("1-2", "3-4", "5-7", "8+"), Source="data")

df_hh <- load_hh_data(couple_probs=c("0.09"), 
                      leaving_probs=c("0.0025"),
                      divorce_probs=c("0.005"),
                      no_runs = 10, duration=200) |>
  group_by(cprob, lprob, dprob, size_ed) |>
  summarise(hh_per=mean(hh_per), Source="synthetic", .groups='drop')

cal_mse(df_hh, real_hh_df)

plot_df <- df_hh |>
  subset(select=c("hh_per", "size_ed", "Source")) |>
  rbind(real_hh_df)

ggplot(plot_df, aes(x=size_ed, y=hh_per*100, group=Source, fill=Source)) +
  geom_col(position="identity", width=0.8, alpha=0.5) +
  # facet_grid(vars(lprob), vars(dprob)) +
  scale_fill_manual(values=c("red", "grey35")) +
  ylab("%") + xlab("Household size") +
  theme_minimal()
ggplot(filter(df_hh, dprob=="0.0005"), aes(x=size, y=hh_per)) +
  geom_bar(stat="identity", position="dodge", width=0.8) +
  scale_x_continuous(limits=c(0.5,14), breaks=seq(1,14,1)) +
  facet_grid(vars(lprob), vars(cprob)) +
  geom_line(data=real_hh_df, aes(x=size, y=hh_per), color="red", size=0.8)
ggplot(filter(df_hh, cprob=="0.08"), aes(x=size, y=hh_per)) +
  geom_bar(stat="identity", position="dodge", width=0.8) +
  scale_x_continuous(limits=c(0.5,14), breaks=seq(1,14,1)) +
  facet_grid(vars(lprob), vars(dprob)) +
  geom_line(data=real_hh_df, aes(x=size, y=hh_per), color="red", size=0.8)


0.01: 0.005705558
0.075: 0.005705558
0.005: 0.004261089
0.0025: 0.003745242
0.001: 0.003633494
0.0005: 0.003672972