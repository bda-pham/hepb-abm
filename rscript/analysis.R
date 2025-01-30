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

colors <- diverge_hsl(3)
colors_light <- c("indianred1", "khaki1", "palegreen")

load_data <- function(duration=30, epi_burn_in=150, t_per_year=13, t_offset = 0) {
  df <- NULL
  burnin = epi_burn_in
  path = sprintf("main/output/disease_%s_%s.hd5",
                               burnin, burnin+duration)
  
  print(here(path))
  cur_df <- h5read(here(path), "state/base_data", compoundAsDataFrame = TRUE)
  cur_df <- cur_df |>
    mutate(total=(S+A+C+T+R+V), year=(t-t_offset)/t_per_year-burnin+1983,
           Acute=A, Chronic=C, Recovered=R, Treated=T, Vaccinated=V)
  cur_df <- cur_df %>% 
    pivot_longer(
      cols = `Acute`:`Vaccinated`, 
      names_to = "state",
      values_to = "count"
    ) 
  print(typeof(cur_df))
  cur_df
}

load_data_from_name <- function(file_name, epi_burn_in=150, t_per_year=13, t_offset = 0) {
  df <- NULL
  burnin = epi_burn_in
  path = sprintf("main/output/%s.hd5", file_name)
  
  print(here(path))
  cur_df <- h5read(here(path), "state/base_data", compoundAsDataFrame = TRUE)
  cur_df <- cur_df |>
    mutate(total=(S+A+C+T+R+V), year=(t-t_offset)/t_per_year-burnin+1983,
           Acute=A, Chronic=C, Recovered=R, Treated=T, Vaccinated=V)
  cur_df <- cur_df %>% 
    pivot_longer(
      cols = `Acute`:`Vaccinated`, 
      names_to = "state",
      values_to = "count"
    ) 
  print(typeof(cur_df))
  cur_df
}


plot_prevalence <- function(prevalence_df, title="", duration=40) {
  plot <- ggplot(filter(prevalence_df, state!="S"), aes(x=year, y=100*count/total, color=state)) +
    geom_line(size=1) +
    scale_y_continuous("Prevalence (%)") +
    xlab("Year") +
    theme_minimal()
  plot
}
setwd(here())
real_prev <- data.frame(prevalence=c(0.0494, 0.0643, 0.0815, 0.0857, 0.1471),
                        year=c(2012, 2007, 2003, 1989, 1985), Source="data")
df <- load_data(duration=40, epi_burn_in = 160, t_offset = 0)
ggplot(filter(df, state!="S")) +
  geom_line(aes(x=year, y=100*count/total, color=state), size=1) +
  geom_point(data=real_prev, aes(x=year, y=prevalence*100), size=2) +
  scale_y_continuous("Prevalence (%)") +
  xlab("Year") +
  theme_minimal()

plot_prevalence(df, duration=40)

df_cover <- load_data_from_name("disease_cover", epi_burn_in = 160) |>
  mutate(scenario="cover")
df_base <- load_data_from_name("disease_base", epi_burn_in = 160) |>
  mutate(scenario="base")
df_access <- load_data_from_name("disease_access", epi_burn_in = 160) |>
  mutate(scenario="access")
df_noaccess <- load_data_from_name("disease_noaccess", epi_burn_in = 160) |>
  mutate(scenario="noaccess")

all_df <- rbind(df_base, df_cover, df_access, df_noaccess) |>
  filter(state=="Chronic")
ggplot(all_df) +
  geom_line(aes(x=year, y=100*count/total, color=scenario), size=1) +
  geom_point(data=real_prev, aes(x=year, y=prevalence*100), size=2) +
  scale_y_continuous("Prevalence (%)") +
  xlab("Year") +
  theme_minimal()

h### age prevalence
load_age_hh_prev_data <- function(mode="prev_age", duration=100, epi_burn_in=200, t_per_year=13, t_offset=0) {
  df <- NULL
  path = sprintf("main/output/disease_%s_%s.hd5", epi_burn_in, epi_burn_in+duration)
  if ( file.exists( path ) ) {
    cur_df <- h5read(path, "prevalence/prev_group", compoundAsDataFrame=FALSE)
    cur_df <- cur_df[[mode]][,t_per_year*(duration)+1]
    
    cur_df <- data.frame(prevalence=cur_df) |>
      mutate(group = row_number())
    if (is.null(df)) {
      df <- cur_df
    } else {
      df <- rbind(df, cur_df)
    }
  } else {
    printf("%s does not exist!", path)
  }
  df
}

plot_age_prevalence <- function(age_prev, title="") {
  x_labs <- c("0-9", "10-24", "25-34", "35-44", "45-54", "55+")
  hospital_labeller <- function(variable,value){
    return(x_labs[value])
  }
  plot <- ggplot(age_prev, aes(x=factor(group), y=mean*100, colour=Source)) +
    geom_point(colour="royalblue", size=2) +
    geom_boxplot(colour="black", size=0.5) +
    scale_y_continuous("Prevalence (%)") +
    scale_x_discrete(labels=x_labs) +
    scale_colour_manual(values=c("red", "grey35"))
    labs(x="Age group") +
    # geom_hline(yintercept=, linetype = 'dotted') +
    theme_minimal()
  plot
}


real_age_prev <- data.frame(prevalence=c(0.003, 0.011, 0.0326, 0.0593, 0.0714, 0.0662),
                            
                            Source="data") |>
  mutate(group = row_number())
age_prev <- load_age_hh_prev_data(duration=40, epi_burn_in=160) |>
  mutate(Source="synthetic")
ggplot(rbind(age_prev, real_age_prev), aes(x=factor(group), y=prevalence*100, fill=Source, group=Source)) +
  geom_col(position="identity", width=0.6, alpha=0.5) +
  # geom_boxplot(colour="black", size=0.5) +
  scale_y_continuous("Prevalence (%)") +
  scale_x_discrete(labels=c("0-9", "10-24", "25-34", "35-44", "45-54", "55+")) +
  scale_fill_manual(values=c("red", "grey35")) +
  labs(x="Age group") +
  # geom_hline(yintercept=, linetype = 'dotted') +
  theme_minimal()


plot_age_prevalence(rbind(age_prev, real_age_prev))
ggsave("agegroup_prev.pdf", width=10, height=8, units="cm", path="./figures")

  age_test <- h5read("main/output/params_mda_rounds=0_q_red=0/seed_0/disease_1_101.hd5", "prevalence/prev_group", compoundAsDataFrame=FALSE)

hh_prev <- load_age_hh_prev_data(mode="prev_household", no_runs=2) |>
  group_by(group) |>
  summarise(mean=mean(prevalence), median=median(prevalence), q05=quantile(prevalence, 0.05), q95=quantile(prevalence, 0.95), .groups='drop')

plot_hh_prevalence <- function(hh_prev, title="") {
  x_labs <- c("1", "2", "3", "4", "5", "6", "7", "8", "9+")
  hospital_labeller <- function(variable,value){
    return(x_labs[value])
  }
  plot <- ggplot(hh_prev, aes(x=factor(group), y=mean*100), colour="royalblue") +
    geom_point(colour="royalblue", size=2) +
    geom_errorbar(aes(ymin=100*q05, ymax=100*q95), colour="royalblue", width=0.5) +
    scale_y_continuous("Prevalence (%)", breaks=seq(0,50,10), limits = c(0, 50)) +
    labs(x="Household size") +
    scale_x_discrete(labels=x_labs) +
    # geom_hline(yintercept=, linetype = 'dotted') +
    theme_minimal()
  plot
}

plot_hh_prevalence(hh_prev)
ggsave("household_prev.pdf", width=10.5, height=8, units="cm", path="./figures")

