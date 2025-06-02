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

load_data <- function(duration=40, epi_burn_in=160, t_per_year=13, t_offset = 0) {
  df <- NULL
  burnin = epi_burn_in
  path = sprintf("main/output/disease_%s_%s.hd5",
                               burnin, burnin+duration)
  
  print(here(path))
  cur_df <- h5read(here(path), "community/base_data", compoundAsDataFrame = TRUE)
  cur_df <- cur_df |>
    mutate(total=(S+A+C+T+R+V), year=(t-t_offset)/t_per_year-burnin+1983,
           Acute=A, Chronic=C, Recovered=R, Treated=T, Vaccinated=V, Site=ifelse(community==0, "Migrant", "Thai"))
  cur_df <- cur_df %>% 
    pivot_longer(
      cols = `Acute`:`Vaccinated`, 
      names_to = "state",
      values_to = "count"
    ) 
  print(typeof(cur_df))
  cur_df
}

load_data_from_name <- function(file_name, burnin=160, t_per_year=13, t_offset = 0, obs="community") {
  df <- NULL
  path = sprintf("main/output/%s.hd5", file_name)
  
  print(here(path))
  cur_df <- h5read(here(path), sprintf("%s/base_data", obs), compoundAsDataFrame = TRUE)
  cur_df <- cur_df |>
    mutate(total=(S+A+C+T+R+V), year=(t-t_offset)/t_per_year-burnin+1983,
           Acute=A, Chronic=C, Recovered=R, Treated=T, Vaccinated=V,
           Site=case_when(community==0 ~ "Surrounding area",
                          community==1 ~ "Border village",
                          TRUE ~ "Thai city"),
           origin=ifelse(as.numeric(origin)==0,"Thai","Migrant"))
  cur_df <- cur_df %>% 
    pivot_longer(
      cols = `Acute`:`Vaccinated`, 
      names_to = "state",
      values_to = "count"
    ) 
  cur_df <- cur_df |>
    group_by(year, state, origin, Site) |>
    summarise(count = count, total = total)
  print(typeof(cur_df))
  cur_df
}

load_all_data <- function(m_access=c(0.8,1), r_access=c(0.1,0.5,0.9), m_mobility=c(1,0.5,0.1), no_runs=5, duration=60, burn_in=160, obs="community", chronic_only=TRUE) {
  data <- NULL
  for (m_a in m_access) {
    for (r_a in r_access) {
      for (m_m in m_mobility) {
        for (seed in 0:(no_runs-1)) {
          path = sprintf("params_new_migrant_access=%s_new_remote_access=%s_new_migrant_mobility=%s/seed_%s/disease_%s_%s",
                         m_a, r_a, m_m, seed, burn_in, burn_in+duration)
          if (no_runs == 1) {
            path = sprintf("params_new_migrant_access=%s_new_remote_access=%s_new_migrant_mobility=%s/disease_%s_%s",
                           m_a, r_a, m_m, burn_in, burn_in+duration)
          }
          df <- load_data_from_name(path, burnin = burn_in, obs=obs)
          if (chronic_only) {
            df <- filter(df, state=="Chronic")
          }
          df <- df |>
            mutate(m_a=m_a, r_a=r_a, m_m=m_m,seed=seed, year_sum = 1984+as.numeric(cut(year, breaks=c(0, seq(1984,2044,1)),
                                                                       right=FALSE, labels=seq(1984,2044,1)))) |>
            group_by(m_a, r_a, m_m, year_sum, seed, state, origin, Site) |>
            summarise(count=mean(count), total=mean(total))
          if (is.null(data)) {
            data <- df
          } else {
            data <- rbind(data, df)
          }
        }
      }
    }
  }
  data
}

load_cf_data <- function(v_access=c(0.2), no_runs=5, duration=60, burn_in=160, obs="community", chronic_only=TRUE) {
  data <- NULL
  for (v_a in v_access) {
    for (seed in 0:(no_runs-1)) {
      path = sprintf("params_new_village_access=%s/seed_%s/disease_%s_%s",
                     v_a, seed, burn_in, burn_in+duration)
      if (no_runs == 1) {
        path = sprintf("params_new_village_access=%s/disease_%s_%s",
                       v_a, burn_in, burn_in+duration)
      }
      df <- load_data_from_name(path, burnin = burn_in, obs=obs)
      if (chronic_only) {
        df <- filter(df, state=="Chronic")
      }
      df <- df |>
        mutate(v_a=v_a,seed=seed, year_sum = 1984+as.numeric(cut(year, breaks=c(0, seq(1984,2044,1)),
                                                                                   right=FALSE, labels=seq(1984,2044,1)))) |>
        group_by(v_a, year_sum, seed, state, origin, Site) |>
        summarise(count=mean(count), total=mean(total))
      if (is.null(data)) {
        data <- df
      } else {
        data <- rbind(data, df)
      }
    }
  }
  data
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
cali_df <- load_all_data(no_runs=50, m_access=c(0.6), r_access=c(0.2), m_mobility = c(1), chronic_only=FALSE)

real_prev <- data.frame(matrix(ncol = 6, nrow = 0))
x <- c("year", "prevalence", "q05", "q95", "Source", "author")
colnames(real_prev) <- x
real_prev[nrow(real_prev)+1,] = list(year=1983, prevalence=0.1418, q05=0.1092, q95=0.1798, Source="data", author="Brown et al")
real_prev[nrow(real_prev)+1,] = list(year=1989, prevalence=0.0857, q05=0.0538, q95=0.1280, Source="data", author="Phuapradit et al")
real_prev[nrow(real_prev)+1,] = list(year=1994, prevalence=0.0471, q05=0.0318, q95=0.0669, Source="data", author="Suwanagool et al")
real_prev[nrow(real_prev)+1,] = list(year=1990, prevalence=0.1163, q05=0.1007, q95=0.1333, Source="data", author="Kozik et al")
real_prev[nrow(real_prev)+1,] = list(year=1998, prevalence=0.1018, q05=0.0798, q95=0.1217, Source="data", author="Ishida et al")
real_prev[nrow(real_prev)+1,] = list(year=1999, prevalence=0.0815, q05=0.0598, q95=0.1078, Source="data", author="Pichainarong et al")
real_prev[nrow(real_prev)+1,] = list(year=2000, prevalence=0.0850, q05=0.0503, q95=0.1326, Source="data", author="Wiwanitkit")
real_prev[nrow(real_prev)+1,] = list(year=2003, prevalence=0.0860, q05=0.0650, q95=0.1070, Source="data", author="Jutavijittum et al")
real_prev[nrow(real_prev)+1,] = list(year=2003, prevalence=0.0643, q05=0.0577, q95=0.0714, Source="data", author="Thanachartwet et al")
real_prev[nrow(real_prev)+1,] = list(year=2004, prevalence=0.0396, q05=0.0349, q95=0.0448, Source="data", author="Chongsrisawat et al")
real_prev[nrow(real_prev)+1,] = list(year=2009, prevalence=0.0494, q05=0.0472, q95=0.0516, Source="data", author="Sangrajrang et al")
real_prev[nrow(real_prev)+1,] = list(year=2014, prevalence=0.0348, q05=0.0301, q95=0.0395, Source="data", author="Posuwan et al")
real_prev[nrow(real_prev)+1,] = list(year=2024, prevalence=0.0168, q05=0.0135, q95=0.0201, Source="data", author="Nilyanimit et al")

cali_plot_df <-  filter(cali_df, state=="Chronic", year_sum<=2024, origin=="Thai") |>
  group_by(year_sum, origin, seed) |>
  summarise(count=sum(count), total=sum(total)) |>
  mutate(prev=count/total) |>
  group_by(year_sum, origin) |>
  summarise(count=mean(count), total=mean(total), prev_mean=mean(prev), q05=quantile(prev, 0.025), q95=quantile(prev, 0.975))
ggplot(cali_plot_df) +
  geom_line(aes(x=year_sum, y=100*prev_mean), linewidth=0.8, colour="black") +
  geom_point(data=real_prev, aes(x=year, y=prevalence*100), size=1.5, colour="salmon1") +
  geom_errorbar(data=real_prev, aes(x=year, ymin=q05*100, ymax=q95*100), width=.2,
                position=position_dodge(.9), colour="salmon2") +
  geom_ribbon(aes(x=year_sum, ymin=100*q05, ymax=100*q95), alpha=0.3, fill="black", colour=NA) +
  scale_y_continuous("Prevalence (%)", expand = c(0, 0)) +
  coord_cartesian(ylim=c(0,18)) +
  xlab("Year") +
  geom_vline(xintercept=2024, linetype=3) +
  geom_text(aes(x=2024, label="present \n", y=12), colour="black", angle=90) +
  theme_minimal() +
  theme(legend.position="top")

cali_plot_df <-  filter(cali_df, state=="Chronic", year_sum<=2024) |>
  group_by(year_sum, Site, seed) |>
  summarise(count=sum(count), total=sum(total)) |>
  mutate(prev=count/total) |>
  group_by(year_sum, Site) |>
  summarise(count=mean(count), total=mean(total), prev_mean=mean(prev), q05=quantile(prev, 0.25), q95=quantile(prev, 0.75))
ggplot(cali_plot_df) +
  geom_line(aes(x=year_sum, y=100*prev_mean, colour=Site), linewidth=1) +
  geom_point(data=real_prev, aes(x=year, y=prevalence*100), size=2) +
  geom_ribbon(aes(x=year_sum, ymin=100*q05, ymax=100*q95, fill=Site), alpha=0.2) +
  scale_y_continuous("Prevalence (%)", expand = c(0, 0)) +
  coord_cartesian(ylim=c(0,15)) +
  xlab("Year") +
  geom_vline(xintercept=2024, linetype=3) +
  geom_text(aes(x=2024, label="present \n", y=12), colour="black", angle=90) +
  theme_minimal()

pop_cali_df <- load_cf_data(no_runs=1, v_access=c(0), chronic_only=FALSE, burn_in = 140)

pop_cali_plot_df <- filter(pop_cali_df, state=="Chronic") |>
  group_by(year_sum, Site, origin) |>
  summarise(total=mean(total), q05=quantile(total, 0.25), q95=quantile(total, 0.75))
ggplot(pop_cali_plot_df, aes(y=total*2, x=year_sum, fill=origin)) +
  geom_area() +
  facet_wrap(~ Site, ncol=3) +
  xlab("Year") +
  ylab("Population") +
  scale_x_discrete(expand = c(0, 0), limits = c(1985,2005, 2025, 2045)) +
  coord_cartesian(xlim=c(1985,2045)) +
  theme_minimal() +
  theme(legend.position="top",
        axis.text.x = element_text(size=8, angle=45, hjust=1, vjust=1),
        panel.spacing = unit(1.2, "lines"))
  



jobs <- "community"
test <- load_data_from_name("params_new_migrant_access=0.6_new_remote_access=0.2_new_migrant_mobility=1/disease_160_220", burnin = 160, obs=obs)



base_df <- load_all_data(obs=obs, no_runs=50, m_access=c(0.6), r_access=c(0.2), m_mobility = c(1)) |>
  mutate(Scenario="Base")

ima_df <- load_all_data(obs=obs, no_runs=50, m_access=c(1), r_access=c(0.2), m_mobility = c(1)) |>
  mutate(Scenario="Improved migrant access")

ira_df <- load_all_data(obs=obs, no_runs=50, m_access=c(0.6), r_access=c(0.6), m_mobility = c(1)) |>
  mutate(Scenario="Improved remote access")

imra_df <- load_all_data(obs=obs, no_runs=50, m_access=c(1), r_access=c(0.6), m_mobility = c(1)) |>
  mutate(Scenario="Improved both access")

imra_df <- load_all_data(obs=obs, no_runs=50, m_access=c(1), r_access=c(0.6), m_mobility = c(1)) |>
  mutate(Scenario="Improved both access")

cf_df <- load_cf_data(obs=obs, no_runs=50, v_access=c(0.2)) |>
  mutate(Scenario="Reduced village access")

all_df <- rbind(base_df, ima_df, ira_df, imra_df, cf_df)
#all_child_df <- all_df
all_com_df <- all_df

Thai_df <- all_child_df |>
  group_by(m_a, r_a, m_m, year_sum, state, Scenario, seed) |>
  summarise(total=sum(total), count=sum(count)) |>
  mutate(Site="All sites", origin="All origins")
all_child_df <- rbind(all_child_df, Thai_df)

site_df <- all_child_df |>
  mutate(Site=ifelse(Site=="Border village" | Site=="Surrounding area", "Border village & region", Site)) |>
  group_by(m_a, r_a, m_m, year_sum, state, Site, Scenario, seed) |>
  summarise(total=sum(total), count=sum(count)) |>
  mutate(prev=count/total) |>
  group_by(m_a, r_a, m_m, year_sum, state, Scenario, Site) |>
  summarise(total=mean(total), count=mean(count), prev_mean=mean(prev), q05=quantile(prev, 0.25), q95=quantile(prev, 0.75)) |>
  mutate(linetype=ifelse(Scenario=="Reduced village access", "dashed", "solid"))

ggplot(filter(site_df, Scenario=="Base" | Scenario=="Improved both access" | Scenario=="Reduced village access"), aes(x=as.numeric(year_sum), y=100*prev_mean, 
                    colour=Scenario, linetype=Scenario)) +
  geom_line(linewidth=0.8) +
  geom_ribbon(aes(ymin=100*q05, ymax=100*q95, fill=Scenario), alpha=0.25, colour=NA) +
  scale_y_continuous("Prevalence (%)", expand = c(0, 0)) +
  coord_cartesian(ylim=c(0,2)) +
  scale_x_continuous("Year", limits = c(2024, NA)) +
  scale_colour_manual(values=c("grey30", "palegreen3", "orangered1"))+
  scale_fill_manual(values=c("grey30", "palegreen3", "orangered1"))+
  scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dotted"))+
  geom_hline(yintercept=0.1, linetype=3) +
  geom_text(aes(x=2026, label="0.1%", y=0.3), colour="black", size=2.8) +
  theme_minimal() +
  theme(legend.position="top",
        axis.text.x = element_text(size=8, angle=45, hjust=1.2, vjust=1.5)) +
  guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
  facet_wrap(~ Site, ncol=4)

origin_df <- all_child_df |>
  group_by(m_a, r_a, m_m, year_sum, state, origin, Scenario, seed) |>
  summarise(total=sum(total), count=sum(count)) |>
  mutate(prev=count/total) |>
  group_by(m_a, r_a, m_m, year_sum, state, Scenario, origin) |>
  summarise(total=mean(total), count=mean(count), prev_mean=mean(prev), q05=quantile(prev, 0.25), q95=quantile(prev, 0.75))
ggplot(filter(origin_df, Scenario=="Base" | Scenario=="Improved both access" | Scenario=="Reduced village access"),
       aes(x=as.numeric(year_sum), y=100*prev_mean, colour=Scenario, linetype=Scenario)) +
  geom_line(linewidth=1) +
  geom_ribbon(aes(ymin=100*q05, ymax=100*q95, fill=Scenario), alpha=0.25, colour=NA) +
  scale_y_continuous("Prevalence (%)", expand = c(0, 0)) +
  coord_cartesian(ylim=c(0,2.5)) +
  scale_x_continuous("Year", limits = c(2024, NA)) +
  scale_colour_manual(values=c("grey30", "palegreen3", "orangered1"))+
  scale_fill_manual(values=c("grey30", "palegreen3", "orangered1"))+
  scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dotted"))+
  geom_hline(yintercept=0.1, linetype=3) +
  geom_text(aes(x=2026, label="0.1%", y=0.3), colour="black", size=2.8) +
  theme_minimal() +
  theme(legend.position="top",
        axis.text.x = element_text(size=8, angle=45, hjust=1, vjust=1)) +
  guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
  facet_wrap(~ origin)

Thai_df <- all_com_df |>
  group_by(m_a, r_a, m_m, year_sum, state, Scenario, seed) |>
  summarise(total=sum(total), count=sum(count)) |>
  mutate(Site="All sites", origin="All origins")
all_com_df <- rbind(all_com_df, Thai_df)
site_df <- all_com_df |>
  mutate(Site=ifelse(Site=="Border village" | Site=="Surrounding area", "Border village & region", Site)) |>
  group_by(m_a, r_a, m_m, year_sum, state, Site, Scenario, seed) |>
  summarise(total=sum(total), count=sum(count)) |>
  mutate(prev=count/total) |>
  group_by(m_a, r_a, m_m, year_sum, state, Scenario, Site) |>
  summarise(total=mean(total), count=mean(count), prev_mean=mean(prev), q05=quantile(prev, 0.05), q95=quantile(prev, 0.95)) |>
  mutate(linetype=ifelse(Scenario=="Reduced village access", "dashed", "solid"))

ggplot(filter(site_df, Scenario=="Base" | Scenario=="Improved both access" | Scenario=="Reduced village access"), aes(x=as.numeric(year_sum), y=100*prev_mean, 
                                                                                                                      colour=Scenario, linetype=Scenario)) +
  geom_line(linewidth=0.8) +
  geom_ribbon(aes(ymin=100*q05, ymax=100*q95, fill=Scenario), alpha=0.25, colour=NA) +
  scale_y_continuous("Prevalence (%)", expand = c(0, 0)) +
  coord_cartesian(ylim=c(0,5)) +
  scale_x_continuous("Year", limits = c(2024, NA)) +
  scale_colour_manual(values=c("grey30", "palegreen3", "orangered1"))+
  scale_fill_manual(values=c("grey30", "palegreen3", "orangered1"))+
  scale_linetype_manual(values=c("solid", "solid", "solid", "solid", "dotted"))+
  theme_minimal() +
  theme(legend.position="top",
        axis.text.x = element_text(size=8, angle=45, hjust=1.2, vjust=1.5)) +
  guides(colour=guide_legend(nrow=2,byrow=TRUE)) +
  facet_wrap(~ Site, ncol=4)


heatmap_df <- load_all_data(obs="children", no_runs=50, m_access=c(0.6, 0.8, 1), r_access=c(0.4, 0.6, 0.8), m_mobility = c(1))
heatmap_all_df <- filter(heatmap_df, year_sum==2045) |>
  group_by(m_a, r_a, m_m, year_sum, state, origin, seed) |>
  summarise(total=sum(total), count=sum(count)) |>
  mutate(Site="All")
heatmap_plot_df <- rbind(heatmap_df, heatmap_all_df)
heatmap_plot_df <- filter(heatmap_plot_df, year_sum==2045) |>
  group_by(m_a, r_a, m_m, year_sum, state, Site, seed) |>
  summarise(total=sum(total), count=sum(count)) |>
  mutate(prev=count/total) |>
  group_by(m_a, r_a, m_m, year_sum, state, Site) |>
  summarise(total=mean(total), count=mean(count), prev_mean=mean(prev), q05=quantile(prev, 0.25), q95=quantile(prev, 0.75))
ggplot(heatmap_plot_df) +
  geom_tile(aes(x=m_a, y=r_a, fill=prev_mean*100)) +
  geom_text(aes(label=round(prev_mean*100,2), x=m_a, y=r_a, fontface = ifelse(round(prev_mean*100,2)<=0.1, "bold", "plain")),size=2.5) +
  scale_fill_gradientn(colors=c("green", "white"), na.value="white") +
  labs(x="Migrant access", y="Remote access", fill="Prevalence\n0-5 years old (%)") +
  theme_minimal() +
  facet_wrap(~ Site, ncol=2)

### age prevalence
load_age_hh_prev_data_from_name <- function(path=NULL, mode="prev_age", duration=60, burn_in=160, year=40, t_per_year=13, t_offset=0) {
  df <- NULL
  if (is.null(path)) {
    path = sprintf("main/output/disease_%s_%s.hd5", burn_in, burn_in+duration)
  }
  if ( file.exists( path ) ) {
    cur_df <- h5read(path, "prevalence/prev_group", compoundAsDataFrame=FALSE)
    
    for (t in range(t_per_year*(year-1),t_per_year*(year))) {
      cur_df <- cur_df[[mode]][,t]
      cur_df <- data.frame(prevalence=cur_df) |>
        mutate(group = row_number(), t=t)
      if (is.null(df)) {
        df <- cur_df
      } else {
        df <- rbind(df, cur_df)
      }
    }
    
    df <- df |>
      group_by(group) |>
      summarise(prevalence = mean(prevalence))
  } else {
    sprintf("%s does not exist!", path)
  }
  df
}

load_all_age_prev_data <- function(m_access=c(0.8,1), r_access=c(0.1,0.5,0.9), m_mobility=c(1,0.5,0.1), no_runs=5, duration=60, burn_in=160, year=40) {
  data <- NULL
  for (m_a in m_access) {
    for (r_a in r_access) {
      for (m_m in m_mobility) {
        for (seed in 0:(no_runs-1)) {
          path = sprintf("main/output/params_new_migrant_access=%s_new_remote_access=%s_new_migrant_mobility=%s/seed_%s/disease_%s_%s.hd5",
                         m_a, r_a, m_m, seed, burn_in, burn_in+duration)
          if (no_runs == 1) {
            path = sprintf("main/output/params_new_migrant_access=%s_new_remote_access=%s_new_migrant_mobility=%s/disease_%s_%s.hd5",
                           m_a, r_a, m_m, burn_in, burn_in+duration)
          }
          df <- load_age_hh_prev_data_from_name(path, year=year) |>
            mutate(m_a=m_a, r_a=r_a, m_m=m_m,seed=seed)
          if (is.null(data)) {
            data <- df
          } else {
            data <- rbind(data, df)
          }
        }
      }
    }
  }
  data
}

plot_age_prevalence <- function(age_prev, title="") {
  x_labs <- c("0-9", "10-24", "25-34", "35-44", "45-54", "55+")
  hospital_labeller <- function(variable,value){
    return(x_labs[value])
  }
  plot <- ggplot(age_prev, aes(x=factor(group), y=prevalence*100, colour=Source)) +
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
                            
                            source="data", year = "2015 Leroi") |>
  mutate(group = row_number())

real_age_prev_posu <- data.frame(prevalence=c(0.002, 0.015, 0.0337, 0.0422, 0.0533, 0.0662),
                            
                            source="data", year = "2014 Posuwan") |>
  mutate(group = row_number())

real_age_prev_2024 <- data.frame(prevalence=c(0.0, 0.001, 0.0075, 0.0163, 0.02595, 0.0373),
                                 
                                 source="data", year = "2024 Nilyanimit") |>
  mutate(group = row_number())
age_prev <- load_all_age_prev_data(duration=60, burn_in=160, year=31, 
                                   m_access=c(0.6), r_access = c(0.2), m_mobility=c(1),
                                   no_runs=5) |>
  mutate(source="synthetic", year="2015 Leroi") |>
  group_by(source, year, group) |>
  summarise(prevalence=mean(prevalence))

age_prev_posu <- load_all_age_prev_data(duration=60, burn_in=160, year=30, 
                                   m_access=c(0.6), r_access = c(0.2), m_mobility=c(1),
                                   no_runs=5) |>
  mutate(source="synthetic", year="2014 Posuwan") |>
  group_by(source, year, group) |>
  summarise(prevalence=mean(prevalence))

age_prev_2024 <- load_all_age_prev_data(duration=60, burn_in=160, year=40, 
                                   m_access=c(0.6), r_access = c(0.2), m_mobility=c(1),
                                   no_runs=5) |>
  mutate(source="synthetic", year="2024 Nilyanimit") |>
  group_by(source, year, group) |>
  summarise(prevalence=mean(prevalence))

age_plot_df <- rbind(rbind(age_prev, age_prev_posu, age_prev_2024)[c("prevalence", "source", "group", "year")], 
                      real_age_prev, real_age_prev_posu, real_age_prev_2024)
ggplot(age_plot_df, aes(x=factor(group), y=prevalence*100, fill=source)) +
  geom_col(aes( group=source), position="dodge", width=0.6, alpha=0.6) +
  #geom_boxplot(aes(color=source), size=0.5) +
  scale_y_continuous("Prevalence (%)") +
  scale_x_discrete(labels=c("0-9", "10-24", "25-34", "35-44", "45-54", "55+")) +
  scale_fill_manual(values=c("red", "grey35")) +
  labs(x="Age group") +
  # geom_hline(yintercept=, linetype = 'dotted') +
  theme_minimal() +
  theme(legend.position="top",
        axis.text.x = element_text(size=8, angle=45, hjust=1, vjust=1.5)) +
  facet_wrap(~ year)

sdf <- load_age_hh_prev_data_from_name("main/output/params_new_migrant_access=0.6_new_remote_access=0.1_new_migrant_mobility=1/seed_0/disease_160_220.hd5", duration=60, burn_in =160)
  
base_age_df <- load_all_age_prev_data(no_runs=5, m_access=c(0.6), r_access=c(0.1), m_mobility = c(1)) |>
  mutate(Scenario="Base")

ima_age_df <- load_all_age_prev_data(no_runs=5, m_access=c(1), r_access=c(0.1), m_mobility = c(1)) |>
  mutate(Scenario="Improved migrant access")

ira_age_df <- load_all_age_prev_data(no_runs=5, m_access=c(0.6), r_access=c(0.5), m_mobility = c(1)) |>
  mutate(Scenario="Improved remote access")

rm_age_df <- load_all_age_prev_data(no_runs=5, m_access=c(0.6), r_access=c(0.1), m_mobility = c(0.5)) |>
  mutate(Scenario="Reduced mobility")

ggplot(rbind(base_age_df, ima_age_df, ira_age_df, rm_age_df), aes(x=Scenario, y=prevalence*100)) +
  geom_boxplot() +
  # geom_boxplot(colour="black", size=0.5) +
  scale_y_continuous("Prevalence (%)") +
  scale_fill_manual(values=c("red", "grey35")) +
  # geom_hline(yintercept=, linetype = 'dotted') +
  theme_minimal() +
  theme(legend.position="top")



plot_age_prevalence(rbind(age_prev,base_age_df))
ggsave("agegroup_prev.pdf", width=10, height=8, units="cm", path="./figures")

age_test <- h5read("main/output/params_new_migrant_access=0.6_new_remote_access=0.2_new_migrant_mobility=1/disease_160_220.hd5", "prevalence/prev_group", compoundAsDataFrame=FALSE)

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

h5read(path, "prevalence/prev_group", compoundAsDataFrame=FALSE)