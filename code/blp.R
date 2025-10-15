# evaluate the performance of using best linear prediction to identify the effect modifer
# with causal forest and DR learner
library(tidyverse)
library(here)

## compare blp output for c15_n1000_sims1000
# read in the blp output files
drl <- read_csv(here("output/sim1000", "Blp_drl_c15_n1000_sims1000_2024-01-06.csv"))
cf <- read_csv(here("output/sim1000", "Blp_cf_c15_n1000_sims1000_2024-01-06.csv"))
# rank using the statistics by group of effect parametrizations
# and rank using the magnitude of the estimate by group of effect parametrizations
drl_rank <- drl %>%
  group_by(method, confounder_number, sample_size, group, seed_number) %>%
  mutate(rank_sta = order(abs(statistic), decreasing = TRUE),
         rank_est = order(abs(estimate), decreasing = TRUE)) %>%
  filter(term=="m") %>%
  group_by(method, confounder_number, sample_size, group) %>%
  summarise(avg_rank_sta = mean(rank_sta), avg_rank_est = mean(rank_est), # get the average rank across simulations
            pct_rank_sta = mean(rank_sta == 1), pct_rank_est = mean(rank_est == 1),
            pct_p_value = mean(p.value < 0.05)) 

cf_rank <- cf %>%
  group_by(method, confounder_number, sample_size, group, seed_number) %>%
  mutate(rank_sta = order(abs(statistic), decreasing = TRUE),
         rank_est = order(abs(estimate), decreasing = TRUE)) %>%
  filter(term=="m") %>%
  group_by(method, confounder_number, sample_size, group) %>%
  summarise(avg_rank_sta = mean(rank_sta), avg_rank_est = mean(rank_est), 
            pct_rank_sta = mean(rank_sta == 1), pct_rank_est = mean(rank_est == 1),
            pct_p_value = mean(p.value < 0.05)) 

compare_rank <- rbind(drl_rank, cf_rank) %>%
  arrange( confounder_number, sample_size, group)

## compare for all scenarios
setwd("H:/RA/HTE/output/sim100_cf")
blp_all <-list.files(path = getwd(), pattern = "Blp") %>% 
  map_dfr(read_csv)

blp_rank <- blp_all %>%
  group_by(method, confounder_number, sample_size, group, seed_number) %>%
  mutate(rank_sta = order(abs(statistic), decreasing = TRUE),
         rank_est = order(abs(estimate), decreasing = TRUE)) %>%
  filter(term=="m") %>%
  group_by(method, confounder_number, sample_size, group) %>%
  summarise(avg_rank_sta = mean(rank_sta), avg_rank_est = mean(rank_est), # get the average rank across simulations
            pct_rank_sta = mean(rank_sta == 1), pct_rank_est = mean(rank_est == 1),
            pct_p_value = mean(p.value < 0.05)) %>%
  arrange(confounder_number, sample_size, group)

write_csv(blp_rank, here("output/sim100_cf", paste0("BLP_compare", "_", Sys.Date(), ".csv")))

#----------------------------------------------------------------------------------
## distribution of rank for nsims = 1000
setwd("H:/RA/HTE/output/sim1000_n10000")
blp_big <-list.files(path = getwd(), pattern = "Blp") %>% 
  map_dfr(read_csv)
setwd("H:/RA/HTE/output/sim1000_n1000")
blp_small <-list.files(path = getwd(), pattern = "Blp") %>% 
  map_dfr(read_csv)
blp_all <- rbind(blp_big, blp_small)

# read in the truth file
truth_data <- read_csv(
  here::here("data", "truth_data_simulation_all_update.csv")
)

# get the truth data for each group
truth <- truth_data %>%
  filter(c_dim == 15) %>%
  mutate(group = row_number(), RD_m0 = round(risk_diff_m0, 3), RD_m1 = round(risk_diff_m1, 3)) %>%
  select(group, RD_m0, RD_m1)

# merge the truth data with the result
blp_truth <- blp_all %>%
  left_join(truth, by = "group") %>%
  select(-group)

blp_rank <- blp_truth %>%
  group_by(method, confounder_number, sample_size, RD_m0, RD_m1, seed_number) %>%
  mutate(rank_sta = order(abs(statistic), decreasing = TRUE),
         rank_est = order(abs(estimate), decreasing = TRUE)) %>%
  filter(term=="m") 

  # group_by(method, confounder_number, sample_size, group) %>%
  # group_walk(~ hist(.x$rank_sta))

  # by group; by sample size and confounder number; by all three
library(viridis)
#install.packages("hrbrthemes")
library(hrbrthemes)
library(gridExtra)

## method 1: overlapped histograms 
plot <- blp_rank %>%
  filter(confounder_number == 15 & sample_size == 10000) %>%
  ggplot( aes(x=rank_sta, color=method, fill=method)) +
  geom_histogram( position = "identity",alpha=0.6) +
  facet_grid(RD_m0~RD_m1) +
  #facet_wrap(~sample_size)
  theme_ipsum() +
  theme(
    legend.position="right",
    panel.spacing = unit(0.1, "lines")
  ) +
  labs(title = "Distribution of BLP rank by subgroup risk difference",
       subtitle = "(Sample size = 10,000, Confounder number = 15)",
       y = "Count",
       x = "Blp Rank")
  
ggsave("blp_n10k_c15.png", plot, path = here::here("output"), width = 12, height = 10)

plot2 <- blp_rank %>%
  filter(confounder_number == 15 & sample_size == 500) %>%
  ggplot( aes(x=rank_sta, color=method, fill=method)) +
  geom_histogram( position = "identity",alpha=0.6) +
  facet_grid(RD_m0~RD_m1) +
  #facet_wrap(~sample_size)
  theme_ipsum() +
  theme(
    legend.position="right",
    panel.spacing = unit(0.1, "lines")
  ) +
  labs(title = "Distribution of BLP rank by subgroup risk difference",
       subtitle = "(Sample size = 500, Confounder number = 15)",
       y = "Count",
       x = "Blp Rank")  

ggsave("blp_n500_c15.png", plot2, path = here::here("output"), width = 12, height = 10)

plot3 <- blp_rank %>%
  filter(confounder_number == 5 & sample_size == 10000) %>%
  ggplot( aes(x=rank_sta, color=method, fill=method)) +
  geom_histogram( position = "identity",alpha=0.6) +
  facet_grid(RD_m0~RD_m1) +
  #facet_wrap(~sample_size)
  theme_ipsum() +
  theme(
    legend.position="right",
    panel.spacing = unit(0.1, "lines")
  ) +
  labs(title = "Distribution of BLP rank by subgroup risk difference",
       subtitle = "(Sample size = 10,000, Confounder number = 5)",
       y = "Count",
       x = "Blp Rank")  
ggsave("blp_n10K_c5.png", plot3, path = here::here("output"), width = 12, height = 10)

## method 2: ditribution side by side (https://ggplot2.tidyverse.org/reference/position_dodge.html)
## use histograms
# plot <- blp_rank %>%
#   filter(confounder_number == 15 & sample_size == 10000) %>%
#   ggplot( aes(x=rank_sta, color=method, fill=method)) +
#   geom_histogram(bins = 31, position = "dodge") +
#   facet_grid(RD_m0~RD_m1) +
#   #facet_wrap(~sample_size)
#   #theme_ipsum() +
#   theme(
#     legend.position="right",
#     panel.spacing = unit(0.1, "lines")
#   ) +
#   #scale_color_grey()+scale_fill_grey() + # change color: https://learn.saylor.org/mod/book/view.php?id=58485&chapterid=45033
#   scale_color_manual(values = c("darkgrey", "black")) + scale_fill_manual(values = c("darkgrey", "black"))+
#   theme_classic() +
#   scale_y_continuous(sec.axis = sec_axis(~ . , name = "RD_m0", breaks = NULL, labels = NULL)) +  # add a secondary axis: https://ggplot2.tidyverse.org/reference/sec_axis.html
#   scale_x_continuous(sec.axis = sec_axis(~ . , name = "RD_m1", breaks = NULL, labels = NULL)) +
#   labs(title = "Distribution of BLP rank by subgroup risk difference",
#        subtitle = "(Sample size = 10,000, Confounder number = 15)",
#        y = "Count",
#        x = "Blp Rank")
# ggsave("blp_n10k_c15_histogram.png", plot, path = here::here("output"), width = 14, height = 10)
# update Oct 8, 2025: revision on format
## use barplots: specify margins = true inside facet_grid() if "all" is needed, but percentage for "all" is incorrect (should divided by 3)
plot <- blp_rank %>%
  filter(confounder_number == 15 & sample_size == 10000) %>%
  ggplot( aes(x=rank_sta, y = 18* after_stat(count)/sum(after_stat(count)), color=method, fill=method)) +
  geom_bar(width = 0.55, position = "dodge") +  # https://www.learnbyexample.org/r-bar-plot-ggplot2/
  facet_grid(
    RD_m0~RD_m1
    ) +
  #facet_wrap(~sample_size)
  #theme_ipsum() +
  # theme(
  #   legend.position="right",
  #   panel.spacing = unit(0.1, "lines"),
  #   legend.text = element_text(size = 30)
  # ) +
  #scale_color_grey()+scale_fill_grey() + # change color: https://learn.saylor.org/mod/book/view.php?id=58485&chapterid=45033
  scale_color_manual(
    values = c("darkgrey", "black"),
    labels = c(
      "BLP_causal_forest" = "Causal Forest",
      "BLP_dr_learner" = "DR Learner"
    )
  ) + 
  scale_fill_manual(
    values = c("darkgrey", "black"),
    labels = c(
      "BLP_causal_forest" = "Causal Forest",
      "BLP_dr_learner" = "DR Learner"
    )
  )+
  theme_classic(base_size = 28) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "RD (M = 0)", breaks = NULL), labels = scales::percent) +  # add a secondary axis: https://ggplot2.tidyverse.org/reference/sec_axis.html
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "RD (M = 1)", breaks = NULL, labels = NULL)) +
  # labs(title = "Distribution of BLP rank by subgroup risk difference",
  #      subtitle = "(Sample size = 10,000, Confounder number = 15)",
  #      y = "Percentage",
  #      x = "Blp Rank")
  labs(y = "Percentage",
       x = "BLP Rank")
#grid.arrange(plot,top='RD_m1', right='RD_m0') # another way to add labels to facets

ggsave("blp_n10k_c15.png", plot, path = here::here("output"), width = 14, height = 10)


plot2 <- blp_rank %>%
  filter(confounder_number == 15 & sample_size == 500) %>%
  ggplot( aes(x=rank_sta, y = 18*after_stat(count)/sum(after_stat(count)), color=method, fill=method)) +
  geom_bar(width = 0.55, position = "dodge") +  # https://www.learnbyexample.org/r-bar-plot-ggplot2/
  facet_grid(
    RD_m0~RD_m1
  ) +
  #facet_wrap(~sample_size)
  #theme_ipsum() +
  # theme(
  #   legend.position="right",
  #   panel.spacing = unit(0.1, "lines")
  # ) +
  #scale_color_grey()+scale_fill_grey() + # change color: https://learn.saylor.org/mod/book/view.php?id=58485&chapterid=45033
  scale_color_manual(
    values = c("darkgrey", "black"),
    labels = c(
      "BLP_causal_forest" = "Causal Forest",
      "BLP_dr_learner" = "DR Learner"
    )
  ) + scale_fill_manual(
    values = c("darkgrey", "black"),
    labels = c(
      "BLP_causal_forest" = "Causal Forest",
      "BLP_dr_learner" = "DR Learner"
    ))+
  theme_classic(base_size = 28) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "RD (M = 0)", breaks = NULL),labels = scales::percent) +  # add a secondary axis: https://ggplot2.tidyverse.org/reference/sec_axis.html
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "RD (M = 1)", breaks = NULL, labels = NULL)) +
  # labs(title = "Distribution of BLP rank by subgroup risk difference",
  #      subtitle = "(Sample size = 500, Confounder number = 15)",
  #      y = "Percentage",
  #      x = "Blp Rank")
  labs(y = "Percentage",
       x = "BLP Rank")

ggsave("blp_n500_c15.png", plot2, path = here::here("output"), width = 14, height = 10)


plot3 <- blp_rank %>%
  filter(confounder_number == 5 & sample_size == 10000) %>%
  ggplot( aes(x=rank_sta, y = 18*after_stat(count)/sum(after_stat(count)), color=method, fill=method)) +
  geom_bar(width = 0.55, position = "dodge") +  # https://www.learnbyexample.org/r-bar-plot-ggplot2/
  facet_grid(
    RD_m0~RD_m1
  ) +
  #facet_wrap(~sample_size)
  #theme_ipsum() +
  # theme(
  #   legend.position="right",
  #   panel.spacing = unit(0.1, "lines")
  # ) +
  #scale_color_grey()+scale_fill_grey() + # change color: https://learn.saylor.org/mod/book/view.php?id=58485&chapterid=45033
  scale_color_manual(
    values = c("darkgrey", "black"),
    labels = c(
      "BLP_causal_forest" = "Causal Forest",
      "BLP_dr_learner" = "DR Learner"
    )) + 
  scale_fill_manual(
      values = c("darkgrey", "black"),
      labels = c(
        "BLP_causal_forest" = "Causal Forest",
        "BLP_dr_learner" = "DR Learner"
    ))+
  theme_classic(base_size = 28) +
  scale_y_continuous(sec.axis = sec_axis(~ . , name = "RD (M = 0)", breaks = NULL), labels = scales::percent) +  # add a secondary axis: https://ggplot2.tidyverse.org/reference/sec_axis.html
  scale_x_continuous(sec.axis = sec_axis(~ . , name = "RD (M = 1)", breaks = NULL, labels = NULL)) +
  # labs(title = "Distribution of BLP rank by subgroup risk difference",
  #      subtitle = "(Sample size = 10,000, Confounder number = 5)",
  #      y = "Percentage",
  #      x = "Blp Rank")
  labs(y = "Percentage",
       x = "BLP Rank")

ggsave("blp_n10K_c5.png", plot3, path = here::here("output"), width = 14, height = 10)

# update: change the count to percentage within each group of method 
#         this change cannot be done to the plot with "all", as it will require multiplying the percentage by 2*3 = 6 for the "all" column to get the within-group percentage

# update Oct 8, 2025: create a panel figure
library(patchwork)

# Assuming your two plots are plot1 and plot2
combined_plot <- plot2 + plot3 + 
  plot_layout(ncol = 2, guides = "collect") +  # collect legends into one
  plot_annotation(tag_levels = 'a') &
  theme(legend.position = "bottom")  # place shared legend at bottom

# Display
combined_plot

ggsave("blp_panel.png", combined_plot, path = here::here("output"), width = 20, height = 10)
