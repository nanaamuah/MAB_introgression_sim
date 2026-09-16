library(dplyr)
library(ggplot2)
library(here)


# 1. Load processed results

replicate_runs <- readRDS(here("data", "processed", "replicate_runs.rds"))

convergence_99 <- readRDS(here("data", "processed", "convergence_99.rds"))

cost_to_99 <- readRDS(here("data", "processed", "cost_to_99.rds"))

pyramid_summary <- readRDS(here("data", "processed", "pyramid_summary.rds"))


# 2. Shared plot settings
# Schemes are shown by both colour and line type because Background and Dense genomic overlap

scheme_colours <- c(
  "No marker" = "#7f7f7f",
  "Foreground" = "#d95f02",
  "Background" = "#1b9e77",
  "Dense genomic" = "#7570b3"
)

save_plot <- function(plot, file) {
  ggsave(here("outputs", "figures", file), plot, width = 7, height = 5, dpi = 300)
}


# 3. Recurrent-parent recovery trajectory
# Mean IBD recovery across 50 replicates; error bars are one standard error

recovery_plot_data <- replicate_runs %>%
  group_by(scheme, generation) %>%
  summarise(mean_recovery = mean(recovery_ibd), se_recovery = sd(recovery_ibd) / sqrt(n()), .groups = "drop")

p_recovery <- recovery_plot_data %>%
  ggplot(aes(x = generation, y = mean_recovery, colour = scheme, linetype = scheme)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_recovery - se_recovery, ymax = mean_recovery + se_recovery), width = 0.1) +
  scale_y_continuous(limits = c(0.7, 1), labels = scales::percent) +
  scale_colour_manual(values = scheme_colours) +
  labs(
    x = "Backcross generation",
    y = "Recurrent-parent genome recovery",
    colour = "Scheme",
    linetype = "Scheme"
  ) +
  theme_minimal(base_size = 12)

p_recovery

save_plot(p_recovery, "recovery_trajectory.png")


# 4. Linkage-drag trajectory
# No marker is left out because most of its plants no longer carry the target allele

drag_plot_data <- replicate_runs %>%
  filter(scheme != "No marker") %>%
  group_by(scheme, generation) %>%
  summarise(mean_drag = mean(linkage_drag_cM), se_drag = sd(linkage_drag_cM) / sqrt(n()), .groups = "drop")

p_drag <- drag_plot_data %>%
  ggplot(aes(x = generation, y = mean_drag, colour = scheme, linetype = scheme)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_drag - se_drag, ymax = mean_drag + se_drag), width = 0.1) +
  scale_colour_manual(values = scheme_colours) +
  labs(x = "Backcross generation", y = "Linkage drag (cM)", colour = "Scheme", linetype = "Scheme") +
  theme_minimal(base_size = 12)

p_drag

save_plot(p_drag, "linkage_drag.png")


# 5. Population size and recovery speed

popsize_plot_data <- cost_to_99 %>%
  group_by(scheme, pop_size) %>%
  summarise(mean_gen_99 = mean(generation), se_gen_99 = sd(generation) / sqrt(n()), .groups = "drop")

p_popsize <- popsize_plot_data %>%
  ggplot(aes(x = pop_size, y = mean_gen_99, colour = scheme, linetype = scheme)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_gen_99 - se_gen_99, ymax = mean_gen_99 + se_gen_99), width = 10) +
  scale_x_continuous(breaks = c(50, 100, 200, 400)) +
  scale_colour_manual(values = scheme_colours) +
  labs(x = "Population size", y = "Generations to 99% recovery", colour = "Scheme", linetype = "Scheme") +
  theme_minimal(base_size = 12)

p_popsize

save_plot(p_popsize, "popsize_vs_recovery.png")


# 6. Cost versus recovery speed
# Each point is one population size; cost is the mean number of genotype data points up to 99% recovery

cost_plot_data <- cost_to_99 %>%
  group_by(scheme, pop_size) %>%
  summarise(mean_generation = mean(generation), mean_cost = mean(cumulative_cost), .groups = "drop")

p_cost <- cost_plot_data %>%
  ggplot(aes(x = mean_cost, y = mean_generation, colour = scheme, linetype = scheme)) +
  geom_line(linewidth = 1) +
  geom_point(aes(shape = factor(pop_size)), size = 3) +
  scale_x_continuous(labels = scales::label_comma()) +
  scale_colour_manual(values = scheme_colours) +
  labs(
    x = "Cumulative genotyping cost (data points)",
    y = "Generations to 99% recovery",
    colour = "Scheme",
    linetype = "Scheme",
    shape = "Population size"
  ) +
  theme_minimal(base_size = 12)

p_cost

save_plot(p_cost, "cost_vs_recovery.png")


# 7. Two-target pyramiding recovery

p_pyramid <- pyramid_summary %>%
  ggplot(aes(x = generation, y = recovery_ibd, colour = scheme, linetype = scheme)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0.7, 1), labels = scales::percent) +
  scale_colour_manual(values = scheme_colours) +
  labs(
    x = "Backcross generation",
    y = "Recurrent-parent genome recovery",
    colour = "Scheme",
    linetype = "Scheme"
  ) +
  theme_minimal(base_size = 12)

p_pyramid

save_plot(p_pyramid, "pyramiding_recovery.png")


# 8. Replicate convergence

p_convergence <- convergence_99 %>%
  ggplot(aes(x = replicate, y = running_mean_99, colour = scheme, linetype = scheme)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = scheme_colours) +
  labs(
    x = "Number of replicates",
    y = "Running mean generations to 99% recovery",
    colour = "Scheme",
    linetype = "Scheme"
  ) +
  theme_minimal(base_size = 12)

p_convergence

save_plot(p_convergence, "convergence.png")


# 9. Marker-estimated versus true IBD recovery
# Points on the dashed line mean the chip estimate equals the true IBD value

p_marker_ibd <- replicate_runs %>%
  ggplot(aes(x = recovery_ibd, y = recovery_marker, colour = scheme, shape = scheme)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_colour_manual(values = scheme_colours) +
  labs(x = "True IBD recovery", y = "Marker-estimated recovery", colour = "Scheme", shape = "Scheme") +
  theme_minimal(base_size = 12)

p_marker_ibd

save_plot(p_marker_ibd, "marker_vs_ibd.png")
