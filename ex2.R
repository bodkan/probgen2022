# WARNING: A model like this is impossible to run using the extremely
# limited resources provided by the free cloud Binder instance. In order
# to even run this model, the population sizes had to be dramatically
# decreased. The result produced by this model will be very different
# to what is shown in the poster. However, the slendr workflow remains
# identical and gives a good idea what an R slendr pipeline can look like.

# import required libraries -----------------------------------------------

library(slendr)

library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(cowplot)
library(forcats)

# define the model,  simulate tree-sequence data --------------------------

map <- world(xrange = c(0, 10), yrange = c(0, 10),
             landscape = region(center = c(5, 5), radius = 5))

# NOTE:
# population sizes had to be scaled down dramatically and the emerging
# spatial clustering pattern will be very different as a result
p1 <- population("pop1", time = 1, N = 500, map = map, competition = 0)
p2 <- population("pop2", time = 1, N = 500, map = map, competition = 6)
p3 <- population("pop3", time = 1, N = 500, map = map, competition = 5)
p4 <- population("pop4", time = 1, N = 500, map = map, competition = 4)
p5 <- population("pop5", time = 1, N = 500, map = map, competition = 3)
p6 <- population("pop6", time = 1, N = 500, map = map, competition = 2)
p7 <- population("pop7", time = 1, N = 500, map = map, competition = 1)

model <- compile_model(
  populations = list(p1, p2, p3, p4, p5, p6, p7),
  generation_time = 1, sim_length = 5000, resolution = 0.1,
  mating = 0.1, dispersal = 0.05
)

# NOTE:
# the amount of sequence to be simulated had to be drastically reduced
# (diversity calculation will be extremely noisy)
slim(model, sequence_length = 100000, recombination_rate = 1e-8)

ts <- ts_load(model) %>% ts_simplify() %>% ts_mutate(mutation_rate = 1e-6)

locations <- ts_data(ts) %>% filter(time == max(time))

heterozygosity <- ts_samples(ts) %>%
  group_by(pop) %>%
  sample_n(50) %>%
  mutate(pi = ts_diversity(ts, name)$diversity)

# plotting code -----------------------------------------------------------

p_ex2_clustering <- ggplot() +
  geom_sf(data = map) +
  geom_sf(data = locations, aes(color = pop), size = 0.05, alpha = 0.25) +
  facet_grid(. ~ pop, switch = "x") +
  xlab("spatial distributions emerged in the simulation") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.background = element_blank()
  ) +
  guides(color = "none")

p_ex2_diversity <- ggplot(heterozygosity, aes(pop, pi, color = pop)) +
  geom_violin(color = "black") +
  geom_jitter(alpha = 0.5) +
  labs(y = "individual heterozygosity") +
  guides(color = "none") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), panel.grid.major.x = element_blank(),
        plot.margin = margin(t = 0.2, r = 0.2, b = -0.1, l = 0.2, "cm"))
plot_grid(
  p_ex2_diversity,
  p_ex2_clustering +
    theme(plot.margin = margin(t = 0, r = 0.4, b = 0, l = 1.8, "cm")),
  nrow = 2,
  rel_heights = c(1, 0.5),
  labels = c("B", "C")
)
