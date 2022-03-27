# Please note that in order to make it possible to run the simulation
# on the free Binder cloud instance, the amount of data to be simulated
# had to be significantly reduced. We also skip the SLiM simulation
# because that would be too computationally intensive for Binder to
# handle.

# import required libraries -----------------------------------------------

library(slendr)

library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(cowplot)
library(forcats)

# define the model,  simulate tree-sequence data --------------------------

o <- population("o", time = 1, N = 100)
c <- population("c", time = 2500, N = 100, parent = o)
a <- population("a", time = 2800, N = 100, parent = c)
b <- population("b", time = 3700, N = 100, parent = a)
x1 <- population("x1", time = 4000, N = 5000, parent = c)
x2 <- population("x2", time = 4300, N = 5000, parent = x1)

gf <- gene_flow(from = b, to = x1, start = 5400, end = 5800, rate = 0.1)

model <- compile_model(
  populations = list(o, a, b, c, x1, x2), gene_flow = gf,
  generation_time = 1, sim_length = 6000
)

plot_model(model, sizes = FALSE, proportions = TRUE)

msprime(model, sequence_length = 100e6, recombination_rate = 1e-8)

ts <- ts_load(model) %>% ts_mutate(mutation_rate = 1e-8)

samples <- ts_samples(ts) %>% group_by(pop) %>% sample_n(25)

divergence <- ts_divergence(ts, split(samples$name, samples$pop))

f4ratio <- ts_f4ratio(
  ts, X = filter(samples, pop %in% c("x1", "x2"))$name,
  A = "a_1", B = "b_1", C = "c_1", O = "o_1"
)

# plotting code -----------------------------------------------------------

divergence %>%
  mutate(backend = "msprime", pair = paste(x, "-", y)) %>%
  ggplot(aes(fct_reorder(pair, divergence), color = backend, shape = backend, divergence)) +
  geom_point(size = 3) +
  xlab("population pair") + ylab("pairwise divergence") +
  theme_minimal() +
  scale_alpha_manual(values = c(1, 0.25)) +
  scale_color_manual(values = c("black", "darkgray")) +
  guides(shape = guide_legend("slendr simulation engine used:"),
         color = guide_legend("slendr simulation engine used:",
                              override.aes = list(size = 3))) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10),
        axis.text.x = element_text(hjust = 1, angle = 45, size = 9),
        axis.title.x = element_blank())

f4ratio %>%
  mutate(backend = "msprime",
         population = gsub("^(.*)_.*$", "\\1", X), alpha = alpha * 100) %>%
  ggplot(aes(population, alpha)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_jitter(aes(color = backend, shape = backend), alpha = 0.75,
              position = position_jitterdodge(jitter.width = 0.25)) +
  ylab(base::expression(italic("f")[4]~"-ratio ancestry proportion [%]")) +
  scale_color_manual(values = c("black", "darkgray")) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 20)) +
  theme(axis.text.x = element_text(size = 11),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank())
