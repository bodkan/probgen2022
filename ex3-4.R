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

# example 3 ---------------------------------------------------------------

source("coordinates.R")

map <- world(xrange = c(-15, 60), yrange = c(20, 65), crs = 3035)

ooa <- population(
  "OOA", time = 50000, N = 500, remove = 23000,
  map = map, center = c(33, 30), radius = 400e3
) %>%
  move(trajectory = ooa_trajectory, start = 50000, end = 40000, snapshots = 30)

ehg <- population(
  "EHG", time = 28000, N = 1000, parent = ooa, remove = 6000,
  map = map, polygon = R1
)

eur <- population(
  "EUR", time = 30000, N = 2000, parent = ooa,
  map = map, polygon = R2
) %>%
  resize(N = 10000, time = 5000, end = 0, how = "exponential")

ana <- population(
  "ANA", time = 25000, N = 4000, parent = ooa, remove = 3000,
  map = map, polygon = R3
) %>%
  expand_range(by = 3e6, start = 10000, end = 7000, polygon = R4, snapshots = 15)

yam <- population(
  "YAM", time = 7000, N = 600, parent = ehg, remove = 2500,
  map = m, polygon = R5
) %>%
  move(trajectory = list(c(15, 50)), start = 5000, end = 3000, snapshots = 10)

gf <- list(
  gene_flow(ana, to = yam, rate = 0.5, start = 6500, end = 5000),
  gene_flow(ana, to = eur, rate = 0.6, start = 8000, end = 6000),
  gene_flow(yam, to = eur, rate = 0.7, start = 3500, end = 3000)
)

model <- compile_model(
  populations = list(ooa, ehg, eur, ana, yam), gene_flow = gf,
  generation_time = 30, resolution = 10e3,
  competition = 130e3, mating = 100e3,
  dispersal = 70e3,
)

samples <- schedule_sampling(
  model, times = seq(0, 50000, by = 1000),
  list(ehg, 20), list(ana, 20), list(yam, 20), list(eur, 20)
)

plot_model(model, sizes = FALSE)

plot_map(model)

# WARNING: We had to dramatcally scale down the amount of sequence
# simulated because Binder cannot handle the amount of data necessary
# for the example from the poster.

slim(
  model, burnin = 200000, sampling = samples,
  sequence_length = 10000, recombination_rate = 1e-8, verbose = TRUE
)

# example 4 ---------------------------------------------------------------

ts <- ts_load(model)

ts_small <- ts_simplify(ts, c("EUR_599", "ANA_322", "EHG_7", "EUR_578", "EUR_501", "YAM_30"))

tree <- ts_phylo(ts_small, i = 1)
nodes <- ts_data(tree)
branches <- ts_branches(tree)

ancestors <- ts_ancestors(ts, "EUR_599")

# plotting code -----------------------------------------------------------

# panel B

library(ggtree)

df <- as_tibble(nodes) %>% select(node = phylo_id, pop)

highlight_nodes <- as_tibble(nodes) %>% dplyr::filter(name == "EUR_599") %>% .$phylo_id
p_tree <- ggtree(tree, aes(color = pop, fill = pop)) %<+% df +
  geom_tiplab(align = TRUE, geom = "label", offset = 2000,
              color = "white", fontface = "bold", size = 3) +
  geom_tiplab(align = TRUE, geom = NULL, linetype = "dotted", size = 0) +
  geom_point2(aes(subset = (node %in% highlight_nodes)), color = "black", size = 3) +
  geom_label2(aes(label = label, subset = !isTip),
              color = "black", size = 3) +
  theme_tree2() +
  theme(legend.position = "none") +
  xlab("time before present [years ago]") +
  scale_x_continuous(limits = c(-55000, 22000), labels = abs,
                     breaks = -c(60, 40, 20, 0) * 1000)
p_tree <- revts(p_tree)
p_tree

# panel C

ggplot() +
  geom_sf(data = map) +
  geom_sf(data = branches, aes(color = parent_pop), size = 0.5) +
  geom_sf(data = filter(nodes, is.na(name)),
          aes(color = pop, shape = pop), size = 5) +
  geom_sf_label(data = nodes[!nodes$sampled, ],
                aes(label = node_id, fill = pop), size = 3) +
  geom_sf_label(data = nodes[nodes$sampled, ],
                aes(label = node_id, fill = pop), size = 3,
                fontface = "bold", color = "white") +
  coord_sf(xlim = c(4047066.1, 8688656.9),
           ylim = c(757021.7, 4972983.3), expand = 0) +
  labs(x = "longitude", y = "latitude") +
  guides(fill = guide_legend("", override.aes = aes(label = ""))) +
  guides(color = "none", shape = "none") +
  theme_bw() +
  theme(legend.position = "bottom")

# panel D
chrom_names <- setNames(
  c("EUR_599 (chromosome 1, node 10)", "EUR_599 (chromosome 2, node 11)"),
  unique(ancestors$node_id)
)

ggplot() +
  geom_sf(data = map) +
  geom_sf(data = ancestors, size = 0.5, aes(alpha = parent_time)) +
  geom_sf(data = sf::st_set_geometry(ancestors, "parent_location"),
          aes(shape = parent_pop, color = parent_pop)) +
  geom_sf(data = filter(ts_data(ts), name == "EUR_599"), size = 3) +
  coord_sf(expand = 0) +
  labs(x = "longitude", y = "latitude") +
  theme_bw() +
  facet_grid(. ~ node_id, labeller = labeller(node_id = chrom_names))
