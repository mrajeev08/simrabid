# Figure for hdsi application

# synthetic example
library(raster)
library(data.table)
library(rgdal)
library(rmapshaper)
library(ggraph)
library(tidygraph)
library(igraph)
library(patchwork)
library(dplyr)
devtools::load_all()

# set up Serengeti raster
sd_shape <- readOGR("data/SD_shape/Serengeti_villages UTM_region.shp")
sd_dissolved <- as(ms_dissolve(sd_shape), "SpatialLinesDataFrame")
pop_data <- read.csv("data/SerengetiPop.csv")
census_data <- read.csv("data/SDcompiled.csv")
rast <- setup_space(shapefile = sd_shape,
                    resolution = 1000, id_col = "VILLCODE",
                    use_fasterize = FALSE)

set.seed(1131)
start_pop <- get.demdata(shapefile = sd_shape, res_meters = 1000,
                         pop = pop_data, census = census_data)$start_pop
start_up <- setup_sim(tmax = 52, start_pop, rast)

vacc_dt <- sim_campaigns(locs = 1:max(rast[], na.rm = TRUE),
                         campaign_prob = 0,
                         coverage = 0,
                         sim_years = 8, burn_in_years = 2,
                         steps_in_year = 52,
                         sample_tstep = function(n) sample.int(52, n, replace = TRUE))

# function to return the line list
get_Idt <- function(names = c("I_dt")) {
  mget(names, envir = parent.frame(3))
}

# simulate
set.seed(1233)
out <- simrabid(start_up, start_vacc = 0, I_seeds = 5, vacc_dt,
                params = list(R0 = 1.2, k = 0.5, iota = 4),
                death_prob = get_prob(rate = 0.48, step = 52), # annual death rate to prob
                waning_prob = get_prob(rate = 1/3, step = 52), # annual waning to prob
                birth_prob = get_prob(rate = 0.52, step = 52), # annual birth rate to prob
                days_in_step = 7, # get tmax from start up!
                generation_fun = gen_gamma,
                dispersal_fun = dispersal_gamma,
                secondary_fun = negbinom_constrained,
                incursion_fun = sim_incursions_pois,
                sequential = TRUE, allow_invalid = FALSE,
                leave_bounds = FALSE, max_tries = 10,
                summary_funs = list(get_Idt = get_Idt),
                track = TRUE,
                row_probs = NULL,
                vacc_type = "coverage")

case_data <- out$get_Idt$I_dt

# filter case data to ones that became infectious in time step t
case_data <- case_data[id != 0]
# filter out all singletons
nprogens <- tabulate(case_data$progen_id)
case_data$nprogens <- nprogens[case_data$id]
case_data <- case_data[!is.na(nprogens) & nprogens > 0]
case_data$detected <- rbinom(nrow(case_data), 1, 0.5)

# time series
ts <- ggplot(data = case_data, aes(x = floor(t_infectious))) +
  geom_bar(aes(fill = factor(detected)), position = "stack", color = "grey50") +
  scale_fill_manual(values = c("white", "darkred"), guide = "none") +
  labs(y = "Number of cases",   x = "Week") +
  cowplot::theme_minimal_hgrid()

# assign to lineage & generate substitutions
transmission <- case_data[progen_id > 0]
trees <- data.frame(from = transmission$progen_id,
                    to = transmission$id)
g <- graph_from_data_frame(trees, directed = TRUE)
lineage <- components(g)$membership
case_data$lineage <- lineage[match(case_data$id, names(lineage))]

ggt <- as_tbl_graph(g)
ggt %>%
  mutate(group = group_components(),
         detected = case_data$detected[match(name, case_data$id)],
         x = layout_as_tree(ggt)[, 1],
         y = layout_as_tree(ggt)[, 2]) -> ggt

group_cols <- c("#d95f02", "#7570b3", "#1b9e77")
names(group_cols) <- 1:3
tree <- ggraph(ggt, layout = "dendrogram") +
  geom_edge_bend2() +
  geom_node_point(aes(shape = factor(detected), color = factor(group)),
                  size = 2, fill = "white") +
  scale_shape_manual(values = c(21, 16), guide = "none") +
  scale_color_manual(values = group_cols,
                     guide = "none") +
  cowplot::theme_map()

values(rast) <- start_pop
pop_df <- as.data.frame(rast, xy = TRUE)
pop_df$pop <- start_pop
names(group_cols) <- 3:1
mapped <- ggplot() +
  geom_raster(data = pop_df, aes(x, y, fill = pop), interpolate = TRUE,
              alpha = 0.75) +
  geom_point(
    data = case_data, aes(x = x_coord, y = y_coord,
                          shape = factor(detected),
                          color = factor(lineage)), size = 4, alpha = 0.75) +
  scale_shape_manual(values = c(21, 16), guide = "none") +
  scale_color_manual(values = group_cols,
                     guide = "none") +
  scale_fill_distiller(na.value = "white", direction = 1,
                       palette = "Greys",
                       trans = "log", guide = "none"
  ) +
  cowplot::theme_map() +
  labs(tag = "B")

# Assemble plots ----------------

# Parameters -----
par_theme <- cowplot::theme_minimal_hgrid(font_size = 12, line_size = 1) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())

# R0
set.seed(123)
secondaries <- data.frame(secondaries = rnbinom(1000, mu = 1.2, size = 0.5))
sec <- ggplot(data = secondaries, aes(x = secondaries)) +
  geom_histogram(binwidth = 1, color = "lightgray", fill = "red", alpha = 1) +
  geom_vline(xintercept = 1.2, color = "red", linetype = 2, size = 1.1) +
  xlab("Secondary cases") +
  labs(tag = "A") +
  ylab("Frequency") +
  par_theme

generation <- data.frame(Days = seq(0, 365, 1),
                         Density = dgamma(seq(0, 365, 1),
                                          shape = 1.46, scale = 16.1))
gen <- ggplot(data = generation, aes (x = Days, y = Density)) +
  geom_line(color = "darkred", size = 1.2) +
  xlab("Generation time \n (days)") +
  geom_vline(xintercept = 22.3, color = "darkred", linetype = 2, size = 1.1,
             alpha = 0.5) +
  par_theme


dispersal <- as.data.frame(list(km = seq(0, 10, 0.01),
                                Density = dgamma(seq(0, 10, 0.01), shape = 0.8,
                                                 scale = 1)))
disp <- ggplot(data = dispersal, aes (x = km, y = Density)) +
  geom_line(color = "#500000", size = 1.2) +
  labs(x = "Distance to next bite \n (km)", y = "") +
  geom_vline(xintercept = 0.88, color = "#500000", linetype = 2, size = 1.2, alpha = 0.75) +
  par_theme

toprow <- sec | gen | disp
bottomrow <- mapped | ts | tree

# figure 1
fig1 <- toprow / bottomrow
ggsave("vignettes/hdsi_fig.jpeg", fig1, width = 7, height = 5)
