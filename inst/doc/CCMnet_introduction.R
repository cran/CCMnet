## ----setup, include=FALSE-----------------------------------------------------
# Set global chunk options
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.width = 6,
  fig.height = 4
)

# Load libraries
library(CCMnet)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(1)

## ----edge-example-poisson-----------------------------------------------------

ccm_sample <- sample_ccm(
  network_stats = list("edges"),
  prob_distr = list("poisson"),
  prob_distr_params = list(list(350)), 
  population = 50
)

summary(ccm_sample)

print(ccm_sample)

# Compare MCMC samples to theoretical Poisson distribution
ccm_sample<- CCM_theoretical_check(ccm_sample, n_sim = 1000)

# Plot density with theoretical distribution
plot(ccm_sample, stats = "edges", type = "hist", include_theoretical = TRUE)

## ----edge-example-uniform-----------------------------------------------------

ccm_sample <- sample_ccm(
  network_stats = list("edges"),
  prob_distr = list("uniform"), 
  prob_distr_params = list(list(0)),
  population = 20,
  sample_size = 20000L,
  burnin = 100000L,
  interval = 1000L,
)

ccm_sample<- CCM_theoretical_check(ccm_sample, n_sim = 10000)
plot(ccm_sample, stats = "edges", type = "hist", include_theoretical = TRUE)

## ----edge-example-np----------------------------------------------------------

n_max <- choose(50, 2)
alpha <- dpois(0:n_max, lambda = 50) + dpois(0:n_max, lambda = 100)
prob_distr_params <- alpha / sum(alpha)

ccm_sample <- sample_ccm(
  network_stats = list("edges"),
  prob_distr = list("np"),
  prob_distr_params = list(list(prob_distr_params)),
  population = 50L,
  sample_size = 10000L,
  burnin = 100000L
)

ccm_sample<- CCM_theoretical_check(ccm_sample, n_sim = 50000)
plot(ccm_sample, stats = "edges", type = "hist", include_theoretical = TRUE)

## ----degree-example-dirmult---------------------------------------------------

ccm_sample<- sample_ccm(network_stats='DegreeDist',
                        prob_distr='DirMult',
                        prob_distr_params=list(list(c(2,21,15,12))), 
                        population = 100L, 
                        sample_size = 10000L,
                        burnin=100000L, 
                        interval=1000L) 

ccm_sample <- CCM_theoretical_check(ccm_sample, n_sim = 1000)
plot(ccm_sample, 
     stats = paste0("deg", 0:3), 
     type = "hist", 
     include_theoretical = TRUE
)    

## ----degmix_triangle-example-normal-------------------------------------------
mean_vec = c(23, 66, 44, 20, 120, 80)
var_mat = matrix(data = c(22, -3, -2, -5, -6, -4,
                          -3, 58, -7, -14, -18, -12,
                          -2, -7, 41, -9, -12, -8,
                          -5, -14, -9, 75, -25, -17,
                          -6, -18, -12, -25, 89, -22,
                          -4, -12, -8, -17, -22, 68), ncol = 6)
prob_distr_params = list(list(mean_vec, var_mat),
                         list(10,3))

ccm_sample <- sample_ccm(network_stats=c('DegMixing', 'Triangles'),
                         prob_distr=c('Normal', 'Normal'),
                         prob_distr_params=prob_distr_params, 
                         sample_size = 10000L,
                         burnin=100000L, 
                         interval=1000L,
                         population=500L) 

ccm_sample <- CCM_theoretical_check(ccm_sample, n_sim = 1000)
plot(ccm_sample, 
     stats = c("DM11", "DM12", "DM13", "DM22", "DM23", "DM33", "triangles"), 
     type = "hist", 
     include_theoretical = TRUE)
CCM_traceplot(ccm_sample, stats = "triangles")

## ----data-school--------------------------------------------------------------

utils::data("faux.mesa.high", package = "ergm", envir = environment())
utils::data("faux.desert.high", package = "ergm", envir = environment())
utils::data("faux.dixon.high", package = "ergm", envir = environment())
utils::data("faux.magnolia.high", package = "ergm", envir = environment())

mesa_net = intergraph::asIgraph(faux.mesa.high)
desert_net = intergraph::asIgraph(faux.desert.high)
dixon_net = intergraph::asIgraph(faux.dixon.high)
magnolia_net = intergraph::asIgraph(faux.magnolia.high)

# Create summary table
hs_summary <- data.frame(
  High_School = c("Mesa High", "Desert High", "Dixon High", "Magnolia High"),
  Nodes = c(
    igraph::vcount(mesa_net),
    igraph::vcount(desert_net),
    igraph::vcount(dixon_net),
    igraph::vcount(magnolia_net)
  ),
  Edges = c(
    igraph::gsize(mesa_net),
    igraph::gsize(desert_net),
    igraph::gsize(dixon_net),
    igraph::gsize(magnolia_net)
  )
)

# Render table
kableExtra::kable(
  hs_summary,
  caption = "Summary of high school friendship networks from the ERGM package"
)


## ----new-school-estimation----------------------------------------------------

#Normal
densities <- hs_summary$Edges / choose(hs_summary$Nodes, 2)
sigma <- sd(densities)
prior <- RBesT::mixnorm(c(1, 0.5, 1), sigma = sigma)
post_normal <- RBesT::postmix(prior, m = mean(densities), n = length(densities))


## ----density-example-new_school-normal----------------------------------------
population = 100
n_samples = 10000L

#CCM
ccm_sample <- sample_ccm(
  network_stats = list("Density"),
  prob_distr = list("Normal"),
  prob_distr_params = list(list(post_normal[2], (post_normal[3])^2)), 
  population = population,
  sample_size = n_samples,
  burnin = 100000L,
  interval = 1000L,
)

#ERGM
net <- network::network(population, 
                        directed = FALSE)
ergm_sample <- simulate(net ~ edges,
                        coef = log(post_normal[2] / (1 - post_normal[2])),
                        nsim = n_samples,
                        output = "stats") / choose(population, 2)

#G(n,m)
m_edges = round(post_normal[2]*choose(population, 2))
er_gnm_list <- replicate(n_samples, 
                         igraph::sample_gnm(n = population, m = m_edges, directed = FALSE), 
                         simplify = FALSE)
er_gnm_stats <- sapply(er_gnm_list, igraph::ecount)
er_gnm_sample <- er_gnm_stats / choose(population, 2)

## ----density-example-new_school-normal-compare--------------------------------

Network_samples <- data.frame(
  value = c(ccm_sample$mcmc_stats$density, ergm_sample, er_gnm_sample),
  model = c(rep("CCMnet", n_samples), rep("ERGM", n_samples), rep("ER", n_samples))
)

ggplot2::ggplot( Network_samples %>% filter(model != "ER"), 
                 aes(x = value, fill = model)
) +
  geom_density(alpha = 0.25, linewidth = .1) +
  geom_vline(
    aes(xintercept = er_gnm_sample[1], colour = "ER"),
    linewidth = 1,
    linetype = "solid"
  ) +
  scale_fill_manual(
    values = c(CCMnet = "red", ERGM = "green")
  ) +
  scale_colour_manual(
    values = c(ER = "blue")
  ) +
  labs(x = "Density", y = "Probability density", fill = "Model", colour = "Model") +
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank())


## ----density-example-sampled_school-beta--------------------------------------

res = data.frame(model = NULL,
                 value = NULL)


dixon_deg = igraph::degree(dixon_net)
prior.unif <- RBesT::mixbeta(c(1, 1, 1))
N = igraph::vcount(dixon_net)

n_samples <- 1000L

for (num_sample_nodes in seq(40,240,40)) {
  dixon_deg_sample = sample(x = dixon_deg, size = num_sample_nodes, replace = FALSE)
  
  r=sum(dixon_deg_sample)
  n=sum(num_sample_nodes*(N-1))
  posterior.sum_beta <- RBesT::postmix(prior.unif, 
                                       n=n, 
                                       r=r)
  
  alpha_post <- posterior.sum_beta[2]
  beta_post  <- posterior.sum_beta[3]
  
  # Infinite-population variance
  var_inf <- (alpha_post * beta_post) /
    ((alpha_post + beta_post)^2 * (alpha_post + beta_post + 1))
  
  # Finite population correction
  fpc <- (N - num_sample_nodes) / (N - 1)
  var_fpc <- var_inf * fpc
  
  # Moment-matched Beta parameters
  mu <- alpha_post / (alpha_post + beta_post)
  S <- mu * (1 - mu) / var_fpc - 1
  posterior.sum_beta[2] <- mu * S
  posterior.sum_beta[3] <- (1 - mu) * S
  
  ccm_sample <- sample_ccm(
    network_stats = list("Density"),
    prob_distr = list("Beta"),
    prob_distr_params = list(list(posterior.sum_beta[2], posterior.sum_beta[3])), 
    population = N,
    sample_size = n_samples,
    burnin = 100000L
  )

  res = bind_rows(res, data.frame(
    model = rep("CCMnet", length(ccm_sample$mcmc_stats$density)),
    value = ccm_sample$mcmc_stats$density,
    sample = num_sample_nodes))
  
  net <- network::network(N, directed = FALSE)
  p = posterior.sum_beta[2]/(posterior.sum_beta[2] + posterior.sum_beta[3])
  theta = log(p / (1-p))
  
  ERGM <- simulate(
    net ~ edges,
    coef = theta,
    nsim = n_samples,
    output = "stats"
  ) / choose(N, 2)

  ER = rep(posterior.sum_beta[2]/(posterior.sum_beta[2] + posterior.sum_beta[3]), n_samples)
  
  res <- bind_rows(res, data.frame(
    value = c(ERGM, ER),
    model = c(rep("ERGM", n_samples), rep("ER", n_samples)),
    sample = c(rep(num_sample_nodes, n_samples), rep(num_sample_nodes, n_samples))
  ))

}

## ----density-example-sampled_school-compare-----------------------------------

# ER vertical lines
ER_lines <- res %>%
  filter(model == "ER") %>%
  group_by(sample, model) %>%
  summarise(xintercept = unique(value), .groups = "drop")

ggplot2::ggplot(
  res %>% filter(model != "ER"),
  aes(x = value, fill = model)
) +
  geom_density(alpha = 0.25, linewidth = .1) +
  
  # ER vertical lines
  geom_vline(
    data = ER_lines,
    aes(xintercept = xintercept, colour = model),
    linetype = "solid",
    linewidth = 1
  ) +
  
  scale_fill_manual(
    values = c(
      CCMnet = "red",
      ERGM   = "green"
    )
  ) +
  scale_colour_manual(
    values = c(
      ER = "blue"
    )
  ) +
  
  labs(
    x = "Density",
    y = "Probability density",
    fill = "Model",
    colour = NULL
  ) +
  
  guides(
    colour = guide_legend(override.aes = list(fill = NA)),
    fill   = guide_legend(override.aes = list(colour = NA))
  ) +
  
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),  
        axis.ticks.y = element_blank()) +
  facet_wrap(~sample, scales =  "free")


## ----density-example-sampled_school-ensemble----------------------------------

ccm_sample <- sample_ccm(
  network_stats = list("Density"),
  prob_distr = list("Normal"),
  prob_distr_params = list(list(post_normal[2], (post_normal[3])^2)), 
  population = 100L,
  sample_size = 10000L,
  burnin = 100000L
)

school_ensemble <- sample_ccm(
  network_stats = list("Density"),
  prob_distr = list("Normal"),
  prob_distr_params = list(list(post_normal[2], (post_normal[3])^2)), 
  population = 100L,
  sample_size = 10L,
  burnin = 1L,
  interval = 1000L,
  initial_g = ccm_sample$g[[1]], 
  use_initial_g = TRUE, 
  stats_only = FALSE)

class(school_ensemble$g)
length(school_ensemble$g)
class(school_ensemble$g[[1]])


