

## Null distribution of sctransform-style Pearson residual variance
## ----------------------------------------------------------------
## Model: x_gc ~ NB(mu = p_g * s_c, theta), i.e. counts explained entirely by sequencing depth -> no biological heterogeneity
## Pearson residual: z = (x - mu) / sqrt(mu + mu^2/theta)
## Under this null, Var_c(z) should be ~1, and the *estimate* over n cells should follow chi^2_{n-1} / (n-1).

library(ggplot2)
library(matrixStats)

set.seed(42)

n_cells <- 3000
n_genes <- 5000
theta   <- 100      # fixed overdispersion, cf. Lause/Berens/Kobak offset model
clip    <- sqrt(n_cells)

## simulation
s <- rlnorm(n_cells, meanlog = 0, sdlog = 0.4)   # per-cell depth factor
s <- s / mean(s)

# gene rates spanning rare -> highly expressed (mean counts ~1e-3 to ~50)
p <- 10^runif(n_genes, -3, 1.7)

mu <- outer(p, s)                                 # genes x cells
x  <- matrix(rnbinom(length(mu), mu = mu, size = theta),
             nrow = n_genes)

##  Pearson residuals against the TRUE mu (oracle null) 
z <- (x - mu) / sqrt(mu + mu^2 / theta)
z_clipped <- pmin(pmax(z, -clip), clip)

resvar_raw  <- rowVars(z)
resvar_clip <- rowVars(z_clipped)


## Gaussian residuals: z ~ N(0,1) instead of NB counts 
z_norm      <- matrix(rnorm(n_genes * n_cells), nrow = n_genes)
resvar_norm <- rowVars(z_norm)

comp_df <- rbind(
  data.frame(value = resvar_norm, source = "Simulated N(0,1) residuals"),
  data.frame(value = resvar_raw,  source = "Simulated null (unclipped)"),
  data.frame(value = resvar_clip, source = "Simulated null (clipped)"))
comp_df$source <- factor(comp_df$source, levels = unique(comp_df$source))


## Gaussian fitted to log10(residual variance)
lv <- log10(resvar_raw)
fit_norm <- data.frame(log_value = seq(min(lv), max(lv), length.out = 1000))
fit_norm$density <- dnorm(fit_norm$log_value,
                          mean = mean(lv),
                          sd   = sd(lv))
fit_norm$value <- 10^fit_norm$log_value

# plot
t_plot <- ggplot(comp_df, aes(x = value, fill = source)) +
  geom_histogram(aes(y = after_stat(density)), bins = 500,
                 position = "identity", alpha = 0.5, colour = NA) +
  geom_line(data = fit_norm, aes(x = value, y = density),
            inherit.aes = FALSE, linetype = 2, linewidth = 0.6) +
  geom_vline(xintercept = 1, linetype = 3, linewidth = 0.4) +
  scale_x_log10(breaks = scales::breaks_log(n = 4),
                labels = scales::label_log()) +
  coord_cartesian(xlim = range(10^-1,10^1), expand = FALSE) +
  theme_classic(base_size = 16) +
  labs(x = "Residual variance (log10)", y = "Density", fill = NULL)

ggsave('droplet/plots/regnoise_vs2/theor_dist.png', t_plot, width = 10, height = 5)


