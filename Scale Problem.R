library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# --- Test Statistics Functions ---

savage_test <- function(x, y) {
  n <- length(x); m <- length(y); N <- n + m
  combined <- c(x, y); ranks <- rank(combined)
  z <- c(rep(1, n), rep(0, m))  
  savage_scores <- map_dbl(1:N, ~sum(1/(N + 1 - ranks[.x]:N)))
  S <- sum(savage_scores * z)
  all_scores <- map_dbl(1:N, ~sum(1/(N + 1 - .x:N)))
  E0 <- n * mean(all_scores)
  Var0 <- (n*m)/(N*(N-1)) * sum((all_scores - mean(all_scores))^2)
  (S - E0)/sqrt(Var0)
}

mood_test <- function(x, y) {
  n <- length(x); m <- length(y); N <- n + m
  R <- rank(c(x, y))[1:n]
  M <- sum((R - (N + 1)/2)^2)
  E0 <- n*(N^2 - 1)/12
  Var0 <- n*m*(N + 1)*(N^2 - 4)/180
  (M - E0)/sqrt(Var0)
}

klotz_test <- function(x, y) {
  n <- length(x)
  m <- length(y)
  N <- n + m
  combined <- c(x, y)
  ranks <- rank(combined)
  z <- c(rep(1, n), rep(0, m))  
  inv_norm_quantiles <- qnorm(ranks / (N + 1))
  K <- sum((inv_norm_quantiles^2) * z) 
  
  all_scores <- qnorm((1:N)/(N+1))^2
  E0 <- n * mean(all_scores)
  
  Var0 <- (n*m)/(N*(N-1)) * sum((all_scores - mean(all_scores))^2)
  
  K_std <- (K - E0) / sqrt(Var0)
  return(K_std)
}

# ==============================================
# Klotz Test Simulation and Visualization
# ==============================================

# Simulation for Klotz test
simulate_klotz <- function(distr, n, m, theta, B = 10000) {
  K_stats <- numeric(B)
  
  for (i in 1:B) {
    if (distr == "Normal") {
      x <- rnorm(n, mean = 0, sd = 1)
      y <- rnorm(m, mean = 0, sd = theta)
    } else if (distr == "Exponential") {
      x <- rexp(n, rate = 1)
      y <- rexp(m, rate = 1 / theta)
    } else if (distr == "Cauchy") {
      x <- rcauchy(n, location = 0, scale = 1)
      y <- rcauchy(m, location = 0, scale = theta)
    } else if (distr == "Geometric") {
      x <- rgeom(n, prob = 0.5)
      y <- rgeom(m, prob = 1 / (1 + theta))
    }
    
    K_stats[i] <- klotz_test(x, y)
  }
  
  return(K_stats)
}

# Run Klotz simulations
distributions <- c("Normal", "Exponential", "Cauchy", "Geometric")
sample_sizes <- list(c(10, 12), c(70, 65), c(240, 250))
theta_values <- c(1, 2, 6)

set.seed(123)
klotz_results <- list()

for (distr in distributions) {
  for (sizes in sample_sizes) {
    n <- sizes[1]
    m <- sizes[2]
    key <- paste(distr, "n=", n, "m=", m)
    
    sim_H0 <- simulate_klotz(distr, n, m, theta = 1, B = 10000)
    sim_H1_theta2 <- simulate_klotz(distr, n, m, theta = 2, B = 10000)
    sim_H1_theta6 <- simulate_klotz(distr, n, m, theta = 6, B = 10000)
    
    klotz_results[[key]] <- list(
      H0 = sim_H0,
      H1_theta2 = sim_H1_theta2,
      H1_theta6 = sim_H1_theta6,
      distr = distr,
      n = n,
      m = m
    )
  }
}

# Plotting for Klotz test under Null Distribution
klotz_plot_list_H0 <- list()

for (i in seq_along(klotz_results)) {
  res <- klotz_results[[i]]
  df <- data.frame(K = res$H0)
  
  p <- ggplot(df, aes(x = K)) +
    geom_histogram(
      aes(y = ..density..),
      bins = 30,
      fill = "#69b3a2",
      color = "white",
      alpha = 0.8,
      linewidth = 0.3
    ) +
    geom_density(color = "#404080", linewidth = 1.2) +
    labs(
      title = paste("Klotz Test:", res$distr, "(n =", res$n, ", m =", res$m, ")"),
      x = "Standardized Klotz Statistic",
      y = "Density"
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
  
  klotz_plot_list_H0[[i]] <- p
}

# Plotting for Klotz test under Alternatives
klotz_plot_list_H1_theta2 <- list()

for (i in seq_along(klotz_results)) {
  res <- klotz_results[[i]]
  df <- data.frame(K = res$H1_theta2)
  
  p <- ggplot(df, aes(x = K)) +
    geom_histogram(
      aes(y = ..density..),
      bins = 30,
      fill = "#69b3a2",
      color = "white",
      alpha = 0.8,
      linewidth = 0.3
    ) +
    geom_density(color = "#404080", linewidth = 1.2) +
    labs(
      title = paste("Klotz Test:", res$distr, "(n =", res$n, ", m =", res$m, ")"),
      x = "Standardized Klotz Statistic",
      y = "Density"
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
  
  klotz_plot_list_H1_theta2[[i]] <- p
}

klotz_plot_list_H1_theta6 <- list()

for (i in seq_along(klotz_results)) {
  res <- klotz_results[[i]]
  df <- data.frame(K = res$H1_theta6)
  
  p <- ggplot(df, aes(x = K)) +
    geom_histogram(
      aes(y = ..density..),
      bins = 30,
      fill = "#69b3a2",
      color = "white",
      alpha = 0.8,
      linewidth = 0.3
    ) +
    geom_density(color = "#404080", linewidth = 1.2) +
    labs(
      title = paste("Klotz Test:", res$distr, "(n =", res$n, ", m =", res$m, ")"),
      x = "Standardized Klotz Statistic",
      y = "Density"
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
  
  klotz_plot_list_H1_theta6[[i]] <- p
}

# Q-Q Plots for Klotz test
klotz_qq_plots <- list()

for (i in seq_along(klotz_results)) {
  res <- klotz_results[[i]]
  qq_data <- data.frame(K = res$H0)
  
  qq_plot <- ggplot(qq_data, aes(sample = K)) +
    stat_qq(distribution = qnorm, size = 1) +
    stat_qq_line(distribution = qnorm, color = "red", linewidth = 1) +
    labs(
      title = paste("Klotz Test:", res$distr, "(n =", res$n, ", m =", res$m, ")"),
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 8))
  
  klotz_qq_plots[[i]] <- qq_plot
}

# ==============================================
# Mood Test Simulation and Visualization
# ==============================================

# Simulation for Mood test
simulate_mood <- function(distr, n, m, theta, B = 10000) {
  M_stats <- numeric(B)
  
  for (i in 1:B) {
    if (distr == "Normal") {
      x <- rnorm(n, mean = 0, sd = 1)
      y <- rnorm(m, mean = 0, sd = theta)
    } else if (distr == "Exponential") {
      x <- rexp(n, rate = 1)
      y <- rexp(m, rate = 1 / theta)
    } else if (distr == "Cauchy") {
      x <- rcauchy(n, location = 0, scale = 1)
      y <- rcauchy(m, location = 0, scale = theta)
    } else if (distr == "Geometric") {
      x <- rgeom(n, prob = 0.5)
      y <- rgeom(m, prob = 1 / (1 + theta))
    }
    
    M_stats[i] <- mood_test(x, y)
  }
  
  return(M_stats)
}

# Run Mood simulations
set.seed(123)
mood_results <- list()

for (distr in distributions) {
  for (sizes in sample_sizes) {
    n <- sizes[1]
    m <- sizes[2]
    key <- paste(distr, "n=", n, "m=", m)
    
    sim_H0 <- simulate_mood(distr, n, m, theta = 1, B = 10000)
    sim_H1_theta2 <- simulate_mood(distr, n, m, theta = 2, B = 10000)
    sim_H1_theta6 <- simulate_mood(distr, n, m, theta = 6, B = 10000)
    
    mood_results[[key]] <- list(
      H0 = sim_H0,
      H1_theta2 = sim_H1_theta2,
      H1_theta6 = sim_H1_theta6,
      distr = distr,
      n = n,
      m = m
    )
  }
}

# Plotting for Mood test under Null Distribution
mood_plot_list_H0 <- list()

for (i in seq_along(mood_results)) {
  res <- mood_results[[i]]
  df <- data.frame(M = res$H0)
  
  p <- ggplot(df, aes(x = M)) +
    geom_histogram(
      aes(y = ..density..),
      bins = 30,
      fill = "#69b3a2",
      color = "white",
      alpha = 0.8,
      linewidth = 0.3
    ) +
    geom_density(color = "#404080", linewidth = 1.2) +
    labs(
      title = paste("Mood Test:", res$distr, "(n =", res$n, ", m =", res$m, ")"),
      x = "Standardized Mood Statistic",
      y = "Density"
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
  
  mood_plot_list_H0[[i]] <- p
}

# Plotting for Mood test under Alternatives
mood_plot_list_H1_theta2 <- list()

for (i in seq_along(mood_results)) {
  res <- mood_results[[i]]
  df <- data.frame(M = res$H1_theta2)
  
  p <- ggplot(df, aes(x = M)) +
    geom_histogram(
      aes(y = ..density..),
      bins = 30,
      fill = "#69b3a2",
      color = "white",
      alpha = 0.8,
      linewidth = 0.3
    ) +
    geom_density(color = "#404080", linewidth = 1.2) +
    labs(
      title = paste("Mood Test:", res$distr, "(n =", res$n, ", m =", res$m, ")"),
      x = "Standardized Mood Statistic",
      y = "Density"
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
  
  mood_plot_list_H1_theta2[[i]] <- p
}

mood_plot_list_H1_theta6 <- list()

for (i in seq_along(mood_results)) {
  res <- mood_results[[i]]
  df <- data.frame(M = res$H1_theta6)
  
  p <- ggplot(df, aes(x = M)) +
    geom_histogram(
      aes(y = ..density..),
      bins = 30,
      fill = "#69b3a2",
      color = "white",
      alpha = 0.8,
      linewidth = 0.3
    ) +
    geom_density(color = "#404080", linewidth = 1.2) +
    labs(
      title = paste("Mood Test:", res$distr, "(n =", res$n, ", m =", res$m, ")"),
      x = "Standardized Mood Statistic",
      y = "Density"
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
  
  mood_plot_list_H1_theta6[[i]] <- p
}

# Q-Q Plots for Mood test
mood_qq_plots <- list()

for (i in seq_along(mood_results)) {
  res <- mood_results[[i]]
  qq_data <- data.frame(M = res$H0)
  
  qq_plot <- ggplot(qq_data, aes(sample = M)) +
    stat_qq(distribution = qnorm, size = 1) +
    stat_qq_line(distribution = qnorm, color = "red", linewidth = 1) +
    labs(
      title = paste("Mood Test:", res$distr, "(n =", res$n, ", m =", res$m, ")"),
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 8))
  
  mood_qq_plots[[i]] <- qq_plot
}

# ==============================================
# Savage Test Simulation and Visualization
# ==============================================

# Simulation for Savage test
simulate_savage <- function(distr, n, m, theta, B = 10000) {
  S_stats <- numeric(B)
  
  for (i in 1:B) {
    if (distr == "Normal") {
      x <- rnorm(n, mean = 0, sd = 1)
      y <- rnorm(m, mean = 0, sd = theta)
    } else if (distr == "Exponential") {
      x <- rexp(n, rate = 1)
      y <- rexp(m, rate = 1 / theta)
    } else if (distr == "Cauchy") {
      x <- rcauchy(n, location = 0, scale = 1)
      y <- rcauchy(m, location = 0, scale = theta)
    } else if (distr == "Geometric") {
      x <- rgeom(n, prob = 0.5)
      y <- rgeom(m, prob = 1 / (1 + theta))
    }
    
    S_stats[i] <- savage_test(x, y)
  }
  
  return(S_stats)
}

# Run Savage simulations
set.seed(123)
savage_results <- list()

for (distr in distributions) {
  for (sizes in sample_sizes) {
    n <- sizes[1]
    m <- sizes[2]
    key <- paste(distr, "n=", n, "m=", m)
    
    sim_H0 <- simulate_savage(distr, n, m, theta = 1, B = 10000)
    sim_H1_theta2 <- simulate_savage(distr, n, m, theta = 2, B = 10000)
    sim_H1_theta6 <- simulate_savage(distr, n, m, theta = 6, B = 10000)
    
    savage_results[[key]] <- list(
      H0 = sim_H0,
      H1_theta2 = sim_H1_theta2,
      H1_theta6 = sim_H1_theta6,
      distr = distr,
      n = n,
      m = m
    )
  }
}

# Plotting for Savage test under Null Distribution
savage_plot_list_H0 <- list()

for (i in seq_along(savage_results)) {
  res <- savage_results[[i]]
  df <- data.frame(S = res$H0)
  
  p <- ggplot(df, aes(x = S)) +
    geom_histogram(
      aes(y = ..density..),
      bins = 30,
      fill = "#69b3a2",
      color = "white",
      alpha = 0.8,
      linewidth = 0.3
    ) +
    geom_density(color = "#404080", linewidth = 1.2) +
    labs(
      title = paste("Savage Test:", res$distr, "(n =", res$n, ", m =", res$m, ")"),
      x = "Standardized Savage Statistic",
      y = "Density"
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
  
  savage_plot_list_H0[[i]] <- p
}

# Plotting for Savage test under Alternatives
savage_plot_list_H1_theta2 <- list()

for (i in seq_along(savage_results)) {
  res <- savage_results[[i]]
  df <- data.frame(S = res$H1_theta2)
  
  p <- ggplot(df, aes(x = S)) +
    geom_histogram(
      aes(y = ..density..),
      bins = 30,
      fill = "#69b3a2",
      color = "white",
      alpha = 0.8,
      linewidth = 0.3
    ) +
    geom_density(color = "#404080", linewidth = 1.2) +
    labs(
      title = paste("Savage Test:", res$distr, "(n =", res$n, ", m =", res$m, ")"),
      x = "Standardized Savage Statistic",
      y = "Density"
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
  
  savage_plot_list_H1_theta2[[i]] <- p
}

savage_plot_list_H1_theta6 <- list()

for (i in seq_along(savage_results)) {
  res <- savage_results[[i]]
  df <- data.frame(S = res$H1_theta6)
  
  p <- ggplot(df, aes(x = S)) +
    geom_histogram(
      aes(y = ..density..),
      bins = 30,
      fill = "#69b3a2",
      color = "white",
      alpha = 0.8,
      linewidth = 0.3
    ) +
    geom_density(color = "#404080", linewidth = 1.2) +
    labs(
      title = paste("Savage Test:", res$distr, "(n =", res$n, ", m =", res$m, ")"),
      x = "Standardized Savage Statistic",
      y = "Density"
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
  
  savage_plot_list_H1_theta6[[i]] <- p
}

# Q-Q Plots for Savage test
savage_qq_plots <- list()

for (i in seq_along(savage_results)) {
  res <- savage_results[[i]]
  qq_data <- data.frame(S = res$H0)
  
  qq_plot <- ggplot(qq_data, aes(sample = S)) +
    stat_qq(distribution = qnorm, size = 1) +
    stat_qq_line(distribution = qnorm, color = "red", linewidth = 1) +
    labs(
      title = paste("Savage Test:", res$distr, "(n =", res$n, ", m =", res$m, ")"),
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 8))
  
  savage_qq_plots[[i]] <- qq_plot
}

#Power Curve

library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(gridExtra)


# --- Simulation Function ---
simulate_test <- function(distr, n, m, theta, test_func, B = 1000) {
  stats <- numeric(B)
  for (i in 1:B) {
    if (distr == "Normal") {
      x <- rnorm(n, 0, 1); y <- rnorm(m, 0, theta)
    } else if (distr == "Exponential") {
      x <- rexp(n, 1); y <- rexp(m, 1/theta)
    } else if (distr == "Cauchy") {
      x <- rcauchy(n, 0, 1); y <- rcauchy(m, 0, theta)
    } else if (distr == "Geometric") {
      x <- rgeom(n, 0.5); y <- rgeom(m, 1/(1 + theta))
    }
    stats[i] <- test_func(x, y)
  }
  stats
}

# --- Power Calculation ---
compute_power <- function(distr, n, m, theta_range, B = 1000, alpha = 0.05) {
  tests <- list(
    Klotz = klotz_test,
    Savage = savage_test,
    Mood = mood_test
  )
  
  map_dfr(names(tests), function(test_name) {
    # Get critical value under H0
    H0_stats <- simulate_test(distr, n, m, 1, tests[[test_name]], B)
    crit_val <- quantile(H0_stats, alpha)
    
    # Compute power for each theta
    map_dfr(theta_range, function(theta) {
      H1_stats <- simulate_test(distr, n, m, theta, tests[[test_name]], B)
      data.frame(
        Test = test_name,
        theta = theta,
        Power = mean(H1_stats < crit_val)
      )
    })
  })
}

# --- Main Analysis ---
distributions <- c("Normal", "Exponential", "Cauchy", "Geometric")
sample_sizes <- list(c(10, 12), c(70, 65), c(240, 250))
theta_range <- seq(1, 1.5, length.out = 20)

set.seed(123)
power_results <- list()

for (distr in distributions) {
  for (sizes in sample_sizes) {
    n <- sizes[1]; m <- sizes[2]
    key <- paste(distr, "n=", n, "m=", m)
    
    power_results[[key]] <- compute_power(
      distr = distr,
      n = n,
      m = m,
      theta_range = theta_range,
      B = 1000
    )
  }
}

# --- Plotting Functions ---
create_power_plot <- function(results, title) {
  ggplot(results, aes(x = theta, y = Power, color = Test)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      title = title,
      x = "Scale Parameter (θ)",
      y = "Test Power",
      color = "Test"
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Create and store all power curve plots
power_plots <- map(names(power_results), function(scenario) {
  create_power_plot(
    power_results[[scenario]],
    title = paste("Power Curve:", scenario)
  )
})
names(power_plots) <- names(power_results)








output_dir <- "C:/Users/SAMSUNG/OneDrive/Desktop/Projects/Location Problem" 


# Function to save plots with consistent naming convention
save_test_plots <- function(plot_list, test_name, condition, results_list) {
  # Create the appropriate subdirectory
  subdir <- file.path(output_dir, paste0(test_name, "_", condition))
  if (!dir.exists(subdir)) {
    dir.create(subdir, recursive = TRUE)
  }
  
  # Save each plot with the correct naming
  for (i in seq_along(plot_list)) {
    # Get distribution and sample size info from results
    distr <- results_list[[i]]$distr
    n <- results_list[[i]]$n
    m <- results_list[[i]]$m
    
    # Determine the prefix based on distribution
    prefix <- switch(distr,
                     "Normal" = "N",
                     "Exponential" = "E",
                     "Cauchy" = "C",
                     "Geometric" = "G")
    
    # Determine the suffix based on sample size
    suffix <- case_when(
      n == 10 & m == 12 ~ "1",
      n == 70 & m == 65 ~ "2",
      n == 240 & m == 250 ~ "3",
      TRUE ~ as.character(i)
    )
    
    # Create the filename
    filename <- paste0(prefix, suffix, ".png")
    
    # Save the plot
    ggsave(
      filename = file.path(subdir, filename),
      plot = plot_list[[i]],
      width = 8,
      height = 6,
      dpi = 300,
      bg = "white"
    )
  }
}

# ==============================================
# Save Klotz Test Plots
# ==============================================

# H0 Plots
save_test_plots(klotz_plot_list_H0, "Klotz", "H0", klotz_results)

# H1 theta=2 Plots
save_test_plots(klotz_plot_list_H1_theta2, "Klotz", "H1", klotz_results)

# H1 theta=6 Plots
save_test_plots(klotz_plot_list_H1_theta6, "Klotz", "H2", klotz_results)

# Q-Q Plots
save_test_plots(klotz_qq_plots, "Klotz", "QQ", klotz_results)

# ==============================================
# Save Mood Test Plots
# ==============================================

# H0 Plots
save_test_plots(mood_plot_list_H0, "Mood", "H0", mood_results)

# H1 theta=2 Plots
save_test_plots(mood_plot_list_H1_theta2, "Mood", "H1", mood_results)

# H1 theta=6 Plots
save_test_plots(mood_plot_list_H1_theta6, "Mood", "H2", mood_results)

# Q-Q Plots
save_test_plots(mood_qq_plots, "Mood", "QQ", mood_results)

# ==============================================
# Save Savage Test Plots
# ==============================================

# H0 Plots
save_test_plots(savage_plot_list_H0, "Savage", "H0", savage_results)

# H1 theta=2 Plots
save_test_plots(savage_plot_list_H1_theta2, "Savage", "H1", savage_results)

# H1 theta=6 Plots
save_test_plots(savage_plot_list_H1_theta6, "Savage", "H2", savage_results)

# Q-Q Plots
save_test_plots(savage_qq_plots, "Savage", "QQ", savage_results)


#Innovation

# Modified Discrete Rank Test for Geometric Distribution
mdrt_test <- function(x, y) {
  n <- length(x); m <- length(y); N <- n + m
  
  # Variance-stabilizing transformation
  trans_x <- sqrt(x + 3/8)
  trans_y <- sqrt(y + 3/8)
  
  combined <- c(trans_x, trans_y)
  ranks <- rank(combined)
  
  # Discrete probability scores
  p_scores <- (ranks - 0.5)/N
  
  # Expected values under H0
  E0 <- n * mean(p_scores)
  
  # Variance under H0 with continuity correction
  Var0 <- (n*m)/(12*N) * (1 + 1/(N+1))
  
  # Test statistic
  S <- sum(p_scores[(1):n])
  (S - E0)/sqrt(Var0)
}

library(ggplot2)
library(dplyr)
library(purrr)

# Simulate Null Distribution
# Simulation for MDRT test (Geometric only)
simulate_mdrt_geom <- function(n, m, theta, B = 10000) {
  M_stats <- numeric(B)
  
  for (i in 1:B) {
    x <- rgeom(n, prob = 0.5)
    y <- rgeom(m, prob = 1/ (1 + theta))
    M_stats[i] <- mdrt_test(x, y)
  }
  
  return(M_stats)
}

# Run MDRT simulations (Geometric only)
set.seed(123)
mdrt_geom_results <- list()

for (sizes in sample_sizes) {
  n <- sizes[1]
  m <- sizes[2]
  key <- paste("Geometric", "n=", n, "m=", m)
  
  sim_H0 <- simulate_mdrt_geom(n, m, theta = 1, B = 10000)
  sim_H1_theta2 <- simulate_mdrt_geom(n, m, theta = 2, B = 10000)
  sim_H1_theta6 <- simulate_mdrt_geom(n, m, theta = 6, B = 10000)
  
  mdrt_geom_results[[key]] <- list(
    H0 = sim_H0,
    H1_theta2 = sim_H1_theta2,
    H1_theta6 = sim_H1_theta6,
    distr = "Geometric",
    n = n,
    m = m
  )
}

# Plotting for MDRT test under Null Distribution (Geometric only)
mdrt_geom_plot_list_H0 <- list()

for (i in seq_along(mdrt_geom_results)) {
  res <- mdrt_geom_results[[i]]
  df <- data.frame(M = res$H0)
  
  p <- ggplot(df, aes(x = M)) +
    geom_histogram(
      aes(y = ..density..),
      bins = 30,
      fill = "#69b3a2",
      color = "white",
      alpha = 0.8,
      linewidth = 0.3
    ) +
    geom_density(color = "#404080", linewidth = 1.2) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1), 
                  color = "red", linetype = 2, linewidth = 1) +
    labs(
      title = paste("MDRT Test:", res$distr, "(n =", res$n, ", m =", res$m, ")"),
      x = "Standardized MDRT Statistic",
      y = "Density"
    ) +
    theme_minimal(base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
  
  mdrt_geom_plot_list_H0[[i]] <- p
}

# Power Curve Comparison

# Simulation function (Geometric only)
simulate_geom_test <- function(n, m, theta, test_func, B = 1000) {
  stats <- numeric(B)
  for (i in 1:B) {
    x <- rgeom(n, 0.5)
    y <- rgeom(m, 1/(1 + theta))
    stats[i] <- test_func(x, y)
  }
  stats
}

# Power calculation (Geometric only)
compute_geom_power <- function(n, m, theta_range, B = 1000, alpha = 0.05) {
  tests <- list(
    MDRT = mdrt_test,
    Klotz = klotz_test,
    Savage = savage_test,
    Mood = mood_test
  )
  
  map_dfr(names(tests), function(test_name) {
    # Get critical value under H0
    H0_stats <- simulate_geom_test(n, m, 1, tests[[test_name]], B)
    crit_val <- quantile(H0_stats, alpha)
    
    # Compute power for each theta
    map_dfr(theta_range, function(theta) {
      H1_stats <- simulate_geom_test(n, m, theta, tests[[test_name]], B)
      data.frame(
        Test = test_name,
        theta = theta,
        Power = mean(H1_stats < crit_val))
    })
  })
}

# Main analysis (Geometric only)
sample_sizes <- list(c(10, 12), c(70, 65), c(240, 250))
theta_range <- seq(1, 1.5, length.out = 20)

set.seed(123)
geom_power_results <- list()

for (sizes in sample_sizes) {
  n <- sizes[1]; m <- sizes[2]
  key <- paste("Geometric", "n=", n, "m=", m)
  
  geom_power_results[[key]] <- compute_geom_power(
    n = n,
    m = m,
    theta_range = theta_range,
    B = 1000
  )
}

# Plotting function 

create_power_plot <- function(results, title) {
  ggplot(results, aes(x = theta, y = Power, color = Test)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      title = title,
      x = "Scale Parameter (θ)",
      y = "Test Power",
      color = "Test"
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Create and store all power curve plots
geom_power_plots <- map(names(geom_power_results), function(scenario) {
  create_power_plot(
    geom_power_results[[scenario]],
    title = paste("Power Curve:", scenario)
  )
})
names(geom_power_plots) <- names(geom_power_results)



# Size Analysis (Type I Error)

# Simulation function for asymptotic size (Geometric only)
simulate_mdrt_size <- function(n, m, B = 10000, alpha = 0.05) {
  rej <- 0
  for (b in 1:B) {
    x <- rgeom(n, prob = 0.5)
    y <- rgeom(m, prob = 0.5)
    stat <- mdrt_test(x, y)
    if (abs(stat) > qnorm(1 - alpha/2)) {
      rej <- rej + 1
    }
  }
  return(rej/B)
}

# Simulation parameters (Geometric only)
sample_sizes <- list(
  small = c(20, 20),
  medium = c(50, 50), 
  large = c(100, 100),
  asymp = c(200, 200))
alpha <- 0.05

# Run simulations
set.seed(123)
results <- data.frame()

for (size_name in names(sample_sizes)) {
  n <- sample_sizes[[size_name]][1]
  m <- sample_sizes[[size_name]][2]
  
  rej_rate <- simulate_mdrt_size(n, m, B = 10000, alpha = alpha)
  
  results <- rbind(results, data.frame(
    Distribution = "Geometric",
    n = n,
    m = m,
    Size = size_name,
    RejectionRate = rej_rate,
    Alpha = alpha
  ))
}

# Print results table
library(knitr)
kable(results, digits = 4, caption = "Simulation Results for MDRT (Geometric Only)")




