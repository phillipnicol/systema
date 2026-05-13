set.seed(1)

n_0 <- 5000
n_1 <- 5000
d <- 5000

alphas <- seq(0, 1, by = 0.025)
n_rep <- 10

cosine_similarity <- function(x, y) {
  sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
}

actual_mean <- numeric(length(alphas))
actual_sd <- numeric(length(alphas))
theoretical <- numeric(length(alphas))

c <- d / n_0
r1 <- 1
r2 <- 1

for (i in seq_along(alphas)) {
  alpha <- alphas[i]
  
  sims <- numeric(n_rep)
  
  for (rep in 1:n_rep) {
    C <- rnorm(d, sd = 1 / sqrt(n_0))
    
    P <- rnorm(d, sd = 1 / sqrt(n_1))
    P[1] <- 1 + P[1]
    
    mu_tilde <- rep(0, d)
    mu_tilde[1] <- alpha
    mu_tilde[2] <- sqrt(1 - alpha^2)
    
    sims[rep] <- cosine_similarity(P - C, mu_tilde - C)
  }
  
  actual_mean[i] <- mean(sims)
  actual_sd[i] <- sd(sims)
  
  kappa <- alpha
  theoretical[i] <- (kappa + c) / sqrt((r1^2 + 2 * c) * (r2^2 + c))
}


library(ggplot2)

# Put everything in a data frame
df <- data.frame(
  alpha = alphas,
  actual_mean = actual_mean,
  actual_sd = actual_sd,
  theoretical = theoretical
)

p <- ggplot(df, aes(x = alpha)) +
  
  # SD ribbon
  geom_ribbon(
    aes(ymin = actual_mean - actual_sd,
        ymax = actual_mean + actual_sd),
    fill = "blue",
    alpha = 0.2
  ) +
  
  # Empirical mean
  geom_line(
    aes(y = actual_mean, color = "Empirical"),
    linewidth = 1.2
  ) +
  
  # Theoretical
  geom_line(
    aes(y = theoretical, color = "Theoretical"),
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  
  # y = x line
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "gray40",
    linetype = "dotted",
    linewidth = 1
  ) +
  
  scale_color_manual(
    values = c("Empirical" = "blue", "Theoretical" = "red")
  ) +
  
  labs(
    x = "True cosine similarity",
    y = "Estimated cosine similarity (shared control)",
    color = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top"
  )


p

#Save to fig

ggsave(p, filename="../fig/theory_cosine.png",
       width=4.61, height=4.52)

