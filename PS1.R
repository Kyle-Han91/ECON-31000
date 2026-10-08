#library if not
library(tidyverse)
library(dplyr)

uiaustria <- grossman::load("uiaustria")

df <- uiaustria %>%
  filter(age < 40, window >= 730) %>%
  mutate(Ti = if_else(is.na(days) | days > 728, 105, pmax(ceiling(days / 7), 1)))


#b
hz <- tibble(t = 1:105) %>%
  mutate(
    N_t = map_dbl(t, ~ sum(df$Ti >= .x)),
    D_t = map_dbl(t, ~ sum(df$Ti == .x)),
    H   = D_t / N_t,
    se  = sqrt(H * (1 - H) / N_t),
    lo  = H - qnorm(0.975) * se,
    hi  = H + qnorm(0.975) * se
  )

hz %>%
  filter(t <= 52) %>%
  ggplot(aes(x = t, y = H)) +
  geom_line(alpha = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.6, color = "steelblue") +
  geom_point(size = 0.6) +
  geom_vline(xintercept = 30, linetype = "dashed", color = "red")+
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  labs(x = "Weeks since layoff (t)", y = "Estimated hazard H(t)",
       title = "Weekly re-employment hazard with 95% CIs") +
  theme_minimal()

#d
test_equal <- function(a,b) {
  A <- hz %>% filter (t==a)
  B <- hz %>% filter (t==b)
  diff <- A$H - B$H
  se_diff <- sqrt(A$se^2 + B$se^2)
  z <- diff / se_diff
  tibble(
    null    = paste0("h(", a, ") = h(", b, ")"),
    H_a     = A$H,
    H_b     = B$H,
    diff    = diff,
    se_diff = se_diff,
    z       = z,
    p_value = 2 * (1- pnorm(abs(z)))
  )
}

bind_rows(
  test_equal(31, 30),
  test_equal(31, 32),
  test_equal(32, 33)
)

#e
wald_constant <- function(H, se) {
  k  <- length(H)
  R  <- cbind(0, diag(k - 1)) - cbind(diag(k - 1), 0)
  V  <- diag(se^2)
  RH <- R %*% H
  W  <- as.numeric(t(RH) %*% solve(R %*% V %*% t(R), RH))
  tibble(W = W, df = k - 1,
         p_value = pchisq(W, df = k - 1, lower.tail = FALSE))
}

hz_e <- hz %>% filter(t %in% 5:25)
wald_constant(hz_e$H, hz_e$se)

#f
test_avg <- function(H, se, t) {
  w    <- case_when(t %in% 5:15  ~  1/11,
                    t %in% 16:25 ~ -1/10,
                    TRUE         ~  0)
  B    <- sum(w * H)
  se_B <- sqrt(sum(w^2 * se^2))
  T0   <- B / se_B
  tibble(avg_early = mean(H[t %in% 5:15]),
         avg_late  = mean(H[t %in% 16:25]),
         B = B, se_B = se_B, T = T0,
         cv = qnorm(0.95),
         p_value = pnorm(T0, lower.tail = FALSE))
}

hz_e <- hz %>% filter(t %in% 5:25)
test_avg(hz_e$H, hz_e$se, hz_e$t) %>% print(width = Inf, digits = 6)

#g
test_max <- function(H, se, t, alpha = 0.05, nsim = 1e4) {
  k     <- length(H)
  m     <- k - 1
  Rm    <- cbind(0, diag(m)) - cbind(diag(m), 0)
  VD    <- Rm %*% diag(se^2) %*% t(Rm)             
  D     <- as.numeric(Rm %*% H) 
  Tj    <- D / sqrt(diag(VD))    
  Omega <- cov2cor(VD)
  
  Z     <- matrix(rnorm(nsim * m), nsim, m) %*% chol(Omega)
  maxZ  <- apply(Z, 1, max)
  
  tibble(Tmax    = max(Tj),
         week    = t[-1][which.max(Tj)],
         cv      = quantile(maxZ, 1 - alpha, names = FALSE),
         p_value = mean(maxZ >= max(Tj)))
}

set.seed(211)
hz_e <- hz %>% filter(t %in% 5:25)
test_max(hz_e$H, hz_e$se, hz_e$t)

#h

h_emp <- hz$H

h_null <- h_emp
h_null[5:25] <- mean(h_emp[5:25])

sim_durations <- function(n, h) {
  Ti <- rep(NA_integer_, n) 
  for (t in seq_along(h)) {
    exit     <- is.na(Ti) & (runif(n) < h[t])
    Ti[exit] <- t
  }
  Ti
}

hazard_stats <- function(Ti, weeks = 5:25) {
  D  <- tabulate(Ti, nbins = 105)
  N  <- rev(cumsum(rev(D)))
  H  <- D / N
  se <- sqrt(H * (1 - H) / N)
  list(H = H[weeks], se = se[weeks], weeks = weeks)
}


reject_all <- function(H, se, weeks, alpha = 0.10, nsim = 5000) {
  m  <- length(H) - 1
  Rm <- cbind(0, diag(m)) - cbind(diag(m), 0)
  VD <- Rm %*% diag(se^2) %*% t(Rm)
  D  <- as.numeric(Rm %*% H)
  
  # (e)
  rej_e <- sum(D * solve(VD, D)) > qchisq(1 - alpha, df = m)
  
  # (f)
  w     <- ifelse(weeks <= 15, 1/11, -1/10)
  rej_f <- sum(w * H) / sqrt(sum(w^2 * se^2)) > qnorm(1 - alpha)
  
  # (g)
  Z     <- matrix(rnorm(nsim * m), nsim, m) %*% chol(cov2cor(VD))
  cv    <- quantile(Z[cbind(1:nsim, max.col(Z))], 1 - alpha, names = FALSE)
  rej_g <- max(D / sqrt(diag(VD))) > cv
  
  c(e = rej_e, f = rej_f, g = rej_g)
}


run_mc <- function(n, draw_T, reps = 1000) {
  rej <- replicate(reps, {
    s <- hazard_stats(draw_T(n))
    if (any(s$se == 0)) return(c(e = NA, f = NA, g = NA)) 
    reject_all(s$H, s$se, s$weeks)
  })
  rowMeans(rej, na.rm = TRUE)
}

draw_null <- function(n) sim_durations(n, h_null)


draw_emp  <- function(n) sample(df$Ti, n, replace = TRUE)


set.seed(211)
results <- expand_grid(dgp = c("null imposed", "empirical"), n = c(1000, 4000)) %>%
  mutate(rate = map2(dgp, n, ~ run_mc(.y, if (.x == "null imposed") draw_null else draw_emp))) %>%
  unnest_wider(rate)

results %>% mutate(across(c(e, f, g), ~ round(.x, 3)))