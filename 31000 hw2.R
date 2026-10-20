library(tidyverse)

uiaustria <- grossman::load("uiaustria")

df <- uiaustria %>%
  filter(age < 40) %>%
  mutate(
    Ti = ceiling(days / 7),
    Ci = ceiling(window / 7),
    R    = !is.na(Ti) & Ti <= Ci,
    Yi    = if_else(R, Ti, Ci)
  )

tmax <- 260

#Sn
Sn_df <- tibble(t = 1:tmax) %>%
  mutate(Sn = map_dbl(t, ~ mean(df$Yi >= .x)))

#Kn
Hn_df <- tibble(r = 1:(tmax - 1)) %>%
  mutate(
    num   = map_dbl(r, ~ sum(df$Yi == .x & df$R)),  # 分子
    den = map_dbl(r, ~ sum(df$Yi >= .x)),         # 分母
    Hn      = num / den
  )

Kn_df <- tibble(
  t  = 1:tmax,
  Kn = c(1, cumprod(1 - Hn_df$Hn))
)

est <- Sn_df %>% left_join(Kn_df, by = "t")

#Bm
Bn_df <- tibble(
  t  = 1:tmax,
  Bn = c(1, exp(-cumsum(Hn_df$Hn)))   # B_n(1) = exp(0) = 1
)

est <- est %>% left_join(Bn_df, by = "t")


ggplot(est, aes(x = t)) +
  geom_step(aes(y = Sn, color = "S_n")) +
  geom_step(aes(y = Kn, color = "K_n")) +
  geom_step(aes(y = Bn, color = "B_n")) +
  labs(x = "Week t", y = "Estimated survival", color = NULL)

ggplot(est, aes(x = t)) +
  geom_step(aes(y = Bn - Sn, color = "B_n - S_n")) +
  geom_step(aes(y = Kn - Sn, color = "K_n - S_n")) +
  geom_hline(yintercept = 0) +
  labs(
    x = "Week t",
    y = "Difference from S_n",
    color = NULL
  ) +
  theme(
    legend.position = "inside",
    legend.justification = c("right", "top")
  )

#h
pool <- uiaustria %>%
  mutate(
    Ti = ceiling(days / 7),
    Ci = ceiling(window / 7)
  ) %>%
  filter(Ci >= 260, !is.na(Ti), Ti <= 260)

length(pool$Ti)
summary(pool$Ti)

s_true <- tibble(t = 1:260) %>%
  mutate(s = map_dbl(t, ~ mean(pool$Ti >= .x)))

s_true %>% filter(t %in% c(1, 13, 52, 104, 208, 260))

tmax  <- 260
z     <- qnorm(0.975)
s_vec <- s_true$s

sim_once <- function(n) {
  Ti <- sample(pool$Ti, n, replace = TRUE)
  Ci <- sample(1:tmax, n, replace = TRUE)
  R  <- Ti <= Ci
  Yi <- pmin(Ti, Ci)
  
  at_risk <- rev(cumsum(rev(tabulate(Yi, nbins = tmax))))   #{Yi >= r}
  exits   <- tabulate(Yi[R], nbins = tmax)                  #{Yi = r, Ri = 1}
  Hn      <- ifelse(at_risk > 0, exits / at_risk, 0)
  

  Sn <- at_risk / n
  Kn <- c(1, cumprod(1 - Hn))[1:tmax]
  Bn <- c(1, exp(-cumsum(Hn)))[1:tmax]
  
  se_S   <- sqrt(Sn * (1 - Sn) / (n - 1))
  term_K <- ifelse(at_risk > 0, Hn / ((1 - Hn) * at_risk), 0)
  term_B <- ifelse(at_risk > 0, Hn * (1 - Hn) / at_risk, 0)
  se_K   <- Kn * sqrt(c(0, cumsum(term_K)))[1:tmax]
  se_B   <- Bn * sqrt(c(0, cumsum(term_B)))[1:tmax]
  

  covers <- function(est, se) coalesce(abs(est - s_vec) <= z * se, FALSE)
  
  tibble(
    t = 1:tmax, s = s_vec, at_risk = at_risk,
    Sn = Sn, Kn = Kn, Bn = Bn,
    cover_S = covers(Sn, se_S),
    cover_K = covers(Kn, se_K),
    cover_B = covers(Bn, se_B),
    Kn_zero = (Kn == 0)
  )
}

set.seed(13)
n_reps <- 1e4
ns <- c(50, 250, 1000, 4000)

sims <- lapply(ns, function(n) {
  replicate(n_reps, sim_once(n), simplify = FALSE) %>%
    bind_rows() %>%
    mutate(n = n)
}) %>%
  bind_rows()

dim(sims)

summ <- sims %>%
  group_by(n, t) %>%
  summarize(
    cover_S       = mean(cover_S),
    cover_K       = mean(cover_K),
    cover_B       = mean(cover_B),
    avg_at_risk   = mean(at_risk),
    share_Kn_zero = mean(Kn_zero),
    .groups = "drop"
  )

p_cov <- ggplot(summ, aes(x = t)) +
  geom_line(aes(y = cover_S, color = "S_n")) +
  geom_line(aes(y = cover_K, color = "K_n")) +
  geom_line(aes(y = cover_B, color = "B_n")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  facet_wrap(~ n, labeller = label_both) +
  scale_color_discrete(breaks = c("S_n", "K_n", "B_n")) +
  labs(
    x = "Week t",
    y = "Coverage of 95% CI for s(t)",
    color = NULL,
    title = "Coverage of nominal 95% confidence intervals"
  )

print(p_cov)

report <- summ %>%
  filter(t %in% c(13, 52, 104, 208)) %>%
  arrange(t, n)

print(report, n = Inf)


#i

sim_once <- function(n) {
  Ti <- sample(pool$Ti, n, replace = TRUE)
  Ci <- sample(1:tmax, n, replace = TRUE)
  R  <- Ti <= Ci
  Yi <- pmin(Ti, Ci)
  
  at_risk <- rev(cumsum(rev(tabulate(Yi, nbins = tmax))))   # #{Yi >= r}
  exits   <- tabulate(Yi[R], nbins = tmax)                  # #{Yi = r, Ri = 1}
  Hn      <- ifelse(at_risk > 0, exits / at_risk, 0)
  
  Sn <- at_risk / n
  Kn <- c(1, cumprod(1 - Hn))[1:tmax]
  Bn <- c(1, exp(-cumsum(Hn)))[1:tmax]
  
  se_S   <- sqrt(Sn * (1 - Sn) / (n - 1))
  term_K <- ifelse(at_risk > 0, Hn / ((1 - Hn) * at_risk), 0)
  term_B <- ifelse(at_risk > 0, Hn * (1 - Hn) / at_risk, 0)
  se_K   <- Kn * sqrt(c(0, cumsum(term_K)))[1:tmax]
  se_B   <- Bn * sqrt(c(0, cumsum(term_B)))[1:tmax]
  
  #new
  Gn     <- sqrt(c(0, cumsum(term_K)))[1:tmax]
  lo_log <- Kn * exp(-z * Gn)
  hi_log <- Kn * exp( z * Gn)
  
  covers <- function(est, se) coalesce(abs(est - s_vec) <= z * se, FALSE)
  
  tibble(
    t = 1:tmax, s = s_vec, at_risk = at_risk,
    Sn = Sn, Kn = Kn, Bn = Bn,
    cover_S = covers(Sn, se_S),
    cover_K = covers(Kn, se_K),
    cover_B = covers(Bn, se_B),
    cover_L = coalesce(lo_log <= s_vec & s_vec <= hi_log, FALSE),  # new in (i)
    Kn_zero = (Kn == 0)
  )
}

set.seed(13)
n_reps <- 1e4   # same as in (h)
ns <- c(50, 250, 1000, 4000)

sims <- lapply(ns, function(n) {
  replicate(n_reps, sim_once(n), simplify = FALSE) %>%
    bind_rows() %>%
    mutate(n = n)
}) %>%
  bind_rows()


summ_i <- sims %>%
  group_by(n, t) %>%
  summarize(
    cover_K       = mean(cover_K),
    cover_L       = mean(cover_L),
    avg_at_risk   = mean(at_risk),
    share_Kn_zero = mean(Kn_zero),
    .groups = "drop"
  )

p_cov_i <- ggplot(summ_i, aes(x = t)) +
  geom_line(aes(y = cover_K, color = "K_n (symmetric)")) +
  geom_line(aes(y = cover_L, color = "K_n (log-transformed)")) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  facet_wrap(~ n, labeller = label_both) +
  labs(
    x = "Week t",
    y = "Coverage of 95% CI for s(t)",
    color = NULL,
    title = "Symmetric vs. log-transformed Kaplan-Meier intervals"
  )

print(p_cov_i)


report_i <- summ_i %>%
  filter(t %in% c(13, 52, 104, 208)) %>%
  arrange(t, n)

print(report_i, n = Inf)
