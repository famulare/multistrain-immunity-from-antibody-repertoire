# twisted copula

# Smooth tail-flip generator with optional two-tail extension
# See inline comments for usage.

smooth_tail_flip <- function(
    x,
    rho_m = 0.9,
    rho_t_hi = -0.9,     # upper-tail correlation target
    rho_t_lo = -0.9,     # lower-tail correlation target; if NULL, uses same as rho_m
    alpha = 0.01,         # tail mass for each tail
    h = 0.01,             # smoothness in z-space (larger = smoother)
    seed = NULL
) {
  stopifnot(is.numeric(x), length(x) >= 5)
  if (!is.null(seed)) set.seed(seed)
  
  n <- length(x)
  eps <- 1e-6
  clamp <- function(v, lo = -0.999, hi = 0.999) pmax(lo, pmin(hi, v))
  
  # empirical marginal
  rank_u <- (rank(x, ties.method = "average")) / (n + 1)
  qF <- function(p) as.numeric(stats::quantile(x, probs = p, type = 8))
  z  <- stats::qnorm(rank_u)
  
  # Cutoffs
  z_hi <- stats::qnorm(1 - alpha)
  z_lo <- stats::qnorm(alpha)
  if (is.null(rho_t_lo)) rho_t_lo <- rho_m
  
  # Smooth weights
  w_hi <- stats::plogis((z - z_hi) / max(h, 1e-6))   # upper tail
  w_lo <- stats::plogis(-(z - z_lo) / max(h, 1e-6))  # lower tail
  
  # Combined tail weight (0 in middle, ~1 in each tail)
  w_total <- pmax(w_hi, w_lo)
  
  # Which target correlation applies
  rho_target <- (w_hi * rho_t_hi + w_lo * rho_t_lo) / pmax(w_total, 1e-8)
  rho_target[!is.finite(rho_target)] <- rho_m
  
  # Generate y given bulk correlation rho_bulk
  generate_y_from_rho_bulk <- function(rho_bulk) {
    e  <- stats::rnorm(n)
    rho_i <- clamp((1 - w_total) * rho_bulk + w_total * rho_target)
    z_y <- rho_i * z + sqrt(pmax(1 - rho_i^2, eps)) * e
    u_y <- stats::pnorm(z_y)
    qF(u_y)
  }
  
  return(generate_y_from_rho_bulk(rho_m))
}


# Plot the quantile-dependent correlation blending function
plot_rho_profile <- function(
    rho_m = 0.4,
    rho_t_hi = -0.9,
    rho_t_lo = NULL,
    alpha = 0.1,
    h = 0.5
) {
  if (is.null(rho_t_lo)) rho_t_lo <- rho_m
  z <- seq(-4, 4, length.out = 500)
  z_hi <- qnorm(1 - alpha)
  z_lo <- qnorm(alpha)
  w_hi <- plogis((z - z_hi) / h)
  w_lo <- plogis(-(z - z_lo) / h)
  w_total <- pmax(w_hi, w_lo)
  rho_target <- (w_hi * rho_t_hi + w_lo * rho_t_lo) / pmax(w_total, 1e-8)
  rho_target[!is.finite(rho_target)] <- rho_m
  rho_i <- (1 - w_total) * rho_m + w_total * rho_target
  q <- pnorm(z)
  plot(q, rho_i, type = "l", lwd = 2, col = "steelblue",
       xlab = "Quantile (of x)", ylab = expression(rho(q)),
       main = "Quantile-dependent correlation profile")
  abline(h = rho_m, col = "gray70", lty = 2)
  abline(v = c(alpha, 1 - alpha), col = "gray80", lty = 3)
  if (!is.null(rho_t_lo))
    legend("topright", legend = c("rho(q)", "rho_m"), col = c("steelblue","gray70"), lty = c(1,2), bty = "n")
}

set.seed(42)
n <- 5000
x <- rexp(n)
x = x/max(x)

# x=runif(n)

# Generate with two-tail flip: both upper and lower tails anticorrelated
y <- smooth_tail_flip(
  x,
  rho_m = 0.9,
  rho_t_hi = -0.9,
  rho_t_lo = -0.9,
  alpha = 0.01,
  h = 0.01
)

# check correlations
cor(x, y)  # overall
plot(x,y)
qx <- quantile(x, 0.9); qy <- quantile(y, 0.9)
qlx <- quantile(x, 0.1); qly <- quantile(y, 0.1)
c(
  overall = cor(x, y),
  center = cor(x[x>qlx & x<qx & y>qly & y<qy], y[x>qlx & x<qx & y>qly & y<qy]),
  upper_tail = cor(x[x>qx & y<qly], y[x>qx & y<qly]),
  lower_tail = cor(x[x<qlx & y>qy], y[x<qlx & y>qy])
)

# visualize the correlation profile
plot_rho_profile(rho_m = 0.4, rho_t_hi = -0.9, rho_t_lo = -0.9, alpha = 0.1, h = 0.1)

