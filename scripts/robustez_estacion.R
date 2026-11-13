# prueba 7, estacionalidad: la b del conjunto con dos armónicos del mes en la tasa.
#
# Regla en docs/pre_specified_rules.md (13-nov-2026, noche), escrita antes de correr:
#   λ(x, m) = exp(a_g + b·x + Σ_{k=1,2} [s_k sin(2πkm/12) + c_k cos(2πkm/12)]); a por hidrómetro (cerrada),
#   b y armónicos comunes a los hidrómetros con 15 eventos o más; bootstrap por años hidrológicos (500, semilla 2026,
#   5 eventos o más por réplica). Pasa si el IC90 de b excluye el 0; se informa si cae en el IC90 principal (0,30–0,70).
# No toca modelo.R ni las otras tablas. b_conj, con_eventos y boot_anios repiten los de robustez_indice.R; el control
# sin armónicos tiene que reproducir la b principal (0,51702835).
# Uso (desde la raíz del repositorio):  Rscript scripts/robustez_estacion.R [--boot 500]

args <- commandArgs(trailingOnly = TRUE)
N_BOOT <- if ("--boot" %in% args) as.integer(args[which(args == "--boot") + 1]) else 500L
source("scripts/modelo.R", encoding = "UTF-8")
set.seed(2026)
M0 <- 3
q5 <- function(v) unname(quantile(v, 0.05)); q95 <- function(v) unname(quantile(v, 0.95))
B_PRINCIPAL <- 0.51702835346803; IC_PRINCIPAL <- c(0.30, 0.70)

mes_num <- function(f) as.integer(format(f, "%m"))
armonicos <- function(f) { m <- mes_num(f); cbind(sin(2 * pi * m / 12), cos(2 * pi * m / 12),
                                                   sin(4 * pi * m / 12), cos(4 * pi * m / 12)) }

# −log verosimilitud perfilada en a_g (cerrada): para cada hidrómetro, exp(a_g) = n_g / Σ_meses exp(b·x + h)/12
pnll_est <- function(par, Ds, con_h = TRUE) {
  b <- par[1]; g <- if (con_h) par[2:5] else rep(0, 4)
  sum(vapply(Ds, function(D) {
    he <- if (length(D$xe)) as.vector(armonicos(D$ye) %*% g) else numeric(0)
    hm <- as.vector(armonicos(D$mm) %*% g)
    a <- log(max(length(D$xe), 0.5)) - log(sum(exp(b * D$xm + hm)) / 12)
    -(sum(a + b * D$xe + he) - sum(exp(a + b * D$xm + hm)) / 12)
  }, 0))
}
ajuste_est <- function(Ds) {
  f0 <- optimize(function(b) pnll_est(c(b, 0, 0, 0, 0), Ds, FALSE), c(-6, 6), tol = 1e-10)
  f <- optim(c(f0$minimum, 0, 0, 0, 0), pnll_est, Ds = Ds, method = "BFGS", control = list(reltol = 1e-12, maxit = 1000))
  f <- optim(f$par, pnll_est, Ds = Ds, method = "Nelder-Mead", control = list(reltol = 1e-12, maxit = 5000))
  list(b0 = f0$minimum, nll0 = f0$objective, par = f$par, nll = f$value, conv = f$convergence)
}

con_eventos <- function(Ds, k) Ds[vapply(Ds, function(D) length(D$xe) >= k, TRUE)]
boot_anios <- function(Ds, estad, n = N_BOOT) {
  ye <- lapply(Ds, function(D) anio_h(D$ye, M0)); ym <- lapply(Ds, function(D) anio_h(D$mm, M0))
  anios <- sort(unique(unlist(ym)))
  idx <- lapply(seq_along(Ds), function(i) list(e = split(seq_along(ye[[i]]), ye[[i]]),
                                                 m = split(seq_along(ym[[i]]), ym[[i]])))
  out <- vector("list", n)
  for (r in seq_len(n)) {
    smp <- sample(anios, length(anios), replace = TRUE)
    Db <- lapply(seq_along(Ds), function(i) {
      s <- remuestrear(idx[[i]]$e, idx[[i]]$m, smp)
      list(xe = Ds[[i]]$xe[s$ie], xm = Ds[[i]]$xm[s$im], ye = Ds[[i]]$ye[s$ie], mm = Ds[[i]]$mm[s$im])
    })
    out[[r]] <- estad(Db)
  }
  unlist(out)
}

Ds <- con_eventos(lapply(URUGUAY, function(s) datos(s, "x6")), 15)
f <- ajuste_est(Ds)
cat(sprintf("control sin armónicos: b = %.8f (principal %.8f)\n", f$b0, B_PRINCIPAL))
if (abs(f$b0 - B_PRINCIPAL) > 1e-6) { cat("el control no reproduce la b principal: se para acá\n"); quit(status = 1) }
lr <- 2 * (f$nll0 - f$nll)
bt <- boot_anios(Ds, function(Db) ajuste_est(con_eventos(Db, 5))$par[1])
r <- data.frame(caso = "with annual cycle (two harmonics)", hidrometros = length(Ds),
                n_eventos = sum(vapply(Ds, function(D) length(D$xe), 0)), b = f$par[1], b_boot_lo = q5(bt),
                b_boot_hi = q95(bt), b_sin_estacion = f$b0, s1 = f$par[2], c1 = f$par[3], s2 = f$par[4], c2 = f$par[5],
                lr_armonicos = lr, p_armonicos = pchisq(lr, 4, lower.tail = FALSE), convergencia = f$conv,
                pasa = q5(bt) > 0, dentro_ic_principal = f$par[1] >= IC_PRINCIPAL[1] & f$par[1] <= IC_PRINCIPAL[2])
write.csv(r, "tables/robustez_estacion.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("prueba 7: b = %.3f [%.3f, %.3f] con armónicos (sin: %.3f); LR armónicos = %.1f (p = %.2g) | %s | %s\n",
            r$b, r$b_boot_lo, r$b_boot_hi, r$b_sin_estacion, lr, r$p_armonicos,
            if (r$pasa) "PASA" else "NO PASA", if (r$dentro_ic_principal) "dentro del IC90 principal" else "FUERA del IC90 principal"))
