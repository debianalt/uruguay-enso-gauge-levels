# modelo de frecuencia de crecidas: datos y funciones.
#
# Lo cargan scripts/frecuencia.R (ajuste principal y validación cruzada) y scripts/robustez.R. Es el código que
# estaba en frecuencia.R hasta el 10-nov-2026, movido sin cambiar una línea: tables/ajuste.csv y
# tables/prob_condicional.csv salen idénticos antes y después (el registro de decisiones del autor).
# Se carga desde la raíz del repositorio, porque lee data/ con ruta relativa.

URUGUAY <- c(60, 61, 63, 64, 65, 67, 68)
XI_MIN <- -0.5; XI_MAX <- 1; FUERA <- 1e10        # ξ acotada como en el Python (Coles 2001, § 2.6 y 4.3)
CHI1 <- qchisq(0.95, 1)
CLASES <- c("Niña (x = −1)" = -1, "neutro (x = 0)" = 0, "Niño (x = +1,5)" = 1.5)

leer <- function(t) read.csv(file.path("data", paste0("crecidas_", t, ".csv")), fileEncoding = "UTF-8",
                             stringsAsFactors = FALSE)
como_logico <- function(v) if (is.logical(v)) v else v %in% c("True", "TRUE", "true", "t")

# ---------- datos ------------------------------------------------------------------------------------
enso <- leer("enso"); enso$mes <- as.Date(enso$mes); enso <- enso[order(enso$mes), ]
media_movil <- function(v, k) as.numeric(stats::filter(v, rep(1 / k, k), sides = 1))  # = pandas rolling(k)
enso$x3 <- media_movil(enso$nino34, 3)
enso$x12 <- media_movil(enso$nino34, 12)
enso$x6_ersst <- media_movil(enso$ersstv6, 6)

serie <- leer("serie"); rownames(serie) <- serie$serie
ev <- leer("evento"); ev$fecha_pico <- as.Date(ev$fecha_pico); ev$mes_valido <- como_logico(ev$mes_valido)
mes <- leer("mes"); mes$mes <- as.Date(mes$mes); mes$valido <- como_logico(mes$valido)
ma <- leer("max_anual"); ma$fecha_max <- as.Date(ma$fecha_max); ma <- ma[como_logico(ma$valido), ]
homo <- leer("homogeneidad"); homo <- homo[!is.na(homo$marca) & homo$marca != "", ]

primer_dia <- function(f) as.Date(format(f, "%Y-%m-01"))
mes_siguiente <- function(f) seq(f, by = "month", length.out = 2)[2]

# quiebre de cada hidrómetro: el mes siguiente al último antes del salto (el más tardío si hay varios)
quiebre <- list()
for (i in seq_len(nrow(homo))) {
  t <- mes_siguiente(as.Date(paste0(homo$ultimo_antes[i], "-01")))
  for (s in strsplit(homo$objeto[i], "−", fixed = TRUE)[[1]]) {
    s <- trimws(s)
    quiebre[[s]] <- if (is.null(quiebre[[s]])) t else max(quiebre[[s]], t)
  }
}

datos <- function(sid, xcol, desde = NULL, hasta = NULL, extra = 0) {
  e <- ev[ev$serie == sid & ev$mes_valido, ]
  m <- mes[mes$serie == sid & mes$valido, ]
  if (!is.null(desde)) { e <- e[e$fecha_pico >= desde, ]; m <- m[m$mes >= desde, ] }
  if (!is.null(hasta)) { e <- e[e$fecha_pico < hasta, ]; m <- m[m$mes < hasta, ] }
  e$mes <- primer_dia(e$fecha_pico)
  e <- merge(e, enso[, c("mes", xcol)], by = "mes"); e <- e[!is.na(e[[xcol]]), ]
  m <- merge(m, enso[, c("mes", xcol)], by = "mes"); m <- m[!is.na(m[[xcol]]), ]
  list(xe = e[[xcol]], y = e$pico_m - e$u_m + extra, xm = m[[xcol]], ye = e$fecha_pico, mm = m$mes,
       u = if (nrow(e)) e$u_m[1] else NA_real_)
}

# año hidrológico: empieza en el mes m0 (marzo en el Uruguay, el mes de aguas más bajas)
anio_h <- function(f, m0) as.integer(format(f, "%Y")) - (as.integer(format(f, "%m")) < m0)

# ---------- verosimilitudes ----------------------------------------------------------------------------
# Poisson: cada evento suma a + b·x; cada mes válido resta los eventos esperados en ese mes, exp(a + b·x)/12.
nll_pois <- function(a, b, xe, xm) -(sum(a + b * xe) - sum(exp(a + b * xm)) / 12)
# Para una b dada, la mejor a sale cerrada: exp(a) = n / Σ_meses exp(b·x)/12.
a_hat <- function(b, xe, xm) log(max(length(xe), 0.5)) - log(sum(exp(b * xm)) / 12)
pnll_pois <- function(b, xe, xm) nll_pois(a_hat(b, xe, xm), b, xe, xm)

# GPD de los excedentes y, con escala exp(c + d·x) y forma ξ.
nll_gpd <- function(c, d, xi, y, xe) {
  if (xi < XI_MIN || xi > XI_MAX) return(FUERA)
  s <- exp(c + d * xe)
  z <- 1 + xi * y / s
  if (any(z <= 0)) return(FUERA)          # un excedente fuera del soporte (ξ < 0 pone un techo)
  if (abs(xi) < 1e-6) return(-sum(-log(s) - y / s))
  -sum(-log(s) - (1 + 1 / xi) * log(z))
}
nm <- function(par, f, maxit) optim(par, f, method = "Nelder-Mead", control = list(reltol = 1e-12, maxit = maxit))
fit_gpd_d <- function(d, y, xe) {          # c y ξ con d fija; tres arranques de ξ, como en el Python
  c0 <- log(mean(y)); best <- NULL
  for (x0 in c(0.1, -0.1, 0.3)) {
    r <- nm(c(c0, x0), function(p) nll_gpd(p[1], d, p[2], y, xe), 4000)
    if (is.null(best) || r$value < best$value) best <- r
  }
  list(c = best$par[1], d = d, xi = best$par[2], nll = best$value)
}
fit_gpd <- function(y, xe) {
  c0 <- log(mean(y)); best <- NULL
  for (x0 in c(0.1, -0.1, 0.3)) {
    r <- nm(c(c0, 0, x0), function(p) nll_gpd(p[1], p[2], p[3], y, xe), 6000)
    if (is.null(best) || r$value < best$value) best <- r
  }
  list(c = best$par[1], d = best$par[2], xi = best$par[3], nll = best$value)
}

# Intervalo del 95 % por verosimilitud perfilada: los valores v del parámetro donde 2·(nll(v) − nll mín.)
# cruza el χ² de 1 grado de libertad (3,84). Misma búsqueda que el Python.
perfil <- function(fun_fijo, est, nll_min, paso = 0.02, tope = 4) {
  g <- function(v) 2 * (fun_fijo(v) - nll_min) - CHI1
  sapply(c(-1, 1), function(sgn) {
    v0 <- est; v <- est + sgn * paso
    while (abs(v - est) < tope && g(v) < 0) { v0 <- v; v <- v + sgn * paso * 2 }
    if (g(v) >= 0 && g(v0) < 0) uniroot(g, sort(c(v0, v)), tol = 1e-9)$root else NA_real_
  })
}

ajustar <- function(D) {
  fp <- optimize(function(b) pnll_pois(b, D$xe, D$xm), c(-6, 6), tol = 1e-10)
  b <- fp$minimum; l1p <- fp$objective
  a0 <- a_hat(0, D$xe, D$xm); l0p <- nll_pois(a0, 0, D$xe, D$xm)
  g <- fit_gpd(D$y, D$xe); g0 <- fit_gpd_d(0, D$y, D$xe)
  list(a = a_hat(b, D$xe, D$xm), b = b, c = g$c, d = g$d, xi = g$xi, a0 = a0, c0 = g0$c, xi0 = g0$xi,
       l1 = -(l1p + g$nll), l0 = -(l0p + g0$nll), lb = -(l1p + g0$nll), ld = -(l0p + g$nll),
       nll_p = l1p, nll_g = g$nll)
}

prob <- function(z, u, a, b, c, d, xi, x) {
  lam <- exp(a + b * x); s <- exp(c + d * x); t <- 1 + xi * (z - u) / s
  S <- if (abs(xi) > 1e-6) ifelse(t > 0, pmax(t, 0)^(-1 / xi), 0) else exp(-(z - u) / s)
  1 - exp(-lam * S)
}

# Dejar una década afuera: ajustar sin ella y medir cuánto mejor predice la década excluida el modelo con El Niño
# que el modelo sin El Niño (Δ en log-verosimilitud). Una década con un pico sobre el techo de la GPD queda aparte.
lodo <- function(D, m0) {
  ye <- anio_h(D$ye, m0) %/% 10; ym <- anio_h(D$mm, m0) %/% 10
  out <- numeric(0)
  for (dec in sort(unique(ym))) {
    te <- ye != dec; tm <- ym != dec
    if (sum(te) < 10 || sum(!te) == 0) next
    f <- ajustar(list(xe = D$xe[te], y = D$y[te], xm = D$xm[tm]))
    xe <- D$xe[!te]; y <- D$y[!te]; xm <- D$xm[!tm]
    g1 <- nll_gpd(f$c, f$d, f$xi, y, xe); g0 <- nll_gpd(f$c0, 0, f$xi0, y, xe)
    if (g1 >= FUERA || g0 >= FUERA) { out <- c(out, NA); next }
    out <- c(out, (-nll_pois(f$a, f$b, xe, xm) - g1) - (-nll_pois(f$a0, 0, xe, xm) - g0))
  }
  out
}

# Bootstrap por años hidrológicos enteros: se sortean años con reposición y entran todos sus eventos y meses.
remuestrear <- function(idx_e, idx_m, smp) {
  k <- as.character(smp)
  list(ie = unlist(idx_e[k], use.names = FALSE), im = unlist(idx_m[k], use.names = FALSE))
}
bootstrap <- function(D, m0, n, niveles, xs) {
  ye <- anio_h(D$ye, m0); ym <- anio_h(D$mm, m0); anios <- sort(unique(ym))
  idx_e <- split(seq_along(ye), ye); idx_m <- split(seq_along(ym), ym)
  bs <- ds <- numeric(0); ps <- list()
  for (r in seq_len(n)) {
    s <- remuestrear(idx_e, idx_m, sample(anios, length(anios), replace = TRUE))
    if (length(s$ie) < 10) next
    f <- ajustar(list(xe = D$xe[s$ie], y = D$y[s$ie], xm = D$xm[s$im]))
    bs <- c(bs, f$b); ds <- c(ds, f$d)
    ps[[length(ps) + 1]] <- outer(niveles, xs, Vectorize(function(z, x) prob(z, D$u, f$a, f$b, f$c, f$d, f$xi, x)))
  }
  list(b = bs, d = ds, p = simplify2array(ps))     # p: niveles × clases × réplicas
}

# ---------- GEV de control (máximo anual con μ = μ0 + μ1·x) ----------------------------------------------
nll_gev <- function(m0, m1, ls, xi, y, x) {
  s <- exp(ls); t <- (y - (m0 + m1 * x)) / s
  if (abs(xi) < 1e-6) return(-sum(-ls - t - exp(-t)))
  z <- 1 + xi * t
  if (any(z <= 0)) return(FUERA)
  -sum(-ls - (1 + 1 / xi) * log(z) - z^(-1 / xi))
}
fit_gev <- function(y, x, m1_fijo = NULL) {
  sc <- sd(y) * sqrt(6) / pi; loc <- mean(y) - 0.5772 * sc; best <- NULL
  for (x0 in c(0.1, -0.1, -0.3, 0.3)) {
    if (is.null(m1_fijo)) {
      r <- nm(c(loc, 0, log(sc), x0), function(p) nll_gev(p[1], p[2], p[3], p[4], y, x), 8000)
      par <- r$par
    } else {
      r <- nm(c(loc, log(sc), x0), function(p) nll_gev(p[1], m1_fijo, p[2], p[3], y, x), 8000)
      par <- c(r$par[1], m1_fijo, r$par[2], r$par[3])
    }
    if (is.null(best) || r$value < best$nll) best <- list(par = par, nll = r$value)
  }
  best
}

# ---------- Uruguay en conjunto: a, c, ξ propios de cada estación; b y d compartidos --------------------------
ajustar_conjunto <- function(Ds) {
  pn <- function(b) sum(vapply(Ds, function(D) pnll_pois(b, D$xe, D$xm), 0))
  gn <- function(d) sum(vapply(Ds, function(D) fit_gpd_d(d, D$y, D$xe)$nll, 0))
  fp <- optimize(pn, c(-6, 6), tol = 1e-10)
  fg <- optimize(gn, c(-2, 2), tol = 1e-6)
  l0p <- pn(0); l0g <- gn(0)
  list(b = fp$minimum, d = fg$minimum, l1 = -(fp$objective + fg$objective), l0 = -(l0p + l0g),
       lb = -(fp$objective + l0g), ld = -(l0p + fg$objective))
}

# ---------- prueba sintética ------------------------------------------------------------------------------
prueba_sintetica <- function(n_rep = 200) {
  set.seed(7)
  x <- enso$x6[enso$mes >= as.Date("1950-01-01") & enso$mes < as.Date("2010-01-01")]
  x <- x[!is.na(x)]
  a <- 0; b <- 0.3; c <- log(2); d <- 0.15; xi <- 0.1
  cub_b <- cub_d <- 0
  for (r in seq_len(n_rep)) {
    k <- rpois(length(x), exp(a + b * x) / 12)
    xe <- rep(x, k)
    s <- exp(c + d * xe)
    y <- s * (runif(length(xe))^(-xi) - 1) / xi           # GPD por inversión
    fp <- optimize(function(v) pnll_pois(v, xe, x), c(-6, 6), tol = 1e-10)
    ib <- perfil(function(v) pnll_pois(v, xe, x), fp$minimum, fp$objective)
    g <- fit_gpd(y, xe)
    id <- perfil(function(v) fit_gpd_d(v, y, xe)$nll, g$d, g$nll)
    cub_b <- cub_b + isTRUE(ib[1] <= b && b <= ib[2])
    cub_d <- cub_d + isTRUE(id[1] <= d && d <= id[2])
  }
  cat(sprintf("prueba sintética: cobertura de b %.1f %%, de d %.1f %% (pasa si ≥ 90 %%)\n",
              100 * cub_b / n_rep, 100 * cub_d / n_rep))
  writeLines(sprintf("%s\ncobertura_b=%.3f\ncobertura_d=%.3f\npasa=%s", Sys.time(), cub_b / n_rep, cub_d / n_rep,
                     cub_b / n_rep >= 0.9 && cub_d / n_rep >= 0.9), "validation/prueba_sintetica.txt")
  cub_b / n_rep >= 0.9 && cub_d / n_rep >= 0.9
}
