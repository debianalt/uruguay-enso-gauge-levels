# ¿El Niño cambia la frecuencia y el tamaño de las crecidas del Uruguay?
#
# Reproducción en R de una implementación previa en Python (data/crecidas_ajuste.csv y
# data/crecidas_prob_condicional.csv son sus salidas). Las reglas de esta reproducción (tolerancias, alcance,
# implementación) se escribieron antes de correr.
#
# Modelo, en palabras:
#   - Un "evento" es una crecida que pasa el umbral u (el percentil 98 de todas las lecturas diarias del
#     hidrómetro). Los eventos ya están armados en data/crecidas_evento.csv; acá no se rearman.
#   - Cuántos eventos hay por año: un proceso de Poisson con tasa λ(x) = exp(a + b·x), donde x es el Niño 3.4
#     promediado en los 6 meses que terminan en el mes del pico (x6). b > 0 quiere decir más crecidas con El Niño:
#     exp(b) es cuántas veces se multiplica la tasa por cada °C.
#   - Cuánto pasa cada pico por encima de u (el excedente y): una distribución generalizada de Pareto (GPD) con
#     escala σ(x) = exp(c + d·x) y forma ξ. d > 0 quiere decir crecidas más altas con El Niño.
#   - Probabilidad de que en un año el río llegue al menos una vez a la cota z:
#     P = 1 − exp(−λ(x)·S(z − u)), con S la probabilidad de que un excedente pase z − u.
#
# Uso (desde la raíz del repositorio):
#   Rscript scripts/frecuencia.R --sintetica     prueba del código con datos simulados (debe dar ≥ 90 %)
#   Rscript scripts/frecuencia.R [--boot 500]    ajuste, tablas y validación cruzada contra el Python

args <- commandArgs(trailingOnly = TRUE)
N_BOOT <- if ("--boot" %in% args) as.integer(args[which(args == "--boot") + 1]) else 500L
SINTETICA <- "--sintetica" %in% args
SEMILLA <- if ("--semilla" %in% args) as.integer(args[which(args == "--semilla") + 1]) else 2026L
SALIDA <- if ("--semilla" %in% args) paste0("_semilla", SEMILLA) else ""   # corridas de control, aparte

source("scripts/modelo.R", encoding = "UTF-8")   # datos, verosimilitudes, ajuste, bootstrap

dir.create("validation", showWarnings = FALSE); dir.create("tables", showWarnings = FALSE)
if (SINTETICA) quit(status = if (prueba_sintetica()) 0 else 1)

# ---------- principal -----------------------------------------------------------------------------------
x_ult <- tail(enso[!is.na(enso$x6), ], 1)
clase_ult <- sprintf("último dato (%s, x = %+.2f)", format(x_ult$mes, "%m-%Y"), x_ult$x6)
ERSST_DESDE <- as.Date("1948-07-01")

variantes <- function(sid) {
  sesgo <- serie[as.character(sid), "sesgo_pico_m"]; if (is.na(sesgo)) sesgo <- 1
  v <- list(list("principal", "x6", NULL, 0), list("mas1m", "x6", NULL, sesgo),
            list("ersstv6", "x6_ersst", ERSST_DESDE, 0),
            list("x3 (exploratoria)", "x3", NULL, 0), list("x12 (exploratoria)", "x12", NULL, 0))
  q <- quiebre[[as.character(sid)]]
  if (!is.null(q)) v[[length(v) + 1]] <- list("post_quiebre", "x6", q, 0)
  v
}

set.seed(SEMILLA)
filas <- list(); probs <- list()
for (sid in URUGUAY) {
  m0 <- serie[as.character(sid), "mes_inicio_hidro"]; nombre_s <- serie[as.character(sid), "nombre"]
  for (v in variantes(sid)) {
    var <- v[[1]]; xcol <- v[[2]]; desde <- v[[3]]; extra <- v[[4]]
    D <- datos(sid, xcol, desde, NULL, extra)
    n <- length(D$xe)
    if (n < 15) {
      filas[[length(filas) + 1]] <- data.frame(serie = as.character(sid), nombre = nombre_s, variante = var,
                                               modelo = "pot", n_eventos = n,
                                               veredicto = "menos de 15 eventos: no se ajusta")
      next
    }
    f <- ajustar(D)
    lrt <- 2 * (f$l1 - f$l0); p <- pchisq(lrt, 2, lower.tail = FALSE)
    ib <- perfil(function(b) pnll_pois(b, D$xe, D$xm), f$b, f$nll_p)
    id <- perfil(function(d) fit_gpd_d(d, D$y, D$xe)$nll, f$d, f$nll_g)
    dl_todo <- lodo(D, m0); dl <- dl_todo[!is.na(dl_todo)]
    explor <- grepl("exploratoria", var)
    gana <- p < 0.05 && length(dl) > 0 && sum(dl) > 0 && mean(dl > 0) >= 0.6
    fila <- data.frame(serie = as.character(sid), nombre = nombre_s, variante = var, modelo = "pot", n_eventos = n,
                       anios = length(unique(anio_h(D$mm, m0))), desde = as.character(min(D$mm)),
                       hasta = as.character(max(D$mm)), u_m = D$u, a = f$a, b = f$b, b_lo = ib[1], b_hi = ib[2],
                       b_boot_lo = NA, b_boot_hi = NA, c = f$c, d = f$d, d_lo = id[1], d_hi = id[2],
                       d_boot_lo = NA, d_boot_hi = NA, xi = f$xi, ll0 = f$l0, ll1 = f$l1, lrt = lrt, p = p,
                       p_b = pchisq(2 * (f$lb - f$l0), 1, lower.tail = FALSE),
                       p_d = pchisq(2 * (f$ld - f$l0), 1, lower.tail = FALSE),
                       aic0 = -2 * f$l0 + 6, aic1 = -2 * f$l1 + 10, lodo_suma = sum(dl), lodo_pos = sum(dl > 0),
                       lodo_n = length(dl), lodo_fuera = sum(is.na(dl_todo)),
                       veredicto = if (explor) NA else if (gana) "se gana su lugar" else "no se gana su lugar")
    # GEV de control, mismos años
    g <- ma[ma$serie == sid, ]
    if (!is.null(desde)) g <- g[g$fecha_max >= desde, ]
    g$mes <- primer_dia(g$fecha_max)
    g <- merge(g, enso[, c("mes", xcol)], by = "mes"); g <- g[!is.na(g[[xcol]]), ]
    if (nrow(g) >= 20) {
      yy <- g$max_m + extra; xx <- g[[xcol]]
      g1 <- fit_gev(yy, xx); g0 <- fit_gev(yy, xx, m1_fijo = 0)
      filas[[length(filas) + 1]] <- data.frame(serie = as.character(sid), nombre = nombre_s, variante = var,
        modelo = "gev", n_eventos = nrow(g), anios = nrow(g), a = g1$par[1], b = g1$par[2], c = g1$par[3],
        xi = g1$par[4], ll0 = -g0$nll, ll1 = -g1$nll, lrt = 2 * (g0$nll - g1$nll),
        p = pchisq(2 * (g0$nll - g1$nll), 1, lower.tail = FALSE), veredicto = "control")
    }
    if (var == "principal") {
      niv <- c(alerta = serie[as.character(sid), "alerta_m"], "evacuación" = serie[as.character(sid), "evacuacion_m"])
      niv <- niv[!is.na(niv) & niv > D$u]
      if (length(niv)) {
        xs <- c(CLASES, x_ult$x6); names(xs)[length(xs)] <- clase_ult
        bt <- bootstrap(D, m0, N_BOOT, niv, xs)
        fila$b_boot_lo <- quantile(bt$b, 0.05); fila$b_boot_hi <- quantile(bt$b, 0.95)
        fila$d_boot_lo <- quantile(bt$d, 0.05); fila$d_boot_hi <- quantile(bt$d, 0.95)
        for (i in seq_along(niv)) for (j in seq_along(xs)) {
          pij <- if (length(dim(bt$p)) == 3) bt$p[i, j, ] else bt$p[j, ]
          probs[[length(probs) + 1]] <- data.frame(serie = sid, nombre = nombre_s, nivel = names(niv)[i],
            z_m = niv[[i]], clase = names(xs)[j], x = xs[[j]],
            p = prob(niv[[i]], D$u, f$a, f$b, f$c, f$d, f$xi, xs[[j]]),
            p_lo = quantile(pij, 0.05), p_hi = quantile(pij, 0.95),
            p_estacionario = prob(niv[[i]], D$u, f$a0, 0, f$c0, 0, f$xi0, xs[[j]]),
            fuera_de_rango = xs[[j]] > max(D$xm))
        }
      }
    }
    filas[[length(filas) + 1]] <- fila
  }
  cat(sid, nombre_s, ": listo\n")
}

# conjunto
for (vc in list(list("principal", 0), list("mas1m", 1))) {
  Ds <- lapply(URUGUAY, function(s) datos(s, "x6", extra = vc[[2]]))
  Ds <- Ds[vapply(Ds, function(D) length(D$xe) >= 15, TRUE)]
  f <- ajustar_conjunto(Ds)
  lrt <- 2 * (f$l1 - f$l0)
  ye <- lapply(Ds, function(D) anio_h(D$ye, 3)); ym <- lapply(Ds, function(D) anio_h(D$mm, 3))
  anios <- sort(unique(unlist(ym)))
  idx <- lapply(seq_along(Ds), function(i) list(e = split(seq_along(ye[[i]]), ye[[i]]),
                                                 m = split(seq_along(ym[[i]]), ym[[i]])))
  nb <- if (vc[[1]] == "principal") N_BOOT else max(100L, N_BOOT %/% 5L)
  bs <- ds <- numeric(nb)
  for (r in seq_len(nb)) {
    smp <- sample(anios, length(anios), replace = TRUE)
    Db <- list()
    for (i in seq_along(Ds)) {
      s <- remuestrear(idx[[i]]$e, idx[[i]]$m, smp)
      if (length(s$ie) >= 5) Db[[length(Db) + 1]] <- list(xe = Ds[[i]]$xe[s$ie], y = Ds[[i]]$y[s$ie],
                                                          xm = Ds[[i]]$xm[s$im])
    }
    fb <- ajustar_conjunto(Db); bs[r] <- fb$b; ds[r] <- fb$d
  }
  q <- function(v, p) unname(quantile(v, p))
  gana <- pchisq(lrt, 2, lower.tail = FALSE) < 0.05 &&
    (q(bs, .05) > 0 || q(bs, .95) < 0 || q(ds, .05) > 0 || q(ds, .95) < 0)
  filas[[length(filas) + 1]] <- data.frame(serie = "uruguay_conjunto", nombre = sprintf("Uruguay, %d hidrómetros", length(Ds)),
    variante = vc[[1]], modelo = "pot", n_eventos = sum(vapply(Ds, function(D) length(D$xe), 0)), anios = length(anios),
    b = f$b, d = f$d, b_boot_lo = q(bs, .05), b_boot_hi = q(bs, .95), d_boot_lo = q(ds, .05), d_boot_hi = q(ds, .95),
    ll0 = f$l0, ll1 = f$l1, lrt = lrt, p = pchisq(lrt, 2, lower.tail = FALSE),
    p_b = pchisq(2 * (f$lb - f$l0), 1, lower.tail = FALSE), p_d = pchisq(2 * (f$ld - f$l0), 1, lower.tail = FALSE),
    veredicto = if (gana) "se gana su lugar" else "no se gana su lugar")
  cat(sprintf("Uruguay en conjunto (%s): b = %+.3f [%+.3f, %+.3f], d = %+.3f [%+.3f, %+.3f]\n", vc[[1]], f$b,
              q(bs, .05), q(bs, .95), f$d, q(ds, .05), q(ds, .95)))
}

rbind_todo <- function(l) {
  cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(d) { for (c in setdiff(cols, names(d))) d[[c]] <- NA; d[cols] }))
}
aj <- rbind_todo(filas); pr <- rbind_todo(probs)
rownames(aj) <- NULL; rownames(pr) <- NULL
if (SALIDA == "") {
  write.csv(aj, "tables/ajuste.csv", row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(pr, "tables/prob_condicional.csv", row.names = FALSE, fileEncoding = "UTF-8")
} else {
  write.csv(aj, paste0("validation/ajuste", SALIDA, ".csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(pr, paste0("validation/prob_condicional", SALIDA, ".csv"), row.names = FALSE, fileEncoding = "UTF-8")
}

# ---------- validación cruzada contra el Python (tolerancias del log del 9-nov-2026) -----------------------
py <- leer("ajuste"); py$serie <- as.character(py$serie)
pyp <- leer("prob_condicional")
chk <- list()
agregar <- function(clave, campo, r, p, tol) {
  dif <- if (is.numeric(r) && is.numeric(p)) r - p else NA
  vacio <- function(v) is.na(v) || identical(as.character(v), "")
  ok <- if (is.na(tol)) (identical(as.character(r), as.character(p)) || (vacio(r) && vacio(p))) else
    (is.na(r) && is.na(p)) || (!is.na(dif) && abs(dif) <= tol)
  fmt <- function(v) if (is.numeric(v)) as.character(signif(v, 6)) else as.character(v)
  chk[[length(chk) + 1]] <<- data.frame(clave = clave, campo = campo, r = fmt(r), python = fmt(p),
                                        dif = dif, tol = tol, ok = isTRUE(ok))
}
for (i in seq_len(nrow(aj))) {
  r <- aj[i, ]
  p <- py[py$serie == r$serie & py$variante == r$variante & py$modelo == r$modelo, ]
  clave <- paste(r$serie, r$variante, r$modelo, sep = " | ")
  if (nrow(p) != 1) { chk[[length(chk) + 1]] <- data.frame(clave = clave, campo = "fila", r = "", python = "falta",
                                                          dif = NA, tol = NA, ok = FALSE); next }
  agregar(clave, "n_eventos", r$n_eventos, p$n_eventos, 0)
  agregar(clave, "veredicto", r$veredicto, p$veredicto, NA)
  if (r$modelo == "gev") { agregar(clave, "mu1", r$b, p$b, 0.05); next }
  for (cc in c("b", "d")) agregar(clave, cc, r[[cc]], p[[cc]], 0.02)
  for (cc in c("b_lo", "b_hi", "d_lo", "d_hi", "b_boot_lo", "b_boot_hi", "d_boot_lo", "d_boot_hi"))
    agregar(clave, cc, r[[cc]], p[[cc]], 0.05)
  agregar(clave, "lodo_n", r$lodo_n, p$lodo_n, 0)
}
for (i in seq_len(nrow(pr))) {
  r <- pr[i, ]
  p <- pyp[pyp$serie == r$serie & pyp$nivel == r$nivel & pyp$clase == r$clase, ]
  clave <- paste(r$serie, r$nivel, r$clase, sep = " | ")
  if (nrow(p) != 1) { chk[[length(chk) + 1]] <- data.frame(clave = clave, campo = "fila", r = "", python = "falta",
                                                          dif = NA, tol = NA, ok = FALSE); next }
  agregar(clave, "p", r$p, p$p, 0.02); agregar(clave, "p_estacionario", r$p_estacionario, p$p_estacionario, 0.02)
  agregar(clave, "p_lo", r$p_lo, p$p_lo, 0.04); agregar(clave, "p_hi", r$p_hi, p$p_hi, 0.04)
}
cv <- do.call(rbind, chk)
write.csv(cv, paste0("validation/crossval_report", SALIDA, ".csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("validación cruzada: %d comparaciones, %d fuera de tolerancia\n", nrow(cv), sum(!cv$ok)))
if (any(!cv$ok)) print(cv[!cv$ok, c("clave", "campo", "r", "python", "dif")], row.names = FALSE)
